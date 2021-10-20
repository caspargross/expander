import csv
import os
import pathlib
from glob import glob

# Load config
configfile: srcdir("config.yml")

# Parse sample targets FOFN
samples = set()
targets = dict()
if  os.path.isfile(config['sample_targets']):
    with open(config['sample_targets']) as f:
        reader=csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) == 2:
                samples.add(row[0])
                targets[row[0]] = row[1]
            else :
                print("Wrong number of columns in row: ", row, len(row))
                pass

# If no FOFN available, check dir for folders with "Sample_" prefix
else:
    wc, = glob_wildcards('Sample_{sample}/')
    samples =  [sample for sample in wc if '/' not in sample]
    for sample in samples:
        targets[sample]=False

print(samples)
print(targets)

rule all:
    input:
        expand("Sample_{sample}/{sample}.aligned.bam", sample = samples),
        expand("Sample_{sample}/{sample}_coverage.pdf", sample = samples),
        expand("Sample_{sample}/{sample}_on_target_count.tsv", sample = samples),
        expand("Sample_{sample}/{sample}.phased.csv", sample = samples)

def get_demuxed_files(wc):
    pb_xml = glob('Sample_'+wc.sample+'/*.consensusreadset.xml')
    return(pb_xml)
    
rule mapping:
    input:
        get_demuxed_files
    output:
        "Sample_{sample}/{sample}.aligned.bam"
    conda:
        "env/pb_tools.yml"
    threads: 8
    params:
        ref = config['ref_genome'],
        # Following parameterset is taken from PacBio Repeat Expansion Scripts
        r = config['mm2.r'],     # Bandwidth for inital chaining and base alignment extension [500]
        A = config['mm2.A'],     # Matching score [2]
        B = config['mm2.B'],     # Mismatching penalty [4]
        z = config['mm2.z'],     # Truncate alignment if diagonal score in DP matrix drops below z [400]
        Z = config['mm2.Z'],     # Not documentend in minimap2 man
        O = config['mm2.O'],     # Gap open penalty [4]
        e = config['mm2.e'],     # Sample a high-frequency minimizer every INT basepairs [500]
        E = config['mm2.E'],     # Gap extension penalty [2]
        Y = "-Y" if config['mm2.Y'] else ""  # -Y Enable softclipping  REMOVED
    shell:
        """
        pbmm2 align  \
            -j {threads} \
            -r {params.r} \
            -A {params.A} \
            -B {params.B} \
            -z {params.z} \
            -Z {params.Z} \
            -O {params.O} \
            -e {params.e} \
            -E {params.E} \
            --sort \
            --sample {wildcards.sample} \
            {input} {params.ref} {output}
        samtools index {output}
        """
    
rule mosdepth_reads:
    input:
        rules.mapping.output
    output:
        "Sample_{sample}/coverage/{sample}.per-base.bed.gz"
    conda:
        "env/pb_tools.yml"
    params:
        target = lambda wildcards: "" if not targets[wildcards.sample] else "--by " + targets[wildcards.sample]
    threads:
        3
    shell:
        """
        mosdepth \
            {params.target} \
            --threads {threads} \
            "Sample_{wildcards.sample}/coverage/{wildcards.sample}" \
            {input}
        """

rule plot_coverage:
    input:
        cov = rules.mosdepth_reads.output
    output:
        plot = "Sample_{sample}/{sample}_coverage.pdf"
    conda:
        "env/R.yml"
    params:
        target = lambda wildcards: targets[wildcards.sample]
    threads:
        1
    script:
        "scripts/plotCoverage.R"

rule count_on_target:
    input:
        bam = rules.mapping.output
    output:
        out = "Sample_{sample}/{sample}_on_target_count.tsv"
    conda:
        "env/R.yml"
    params:
        target = lambda wildcards: targets[wildcards.sample],
        filter = ""
    threads:
        1
    script:
        "scripts/countOnTarget.R"

rule bam_to_fasta:
# Basic Bam to Fasta conversion. NO reverse complement for reads on negative Strand!
    input:
        bam = rules.mapping.output
    output:
        "Sample_{sample}/{sample}.onTarget.fasta"
    params:
        target = lambda wildcards: "" if not targets[wildcards.sample] else "-L" + targets[wildcards.sample]
    threads:
        1
    conda:
        "env/pb_tools.yml"
    shell:
        """
        samtools view {params.target} {input} | \
        awk '{{OFS="\\t"; print ">"$1"\\n"$10}}' > {output} 
        """

rule trf:
    input: 
        "Sample_{sample}/{sample}.onTarget.fasta"
    output:
        "Sample_{sample}/{sample}.trf.dat"
    params:
        match = config['trf.match'],
        mismatch = config['trf.mismatch'],
        delta = config['trf.delta'],
        pm = config['trf.pm'],
        pi = config['trf.pi'],
        minscore = config['trf.minscore'],
        maxperiod = config['trf.maxperiod']
    threads: 1
    conda:
        "env/trf.yml"
    shell:
        """
        trf {input} \
            {params.match} \
            {params.mismatch} \
            {params.delta} \
            {params.pm} \
            {params.pi} \
            {params.minscore} \
            {params.maxperiod} \
            -d -h -ngs > {output}
        """

rule parse_trf:
    input:
        trf = rules.trf.output,
        bam = "Sample_{sample}/{sample}.aligned.bam"
    output:
        csv = "Sample_{sample}/{sample}.repeats_raw.csv"
    conda:
        "env/R.yml"
    params:
        target = lambda wildcards: targets[wildcards.sample],
        min_period = config['f.min_period'],
        max_period = config['f.max_period'],
        min_copies = config['f.min_copies'],
        min_match  = config['f.min_match'],
        max_indels = config['f.max_indels'],
        min_score  = config['f.min_score'],
        min_entropy= config['f.min_entropy']
    script:
        "scripts/trfOutputParser.R"

rule phaseRepeats:
    input:
        csv = rules.parse_trf.output.csv
    output:
        csv_phased = "Sample_{sample}/{sample}.phased.csv"
    conda:
        "env/R.yml"
    params:
        out_prefix = "Sample_{sample}/{sample}",
        ch = config['clust.height'],
        ach = config['clust.allele_height'],
        plot_clusters = False
    script:
        "scripts/phaseRepeats.R"
    
rule plotRepeats:
    input:
        csv_phased = rules.phaseRepeats.output.csv_phased
    output:
        plots = directory("Sample_{sample}/plots"),
        csv = "Sample_{sample}/{sample}_repeat_summary.csv"
    conda:
        "env/R.yml"
    params:
        out_prefix = "Sample_{sample}/{sample}"
    script:
        "scripts/plotRepeats.R"