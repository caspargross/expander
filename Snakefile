ref_genome = "/mnt/share/data/genomes/GRCh38.fa"

WC = glob_wildcards('Sample_{sample}/demultiplex.{barcode1}--{barcode2}.bam')

#print(WC)

targets = { '21020Ia049_ataxia_PacBio' : 'targets/ataxia_panel.bed',
    '21020Ia049_ATXN3_PacBio' : 'targets/ATXN3.bed',
    '21020Ia048_ATXN3_PacBio' : 'targets/ATXN3.bed',
    '21020Ia049_ATXN3_IDT' : 'targets/ATXN3.bed',
    '21020Ia044_ataxia_PacBio' : 'targets/ataxia_panel.bed',
    '21020Ia048_ATXN3_IDT' : 'targets/ATXN3.bed',
    '21020Ia044_ATXN3_PacBio' : 'targets/ATXN3.bed',
    '21020Ia044_ATXN3_IDT' : 'targets/ATXN3.bed',
    '21020Ia048_ATXN3_PacBio_low' : 'targets/ATXN3.bed'}

rule all:
    input:
        expand("Sample_{sample}/{sample}.aligned.bam", sample = WC.sample),
        expand("Sample_{sample}/{sample}_coverage.pdf", sample = WC.sample),
        expand("Sample_{sample}/RepeatExpansionTool", sample = WC.sample),
        expand("Sample_{sample}/{sample}_on_target_count.csv", sample = WC.sample),
        expand("Sample_{sample}/{sample}.repeats.csv", sample = WC.sample)

def get_demuxed_files(wc):
    i = WC.sample.index(wc.sample)
    bc1 = WC.barcode1[i]
    bc2 = WC.barcode2[i]
    return("Sample_"+ wc.sample + "/demultiplex."+bc1+"--"+bc2+".consensusreadset.xml")
    
rule mapping:
    input:
        get_demuxed_files
    output:
        "Sample_{sample}/{sample}.aligned.bam"
    conda:
        "environment.yml"
    threads: 8
    params:
        ref = ref_genome,
        # Following parameterset is taken from PacBio Repeat Expansion Scripts
        r = 10000,    # Bandwidth for inital chaining and base alignment extension [500]
        A = 2,        # Matching score [2]
        B = 5,        # Mismatching penalty [4]
        z = 400,      # Truncate alignment if diagonal score in DP matrix drops below z [400]
        Z = 50,       # Not documentend in minimap2 man
        O = 56,       # Gap open penalty [4]
        e = 4,        # Sample a high-frequency minimizer every INT basepairs [500]
        E = 0         # Gap extension penalty [2]
                      # -Y Enable softclipping  
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
            -Y \
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
        "environment.yml"
    params:
        target = lambda wildcards: targets[wildcards.sample]
    threads:
        3
    shell:
        """
        mosdepth \
            --by {params.target} \
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
        "environment_R.yml"
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
        out = "Sample_{sample}/{sample}_on_target_count.csv"
    conda:
        "environment_R.yml"
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
        target = lambda wildcards: targets[wildcards.sample]
    threads:
        1
    shell:
        """
        samtools view -L {params.target} {input} | \
        awk '{{OFS="\\t"; print ">"$1"\\n"$10}}' > {output} 
        """

rule trf:
    input: 
        "Sample_{sample}/{sample}.onTarget.fasta"
    output:
        "Sample_{sample}/{sample}.trf.dat"
    params:
        match = 2,
        mismatch = 7,
        delta = 7,
        pm = 80,
        pi = 10,
        minscore = 50,
        maxperiod = 10
    threads: 1
    conda:
        "environment_trf.yml"
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
        "environment_R.yml"
    params:
        target = lambda wildcards: targets[wildcards.sample]
    script:
        "scripts/trfOutputParser.R"

rule phaseRepeats:
    input:
        csv = rules.parse_trf.output.csv
    output:
        csv_phased = "Sample_{sample}/{sample}.repeats.csv"
    conda:
        "environment_R.yml"
    params:
        out_prefix = "Sample_%s/%s" % wildcards.sample,
        cluster_height = config['cluster_height'],
        allele_cluster_height  =  config['allele_cluster_height']
    script:
        "scripts/phaseRepeats.R"
    
rule plotRepeaets:
    input:
        csv_phased = rules.phaseRepeats.output.csv_phased
    output:
        directory("Sample_{sample}/plots")
    conda:
        "environment_R.yml"
    params:
        out_prefix = "Sample_%s/plots/" % wildcards.sample,
        pdf = True,
        jpg = False
    script:
        "scripts/plotRepeats.R"
