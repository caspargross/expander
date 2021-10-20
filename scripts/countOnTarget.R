#!/usr/bin/env Rscript

'
Count on target Reads for a Bamfile and a target BED File

Usage: 
    countOnTarget.R [options] <targets.bed> <alignment.bam>

Options:
    -h, --help                          Show this screen
    -o <outfile>, --output <outfile>    Output filename (*.png, *.pdf, *.jpg)
                                        [default: OnTargetCounts.csv]
    -f <BamFilter>                      Filter alignments (Samtools Hex Format)
' -> doc

if (exists('snakemake')) {
    #print(snakemake)
    alignment <- snakemake@input[['bam']]
    targets <- snakemake@params[['target']]
    filter <- snakemake@params[['filter']]
    outfile <- snakemake@output[['out']]
} else {
    library(docopt)
    args <- docopt(doc, version = 'Naval Fate 2.0')
    #print(args)
    alignment <- args$'<alignment.bam>'
    targets <- args$'<targets.bed>'
    filter <- args$filter
    outfile <- args$output
}

#options(warn = -1)
suppressMessages(library(tidyverse))
suppressMessages(library(GenomicRanges))
suppressMessages(library(Rsamtools))
suppressMessages(library(bamsignals))

# Load Bamfile
bam_alignment <- BamFile(alignment)

# Load Targets
targ <- read_tsv(targets,
    col_names = c('chr', 'start', 'end', 'strand', 'gene', 'disease', 'motif', 
    'repeat_start', 'repeat_end', 'normal_repeat_count', 'pathological_repeat_count', 'health_condition'),
    col_types = 'fiic') %>%
    filter(chr %in% c(paste0('chr', 1:22), 'chrX', 'chrY'))

gr_targ <- makeGRangesFromDataFrame(targ)
target_count <- bamCount(alignment, gr_targ)

bam_stats <- countBam(bam_alignment)
genome_length <- sum(seqinfo(bam_alignment)@seqlengths)

# Create table with enrichment output
out_table <- targ %>%
    mutate(location = paste0(chr, ':', start, '-', end, ':', strand)) %>%
    mutate(length = end - start) %>%
    mutate(read_count = bamCount(alignment, gr_targ)) %>%
    select(c(gene, location, length, read_count))

#out_table <- out_table %>%    
#    mutate(expected = ((read_count/length)/(bam_stats$records/genome_length))) %>%
#    mutate(enrichment = read_count/expected) 

out_table <- out_table %>%    
    mutate(target_read_density  = (read_count/length)) %>%
    mutate(all_read_density= bam_stats$records/genome_length) %>%
    mutate(enrichment_factor = target_read_density/all_read_density)

out_table <- out_table %>%
    add_row(gene = "#OFFTARGET", 
        location = "whole genome excluding targets", 
        length = genome_length - sum(out_table$length),
        read_count = bam_stats$records - sum(target_count)) %>%
    mutate(read_percentage = read_count *100 / bam_stats$records)

write_tsv(out_table %>% select("gene", "location", "length", "read_count", "enrichment_factor", "read_percentage"), outfile)
