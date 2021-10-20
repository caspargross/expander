#!/usr/bin/env Rscript

'
Parses the .dat output  files of TRF.
Repeats are then filtered and annotated with genomic coordinates from alignment.

Usage: 
    trfOutputParser.R [options] <trf.dat> <alignment.bam>

Options:
    -h, --help                              Show this screen
    -t <targ.bed>, --targets <targ.bed>     Target file with repeat expansion locations
    -o <out.csv>, --output <out.csv>        Output filename [default: trf_repeats.csv]
' -> doc

if (exists('snakemake')) {
    #print(snakemake)
    bam <- snakemake@input[['bam']]
    targets <- snakemake@params[['target']]
    trf <- snakemake@input[['trf']]
    filter <- snakemake@params[['filter']]
    outfile <- snakemake@output[['csv']]
    MIN_PERIOD_SIZE = snakemake@params[['min_period']]
    MAX_PERIOD_SIZE = snakemake@params[['max_period']]
    MIN_COPIES_ALIGNED = snakemake@params[['min_copies']]
    MIN_MATCHES_ADJACENT = snakemake@params[['min_match']]
    MAX_INDELS = snakemake@params[['max_indels']]
    MIN_ALIGNMENT_SCORE = snakemake@params[['min_score']]
    MIN_ENTROPY = snakemake@params[['min_entropy']]
} else {
    library(docopt)
    args <- docopt(doc, version = 'Naval Fate 2.0')
    bam <- args$'<alignment.bam>'
    targets <- args$targets
    trf <- args$trf
    filter <- args$filter
    outfile <- args$output
    MIN_PERIOD_SIZE = 3
    MAX_PERIOD_SIZE = 7
    MIN_COPIES_ALIGNED = 5
    MIN_MATCHES_ADJACENT = 50
    MAX_INDELS = 30
    MIN_ALIGNMENT_SCORE = 50
    MIN_ENTROPY = 0.5
}

suppressMessages(library(tidyverse))
suppressMessages(library(GenomicRanges))
suppressMessages(library(Rsamtools))
suppressMessages(library(bamsignals))

trf_colnames <- c(
    "repeat_start",
    "repeat_end",
    "period_size",
    "n_copies_aligned",
    "size_consensus_pattern",
    "percent_matches_adjacent_copies",
    "percent_indels_adjacent_copies",
    "alignment_score",
    "percent_C",
    "percent_G",
    "percent_T",
    "percent_A",
    "entropy_measure",
    "consensus_sequence",
    "repeat_sequence",
    "flanking_downstream",
    "flanking_upstream"
    )

parse_trf <- function(dat_file) {
# Parses the TRF .dat File line by line
    print(paste("Parsing TRF output:", dat_file))

    i <- 0
    current_part <- ""
    seq_name <- ""
    tbl_l <- list()
    con <- file(dat_file, "r")

    while (TRUE) {
        line <- readLines(con, n = 1)
        i <- i + 1
        if (i %% 1000 == 0) print(paste(i , "repeats analysed"))
        if (length(line) == 0) {
            read_part(current_part, seq_name)
            break
        }
        if (stringr::str_starts(line, "@")) {
            print(line)
            if (seq_name != "") tbl_l[[seq_name]] <- read_part(current_part, seq_name)
            seq_name <- gsub("@", "", line)
            current_part <- ""
            next
            }
        current_part <- paste(current_part, line, sep = "\n")
    }
    close(con)
    return(bind_rows(tbl_l))
}

read_part <- function(txt, read_name){
# Creates a tibble from a string block of repeats and applys filters
    if (txt == "\n") return(NULL)
    tbl <- read_delim(txt,
        delim = " ",
        col_names = trf_colnames,
        col_types = "nnnnnnnnnnnnncc") %>%
        mutate("seq_name" = read_name) %>%
        filter(period_size >= MIN_PERIOD_SIZE) %>%
        filter(period_size <= MAX_PERIOD_SIZE) %>%
        filter(n_copies_aligned >= MIN_COPIES_ALIGNED) %>%
        filter(percent_matches_adjacent_copies >= MIN_MATCHES_ADJACENT) %>%
        filter(percent_indels_adjacent_copies <= MAX_INDELS) %>%
        filter(alignment_score >= MIN_ALIGNMENT_SCORE) %>%
        filter(entropy_measure >= MIN_ENTROPY)

    return(tbl)
}

parse_bam <- function(bam, target=FALSE){
# Parses bam file to extract aligned positions
    if (target != FALSE){
        targ <- read_tsv(targets,
        col_names = c('chr', 'start', 'end', 'strand', 'gene', 'disease', 'repeat_start', 'repeat_end', 'location', 'normal', 'expanded', 'health_condition'),
        col_types = 'fiicccciiccc',
        skip = 1) %>%
        filter(chr %in% c(paste0('chr', 1:22), 'chrX', 'chrY'))
        gr_targ <- makeGRangesFromDataFrame(targ)
    } else {
       gr_targ = GRanges()
    }

    param <- ScanBamParam(which=gr_targ, what=c("qname", "rname", "strand", "pos", "qwidth", "cigar", "seq"))
    aln <- scanBam(bam, param=param)
    return(aln[[1]])
}

dat_trf <- parse_trf(trf)

b <- parse_bam(bam)
string_set <- b[7]$seq
dt_bam <- as_tibble(b[1:6])

# Merge datasets (Inner) 
dt <- inner_join(dat_trf, dt_bam, by = c("seq_name" = "qname"))

print(paste("Rows TRF: ", nrow(dat_trf)))
print(paste("Rows Reads: ", nrow(dt_bam)))
print(paste("Rows Merge: ", nrow(dt)))

write_csv(dt, file = outfile)