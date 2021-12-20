suppressMessages(library(tidyverse))
suppressMessages(library(GenomicRanges))
suppressMessages(library(Rsamtools))
suppressMessages(library(ggbio))
source("scripts/cigar_functions.R")

bam <- "/mnt/projects/research/21073N_SEQ418081722_mitopathy/Sample_21073LRa009/21073LRa009.aligned.bam"

target <- "#Chrom	Start	End	Strand	Gene	Disease	Repeat	Rpt_Start	Rpt_End Location	Normal	Expanded	Health condition
chr4	39285456	39368381	-	RFC1	ataxia	AAAAG	39348424	39348479	<5x	>10x	TLD/ALS
"

FLANKING_WIDTH = 50

# Load target file
t <- read_tsv(target,
    col_names = c('chr', 'start', 'end', 'strand', 'gene', 'disease','repeat_sequence','repeat_start', 'repeat_end', 'location', 'normal', 'expanded', 'health_condition'),
    col_types = 'fiifccciiccc',
    skip = 1) %>%
    filter(chr %in% c(paste0('chr', 1:22), 'chrX', 'chrY'))


gr_targ <- makeGRangesFromDataFrame(t)

gr_flanking <- t %>%
    mutate(flanking_start = repeat_start - FLANKING_WIDTH) %>%
    mutate(flanking_end = repeat_end + FLANKING_WIDTH) %>%
    select(chr, flanking_start, flanking_end) %>%
    makeGRangesFromDataFrame


param <- ScanBamParam(
        which=gr_flanking,
        what=c("qname", "rname", "strand", "pos", "qwidth", "cigar", "seq"),
        flag = scanBamFlag(isSupplementaryAlignment = FALSE),
        mapqFilter = 1)

b <- BamFile(bam, paste0(bam, ".bai"))
aln <- scanBam(b, param=param)[[1]]

# Get downstream seqs
classified_reads <- tibble("Read name" = character(), "Alignment Type" = factor(), "q_start" = numeric(), "q_end" = numeric(), seq = character())

for (i in seq_along(aln$rname)){
    f_start <- t$repeat_start - aln$pos[i] - FLANKING_WIDTH
    f_end <- t$repeat_end - aln$pos[i] + FLANKING_WIDTH
    end_pos <- get_end_pos(aln$cigar[i])
    
    print(paste(aln$rname[i], "flank_DS_start", f_start,  "flank_US_end", f_end))
    print(paste("Read pos end", end_pos$rpos_end, "Query pos end", end_pos$qpos_end))
    
    # Case 1: Full alignment, downstream and upstream flankin regions are present
    if (f_start > 0 & end_pos$rpos_end > f_end) {
        print(paste("Full alignment:", aln$qname[i]))
        q_start <- rpos_to_qpos(aln$cigar[i], f_start )
        q_end <- rpos_to_qpos(aln$cigar[i], f_end )
        seq = subseq(aln$seq[i], start = q_start, end = q_end)
        classified_reads <- add_row(classified_reads, "Read name" = aln$qname[i], "Alignment Type" = "Full", "q_start" = q_start, "q_end"= q_end, "seq"= as.character(seq))
    } 

    # Case 2: Only downstream flanking region present
    else if (f_start > 0 & f_end < end_pos$qpos_end) {
        print(paste("Downstram flanking alignment:", aln$qname[i]))
        q_start <- rpos_to_qpos(aln$cigar[i], f_start)
        q_end <- end_pos$qpos_end
        seq = subseq(aln$seq[i], start = q_start, end = q_end)
    classified_reads <- add_row(classified_reads, "Read name" = aln$qname[i], "Alignment Type" = "Downstream", "q_start" = q_start, "q_end"= q_end, "seq"= as.character(seq))
    }

    # Case 3: Only upstream flankin region is present
    else if (f_start < 0 & end_pos$qpos_end > f_end) {
        print(paste("Upstream flanking alignment:", aln$qname[i]))
        q_start <- 1
        q_end <- rpos_to_qpos(aln$cigar[i], f_end)
        seq = subseq(aln$seq[i], start = q_start, end = q_end)
    classified_reads <- add_row(classified_reads, "Read name" = aln$qname[i], "Alignment Type" = "Upstream", "q_start" = q_start, "q_end"= q_end, "seq"= as.character(seq))
    }

    # Case 4: other
    else {
        print(paste("Unrecognized flanking regions:", aln$qname[i]))
    }

}

# Export as text file
write_delim(classified_reads, "sample_classified_reads.tsv", delim="\t")

# Create multiple sequence alignment
library(msa)
repeat_msa <- msa(classified_reads$seq, method="ClustalOmega", type="dna")
msaPrettyPrint(repeat_msa, output="asis")
repeat_msa
View(classified_reads)

# Create length histogram

