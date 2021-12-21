'
Todo: Create docstring

Usage: 
    flankingStats.R [options] <alignment.bam> <target.bed>

Options:
    -o <out_p>, --output <out_p>     Output prefix [default: re]
    -f <int>, --flanking_width <int> Width of flanking sequence [default: 200]
    -h, --help                       Show this screen
' -> doc

if (exists('snakemake')) {
    bam <- snakemake@input[['bam']]
    target <- snakemake@params[['target']]
    output <- snakemake@output[['tsv']]
    FLANKING_WIDTH <- snakemake@params[['flanking_width']]

} else {
    library(docopt)
    args <- docopt(doc, version = 'Naval Fate 2.0')
    bam <- args$'<alignment.bam>'
    target <- args$'<target.bed>'
    out_p <- args$output
    FLANKING_WIDTH <- args$flanking_width
}

suppressMessages(library(tidyverse))
suppressMessages(library(GenomicRanges))
suppressMessages(library(Rsamtools))
suppressMessages(library(ggplot2))
suppressMessages(library(msa))


source("/mnt/users/ahgrosc1/pipeline/expander/scripts/cigarFunctions.R")

#bam <- "/mnt/projects/research/21073N_SEQ418081722_mitopathy/Sample_21073LRa009/21073LRa009.aligned.bam"
#target <- "#Chrom	Start	End	Strand	Gene	Disease	Repeat	Rpt_Start	Rpt_End Location	Normal	Expanded	Health condition
#chr4	39285456	39368381	-	RFC1	ataxia	AAAAG	39348424	39348479	<5x	>10x	TLD/ALS
#"
#FLANKING_WIDTH = 200
out_p <- str_replace(output, ".reads.tsv", "")

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

# Load Bam File
param <- ScanBamParam(
        which=gr_flanking,
        what=c("qname", "rname", "strand", "pos", "qwidth", "cigar", "seq"),
        flag = scanBamFlag(isSupplementaryAlignment = FALSE),
        mapqFilter = 1)

b <- BamFile(bam, paste0(bam, ".bai"))
aln <- scanBam(b, param=param)[[1]]

# TODO Loop over all target regions

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
            if (q_end - q_start > FLANKING_WIDTH){
                        classified_reads <- add_row(classified_reads, "Read name" = aln$qname[i], "Alignment Type" = "Full", "q_start" = q_start, "q_end"= q_end, "seq"= as.character(seq))
            }    } 

        # Case 2: Only downstream flanking region present
        else if (f_start > 0 & f_end < end_pos$qpos_end) {
            print(paste("Downstram flanking alignment:", aln$qname[i]))
            q_start <- rpos_to_qpos(aln$cigar[i], f_start)
            q_end <- end_pos$qpos_end
            seq = subseq(aln$seq[i], start = q_start, end = q_end)
            if (q_end - q_start > FLANKING_WIDTH){
                        classified_reads <- add_row(classified_reads, "Read name" = aln$qname[i], "Alignment Type" = "Downstream", "q_start" = q_start, "q_end"= q_end, "seq"= as.character(seq))
            }    }

        # Case 3: Only upstream flankin region is present
        else if (f_start < 0 & end_pos$qpos_end > f_end ) {
            print(paste("Upstream flanking alignment:", aln$qname[i]))
            q_start <- 1
            q_end <- rpos_to_qpos(aln$cigar[i], f_end)
            seq = subseq(aln$seq[i], start = q_start, end = q_end)

            if (q_end - q_start > FLANKING_WIDTH){
                        classified_reads <- add_row(classified_reads, "Read name" = aln$qname[i], "Alignment Type" = "Upstream", "q_start" = q_start, "q_end"= q_end, "seq"= as.character(seq))
            }
        }

        # Case 4: other
        else {
            print(paste("Unrecognized flanking regions:", aln$qname[i]))
        }

    }
    classified_reads <- classified_reads %>% 
        mutate(repeat_length = nchar(seq) - 2*FLANKING_WIDTH)

    # Export as text file
    write_delim(classified_reads, paste0(out_p, ".reads.tsv"), delim="\t")

    # Create multiple sequence alignment
    read_set <- as(classified_reads$seq, "BStringSet")
    names(read_set) <- classified_reads$'Read name'
    repeat_msa <- msa(classified_reads$seq, method="ClustalOmega", type="dna")
    repeat_msa_set <- as(repeat_msa, "BStringSet")

    # Export multiple sequence alignments
    saveWidth <- getOption("width")
    options(width=200)
    sink(paste0(out_p, ".msa.txt"))
    print(repeat_msa, show ="complete", halfNrow = -1)
    sink()
    options(width=saveWidth)
    writeXStringSet(repeat_msa_set, file=paste0(out_p, ".msa.fasta"))

    # Create length histogram
    p1 <- ggplot(classified_reads) +
        geom_histogram(aes(x = repeat_length, fill = `Alignment Type`), col = "grey20") + 
        theme_classic() +
        scale_fill_brewer(palette = "Pastel2") 

    ggsave(paste0(out_p, ".lengths.png"), plot=p1)


    # Plot schematic representation
    # Median length of full alignments
    median_full_length = classified_reads %>%
        filter(`Alignment Type` == "Full")%>%
        summarise(m =median(repeat_length)) %>%
        pull(m)

    dt <- classified_reads %>%
        mutate(repeat_sequence = strsplit(seq, split="")) %>%
        unnest(repeat_sequence) %>%
        group_by(`Read name`) %>%
        mutate(x = row_number()) %>%
        ungroup() %>%
        rowwise %>%
        mutate(x = ifelse(`Alignment Type` == "Upstream", median_full_length - (repeat_length  - x),x ))

    dt$'Read name' <- as.factor(dt$'Read name')

    p2 <- ggplot(dt, aes(x = x, y = `Read name`, fill = repeat_sequence)) +
        geom_hline(aes(yintercept=`Read name`, col = `Alignment Type`)) +
        geom_raster() + 
        geom_vline(xintercept = 1, col = "firebrick") + 
        geom_vline(xintercept = FLANKING_WIDTH, col = "firebrick") + 
        geom_vline(xintercept = median_full_length + FLANKING_WIDTH, col = "firebrick") + 
        geom_vline(xintercept = median_full_length + 2*FLANKING_WIDTH, col = "firebrick") + 
        scale_fill_brewer(palette = "Set2") +
        scale_colour_brewer(palette = "Pastel2") +
        theme_classic()
    ggsave(paste0(out_p, ".sequences.png"), plot=p2, width = 25)

