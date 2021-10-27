#!/usr/bin/env Rscript
'
Create repeat expansion plots. Takes data from phaseRepeats script as input.
Default output in .PDF format.

Usage: 
    plotRepeats.R [options] <phased_repeats.csv>

Options:
    -o <out_p>, --output <out_p>    Output prefix [default: re]
    -p, --png                       Create .png plots
    -j, --jpg                       Create .jpg plots
    -h, --help                      Show this screen
' -> doc

if (exists('snakemake')) {
    repeats <- snakemake@input[['csv_phased']]
    out_p <- snakemake@params[['out_prefix']]
    pdf <- as.logical(snakemake@output[['pdf']])
    jpg <- as.logical(snakemake@output[['jpg']])
} else {
    library(docopt)
    args <- docopt(doc, version = 'Naval Fate 2.0')
    repeats <- args$'<phased_repeats.csv>'
    out_p <- args$output
    pdf <- args$pdf
    jpg <- args$jpg
}

ext <- "pdf"

suppressMessages(library(tidyverse))
suppressMessages(library(dendextend))

dir.create(dirname(out_p))


plot_repeat_lengths <- function(dt) {
# Plot repeat length histogram
title = dt$locus_id[[1]]
p_rl <- ggplot(dt, aes(x = repeat_length, fill = allele )) +
    geom_histogram(aes(color = allele), color = "grey20", binwidth = 1) +
    ggtitle(title) +
    theme_bw() +
    labs(x="Motif repeats")

ggsave(paste(out_p,dt$locus_id[[1]],"repeat_lengths",ext, sep="."), plot = p_rl)
#if (png) ggsave(paste0(out_p,".repeat_lengths.png"), p = p_rl)
#if (jpg) ggsave(paste0(out_p,".repeat_lengths.jpg"), p = p_rl)
return(dt)
} 

plot_waterfall <- function(dts) {
# Create waterfall plot
dts <- dts %>% 
    arrange(desc(repeat_start_genomic + repeat_length)) %>%
    mutate(seq_name = factor(seq_name, unique(seq_name))) %>%
    mutate(repeat_matches = str_locate_all(repeat_sequence, consensus_sequence)) %>%
    unnest(repeat_matches) %>%
    mutate(x = ((repeat_matches[,1] + repeat_matches[,2])/2) + repeat_start_genomic) %>%
    rowwise() %>%
    mutate(w = nchar(consensus_sequence)) 

p_wf <- ggplot() +
    geom_tile(data = dts, aes(
            y = seq_name,
            height=0.9 ,
            x= x ,
            width = w,
            fill = consensus_sequence,
            ),
        color = "grey20",
        size = 0.2) +
    theme_bw() +
    ggtitle(dts$locus_id[[1]]) + 
    theme(axis.text.y = element_blank())

ggsave(paste(out_p,dts$locus_id[[1]],"waterfall",ext, sep="."), width = 15, plot = p_wf)
return(dts)
}

dt_phased <- read_csv(repeats)

dt_phased %>%
    group_by(locus) %>%
    group_modify(~plot_repeat_lengths(.x)) %>%
    group_modify(~plot_waterfall(.x))