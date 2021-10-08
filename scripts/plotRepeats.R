#!/usr/bin/env Rscript
'
Create repeat expansion plots. Takes data from phaseRepeats script as input.
Default output in .PNG format.

Usage: 
    plotRepeats.R [options] <phased_repeats.csv>

Options:
    -o <out_p>, --output <out_p>    Output prefix [default: re]
    -p, --pdf                       Create >pdf plots
    -j, --jpg                       Create .jpg plots
    -h, --help                      Show this screen
' -> doc

if (exists('snakemake')) {
    #print(snakemake)
    repeats <- snakemake@input[['csv']]
    outfile <- snakemake@output[['csv']]
    outfile <- snakemake@output[['png']]
} else {
    library(docopt)
    args <- docopt(doc, version = 'Naval Fate 2.0')
    targets <- args$targets
    outfile <- args$output
}

suppressMessages(library(tidyverse))
suppressMessages(library(dendextend))

plot_repeat_lengths <- function(dt) {
# Plot repeat length histogram
title = paste(
    "Repeat: ", names(sort(table(dt$rname), decreasing = T)[1]),
    median(dt$pos+dt$repeat_start),
    names(sort(table(dt$consensus_sequence), decreasing = T)[1])
    ) 

p_rl <- ggplot(dt, aes(x = floor((repeat_end - repeat_start)/3), fill = allele )) +
    geom_histogram(aes(color = allele), color = "grey20", binwidth = 1) +
    ggtitle(title) 
    theme_bw() +
    labs(x="Motif repeats")

ggsave(paste0(SAMPLE,"_repeat_lengths.png"), p = p_rl)
if (pdf) ggsave(paste0(SAMPLE,"_repeat_lengths.pdf"), p = p_rl)
if (jpg) ggsave(paste0(SAMPLE,"_repeat_lengths.jpg"), p = p_rl)
} 

plot_waterfall <- function(dts) {
# Create waterfall plot
dts <- dt_phased %>% select(c("seq_name", "allele", "repeat_start", "repeat_end", "consensus_sequence", "size_consensus_pattern","repeat_sequence")) %>% 
    arrange(desc(repeat_end)) %>%
    mutate(seq_name = factor(seq_name, unique(seq_name))) %>%
    mutate(repeat_matches = str_locate_all(repeat_sequence, consensus_sequence)) %>%
    unnest(repeat_matches)

p_wf <- ggplot(dts) +
    geom_tile(aes(
        y = seq_name,
        height=0.9 ,
        x= (repeat_matches[,1] + repeat_matches[,2])/2,
        width = size_consensus_pattern,
        fill = consensus_sequence,),
        color = "grey20",
        size = 0.2) +
    theme_bw() +
    theme(
        axis.text.y = element_blank())

ggsave(paste0(SAMPLE,"_repeat_lengths.png"), width = 15, plot = p_wf)
if (pdf) ggsave(paste0(SAMPLE,"_repeat_lengths.png"), width = 15, plot = p_wf)
if (jpg) ggsave(paste0(SAMPLE,"_repeat_lengths.png"), width = 15, plot = p_wf)
}

dt_phased <- read_csv(
plot_repeat_lengths(dt_phased)
plot_waterfall(dt_phased)