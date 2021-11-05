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

# Calculate plot limits for nice and uniform display
WMIN <- 50
PADDING <- 5
xmin <- min(dt$repeat_length) - PADDING
xmax <- max(dt$repeat_length) + PADDING
if ((xmax - xmin) < WMIN) {
    delta <- floor(WMIN - (xmax - xmin))
    xmin <- xmin - floor(delta/2)
    xmax <- xmax + floor(delta/2)
}

# Plot repeat length histogram
title = dt$locus_id[[1]]
p_rl <- 
    ggplot(dt, aes(x = repeat_length, fill = as.factor(allele))) +
    geom_histogram(color = "grey20", binwidth = 1) +
    ggtitle(title) +
    theme_bw() +
    xlim(xmin, xmax) + 
    scale_fill_brewer(palette = "Set3") +
    labs(x="Motif repeats", y="Read count", fill = "Allele")

ggsave(paste(out_p,dt$locus_id[[1]],"repeat_lengths",ext, sep="."), plot = p_rl)
#if (png) ggsave(paste0(out_p,".repeat_lengths.png"), p = p_rl)
#if (jpg) ggsave(paste0(out_p,".repeat_lengths.jpg"), p = p_rl)
return(dt)
} 

plot_waterfall <- function(dt) {
# Create waterfall plot
dts <- 
    dt %>% 
    arrange(desc(repeat_start_genomic + repeat_length)) %>%
    mutate(seq_name = factor(seq_name, unique(seq_name))) %>%
    mutate(repeat_matches = str_locate_all(repeat_sequence, as.character(consensus_sequence))) %>%
    unnest(repeat_matches) %>%
    mutate(x = ((repeat_matches[,1] + repeat_matches[,2])/2) + repeat_start_genomic) %>%
    rowwise() %>%
    mutate(w = nchar(as.character(consensus_sequence)))

# Add sequence chars to the plot
dtt <- 
    dt %>% 
    mutate(repeat_sequence = strsplit(repeat_sequence, split="")) %>%
    unnest(repeat_sequence) %>%
    group_by(seq_name) %>%
    mutate(x = row_number() + repeat_start_genomic) %>%
    summarise(x, repeat_sequence, allele, consensus_sequence)  %>% as.data.frame()

set.seed(123)
mypal <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set3"))

p_wf <- 
    ggplot() +
    geom_text(data = dtt, 
        aes(
            y = seq_name,
            x=x,
            label = repeat_sequence)
        ) +
    geom_tile(data=dts, 
        aes(
            y = seq_name,
            x=x,
            fill = consensus_sequence,
            width =w
        ),
        color = "grey20",
        height = 0.9,
        size = 0.2) +
    labs(x = "genomic position", y = "reads") +
    ggtitle(dts$locus_id[[1]]) + 
    facet_grid(rows = vars(allele), space = "free_y", scales = "free_y") +
    scale_fill_discrete(type = mypal(length(levels(dts$consensus_sequence))+2)) +
    theme_bw() +
    theme(
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position="bottom")

ggsave(paste(out_p,dt$locus_id[[1]],"waterfall",ext, sep="."),
    width = 2500,
    height = min(300 + nrow(dt)*70, 5000),
    units = "px",
    scale = 1.2,
    plot = p_wf)
return(dt)
}

dt_phased <- read_csv(repeats, col_types = "fcfnfcnnfcc")

dt_phased %>%
    group_by(locus) %>%
    group_modify(~plot_repeat_lengths(.x)) %>%
    group_modify(~plot_waterfall(.x))
