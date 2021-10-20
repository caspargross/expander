#!/usr/bin/env Rscript

'
Plot genomic coverage profiles from Mosdepth output

Usage: 
    plotCoverage.R [options] <sample.per-base.bed.gz>

Options:
    -h, --help                          Show this screen
    -t <bedfile>, --targets <bedfile>   Annotate target regions from Bedfile
    -o <outfile>, --output <outfile>    Output filename (*.png, *.pdf, *.jpg)
                                        [default: coverage_plot.pdf]
    -s <sample>, --sample <sample>      Samplename
' -> doc

if (exists('snakemake')) {
    #print(snakemake)
    cov.bed.gz <- snakemake@input[['cov']]
    targets <- snakemake@params[['target']]
    sample <- snakemake@wildcards[['sample']]
    outfile <- snakemake@output[['plot']]
} else {
    library(docopt)
    args <- docopt(doc, version = 'Naval Fate 2.0')
    #print(args)
    cov.bed.gz <- args$'<sample.per-base.bed.gz>'
    targets <- args$targets
    sample <- args$sample
    outfile <- args$output
}

#options(warn = -1)
suppressMessages(library(tidyverse))
suppressMessages(library(vroom))
suppressMessages(library(ggrepel))
cov <- vroom(cov.bed.gz,  
    col_types = c('f', 'i', 'i', 'i'), 
    col_names = c('chr', 'start', 'end', 'cov')) %>%
    filter(chr %in% c(paste0('chr', 1:22), 'chrX', 'chrY'))
cov_long <- cov %>%
    pivot_longer(c("start", "end"), names_to = "location", values_to = "position")

p1 <- ggplot() +
    geom_step(data= cov_long, aes(x=position, y = cov)) +
 #   geom_rect(data = targ, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), fill = "blue", alpha = 0.2) +
    facet_grid(cols = vars(chr), scales = "free_x", space = "free_x", switch = "x") + 
    theme_classic() + 
    theme(
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_text(size =7),
#        plot.margin = unit(c(0,0,0,0), "pt"),
        panel.spacing = unit(0, "pt"),
        ) +
    labs(
        x = element_blank(),
        y = paste(sample, "Read Cov") )

if (targets != FALSE) {
    print("Trying to load targets")

    targ <- read_tsv(targets,
        col_names = c('chr', 'start', 'end', 'strand', 'gene', 'disease', 'repeat_start', 'repeat_end', 'location', 'normal', 'expanded', 'health_condition'),
        col_types = 'fiicccciiccc',
        skip = 1) %>%
        filter(chr %in% c(paste0('chr', 1:22), 'chrX', 'chrY'))

    p1 <<- p1 + geom_text_repel(
        data = targ,
        aes(x = (start + end) / 2 ,  y = max(cov_long$cov), label = gene),
        size = 2, angle = 90, hjust = 1)

    print(targ)
}

ggsave(filename = outfile, plot = p1,  width = 15, height = 4)

#Chrom	Start	End	Strand	Gene	Disease	Repeat	Rpt_Start	Rpt_End Location	Normal	Expanded	Health condition