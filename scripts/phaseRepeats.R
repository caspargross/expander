
'
Cluster repeats and phases alleles.

Usage: 
    trfOutputParser.R [options] <phased_repeats.csv>

Options:
    -h, --help                              Show this screen
    -o <out.png>, --output <out.csv>        Output filename [default: trf_repeats.csv]
    -c --cluster_height                     Cutoff height for loci clusterung [default: 300]
    -a --allele_cluster_height              Minimal cutoff height for allele clustering tree [default: 1]
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


# Constant
CLUSTER_HEIGHT = 300
MIN_ALLELE_CLUSTER_HEIGHT = 30

# Debug
repeats <- "/mnt/projects/research/21020I_PB_NoAmp/Sample_21020Ia049_ataxia_PacBio/21020Ia049_ataxia_PacBio.repeats.csv"
repeats <- "/mnt/projects/research/21020I_PB_NoAmp/Sample_21020Ia049_ATXN3_PacBio/21020Ia049_ATXN3_PacBio.repeats.csv"
dt <- read_csv(repeats) 
# Cluster repeats
clust <- dt %>% 
    select(repeat_start, repeat_end, n_copies_aligned, consensus_sequence, rname, pos) %>%
    dist() %>%
    hclust()

# Plot cluster Dendrogram
dend <- as.dendrogram(clust)
dend <- dend %>% color_branches(h= CLUSTER_HEIGHT)
pdf("repeat_clusters.pdf")
plot(dend)
dev.off()

# Assign repeat categories
dt <- dt %>% mutate(category = as.factor(cutree(clust, h = CLUSTER_HEIGHT)))
print(paste(length(levels(dt$category)), "Repeat groups found"))

# Assign allelles
cluster_alleles <- function(dt) {

    clusters <- dt %>% 
        select(repeat_start, repeat_end, n_copies_aligned, consensus_sequence, alignment_score) %>%
        dist %>%
        hclust

    h <- MIN_ALLELE_CLUSTER_HEIGHT
    n_clusters <- max(cutree(clusters, h=h))
    
    # Increase cutoff until a maximum of 2 cluster (= allelles) was found
    while (n_clusters > 2) {
        h <- h+(h*0.5)
        #print(paste("Cutting tree with h = ", h, "N clust", n_clusters))
        n_clusters <- max(cutree(clusters, h = h))
    }
    
    dt <- add_column(dt, allele = as.factor(cutree(clusters, h = h)))
    return(dt)
}


dt_phased <- dt %>% 
    group_by(category) %>%
    group_modify(~cluster_alleles(.x))