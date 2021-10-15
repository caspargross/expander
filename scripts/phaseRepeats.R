
'
Cluster repeats and phase alleles.

Usage: 
    trfOutputParser.R [options] <phased_repeats.csv>

Options:
    -h, --help                              Show this screen
    -t <t.bed>, --targets <t.bed>           Target file used to annotate reads
    -o <out_p>, --output <out_p>            Output prefix [default: repeat_clusters]
    -c --cluster_height                     Cutoff height for loci clusterung [default: 300]
    -a --allele_cluster_height              Minimal cutoff height for allele clustering tree [default: 1]
    -s --sample                             Sample name [default: basename inputfile]
' -> doc

if (exists('snakemake')) {
    #print(snakemake)
    repeats <- snakemake@input[['csv']]
    out_p <- snakemake@params[['out_prefix']]
    targ <- snakemake@params[['targets']]
    CLUSTER_HEIGHT <- snakemake@params[['ch']]
    MIN_ALLELE_CLUSTER_HEIGHT <- snakemake@params[['ach']]
    sample <- snakemake@wildcards[['sample']]

} else {
    library(docopt)
    args <- docopt(doc, version = 'Naval Fate 2.0')
    repeats <- args$'<phased_repeats.csv'
    out_p <- args$'output'
    targ <- args$'targets'
    CLUSTER_HEIGHT <- args$'cluster_height'
    MIN_ALLELE_CLUSTER_HEIGHT <- args$'allele_cluster_height'
    sample <- tools::file_path_sans_ext(basename(repeats))
}

suppressMessages(library(tidyverse))
suppressMessages(library(dendextend))
#suppressMessages(library(stringdist))

dt <- read_csv(repeats) 

# Cluster repeats
clust <- dt %>% 
    dplyr::select(n_copies_aligned, consensus_sequence, rname, pos, repeat_start) %>%
    dist() %>%
    hclust()

# Assign repeat categories
dt <- dt %>% mutate(locus = as.factor(cutree(clust, h = CLUSTER_HEIGHT)))
print(paste(length(levels(dt$locus)), "Repeat groups found"))

# Plot cluster Dendrogram
dend <- as.dendrogram(clust)
pdf(paste0(out_p, ".repeat_clusters.pdf"), width = 12, height = 7)
dend %>%
    color_branches(h= CLUSTER_HEIGHT) %>%
    raise.dendrogram(1) %>%
    set("labels", "") %>% 
    plot(main = paste("Loci", sample))
dend %>% rect.dendrogram(h= CLUSTER_HEIGHT, cluster = dt$locus, border = 8, lty = 8)
colored_bars(colors = as.numeric(factor(dt$rname)), dend = dend, rowLabels = "chr")
dev.off()

# Assign allelles
cluster_alleles <- function(dt) {

    tryCatch( {
        clusters <- dt %>% 
                select(n_copies_aligned, consensus_sequence, repeat_end) %>%
#               stringdistmatrix %>%
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
    
        },
        error=function(cond) {
            message(paste("Cannot cluster for loci", dt$locus[[1]]))
            message("Here's the original error message:")
            message(cond)
            # Choose a return value in case of error
            dt <- add_column(dt, allele=1)
            return(dt)
        },
        finally={
            return(dt)
        }
    )
}

dt_phased <- dt %>% 
    group_by(locus) %>%
    group_modify(~cluster_alleles(.x))

# Write output table
dt_phased %>% 
    select("rname", "pos", "repeat_start", "locus", "n_copies_aligned", "consensus_sequence", "allele", "seq_name", "repeat_sequence", "repeat_end") %>% 
    write_csv(paste0(out_p, ".phased.csv"))

#TODO Annotate Repeats if target file exists
