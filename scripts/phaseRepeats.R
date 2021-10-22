
'
Cluster repeats and phase alleles.

Usage: 
    trfOutputParser.R [options] <phased_repeats.csv>

Options:
    -h, --help                              Show this screen
    -t <targ.bed>, --targets <targ.bed>     Target file used to annotate reads
    -o <out_p>, --output <out_p>            Output prefix [default: repeat_clusters]
    -c --cluster_height                     Cutoff height for loci clusterung [default: 300]
    -a --allele_cluster_height              Minimal cutoff height for allele clustering tree [default: 1]
    -s --sample                             Sample name [default: basename inputfile]
    -p --plot_allele_clusters               Add plot for separate allele clusters
' -> doc

if (exists('snakemake')) {
    #print(snakemake)
    repeats <- snakemake@input[['csv']]
    out_p <- snakemake@params[['out_prefix']]
    target <- snakemake@params[['targets']]
    CLUSTER_HEIGHT <- snakemake@params[['ch']]
    MIN_ALLELE_CLUSTER_HEIGHT <- snakemake@params[['ach']]
    sample <- snakemake@wildcards[['sample']]
    plot_allele_clusters <- snakemake@params[['plot_clusters']]

} else {
    library(docopt)
    args <- docopt(doc, version = 'Naval Fate 2.0')
    repeats <- args$'<phased_repeats.csv'
    out_p <- args$'output'
    target <- args$'targets'
    CLUSTER_HEIGHT <- args$'cluster_height'
    MIN_ALLELE_CLUSTER_HEIGHT <- args$'allele_cluster_height'
    sample <- tools::file_path_sans_ext(basename(repeats))
    plot_allele_clusters <- args$'plot_allele_clusters'
}

# ---------------------- #
#       Functions        #
# -----------------------#

# Consumption table based on
# https://github.com/samtools/htslib/blob/7ecf4e4153d3d3d2ec0adb24611369494a41427d/htslib/sam.h#L129-L143
op_consumes <- list(
    'M' = c("q"=1, "r"=1),
    'I' = c("q"=1, "r"=0),
    'D' = c("q"=0, "r"=1),
    'N' = c("q"=0, "r"=1),
    'S' = c("q"=1, "r"=0),
    'H' = c("q"=0, "r"=0),
    'P' = c("q"=0, "r"=0),
    '=' = c("q"=1, "r"=1),
    'X' = c("q"=1, "r"=1)
);

parse_cigar <- function(c) {
# Parses cigar string into two lists with operators and values
    num <- as.numeric(strsplit(c, "[MIDNSHP=X]")[[1]])
    op <- tail(strsplit(c, "[0-9]+")[[1]], n = -1)
    cigar <- list(op = op, num = num)
    return(cigar)
}

cluster_alleles <- function(dt) {
# Assign allelles
# Cluster by repeat length (n_copies_aligned and genomic repeat_end)
    tryCatch( {
        clusters <- dt %>% 
                select(n_copies_aligned, repeat_start_genomic + (repeat_end - repeat_start)) %>%
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

            if (plot_allele_clusters) {
                d <- as.dendrogram(clusters)
                pdf(paste0(out_p, dt$locus[[1]], ".allele_clusters.pdf"), width = 12, height = 7)
                dend %>%
                color_branches(h= CLUSTER_HEIGHT) %>%
                raise.dendrogram(1) %>%
                set("labels", "") %>% 
                plot(main = paste("Clustered reads for ", dt$locus[[1]]))
                dend %>% rect.dendrogram(h= CLUSTER_HEIGHT, cluster = dt$locus, border = 8, lty = 8)
                dev.off()
            }
    
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

qpos_to_rpos <- function(cigar, qpos){ 
# Calculates referenfe position for a given query position and cigar string
    c <- parse_cigar(cigar)
    i = 1
    iq = 0
    ir = 0
    while (iq <= qpos) {
        num <- c$num[i]
        op  <- c$op[i]
        iq <- iq + op_consumes[[op]]['q']*num
        ir <- ir + op_consumes[[op]]['r']*num
        #cat(paste("i:", i, "iq:", iq, "ir:", ir, sep="\t"), "\n")
        i <- i+1
    }
    # Substract overshoot
    overshoot <- iq - qpos
    rpos <- ir - op_consumes[[op]]['r']*overshoot
    return(unname(rpos))
}


# ---------------------- #
#         Program        #
# -----------------------#

suppressMessages(library(tidyverse))
suppressMessages(library(dendextend))
suppressMessages(library(cluster))

# Load reads
dt <- read_csv(repeats, col_types = 'iidddddddddddcccccffiic')
 
# Load target file
if (target != FALSE){
    targ <- read_tsv(target,
    col_names = c('chr', 'start', 'end', 'strand', 'gene', 'disease','repeat_sequence','repeat_start', 'repeat_end', 'location', 'normal', 'expanded', 'health_condition'),
    col_types = 'fiifccciiccc',
    skip = 1) %>%
    filter(chr %in% c(paste0('chr', 1:22), 'chrX', 'chrY'))
}

# Correct start positions
dt <- dt %>% 
    rowwise %>%
    mutate(`repeat_start_genomic` = (qpos_to_rpos(cigar, repeat_start) + pos))

# Cluster repeats
# Use "gower" clustering for mixed datatypes with factors and numeric values
# Use location information (Chromosome and genomci repeat start) only for loci
clust <- dt %>% 
    dplyr::select(rname, repeat_start_genomic) %>%
    cluster::daisy(., metric = "gower") %>%
    hclust()

# Assign repeat categories
dt <- dt %>%
    ungroup %>%
    mutate(`locus` = as.factor(cutree(clust, h = CLUSTER_HEIGHT)))
print(paste(length(levels(dt$locus)), "Repeat groups found"))

# Annotate with target file
dt %>%
    group_by(locus) %>%
    rowwise %>%
    mutate(gene = get_gene(rname, median(repeat_start_genomic))) %>%
    select(pos, rname, gene)

get_gene <- function (r, p){
# Extracts gene name from target bed for chr and position
    gene <- targ %>% 
        filter(as.character(chr) == as.character(r)) %>%
        filter(start < p) %>%
        filter(end > p) %>%
        pull(gene)

    gene <- ifelse(length(gene)>0, gene, "unassigned")
    return(gene)
}


# Create repeat locus id, either with target file or unique location/repeat name

# Plot cluster Dendrogram
dend <- as.dendrogram(clust)
pdf(paste0(out_p, ".repeat_clusters.pdf"), width = 12, height = 7)
dend %>%
    color_branches(h= CLUSTER_HEIGHT) %>%
#   raise.dendrogram(0.1) %>%
    set("labels", "") %>% 
    plot(main = paste("Loci", sample))
dend %>% rect.dendrogram(h= CLUSTER_HEIGHT, cluster = dt$locus, border = 8, lty = 8)
colored_bars(colors = as.numeric(factor(dt$rname)), dend = dend, rowLabels = "chr")
dev.off()

# Cluster alleles
dt_phased <- dt %>% 
    group_by(locus) %>%
    group_modify(~cluster_alleles(.x))

# Write output table
dt_phased %>% 
    select("rname", "repeat_start_genomic", "repeat_start", "locus", "n_copies_aligned", "consensus_sequence", "size_consensus_pattern", "allele", "seq_name", "repeat_sequence", "repeat_end") %>% 
    write_csv(paste0(out_p, ".phased.csv"))


