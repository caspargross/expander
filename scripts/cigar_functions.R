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

get_pos <- function(cigar, pos){
# Calculates query and reference coordinates for a sequence length pos.
# Result overhangs the correct position, needs to be substracted in next step
    c <- parse_cigar(cigar)
    i = 1
    iq = 0
    ir = 0
    while (ir <= pos) {
        num <- c$num[i]
        op  <- c$op[i]
        iq <- iq + op_consumes[[op]]['q']*num
        ir <- ir + op_consumes[[op]]['r']*num
        cat(paste("i:", i, "iq:", iq, "ir:", ir, sep="\t"), "\n")
        i <- i+1
    }
    return(c(iq, ir))
}

rpos_to_qpos <- function(cigar, rpos){   
# Calculates query positions for a given reference position and cigar string 
    pos <- get_pos(cigar, rpos)
    # Substract overshoot
    qpos <- pos["q"] - (pos["r"] - rpos)
    return(qpos)
}

qpos_to_rpos <- function(cigar, qpos){ 
# Calculates referenfe position for a given query position and cigar string
    pos <- get_pos(cigar, qpos)
    # Substract overshoot
    rpos <- pos["r"] - (pos["q"] - qpos)
    return(rpos)
}

rpos_to_qpos("76=1I99=1X1=1X82=1I71=1X45=2I20=1I21=1I36=1X6=1I25", 300)
qpos_to_rpos("76=1I99=1X1=1X82=1I71=1X45=2I20=1I21=1I36=1X6=1I25", 300)
