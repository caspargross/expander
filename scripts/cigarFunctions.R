# Consumption table based on
# https://github.com/samtools/htslib/blob/7ecf4e4153d3d3d2ec0adb24611369494a41427d/htslib/sam.h#L129-L143
op_consumes <- list(
    'M' = c("q"=1, "r"=1),
    'I' = c("q"=1, "r"=0),
    'D' = c("q"=0, "r"=1),
    'N' = c("q"=0, "r"=1),
    'S' = c("q"=1, "r"=0),
    'H' = c("q"=1, "r"=0),  #Changed!
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


get_end_pos <- function(cigar) {
    c <- parse_cigar(cigar)
    iq = 0
    ir = 0
    for (i in seq_along(c$num)) {
        num <- c$num[i]
        op  <- c$op[i]
        iq <- iq + op_consumes[[op]]['q']*num
        ir <- ir + op_consumes[[op]]['r']*num
        i <- i+1
    }
    return(list(rpos_end = ir, qpos_end = iq))
}


qpos_to_rpos <- function(cigar, qpos){ 
# Calculates reference position for a given query position and cigar string
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

rpos_to_qpos <- function(cigar, rpos){ 
# Calculates reference position for a given query position and cigar string
    c <- parse_cigar(cigar)
    i = 1
    iq = 0
    ir = 0
    while (ir <= rpos) {
        num <- c$num[i]
        op  <- c$op[i]
        iq <- iq + op_consumes[[op]]['q']*num
        ir <- ir + op_consumes[[op]]['r']*num
        #cat(paste("i:", i, "iq:", iq, "ir:", ir, sep="\t"), "\n")
        i <- i+1
    }
    # Substract overshoot
    overshoot <- ir - rpos
    qpos <- iq - op_consumes[[op]]['q']*overshoot
    return(unname(qpos))
}

