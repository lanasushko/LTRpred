#' @title Run LTR_HARVEST_parallel to predict putative LTR Retrotransposons
#' @description This function implements an interface between R and
#' the LTRharvest command line tool to predict putative LTR retrotransposons from R.
#' @param genome.file path to the genome file in \code{fasta} format.



# If we want harvest parameters to be modifiable, I should change the shujun's parallel script
# LTR_HARVEST_parallel_modified2 -seq $seq -threads 10 -size 1000000 -time 300 -try1 1

pathtoltrharves=''

LTRharvest <- function(seq,
                       size  = 1000000,
                       time        = 300,
                       threads     = 1) {



# To run like from bash 
    system(paste0(pathtoltrharvest,"LTR_HARVEST_parallel_modified2 -seq  ",seq,"-threads ",threads,"-size ",size,"-time", time,"-try1",1))
                       }