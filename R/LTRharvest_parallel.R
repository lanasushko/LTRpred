#' @title Run LTR_HARVEST_parallel to predict putative LTR Retrotransposons
#' @description This function implements an interface between R and
#' the LTRharvest command line tool to predict putative LTR retrotransposons from R.
#' @param genome.file path to the genome file in \code{fasta} format.
#' @param index.file specify the name of the enhanced suffix array index file that is computed
#'  by \code{suffixerator}. This opten can be used in case the suffix file was previously 
#'  generated, e.g. during a previous call of this function. In this case the suffix array index
#'  file does not need to be re-computed for new analyses. This is particularly useful when 
#'  running \code{LTRharvest} with different parameter settings.
#' @param range define the genomic interval in which predicted LTR transposons shall be reported
#' . In case \code{range[1] = 1000} and \code{range[2] = 10000} then candidates are only 
#' reported if they start after position 1000 and end before position 10000 in their respective 
#' sequence coordinates. If \code{range[1] = 0} and \code{range[2] = 0}, 
#' so \code{range = c(0,0)} (default) then the entire genome is being scanned.
#' @param seed  the minimum length for the exact maximal repeats. Only repeats with the specified minimum length are considered in all subsequent analyses. Default is \code{seed = 30}.
#' @param minlenltr minimum LTR length. Default is \code{minlenltr = 100}. 
#' @param maxlenltr maximum LTR length. Default is \code{maxlenltr = 3500}.
#' @param mindistltr minimum distance of LTR starting positions. Default is \code{mindistltr = 4000}.
#' @param maxdistltr maximum distance of LTR starting positions. Default is \code{maxdistltr = 25000}.
#' @param similar minimum similarity value between the two LTRs in percent. \code{similar = 70}.
#' @param mintsd minimum target site duplications (TSDs) length. If no search for TSDs
#' shall be performed, then specify \code{mintsd = NULL}. Default is \code{mintsd = 4}.
#' @param maxtsd maximum target site duplications (TSDs) length. If no search for TSDs
#' shall be performed, then specify \code{maxtsd = NULL}. Default is \code{maxtsd = 20}.
#' @param vic number of nucleotide positions left and right (the vicinity) of the predicted
#'  boundary of a LTR that will be searched for TSDs and/or one motif (if specified). 
#'  Default is \code{vic = 60}.
#' @param overlaps specify how overlapping LTR retrotransposon predictions shall be treated. 
#' If \code{overlaps = "no"} is selected, then neither nested nor overlapping predictions will be reported in the output. In case \code{overlaps = "best"} is selected then in the case of two or more nested or overlapping predictions, solely the LTR retrotransposon prediction with
#' the highest similarity between its LTRs will be reported.
#' If \code{overlaps = "all"} is selected then all LTR retrotransposon predictions 
#' will be reported whether there are nested and/or overlapping predictions or not. 
#' Default is \code{overlaps = "best"}.
#' @param xdrop specify the xdrop value (> 0) for extending a seed repeat in both directions
#'  allowing for matches, mismatches, insertions, and deletions. The xdrop extension process
#'   stops as soon as the extension involving matches, mismatches, insersions, and deletions 
#'   has a score smaller than T -X, where T denotes the largest score seen so far. Default is \code{cdrop = 5}.
#' @param mat specify the positive match score for the X-drop extension process. Default is \code{mat = 2}.
#' @param mis specify the negative mismatch score for the X-drop extension process. Default is \code{mis = -2}.
#' @param ins specify the negative insertion score for the X-drop extension process. Default is \code{ins = -3}.
#' @param del specify the negative deletion score for the X-drop extension process. Default is \code{del = -3}.
#' @param motif specify 2 nucleotides for the starting motif and 2 nucleotides for the ending
#'  motif at the beginning and the ending of each LTR, respectively.
#'  Only palindromic motif sequences - where the motif sequence is equal to its complementary
#'  sequence read backwards - are allowed, e.g. \code{motif = "tgca"}. Type the nucleotides without any space
#'  separating them. If this option is not selected by the user, candidate pairs will not be
#'  screened for potential motifs. If this options is set but no allowed number of
#'  mismatches is specified by the argument \code{motifmis} and a search for the exact 
#'  motif will be conducted. If \code{motif = NULL} then no explicit motif is being specified.
#' @param motifmis allowed number of mismatches in the TSD motif specified in \code{motif}. 
#' The number of mismatches needs to be between [0,3].  Default is \code{motifmis = 0}.
#' @param output.path a path/folder to store all results returned by \code{LTRharvest}. 
#' If \code{output.path = NULL} (Default) then a folder with the name of the input genome file
#' will be generated in the current working directory of R and all results are then stored in this folder.
#' @param verbose logical value indicating whether or not detailed information shall be printed on the console.
#' @param threads specify the number of threads to be used in parallel.
#' @param cut.size specify the size in bp of the pieces that the sequence will be cut in (used for running LTRharvest in parallel).
#' @param time set maximum time for a thread to run. After [time] seconds, the child thread is killed.

# If we want harvest parameters to be modifiable, I should change the shujun's parallel script
# LTR_HARVEST_parallel_modified2 -seq $seq -threads 10 -size 1000000 -time 300 -try1 1

pathtoltrharvest='/tmp/global2/ssushko/intactLTR_treeoflife/dev_LTRpred/LTRpred2/bin/LTR_HARVEST_parallel/'

LTRharvest_par <- function(genome.file,
                       index.file  = NULL,
                       range       = c(0,0),
                       seed        = 30,
                       minlenltr   = 100,
                       maxlenltr   = 3500,
                       mindistltr  = 4000,
                       maxdistltr  = 25000,
                       similar     = 70,
                       mintsd      = 4,
                       maxtsd      = 20,
                       vic         = 60,
                       overlaps    = "no",
                       xdrop       = 5,
                       mat         = 2,
                       mis         = -2,
                       ins         = -3,
                       del         = -3,
                       motif       = NULL,
                       motifmis    = 0,
                       output.path = NULL,
                       verbose     = TRUE,
                       threads     = 1,
                       cut.size    = 5000000,
                       time        = 300){


    harvest_para=paste("-minlenltr", minlenltr, "-maxlenltr", maxlenltr, "-mintsd", mintsd, "-maxtsd", maxtsd, "-motif", motif, "-motifmis", motifmis, "-similar", similar, "-vic", vic, "-seed", seed, "-seqids", "no")

    message("Run LTRharvest...")

    system(paste0(pathtoltrharvest,"LTR_HARVEST_parallel_modifiableparams -seq  ",genome.file," -harvest_para ",harvest_para," -threads ",threads," -size ",cut.size," -time ", time," -try1 ","1"," -gt ","/tmp/global2/ssushko/conda/envs/EDTA/bin/gt"))
    # system(paste0(pathtoltrharvest,"LTR_HARVEST_parallel_modifiableparams -seq  ",genome.file," -harvest_para ",harvest_para," -threads ",threads," -size ",cut.size," -time ", time," -try1 ","1"," -gt ","/usr/bin/gt"," > out.out"))
                       }
