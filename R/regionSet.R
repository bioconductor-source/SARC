#' @title regionSet
#'
#' @description Stores .cov and .bed file into a multiassayexperiment object. Will use genomic locations from the .bed file to find the start and end points of the CNVs from the .cov file and store these as two lists.
#'
#' @param bed CNV bed files stored as dataframes. Chromosomes, Start of CNV, End of CNV, Type of CNV and
#' some other information can be kept here.
#' @param cov cov file from WES platform/ sequencer kit or if WGS regular intervals. Stored
#' as a dataframe - genomic locations as rows and samples as columns.
#'
#' @return MultiAssayExperiment Object with two lists: Start sites and End sites.
#'
#' @import MultiAssayExperiment
#' @export
#'
#' @examples
#' data("test_bed")
#' data("test_cov")
#' SARC <- regionSet(bed = test_bed, cov = test_cov)
regionSet <- function(bed, cov){

  if (missing(bed)) stop('bed is missing. Add bed dataframe. Ideally the most recently created bed file should be used.')

  if (missing(cov)) stop('cov is missing. Add a cov dataframe. This is a coverage file consisting of genomic locations and read depths for each sample in the WES/WGS cohort.')

  #make sure cov is not factorized

  cov$CHROM <- as.character(cov$CHROM)

  #Create MA object and other vectors to be stored in MA

  MA <- MultiAssayExperiment()

  cov.mini <- list()

  start.mini <- vector()

  end.mini <- vector()

  for (i in seq_len(nrow(bed))) {

    #Break cov into chromosomes - for speed

    cov.mini[[i]] <- cov[which(cov$CHROM==bed$CHROM[i]),]

    #Make iteration objects

    c <- cov.mini[[i]]

    s <- bed$START[i]

    e <- bed$END[i]

    #Get ith START from cov

    s.min <- which.min(abs(c$START-s))

    start.mini[i] <- c$ID[s.min]

    #Get ith END from cov

    e.min <- which.min(abs(c$END-e))

    end.mini[i] <- c$ID[e.min]

  }

  #store in MA

  metadata(MA)[["REGIONSTART"]] <- as.numeric(start.mini)

  metadata(MA)[["REGIONEND"]] <- as.numeric(end.mini)

  #clean up

  rm(s.min, s, e.min, e, i, start.mini, end.mini)

  #Return a new MA object

  return(MA)

}
