#' regionSplit
#'
#' @description Splits the large cov file into multiple smaller cov files. Each one is specific for each CNV entry in the bam file.
#'
#' @param MA MultiAssayExperiment object used to store all information.
#' @param bed Bed file containing CNVs which the user wishes to generate plots for. It is recommended that the most recently created bed file is used. Check print(MA) to see more bed files created by SARC.
#' @param cov cov file from WES platform/ sequencer kit or if WGS regular intervals. Stored as a dataframe - genomic locations as rows and samples as columns.
#' @param startlist List of start sites created from SARC::regionSet and stored as metadata.
#' @param endlist List of end sites created from SARC::regionSet and stored as metadata.
#'
#' @return A list of dataframes - each one a smaller cov file.
#'
#' @export
#'
#' @examples
#' data("test_bed")
#' data("test_cov")
#' SARC <- regionSet(bed = test_bed, cov = test_cov)
#' SARC <- regionSplit(MA = SARC, bed = test_bed, cov = test_cov,
#'                     startlist = metadata(SARC)[[1]],
#'                     endlist = metadata(SARC)[[2]])
regionSplit <- function(MA, bed, cov, startlist, endlist){

  if (missing(MA)) stop('MA is missing. Add a MultiAssayExperiment object to store data efficiently.')

  if (missing(bed)) stop('bed is missing. Add bed dataframe. Ideally the most recently created bed file should be used.')

  if (missing(cov)) stop('cov is missing. Add a cov dataframe. This is a coverage file consisting of genomic locations and read depths for each sample in the WES/WGS cohort.')

  if (missing(startlist)) stop('start is missing. Add a startlist site list. Should be stored as metadata after SARC::regionSet.')

  if (missing(endlist)) stop('end is missing. Add a endlist site list. Should be stored as metadata after SARC::regionSet.')

  #create split list object

  split.cov <- mapply(splitHelp, y = startlist, x = endlist,

                      MoreArgs=list(cov), SIMPLIFY = FALSE)

  #name each dataframe within list

  names(split.cov) <- paste0("CNV.", seq(1:nrow(bed)))

  suppressWarnings(metadata(MA)[['SPLITCOV']] <- split.cov)

  #return a new MA object

  return(MA)

}
