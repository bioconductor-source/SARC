#' anovaOnCNV
#'
#' @description Apply anova on the CNV data. For each CNV listed in the input cnv file, the read-depths will be used to calculate if this region of the DNA, in this sample, is significantly different from other samples in the cohort of samples, at this same region.
#'
#' @param RE  RaggedExperiment object used to store all information.
#' @param cnv List of CNVs in a dataframe containing CNVs from detection algorithms/ pipelines.
#' @param anovacov List of dataframes which contain the required information to perform anova for each region with a detected CNV. Generated from prepAnova, and stored as metadata.
#' @param nameofnewdf Name of new dataframe to be saved in metadata(RE)[['CNVlist']]. Default is CNVanova.
#'
#' @return A new cnv dataframe with an additional column for the p-value from anova. Will be stored as an experiment in the RE object.
#' @export
#'
#' @examples
#' data("test_cnv")
#' data("test_cov")
#' SARC <- regionSet(cnv = test_cnv, cov = test_cov)
#' SARC <- regionSplit(RE = SARC, cnv = metadata(SARC)[['CNVlist']][[1]],
#'                     startlist = metadata(SARC)[[2]],
#'                     endlist = metadata(SARC)[[3]])
#' SARC <- prepAnova(RE = SARC, cnv = metadata(SARC)[['CNVlist']][[1]],
#'                  start = metadata(SARC)[[2]], end = metadata(SARC)[[3]])
#' SARC <- anovaOnCNV(RE = SARC, cnv = metadata(SARC)[['CNVlist']][[1]],
#'                    anovacov = metadata(SARC)[[4]])
anovaOnCNV <- function(RE, cnv, anovacov, nameofnewdf="CNVanova"){

  if (missing(RE)) stop('RE is missing. Add a RaggedExperiment object to store data efficiently.')

  if (missing(cnv)) stop('cnv is missing. Add cnv dataframe file. Ideally the most recently created cnv file should be used.')

  if (missing(anovacov)) stop('anovacov is missing. Add a list of prepped anovas dataframes. Should be stored as metadata after prepAnova')

  #apply anova on the region of the CNV

  anovaRes <- vector()

  anovaRes <- lapply(anovacov, function(x){

    x <- as.data.frame(x)

    x <- applyAnova(x)

    return(x)

  })

  #store resulting p values in the cnv file
  cnv$ANOVA <- unlist(anovaRes)

  #store new dataframe in the RE object
  metadata(RE)[["CNVlist"]][[nameofnewdf]] <- cnv

  #return MA object

  return(RE)

}
