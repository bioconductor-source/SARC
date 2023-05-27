#' anovaOnCNV
#'
#' @description Apply anova on the CNV data. For each CNV listed in the input bed file, the read-depths will be used to calculate if this region of the DNA, in this sample, is significantly different from other samples in the cohort of samples, at this same region.
#'
#'
#' @param MA  MultiAssayExperiment object used to store all information.
#' @param bed CNV bed files stored as dataframes. Chromosomes, Start of CNV, End of CNV, Type of CNV and some other information can be kept here. It is recomended to use the most recently generated bed file. Check print(MA) to see if bed files from other anlyses have been generated.
#' @param anovacov List of dataframes which contain the required information to perform anova for each region with a detected CNV. Generated from prepAnova, and stored as metadata.
#'
#' @return A new bed file with an additional column for the p-value from anova. Will be stored as an experiment in the MA object.
#' @export
#'
#' @examples
#' data("test_bed")
#' data("test_cov")
#' SARC <- regionSet(bed = test_bed, cov = test_cov)
#' SARC <- regionSplit(MA = SARC, bed = test_bed, cov = test_cov,
#'                     startlist = metadata(SARC)[[1]],
#'                     endlist = metadata(SARC)[[2]])
#' SARC <- prepAnova(MA = SARC, bed = test_bed, cov = test_cov,
#'                  start = metadata(SARC)[[1]], end = metadata(SARC)[[2]])
#' SARC <- anovaOnCNV(MA = SARC, bed = test_bed,
#'                    anovacov = metadata(SARC)[[3]])
anovaOnCNV <- function(MA, bed, anovacov){

  if (missing(MA)) stop('MA is missing. Add a MultiAssayExperiment object to store data efficiently.')

  if (missing(bed)) stop('bed is missing. Add bed dataframe file. Ideally the most recently created bed file should be used.')

  if (missing(anovacov)) stop('anovacov is missing. Add a list of prepped anovas dataframes. Should be stored as metadata after prepAnova')

  #apply anova on the region of the CNV

  anovaRes <- vector()

  anovaRes <- lapply(anovacov, function(x){

    x <- applyAnova(x)

    return(x)

  })

  #store resulting p values in the bed file

  bed$ANOVA <- unlist(anovaRes)

  #store data in the MA object
  MA2 <- suppressWarnings(suppressMessages(MultiAssayExperiment(list(BEDANOVA = bed))))

  MA <- suppressWarnings(suppressMessages(c(MA, MA2)))

  #return MA object

  return(MA)

}
