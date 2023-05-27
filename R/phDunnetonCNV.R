#' phDunnetonCNV
#'
#' @description Applies post-hoc test Dunnet's test on each CNV in the bed file. Performs a pair-wise test between the sample where the CNV was detected and each other sample at the same genomic region. WARNING - very slow, only use if a small number of CNVs / samples are available. Other tests e.g. Quantile Distribution are more effective with larger cohorts.
#'
#' @param MA MultiAssayExperiment object used to store all information.
#' @param bed  CNV bed files stored as dataframes. Chromosomes, Start of CNV, End of CNV, Type of CNV and some other information can be kept here. It is reccomended to use the most recently created bed file. Use print(MA) to check.
#' @param anovacov List of dataframes which contain the required information to perform anova for each region with a detected CNV. Generated from prepAnova.
#'
#' @return A new .bed file with an additional column of the sumlog p-value from each dunnet test performed for the CNV entry.
#' @export
#'
phDunnetonCNV <- function(MA, bed, anovacov){

  if (missing(MA)) stop('MA is missing. Add a MultiAssayExperiment object to store data efficiently.')

  if (missing(bed)) stop('bed is missing. Add bed dataframe. Ideally the most recently created bed file should be used.')

  if (missing(anovacov)) stop('anovacov is missing. Add a list of prepped anovas dataframes. Should be stored as metadata after prepAnova')

  #store list of controls for dunnet analysis

  controls <- as.list(bed$SAMPLE)

  #apply dunnet testing

  res <- mapply(applyDunnet, preppedcov = anovacov, control=controls)

  #store p-values in the bed file

  bed$Dunnet <- unlist(res)

  #store new bed file in the MA object

  MA2 <- suppressWarnings(suppressMessages(MultiAssayExperiment(list(BEDDUNNET = bed))))

  MA <- suppressWarnings(suppressMessages(c(MA, MA2)))

  #return a new MA object

  return(MA)
}
