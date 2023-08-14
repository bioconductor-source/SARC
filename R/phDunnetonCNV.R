#' phDunnetonCNV
#'
#' @description Applies post-hoc test Dunnet test on each CNV in the cnv file. Performs a pair-wise test between the sample where the CNV was detected and each other sample at the same genomic region. WARNING - very slow, only use if a small number of CNVs / samples are available. Other tests e.g. Quantile Distribution are more effective with larger cohorts.
#'
#' @param RE RaggedExperiment object used to store all information.
#' @param cnv  List of CNVs in a dataframe containing CNVs from detection algorithms/ pipelines. It is recommended to use the most recently created cnv file. Use print(RE) to check.
#' @param anovacov List of dataframes which contain the required information to perform anova for each region with a detected CNV. Generated from prepAnova.
#' @param nameofnewdf Name of new dataframe to be saved in metadata(RE)[['CNVlist']]. Default is CNVdunnet.
#'
#' @return A new cnv dataframe with an additional column of the sumlog p-value from each dunnet test performed for the CNV entry.
#' @export
#'
phDunnetonCNV <- function(RE, cnv, anovacov, nameofnewdf="CNVdunnet"){

  if (missing(RE)) stop('RE is missing. Add a RaggedExperiment object to store data efficiently.')

  if (missing(cnv)) stop('cnv is missing. Add cnv dataframe. Ideally the most recently created cnv file should be used.')

  if (missing(anovacov)) stop('anovacov is missing. Add a list of prepped anovas dataframes. Should be stored as metadata after prepAnova')

  #store list of controls for dunnet analysis

  controls <- as.list(cnv$SAMPLE)

  #apply dunnet testing

  res <- mapply(applyDunnet, preppedcov = anovacov, control=controls)

  #store p-values in the cnv file

  cnv$Dunnet <- unlist(res)

  #store new cnv file in the RE object

  metadata(RE)[["CNVlist"]][[nameofnewdf]] <- cnv

  #return a new RE object

  return(RE)
}
