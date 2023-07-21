#' regionSplit
#'
#' @description Splits the large cov file into multiple smaller cov files. Each one is specific for each CNV entry in the bam file.
#'
#' @param RE RaggedExperiment object used to store all information.
#' @param cnv List of CNVs in a dataframe containing CNVs from detection algorithms/ pipelines. It is recommended that the most recently created cnv file is used. Check print(RE) to see more cnv files created by SARC.
#' @param startlist List of start sites created from SARC::regionSet and stored as metadata.
#' @param endlist List of end sites created from SARC::regionSet and stored as metadata.
#'
#' @return A list of dataframes - each one a smaller cov file.
#'
#' @export
#'
#' @examples
#' data("test_cnv")
#' data("test_cov")
#' SARC <- regionSet(cnv = test_cnv, cov = test_cov)
#' SARC <- regionSplit(RE = SARC, cnv =  metadata(SARC)[['CNVlist']][[1]],
#'                     startlist = metadata(SARC)[[2]],
#'                     endlist = metadata(SARC)[[3]])
regionSplit <- function(RE, cnv, startlist, endlist){

  if (missing(RE)) stop('RE is missing. Add a RaggedExperiment object to store data efficiently.')

  if (missing(cnv)) stop('cnv is missing. Add cnv dataframe. Ideally the most recently created cnv file should be used.')

  if (missing(startlist)) stop('start is missing. Add a startlist site list. Should be stored as metadata after SARC::regionSet.')

  if (missing(endlist)) stop('end is missing. Add a endlist site list. Should be stored as metadata after SARC::regionSet.')

  #create split list object

  split.cov <- mapply(splitHelp, x = startlist, y = endlist,

                      MoreArgs=list(RE@assays@unlistData), SIMPLIFY = FALSE)

  #name each dataframe within list

  names(split.cov) <- paste0("CNV.", seq(1:nrow(cnv)))

  suppressWarnings(metadata(RE)[['SPLITCOV']] <- split.cov)

  #return a new RE object

  return(RE)

}
