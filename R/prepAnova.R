#' prepAnova
#'
#' @description Function for annovaOnCNV.This function sets up samples for anova analysis. Can also be used as input for Dunnet analysis.
#'
#' @param RE RaggedExperiment object used to store all information.
#' @param cnv List of CNVs in a dataframe containing CNVs from detection algorithms/ pipelines. It is recommended that the most recently created cnv file is used. Check print(RE) to see more cnv files created by SARC.
#' @param startlist List of start sites created from SARC::regionSet and stored as metadata.
#' @param endlist List of end sites created from SARC::regionSet and stored as metadata.
#'
#' @return Additional list of dataframes which will be needed to calculate anova for each region with a detected CNV.
#' @export
#' @importFrom reshape2 melt
#'
#' @examples
#' data("test_cnv")
#' data("test_cov")
#' SARC <- regionSet(cnv = test_cnv, cov = test_cov)
#' SARC <- regionSplit(RE = SARC, cnv = metadata(SARC)[['CNVlist']][[1]],
#'                     startlist = metadata(SARC)[[2]],
#'                     endlist = metadata(SARC)[[3]])
#' SARC <- regionMean(RE = SARC, cnv = metadata(SARC)[['CNVlist']][[1]],
#'                    splitcov = metadata(SARC)[[4]])
#' SARC <- regionQuantiles(RE = SARC, cnv = metadata(SARC)[['CNVlist']][[2]],
#'                       meancov = metadata(SARC)[[5]], q1=.1, q2=.9)
#' SARC <- prepAnova(RE = SARC, cnv = metadata(SARC)[['CNVlist']][[3]],
#'                  startlist = metadata(SARC)[[2]],
#'                  endlist = metadata(SARC)[[3]])
prepAnova <- function(RE, cnv, startlist, endlist){

  if (missing(RE)) stop('RE is missing. Add a RaggedExperiment object to store data efficiently.')

  if (missing(cnv)) stop('cnv is missing. Add cnv dataframe. Ideally the most recently created cnv file should be used.')

  if (missing(startlist)) stop('start is missing. Add a startlist site list. Should be stored as metadata after SARC::regionSet.')

  if (missing(endlist)) stop('end is missing. Add a endlist site list. Should be stored as metadata after SARC::regionSet.')

  #make sure start and end lists are numerical

  start.mini <- as.numeric(startlist)

  end.mini <- as.numeric(endlist)

  #apply prepBase on each start and end site in the start and end lists

  res <- mapply(prepBase, start=start.mini, end=end.mini)

  #create start and end objects

  starts <- res[1,]
  ends <- res[2,]

  #make Grange to apply anova for each CNV

  anova.cov <- list()

  for (i in seq_len(nrow(cnv))) {

    s <- starts[i]
    e <- ends[i]

    anova.cov[[i]] <- RE@assays@unlistData[s:e]

    names(anova.cov)[i] <- paste0("CNV.", i)

  }

  #store as a list in RE

  metadata(RE)[["ANOVACOV"]] <- anova.cov

  #retun a new RE object

  return(RE)

}
