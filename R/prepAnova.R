#' prepAnova
#'
#' @description Function for annovaOnCNV.This function sets up samples for anova analysis. Can also be used as input for dunnet analysis.
#'
#' @param MA MultiAssayExperiment object used to store all information.
#' @param bed Bed file containing CNVs which the user wishes to generate plots for. It is recommended that the most recently created bed file is used. Check print(MA) to see more bed files created by SARC.
#' @param cov cov file from WES platform/ sequencer kit or if WGS regular intervals. Stored as a dataframe - genomic locations as rows and samples as columns.
#' @param startlist List of start sites created from SARC::regionSet and stored as metadata.
#' @param endlist List of end sites created from SARC::regionSet and stored as metadata.
#'
#' @return Additional list of dataframes which will be needed to calculate anova for each region with a detected CNV.
#' @export
#' @importFrom reshape2 melt
#'
#' @examples
#' data("test_bed")
#' data("test_cov")
#' SARC <- regionSet(bed = test_bed, cov = test_cov)
#' SARC <- regionSplit(MA = SARC, bed = test_bed, cov = test_cov,
#'                     startlist = metadata(SARC)[[1]],
#'                     endlist = metadata(SARC)[[2]])
#' SARC <- regionMean(MA = SARC, bed = test_bed, splitcov = metadata(SARC)[[3]])
#' SARC <- regionQuantiles(MA = SARC, bed = experiments(SARC)[[1]],
#'                        meancov = metadata(SARC)[[4]], q1=.15, q2=.85)
#' SARC <- prepAnova(MA = SARC, bed = experiments(SARC)[[2]], cov = test_cov,
#'                  startlist = metadata(SARC)[[1]], endlist = metadata(SARC)[[2]])
prepAnova <- function(MA, bed, cov, startlist, endlist){

  if (missing(MA)) stop('MA is missing. Add a MultiAssayExperiment object to store data efficiently.')

  if (missing(bed)) stop('bed is missing. Add bed dataframe. Ideally the most recently created bed file should be used.')

  if (missing(cov)) stop('cov is missing. Add a cov dataframe. This is a coverage file consisting of genomic locations and read depths for each sample in the WES/WGS cohort.')

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

  #make dataframes to apply anova for each CNV

  anova.cov <- list()

  for (i in seq_len(nrow(bed))) {

    anova.cov[[i]] <- cov[starts[i]:ends[i],]

    names(anova.cov)[i] <- paste0("CNV.", i)

  }

  #store as a list in MA

  metadata(MA)[["ANOVACOV"]] <- anova.cov

  #retun a new MA object

  return(MA)

}
