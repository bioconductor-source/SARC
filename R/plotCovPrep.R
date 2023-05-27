#' plotCovPrep
#'
#' @param MA MultiAssayExperiment object used to store all information.
#' @param bed Bed file containing CNVs for analysis. Must use one which contains at least MeanScore, Qlow, Qhigh and anova p-values. If Post-hoc tests have also been performed, these results can also be used by this function.
#' @param cov cov file from WES platform/ sequencer kit or if WGS regular intervals. Stored as a dataframe - genomic locations as rows and samples as columns.
#' @param n1 How many rows to pad the start site of each CNV - dafault is 0.
#' @param n2 How many rows to pad the end site of each CNV - dafault is 0.
#' @param startlist List of start sites created from SARC::regionSet and stored as metadata.
#' @param endlist List of end sites created from SARC::regionSet and stored as metadata.
#'
#' @return A new list of dataframes, each df is a region where a CNV was detected.
#' @export
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
#'                  start = metadata(SARC)[[1]], end = metadata(SARC)[[2]])
#' SARC <- anovaOnCNV(MA = SARC, bed = experiments(SARC)[[2]],
#'                    anovacov = metadata(SARC)[[7]])
#' SARC <- plotCovPrep(MA = SARC, bed = experiments(SARC)[[3]], cov = test_cov,
#'                     startlist = metadata(SARC)[[1]],
#'                     endlist = metadata(SARC)[[2]])
plotCovPrep <- function(MA, bed, cov, n1=0, n2=0, startlist, endlist){

  if (missing(MA)) stop('MA is missing. Add a MultiAssayExperiment object to store data efficiently.')

  if (missing(bed)) stop('bed is missing. Add bed dataframe. Ideally the most recently created bed file should be used.')

  if (missing(cov)) stop('cov is missing. Add a cov dataframe. This is a coverage file consisting of genomic locations and read depths for each sample in the WES/WGS cohort.')

  if (missing(startlist)) stop('start is missing. Add a startlist site list. Should be stored as metadata after SARC::regionSet.')

  if (missing(endlist)) stop('end is missing. Add a endlist site list. Should be stored as metadata after SARC::regionSet.')

  rownames(bed) <- NULL

  #add IDs to bed file

  bed$ID <- rownames(bed)

  #turn it into a list

  bedid <- as.list(as.numeric(bed$ID))

  #apply altPos function to each CNV in list

  newpos <- lapply(bedid, function(x){

    altPos(start = startlist[x],

           end = endlist[x], n1=n1, n2=n2)

  })

  #add padding

  altcov <- lapply(newpos, function(z){

    splitHelp(cov = cov, x = z[1], y = z[2])

  })

  #add CNV names

  names(altcov) <- paste0("CNV.", unlist(bedid))

  #add output to MA

  metadata(MA)[["COVPREPPED"]] <- altcov

  #return new MA object

  return(MA)
}
