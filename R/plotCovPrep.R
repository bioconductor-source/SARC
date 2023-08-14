#' plotCovPrep
#'
#' @description Prepares plotting of coverage at specific regions of the DNA - where CNVs have been detected.
#'
#' @param RE MultiAssayExperiment object used to store all information.
#' @param cnv List of CNVs in a dataframe containing CNVs from detection algorithms/ pipelines. Must use one which contains at least MeanScore, Qlow, Qhigh and anova p-values. If Post-hoc tests have also been performed, these results can also be used by this function.
#' @param n1 How many rows to pad the start site of each CNV - dafault is 0.
#' @param n2 How many rows to pad the end site of each CNV - dafault is 0.
#' @param startlist List of start sites created from SARC::regionSet and stored as metadata.
#' @param endlist List of end sites created from SARC::regionSet and stored as metadata.
#'
#' @return A new list of dataframes, each df is a region where a CNV was detected.
#' @export
#'
#' @examples
#' data("test_cnv")
#' test_cnv <- test_cnv[c(1:3),]
#' data("test_cov")
#' SARC <- regionSet(cnv = test_cnv, cov = test_cov)
#' SARC <- regionSplit(RE = SARC, cnv = metadata(SARC)[['CNVlist']][[1]],
#'                     startlist = metadata(SARC)[[2]],
#'                     endlist = metadata(SARC)[[3]])
#' SARC <- plotCovPrep(RE = SARC, cnv = metadata(SARC)[['CNVlist']][[1]],
#'                     startlist = metadata(SARC)[[2]],
#'                     endlist = metadata(SARC)[[3]])
plotCovPrep <- function(RE, cnv, n1=0, n2=0, startlist, endlist){

  if (missing(RE)) stop('RE is missing. Add a RaggedExperiment object to store data efficiently.')

  if (missing(cnv)) stop('cnv is missing. Add cnv dataframe. Ideally the most recently created cnv file should be used.')

  if (missing(startlist)) stop('start is missing. Add a startlist site list. Should be stored as metadata after SARC::regionSet.')

  if (missing(endlist)) stop('end is missing. Add a endlist site list. Should be stored as metadata after SARC::regionSet.')

  rownames(cnv) <- NULL

  #add IDs to cnv file

  cnv$ID <- rownames(cnv)

  #turn it into a list

  cnvid <- as.list(as.numeric(cnv$ID))

  #apply altPos function to each CNV in list

  newpos <- lapply(cnvid, function(x){

    altPos(start = startlist[x],

           end = endlist[x], n1=n1, n2=n2)

  })

  #add padding

  altcov <- lapply(newpos, function(z){

    splitHelp(cov = RE@assays@unlistData, x = z[1], y = z[2])

  })

  #add CNV names

  names(altcov) <- paste0("CNV.", unlist(cnvid))

  #add output to RE

  metadata(RE)[["COVPREPPED"]] <- altcov

  #return new RE object

  return(RE)

}
