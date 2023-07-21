#' regionGrangeMake
#'
#' @description Makes Grange objects for each region with a detected CNV.
#'
#' @param RE RaggedExperiment object used to store all information.
#' @param covprepped  List of dataframes, each dataframe is small cov file which ranges the detected CNV. This should be the padded cov file created with SARC::plotCovPrep and will be in the matadata.
#' @param range Number of rows (genomic ranges) to group together as one entry for grange. As some WES platforms will have very short genomic start-end coordinates, this is useful to extract the exon/ gene data. Default is 10.
#' @param gap The gap between the start of one grange entry to the next entry. Default is 10.
#'
#' @return List of Grange Objects, one for each CNV detected.
#' @export
#'
#' @examples
#' data("test_cnv")
#' test_cnv <- test_cnv[c(1:3),]
#' data("test_cov")
#' SARC <- regionSet(cnv = test_cnv, cov = test_cov)
#' SARC <- regionSplit(RE = SARC, cnv =  metadata(SARC)[['CNVlist']][[1]],
#'                     startlist = metadata(SARC)[[2]],
#'                     endlist = metadata(SARC)[[3]])
#' SARC <- plotCovPrep(RE = SARC, cnv = metadata(SARC)[['CNVlist']][[1]],
#'                     startlist = metadata(SARC)[[2]],
#'                     endlist = metadata(SARC)[[3]])
#' SARC <- regionGrangeMake(RE = SARC, covprepped = metadata(SARC)[[4]])
regionGrangeMake <- function(RE, covprepped, range=10, gap=10){

  if (missing(RE)) stop('RE is missing. Add a RaggedExperiment object to store data efficiently.')

  if (missing(covprepped)) stop('RE is missing. Add a covprepped list from SARC::plotCovPrep.')

  #apply regionGrange function for each prepared cov file

  allRanges <- lapply(covprepped, function(x){

    y <- regionGrange(x, range, gap)

    return(y)

  })

  #store output in RE

  metadata(RE)[["COVGRANGES"]] <- allRanges

  #return new RE object

  return(RE)
}
