#' regionGrangeMake
#'
#' @description Makes Grange objects for each region with a detected CNV.
#'
#' @param MA MultiAssayExperiment object used to store all information.
#' @param covprepped  List of dataframes, each dataframe is small cov file which ranges the detected CNV. This should be the padded cov file created with SARC::plotCovPrep and will be in the matadata.
#' @param range Number of rows (genomic ranges) to group together as one entry for grange. As some WES platforms will have very short genomic start-end coordinates, this is useful to extract the exon/ gene data. Default is 10.
#' @param gap The gap between the start of one grange entry to the next entry. Default is 10.
#'
#' @return List of Grange Objects, one for each CNV detected.
#' @export
#'
#' @examples
#' data("test_bed")
#' data("test_cov")
#' SARC <- regionSet(bed = test_bed, cov = test_cov)
#' SARC <- regionSplit(MA = SARC, bed = test_bed, cov = test_cov,
#'                     startlist = metadata(SARC)[[1]],
#'                     endlist = metadata(SARC)[[2]])
#' SARC <- plotCovPrep(MA = SARC, bed = test_bed, cov = test_cov,
#'                     startlist = metadata(SARC)[[1]],
#'                     endlist = metadata(SARC)[[2]])
#' SARC <- regionGrangeMake(MA = SARC, covprepped = metadata(SARC)[[4]])
regionGrangeMake <- function(MA, covprepped, range=10, gap=10){

  if (missing(MA)) stop('MA is missing. Add a MultiAssayExperiment object to store data efficiently.')

  if (missing(covprepped)) stop('MA is covprepped Add a covprepped list from SARC::plotCovPrep.')

  #apply regionGrange function for each prepared cov file

  allRanges <- lapply(covprepped, function(x){

    y <- regionGrange(x, range, gap)

    return(y)

  })

  #store output in MA

  metadata(MA)[["COVGRANGES"]] <- allRanges

  #return new MA object

  return(MA)
}
