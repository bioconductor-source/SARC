#' setupCNVplot
#'
#' @description Sets up dataframes for plotting the CNVs. If there are grange objects with no genes/ exons, it may be the CNVs are too small. These CNVs should be removed to not lead to errors.
#'
#' @param MA MultiAssayExperiment object used to store all information.
#' @param namedgranges List of grange objects, one for each CNV. These objects should have exons and genes, and will be in metadata after being created by SARC::addExonsGenes.
#' @param covprepped List of dataframes, one for each CNV. This list should be in metadata after bring created by SARC::plotCovPrep.
#'
#' @return List of dataframes, each dataframe can be plotted using plotting functions from this package.
#' @export
#'
#' @examples
#' data("test_bed")
#' data("test_cov")
#' test_bed <- test_bed[c(1),]
#' SARC <- regionSet(bed = test_bed, cov = test_cov)
#' SARC <- plotCovPrep(MA = SARC, bed = test_bed, cov = test_cov,
#'                    startlist = metadata(SARC)[[1]],
#'                    endlist = metadata(SARC)[[2]])
#' SARC <- regionGrangeMake(MA = SARC, covprepped = metadata(SARC)[[3]])
#'
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' library(Homo.sapiens)
#' TxDb(Homo.sapiens) <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' tx <- transcriptsBy(Homo.sapiens, columns = "SYMBOL")
#' txgene <- tx@unlistData
#'
#' SARC <- addExonsGenes(MA = SARC, covgranges = metadata(SARC)[[4]],
#'                       txdb = txdb, txgene = txgene)
#'
#' SARC <- setupCNVplot(MA = SARC, namedgranges =  metadata(SARC)[[5]],
#'                   covprepped = metadata(SARC)[[3]])
setupCNVplot <- function(MA, namedgranges, covprepped){

  if (missing(MA)) stop('MA is missing. Add a MultiAssayExperiment object to store data efficiently.')

  if (missing(namedgranges)) stop('namedgranges is missing. Add a list of named grange objects for plotting. Should be created from SARC::addExonsGenes.')

  if (missing(covprepped)) stop('covprepped is missing. Add a list of small cov dataframes that have been prepared for plotting. Should be created from SARC::plotCovPrep.')

  #apply the setupPlotStep1 function

  sps1 <- setupPlotStep1(namedgranges = namedgranges)

  #remove any null (no genes) CNVs

  ns <- names(namedgranges)

  if (isTRUE(list(NULL) %in% sps1)) {

    if (isTRUE(list(NULL) %in% sps1)) {

      s <- lapply(sps1, is.null)

      n <- names(s[s==TRUE])

      #communicate which samples should be removed if they were not associated with any genes

      paste0(n, " should be removed as it may be too small to evaluate. Removing from further analysis. Please also remove from .bed file.")

      sps1[n] <- NULL

      ns <- names(sps1)

    } else if (isFALSE(list(NULL) %in% sps1)) {

      print("All CNVs are fine for further evaluation")

    }

  }

  covprepped2 <- covprepped[ns]

  #apply setupPlotStep2 on the output of setupPlotStep1

  sps2 <- suppressWarnings(setupPlotStep2(covprepped = covprepped2, sps1 = sps1))

  names(sps2) <- names(sps1)

  #store as list

  metadata(MA)[["CNVPLOTLIST"]] <- sps2

  #return as new MA object

  return(MA)
}
