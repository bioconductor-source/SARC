#' setupCNVplot
#'
#' @description Sets up dataframes for plotting the CNVs. If there are grange objects with no genes/ exons, it may be the CNVs are too small. These CNVs should be removed to not lead to errors.
#'
#' @param RE RaggedExperiment object used to store all information.
#' @param namedgranges List of grange objects, one for each CNV. These objects should have exons and genes, and will be in metadata after being created by SARC::addExonsGenes.
#' @param covprepped List of dataframes, one for each CNV. This list should be in metadata after bring created by SARC::plotCovPrep.
#'
#' @return List of dataframes, each dataframe can be plotted using plotting functions from this package.
#' @export
#'
#' @examples
#' if (requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
#' require("TxDb.Hsapiens.UCSC.hg38.knownGene")
#' } else {}
#'
#' if (requireNamespace("Homo.sapiens", quietly = TRUE)) {
#' require("Homo.sapiens")
#' } else {}
#'
#' data("test_cnv")
#' test_cnv <- test_cnv[c(1),]
#' data("test_cov")
#' SARC <- regionSet(cnv = test_cnv, cov = test_cov)
#' SARC <- regionSplit(RE = SARC, cnv =  metadata(SARC)[['CNVlist']][[1]],
#'                     startlist = metadata(SARC)[[2]],
#'                     endlist = metadata(SARC)[[3]])
#' SARC <- plotCovPrep(RE = SARC, cnv = metadata(SARC)[['CNVlist']][[1]],
#'                     startlist = metadata(SARC)[[2]],
#'                     endlist = metadata(SARC)[[3]])
#' SARC <- regionGrangeMake(RE = SARC, covprepped = metadata(SARC)[[4]])
#'
#' TxDb(Homo.sapiens) <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' tx <- transcriptsBy(Homo.sapiens, columns = "SYMBOL")
#' txgene <- tx@unlistData
#'
#' SARC <- addExonsGenes(RE = SARC, covgranges = metadata(SARC)[[6]],
#'                       txdb = txdb, txgene = txgene)
#' SARC <- setupCNVplot(RE = SARC, namedgranges =  metadata(SARC)[[7]],
#'                   covprepped = metadata(SARC)[[4]])
setupCNVplot <- function(RE, namedgranges, covprepped){

  if (missing(RE)) stop('RE is missing. Add a RaggedExperiment object to store data efficiently.')

  if (missing(namedgranges)) stop('namedgranges is missing. Add a list of named grange objects for plotting. Should be created from SARC::addExonsGenes.')

  if (missing(covprepped)) stop('covprepped is missing. Add a list of small cov dataframes that have been prepared for plotting. Should be created from SARC::plotCovPrep.')

  #apply the setupPlotStep1 function

  sps1 <- setupPlotStep1(namedgranges = namedgranges)

  #remove any null (no genes) CNVs

  ns <- names(namedgranges)

  if (isTRUE(list(NULL) %in% sps1)) {

    s <- lapply(sps1, is.null)

    n <- names(s[s==TRUE])

    #communicate which samples should be removed if they were not associated with any genes

    print(paste0(n, " should be removed as it may be too small to evaluate. Removing from further analysis. Please also remove from .cnv file."))

    sps1[n] <- NULL

    ns <- names(sps1)

  } else if (isFALSE(list(NULL) %in% sps1)) {

    print("All CNVs are fine for further evaluation")

  }

  covprepped2 <- covprepped[ns]

  #apply setupPlotStep2 on the output of setupPlotStep1

  sps2 <- suppressWarnings(setupPlotStep2(covprepped = covprepped2, sps1 = sps1))

  names(sps2) <- names(sps1)

  #store as list

  metadata(RE)[["CNVPLOTLIST"]] <- sps2

  #return as new RE object

  return(RE)
}
