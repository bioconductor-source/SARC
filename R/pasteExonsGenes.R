#' pasteExonsGenes
#'
#' @description Add Genes and Exons to the cnv file in order to show the range the detected CNV extends. This is useful information for clinicians/ diagnostic scientists when interpreting the variants. Note - if padding was performed during plotCovPrep step, the CNV range would be extended by the padding so extra exons may be included in these cases.
#'
#' @param RE RaggedExperiment object used to store all information.
#' @param setup List of dataframes which have been processed for plotting. This will be stored as metadata in the RE object after SARC::setupCNVPlot.
#' @param cnv cnv file containing CNVs which the user wishes to generate plots for. It is recommended that the most recently created cnv file is used. Check print(RE) to see more cnv files created by SARC.
#' @param nameofnewdf Name of new dataframe to be saved in metadata(RE)[['CNVlist']]. Default is cnvExonsGenes.
#'
#' @return A new cnv file with an additional column which lists the genes and exons which the detected variant ranges.
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
#' SARC <- pasteExonsGenes(RE = SARC, setup =  metadata(SARC)[[8]],
#'                        cnv = metadata(SARC)[['CNVlist']][[1]])
#'
pasteExonsGenes <- function(RE, setup, cnv, nameofnewdf="cnvExonsGenes"){

  if (missing(RE)) stop('RE is missing. Add a RaggedExperiment object to store data efficiently.')

  if (missing(setup)) stop('setup is missing. Should be stored as metadata after SARC::setupCNVPlot.')

  if (missing(cnv)) stop('cnv is missing. Add cnv dataframe. Ideally the most recently created cnv file should be used.')

  #separate out the gene.exon strings
  genomics <- lapply(setup, function(x) {

    y <- paste0(unique(x$FEATURE), collapse = "|")

  })

  #organise the dataframe

  genomics <- as.data.frame(genomics)

  genomics <- t(genomics)

  #add this to the cnv file

  cnv$Gene.Exon <- genomics[,1]

  #add the new cnv file into the RE object

  metadata(RE)[["CNVlist"]][[nameofnewdf]] <- cnv

  #return RE object

  return(RE)
}
