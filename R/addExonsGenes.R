#' addExonsGenes
#'
#' @description For the length of the CNV (+padding, if this was done), exons and gene symbol information will be attached in each Grange object.
#'
#' @param RE RaggedExperiment object used to store all information.
#' @param covgranges List of grange objects created by the SARC::regionGrangeMake function. This should be found in the metadata of the RE object.
#' @param txdb Loaded txDB object e.g. TxDb.Hsapiens.UCSC.hg38.knownGene. Standard usage from the GenomicFeatures tutorial https://kasperdanielhansen.github.io/genbioconductor/html/GenomicFeatures.html.
#' @param txgene List of genes and exons associated with each genomic region of interest. Requires specific species package to make e.g. Homo.sapiens for human samples.
#'
#' @return A list of grange objects with exons and genes added. One grange object for each detected CNV. Will be stored in the RE object provided as metadata.
#' @export
#' @import GenomicFeatures
#' @importFrom plyranges join_overlap_inner
#' @importFrom IRanges subsetByOverlaps
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
#' SARC <- plotCovPrep(RE = SARC, cnv = metadata(SARC)[['CNVlist']][[1]],
#'                    startlist = metadata(SARC)[[2]],
#'                    endlist = metadata(SARC)[[3]])
#' SARC <- regionGrangeMake(RE = SARC, covprepped = metadata(SARC)[[4]])
#'
#' TxDb(Homo.sapiens) <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' tx <- transcriptsBy(Homo.sapiens, columns = "SYMBOL")
#' txgene <- tx@unlistData
#' SARC <- addExonsGenes(RE = SARC, covgranges = metadata(SARC)[[5]],
#'                       txdb = txdb, txgene = txgene)
addExonsGenes <- function(RE, covgranges, txdb, txgene){

  if (missing(RE)) stop('RE is missing. Add a RaggedExperiment object to store data efficiently.')

  if (missing(covgranges)) stop('covgranges is missing. Add list of grange objects. This should be stored in the RE object after using rag::regionGrange.')

  if (missing(txdb)) stop('txdb is missing. Follow instructions in the vignette to create the txdb object or read https://kasperdanielhansen.github.io/genbioconductor/html/GenomicFeatures.html.')

  if (missing(txgene)) stop('txgene is missing. Follow instructions in the vignette to create the txgene object or read https://kasperdanielhansen.github.io/genbioconductor/html/GenomicFeatures.html.')

  #add exons

  txexons <- lapply(covgranges, function(x){

    lapply(x, function(y){

      ex <- subsetByOverlaps(exonsBy(txdb, by="tx"), y, ignore.strand = TRUE)

      ex <- ex@unlistData

    })

  })
  #add genes

  newRange <- lapply(txexons, function(x){

    lapply(x, function(y){

      join_overlap_inner(y, txgene)

    })

  })

  #store in RE

  metadata(RE)[["NAMEDGRANGES"]] <- newRange

  #return RE

  return(RE)

}
