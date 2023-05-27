#' addExonsGenes
#'
#' @description For the length of the CNV (+padding, if this was done), exons and gene symbol information will be attached in each Grange object.
#'
#' @param MA MultiAssayExperiment object used to store all information.
#' @param covgranges List of grange objects created by the SARC::regionGrangeMake function. This should be found in the metadata of the MA object.
#' @param txdb Loaded txDB object e.g. TxDb.Hsapiens.UCSC.hg38.knownGene. Standard usage from the GenomicFeatures tutorial https://kasperdanielhansen.github.io/genbioconductor/html/GenomicFeatures.html.
#' @param txgene List of genes and exons associated with each genomic region of interest. Requires specific species package to make e.g. Homo.sapiens for human samples.
#'
#' @return A list of grange objects with exons and genes added. One grange object for each detected CNV. Will be stored in the MA object provided as metadata.
#' @export
#' @import GenomicFeatures
#' @importFrom plyranges join_overlap_inner
#' @importFrom IRanges subsetByOverlaps
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
addExonsGenes <- function(MA, covgranges, txdb, txgene){

  if (missing(MA)) stop('MA is missing. Add a MultiAssayExperiment object to store data efficiently.')

  if (missing(covgranges)) stop('covgranges is missing. Add list of grange objects. This should be stored in the MA object after using SARC::regionGrange.')

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

  #store in MA

  metadata(MA)[["NAMEDGRANGES"]] <- newRange

  #return MA.

  return(MA)

}
