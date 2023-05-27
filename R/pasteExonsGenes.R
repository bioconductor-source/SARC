#' pasteExonsGenes
#'
#' @description Add Genes and Exons to the bed file in order to show the range the detected CNV extends. This is useful information for clinicians/ diagnostic scientists when interpreting the variants. Note - if padding was performed during plotCovPrep step, the CNV range would be extended by the padding so extra exons may be included in these cases.
#'
#' @param MA MultiAssayExperiment object used to store all information.
#' @param setup List of dataframes which have been processed for plotting. This will be stored as metadata in the MA object after SARC::setupCNVPlot.
#' @param bed Bed file containing CNVs which the user wishes to generate plots for. It is recommended that the most recently created bed file is used. Check print(MA) to see more bed files created by SARC.
#'
#' @return A new bed file with an additional column which lists the genes and exons which the detected variant ranges.
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
#' SARC <- setupCNVplot(MA = SARC, namedgranges =  metadata(SARC)[[5]],
#'                   covprepped = metadata(SARC)[[3]])
#' SARC <- pasteExonsGenes(MA = SARC, setup =  metadata(SARC)[[6]],
#'                        bed = test_bed)
#'
pasteExonsGenes <- function(MA, setup, bed){

  if (missing(MA)) stop('MA is missing. Add a MultiAssayExperiment object to store data efficiently.')

  if (missing(setup)) stop('setup is missing. Should be stored as metadata after SARC::setupCNVPlot.')

  if (missing(bed)) stop('bed is missing. Add bed dataframe. Ideally the most recently created bed file should be used.')

  #separate out the gene.exon strings
  genomics <- lapply(setup, function(x) {

    y <- paste0(unique(x$FEATURE), collapse = "|")

  })

  #organise the dataframe

  genomics <- as.data.frame(genomics)

  genomics <- t(genomics)

  #add this to the bed file

  bed$Gene.Exon <- genomics[,1]

  #add the new bed file into the MA object

  MA2 <- suppressWarnings(suppressMessages(MultiAssayExperiment(list(BEDGENEEXON = bed))))

  MA <- suppressWarnings(suppressMessages(c(MA, MA2)))

  #return MA object

  return(MA)
}
