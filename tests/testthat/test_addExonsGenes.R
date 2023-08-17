#this script will test the addition of exons and genes to Grange objects
#this is a vital part of
library(testthat)
data("test_cnv")
data("test_cov")

#for speed we will just keep one CNV
test_cnv <- test_cnv[c(1),]
SARC <- regionSet(cnv = test_cnv, cov = test_cov)
#test cov is attached as assay
expect_equal(start(SARC@assays@unlistData),test_cov$START)
expect_equal(end(SARC@assays@unlistData),test_cov$END)
expect_equal(SARC@assays@unlistData$ID, test_cov$ID)

#and test some samples are the same
expect_equal(SARC@assays@unlistData$SampleA, test_cov$SampleA)
expect_equal(SARC@assays@unlistData$SampleB, test_cov$SampleB)

#test cnv file is embedded in metadata as nested object
expect_equal(metadata(SARC)[[1]][[1]], test_cnv)

#test regionSet makes two lists
expect_equal(length(names(metadata(SARC))), 3)
SARC <- plotCovPrep(RE = SARC, cnv = metadata(SARC)[[1]][[1]],
                   startlist = metadata(SARC)[[2]],
                   endlist = metadata(SARC)[[3]])

SARC <- regionGrangeMake(RE = SARC, covprepped = metadata(SARC)[[4]])

#test the nesting of GRange objects
expect_equal(class(metadata(SARC)[[5]][[1]][[1]])[1], "GRanges")

if (requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
require("TxDb.Hsapiens.UCSC.hg38.knownGene")
} else {}

if (requireNamespace("Homo.sapiens", quietly = TRUE)) {
require("Homo.sapiens")
} else {}

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
TxDb(Homo.sapiens) <- TxDb.Hsapiens.UCSC.hg38.knownGene
tx <- transcriptsBy(Homo.sapiens, columns = "SYMBOL")
txgene <- tx@unlistData

SARC <- addExonsGenes(RE = SARC, covgranges = metadata(SARC)[[5]],
                      txdb = txdb, txgene = txgene)

#same regions of DNA are being used
expect_equal(names(metadata(SARC)[[6]][[1]]), names(metadata(SARC)[[5]][[1]]))
#check FCGR2A is the gene for the first CNV
x <- metadata(SARC)[[6]][[1]][[1]]
expect_equal(x$SYMBOL[[1]],"FCGR2A")

#save for plotting tests - will be hashed out
#saveRDS(SARC, file = "tests/testthat/test1.rds")
