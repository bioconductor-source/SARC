#test addExonsGenes
library(SARC)
library(testthat)
data("test_bed")
data("test_cov")

test_bed <- test_bed[c(1),]

SARC <- regionSet(bed = test_bed, cov = test_cov)

#test regionSet makes two lists
expect_equal(length(names(metadata(SARC))), 2)

SARC <- plotCovPrep(MA = SARC, bed = test_bed, cov = test_cov,
                   startlist = metadata(SARC)[[1]],
                   endlist = metadata(SARC)[[2]])

SARC <- regionGrangeMake(MA = SARC, covprepped = metadata(SARC)[[3]])

#make GRange objects
expect_equal(class(metadata(SARC)[[4]][[1]][[1]])[1], "GRanges")

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Homo.sapiens)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
TxDb(Homo.sapiens) <- TxDb.Hsapiens.UCSC.hg38.knownGene
tx <- transcriptsBy(Homo.sapiens, columns = "SYMBOL")
txgene <- tx@unlistData

SARC <- addExonsGenes(MA = SARC, covgranges = metadata(SARC)[[4]],
                      txdb = txdb, txgene = txgene)

#same regions of DNA are being used
expect_equal(names(metadata(SARC)[[5]][[1]]), names(metadata(SARC)[[4]][[1]]))

#save for more tests
#saveRDS(SARC, file = "tests/testthat/test1.rds")
