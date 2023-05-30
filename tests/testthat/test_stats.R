#test stats

library(SARC)
library(testthat)

data("test_bed")
data("test_cov")

SARC <- regionSet(bed = test_bed, cov = test_cov)
SARC <- regionSplit(MA = SARC, bed = test_bed, cov = test_cov,
                     startlist = metadata(SARC)[[1]],
                      endlist = metadata(SARC)[[2]])
SARC <- regionMean(MA = SARC, bed = test_bed, splitcov = metadata(SARC)[[3]])
SARC <- regionQuantiles(MA = SARC, bed = experiments(SARC)[[1]],
                        meancov = metadata(SARC)[[3]], q1 =.1, q2 = .9)
SARC <- prepAnova(MA = SARC, bed = experiments(SARC)[[2]], cov = test_cov,
                 start = metadata(SARC)[[1]], end=metadata(SARC)[[2]])
SARC <- anovaOnCNV(MA = SARC, bed = experiments(SARC)[[2]],
                  anovacov = metadata(SARC)[[7]])
SARC <- cnvConfidence(MA = SARC, bed = experiments(SARC)[[3]])


x <- experiments(SARC)[[4]]

#expect columns
expect_equal(ncol(x), 13)
expect_equal(names(x)[12], "CNV.SCORE")
expect_equal(names(x)[13], "CNV.TIER")
