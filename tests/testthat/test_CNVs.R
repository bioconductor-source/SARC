#test if the correct regions are split
library(testthat)
data("test_cnv")
data("test_cov")

#for speed we will just keep one CNV
SARC <- regionSet(cnv = test_cnv, cov = test_cov)

S <- metadata(SARC)[[2]]
E <- metadata(SARC)[[3]]

test_cov <- cbind(rownames(test_cov), test_cov)
colnames(test_cov)[1] <- "ID"

#test first CNV
S1 <- test_cov$START[S[1]]
E1 <- test_cov$END[E[1]]

#expect the difference between what the detection algorithm read and the bed file is below 10 nucleotide
#for start it is usually 1 and for end it is usually 0
expect_equal(test_cnv$START[1] -  S1, 1)

expect_equal(test_cnv$END[1] -  E1, 0)

#test final CNV too
S2 <- test_cov$START[S[15]]
E2 <- test_cov$END[E[15]]

expect_equal(test_cnv$START[15] -  S2, 1)

expect_equal(test_cnv$END[15] -  E2, 0)
