#test plotting
library(SARC)
library(testthat)
data("test_bed")
#data("test_cov")

test_bed <- test_bed[c(1),]

SARC <- readRDS("test1.rds")

#setuplot
SARC <- setupCNVplot(MA = SARC, namedgranges = metadata(SARC)[[5]],
                     covprepped = metadata(SARC)[[3]])
#expect dataframes
expect_equal(class(metadata(SARC)[[6]][[1]]), "data.frame")

names <- c("SAMPLE","READ","START","GENE","EXON","FEATURE")

expect_equal(colnames(metadata(SARC)[[6]][[1]]), names)

#plot
p <- plotCNV(bed = test_bed, setup = metadata(SARC)[[6]], FilteredCNV = 1)
expect_equal(class(p)[[1]], "gg")
