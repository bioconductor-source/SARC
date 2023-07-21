#test plotting
#this is a quicker test to just check the plot is being created
library(SARC)
library(testthat)
data("test_cnv")

test_cnv <- test_cnv[c(1),]
SARC <- readRDS("test1.rds")

#setuplot
SARC <- setupCNVplot(RE = SARC, namedgranges = metadata(SARC)[[6]],
                     covprepped = metadata(SARC)[[4]])
#expect dataframes
expect_equal(class(metadata(SARC)[[7]][[1]]), "data.frame")
names <- c("SAMPLE","READ","START","GENE","EXON","FEATURE")
expect_equal(colnames(metadata(SARC)[[7]][[1]]), names)

#plot
p <- plotCNV(cnv = test_cnv, setup = metadata(SARC)[[7]], FilteredCNV = 1)
expect_equal(class(p)[[1]], "gg")
#check the title/ regions are the same
expect_equal(p$labels$title, paste0(test_cnv$CHROM,":",test_cnv$START,":",test_cnv$END))
