#test stats
#this script will test the stats performed on the first CNV are as expected

library(SARC)
library(testthat)

data("test_cnv")
data("test_cov")

SARC <- regionSet(cnv = test_cnv, cov = test_cov)
SARC <- regionSplit(RE = SARC, cnv = metadata(SARC)[[1]][[1]],
                     startlist = metadata(SARC)[[2]],
                      endlist = metadata(SARC)[[3]])
#calculate mean scores
SARC <- regionMean(RE = SARC, cnv = metadata(SARC)[[1]][[1]],
                   splitcov = metadata(SARC)[[4]])
#check calculation of first CNV
x <- metadata(SARC)[[1]][[2]]
y <- as.data.frame(metadata(SARC)[[4]][[1]])
y <- y[,-c(1:6)]
c <- colMeans(y)
c <- c/mean(c)
expect_equal(as.numeric(round(c[1],2)), x$MeanScore[1])

#check quantile distributions
SARC <- regionQuantiles(RE = SARC, cnv = metadata(SARC)[[1]][[2]],
                        meancov = metadata(SARC)[[5]], q1 = 0.1,
                        q2 = 0.9)
#check calculation of first deletion and duplication
#should be 1
x <- metadata(SARC)[[1]][[3]]

y <- metadata(SARC)[[5]][[1]]
y <- as.data.frame(replace(y, y==0, 0.001))
y <- as.matrix(y)
l <- nrow(y)
q <- as.numeric(quantile(y[l,], probs = c(0.1, 0.9)))
r <- as.numeric(y[nrow(y),])
yhigh <- which(r >= q[2])
xhigh <- colnames(y[,yhigh])

expect_equal(xhigh[1], x$SAMPLE[1])
expect_equal("1", x$Qhigh[1])

#should be 0
y <- metadata(SARC)[[5]][[2]]
y <- as.data.frame(replace(y, y==0, 0.001))
y <- as.matrix(y)
l <- nrow(y)
q <- as.numeric(quantile(y[l,], probs = c(0.1, 0.9)))
r <- as.numeric(y[nrow(y),])
ylow <- which(r <= q[1])
xlow <- colnames(y[,ylow])


expect_false(isTRUE(all.equal(xlow[1], x$SAMPLE[2])))
expect_equal("0", x$Qlow[2])

#test anova
SARC <- prepAnova(RE = SARC, cnv = metadata(SARC)[[1]][[3]],
                 start = metadata(SARC)[[2]], end=metadata(SARC)[[3]])
SARC <- anovaOnCNV(RE = SARC, cnv = metadata(SARC)[[1]][[3]],
                  anovacov = metadata(SARC)[[8]])
x <- metadata(SARC)[[1]][[4]]

#check anova on first CNV
y <- as.data.frame(metadata(SARC)[[8]][[1]])
y$LOC <- paste0(y$seqnames, ":", y$start, "-", y$end)
y <- y[,-c(1:6)]
y <- reshape2::melt(y, id.vars = "LOC", variable.name = "SAMPLE",
                    value.name = "READS")
anov <- aov(READS ~ SAMPLE, data = y)
n <- summary(anov)[[1]][,5][[1]]
expect_equal(n, x$ANOVA[1])

#check confidence
SARC <- cnvConfidence(RE = SARC, cnv = metadata(SARC)[[1]][[4]])
x <- metadata(SARC)[[1]][[5]]
#expect columns
expect_equal(ncol(x), 13)
expect_equal(names(x)[12], "CNV.SCORE")
expect_equal(names(x)[13], "CNV.TIER")

#expect stats of first CNV
x$MeanScore[1] #1.44
x$Qhigh[1] #1
x$ANOVA[1] #0.09
#expect confidence == 2
expect_equal(x$CNV.SCORE[1], 2)
