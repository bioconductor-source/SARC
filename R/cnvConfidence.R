#' cnvConfidence
#'
#' @description Flags the confidence of the CNV being a true CNV. The flags range from Very Unconfident - Very Confident. Our clinicians preferred variants to be flagged rather than filtered out - so we simply do this. Confidence scores will be added to the cnv file given. The cnv file should be one which contains all the statistical results from other functions of the SARC package.
#'
#' @param RE RaggedExperiment object used to store all information.
#' @param cnv List of CNVs in a dataframe containing CNVs from detection algorithms/ pipelines. Must use one which contains at least MeanScore, Qlow, Qhigh and anova p-values. If Post-hoc tests have also been performed, these results can also be used by this function.
#' @param ph If Dunnet tests were also performed, the p-value from the tests can be taken into account. Default is FALSE.
#' @param m1 Value used to cut-off the highest mean score allowed for HOMOZYGOUS DELECTIONS, and lowest cut-off for HETEROZYGOUS DELETIONS. Default is 0.2.
#' @param m2 Value used to cut-off the highest mean score allowed for HETEROZYGOUS DELETINS. Default is 0.8.
#' @param m3 Value used to cut-off the lowest mean score allowed for HETEROZYGOUS DUPLICATIONS. Default is 1.2.
#' @param m4 Value used to cut-off the lowest mean score allowed for HETEROZYGOUS DUPLICATIONS. Default is 1.8.
#'
#' @return A new cnv file with additional columns which describe our level of confidence of the detected CNV being a true CNV.
#' @export
#'
#' @examples
#' data("test_cnv")
#' test_cnv <- test_cnv[c(1:3),]
#' data("test_cov")
#' SARC <- regionSet(cnv = test_cnv, cov = test_cov)
#' SARC <- regionSplit(RE = SARC, cnv = metadata(SARC)[['CNVlist']][[1]],
#'                      startlist = metadata(SARC)[[2]],
#'                       endlist = metadata(SARC)[[3]])
#' SARC <- regionMean(RE = SARC, cnv = metadata(SARC)[['CNVlist']][[1]],
#'                   splitcov = metadata(SARC)[[4]])
#' SARC <- regionQuantiles(RE = SARC, cnv = metadata(SARC)[['CNVlist']][[2]],
#'                         meancov = metadata(SARC)[[3]], q1 =.1, q2 = .9)
#' SARC <- prepAnova(RE = SARC, cnv = metadata(SARC)[['CNVlist']][[3]],
#'                  start = metadata(SARC)[[2]], end=metadata(SARC)[[3]])
#' SARC <- anovaOnCNV(RE = SARC, cnv = metadata(SARC)[['CNVlist']][[3]],
#'                   anovacov = metadata(SARC)[[8]])
#' SARC <- cnvConfidence(RE = SARC, cnv = metadata(SARC)[['CNVlist']][[4]])
cnvConfidence <- function(RE, cnv, ph=FALSE, m1=0.2, m2=0.8, m3=1.2, m4=1.8){

  if (missing(RE)) stop('RE is missing. Add a RaggedExperiment object to store data efficiently.')

  if (missing(cnv)) stop('cnv is missing. Add cnv dataframe. Ideally the most recently created cnv file should be used.')

  #change scientific numeric numbers to numerals
  options(scipen=99)

  x <- cnv

  #contrast pipleine value to mean scores

  x$MS.SCORE <- ""

  #HOM DEL
  x$MS.SCORE <- ifelse(x$VALUE == 0 & x$MeanScore <= m1, 1, 0)

  #HET DEL
  x$MS.SCORE <- ifelse(x$VALUE == 1 & x$MeanScore >= m1 & x$MeanScore <= m2,
                       1, x$MS.SCORE)

  #HET DUP
  x$MS.SCORE <- ifelse(x$VALUE == 3 & x$MeanScore >= m3, 1, x$MS.SCORE)

  #HOM DUP
  x$MS.SCORE <- ifelse(x$VALUE >= 4 & x$MeanScore >= m4, 1, x$MS.SCORE)

  #score based on quantile distribution.

  x$QD.SCORE <- ""

  x$QD.SCORE <- ifelse(x$TYPE == "DEL" & x$Qlow == 1, 1, 0)

  x$QD.SCORE <- ifelse(x$TYPE == "DUP" & x$Qhigh == 1, 1, x$QD.SCORE)

  #score based on anova

  x$A.SCORE <- ""

  x$A.SCORE <-  ifelse(x$ANOVA <= 0.05, 1, 0)

  #if no posthoc was performed, tally up and remove columns

  if (ph == FALSE) {

    x$CNV.SCORE <- x$MS.SCORE + x$QD.SCORE + x$A.SCORE

    x$MS.SCORE <- x$QD.SCORE <- x$A.SCORE <- NULL

  }

  #if posthoc was performed add this p-value to scores too

  else if (ph == TRUE) {

    x$D.SCORE <- ""

    x$D.SCORE <-  ifelse(x$Dunnet <= 0.05, 1, 0)

    x$CNV.SCORE <- x$MS.SCORE + x$QD.SCORE + x$A.SCORE + x$D.SCORE

    x$MS.SCORE <- x$QD.SCORE <- x$A.SCORE <- x$D.SCORE <- NULL

  }

  #restore acceptable digits
  options(scipen=3)

  #rank scores of CNVs

  x$CNV.TIER <- ""

  x$CNV.TIER <- ifelse(x$CNV.SCORE == 0, "VERY UNCONFIDENT", x$CNV.TIER)

  x$CNV.TIER <- ifelse(x$CNV.SCORE == 1, "UNCONFIDENT", x$CNV.TIER)

  x$CNV.TIER <- ifelse(x$CNV.SCORE == 2, "CONFIDENT", x$CNV.TIER)

  x$CNV.TIER <- ifelse(x$CNV.SCORE == 3, "VERY CONFIDENT", x$CNV.TIER)

  x$CNV.TIER <- ifelse(x$CNV.SCORE >= 4, "SUPER CONFIDENT", x$CNV.TIER)

  #add df to RE object

  metadata(RE)[["CNVlist"]][["CNVclass"]] <- x

  #return RE object

  return(RE)
}
