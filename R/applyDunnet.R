#' applyDunnet
#'
#' @description Internal function used to calculate the dunnet scores for each CNV found. Dunnet is a powerful but slow post-hoc test. It is recomended when there are fewer than 100 samples in the COV file, and less than 100 CNVs to test.
#'
#' @param preppedcov List of dataframes, each df is a prepped anova df for each region where a CNV was detected. Should be in the MA object as metadata, after the prepAnova function.
#' @param control Single Sample where the CNV was detected. Will be used as control for pairwise tests. This is automatically performed by phDunnetonCNV as a loop.
#'
#' @noRd
#' @importFrom DescTools DunnettTest
#' @importFrom metap sumlog
#' @importFrom stats relevel
#' @import multtest
applyDunnet <- function(preppedcov, control){

  if (missing(preppedcov)) stop('preppedcov is missing. Add dataframe from SARC:prepAnova. Should be stored as metadata..')

  if (missing(control)) stop('control is missing. Add string of sample which had the CNV in question identified in it.')

  df <- preppedcov

  #simplify a catagorical column

  df$LOC <- paste0(df$CHROM, ":", df$START, "-", df$END)

  #remove unneeded columns from cov file
  df <- df[5:(ncol(df))]

  #re-organise to have thre columns : CNV, sample, value(read depth)
  y <- melt(df, id.vars = "LOC", variable.name = "SAMPLE", value.name = "READS")

  #use a specific sample (bed file) as the control to compare against each other sample
  y$Sample <- relevel(y$SAMPLE, ref = control)

  #apply the dunnet test
  D <- DunnettTest(x = y$READS, g = y$SAMPLE, control = control)

  #retreive all p-values
  ddt <- as.data.frame(D[[1]])

  #average p-value
  dun <- sumlog(ddt$pval)[[3]]

  #return averaged p-value

  return(dun)

}
