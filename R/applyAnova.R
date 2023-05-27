#' applyAnova
#'
#' @description Internal function to apply anova using each prepped anova data frame inside the list created by prepAnova.
#'
#' @param df List of prepped prepped anova data frames from the prepAnova function. Will be stored in the MA object after the prepAnova function.
#'
#' @noRd
#' @importFrom stats aov
#' @importFrom reshape2 melt
applyAnova <- function(df){

  if (missing(df)) stop('df is missing. Add dataframe from SARC:prepAnova. Should be stored as metadata.')

  df$LOC <- paste0(df$CHROM, ":", df$START, "-", df$END)

  #Remove unneeded columns

  df <- df[5:(ncol(df))]

  #Reshape to perform the anova formula

  y <- melt(df, id.vars = "LOC", variable.name = "SAMPLE", value.name = "READS")

  #perform anova using the READS and SAMPLE in the formula

  anov <- aov(READS ~ SAMPLE, data = y)

  #store the p values from anova

  n <- summary(anov)[[1]][,5][[1]]

}
