#' setupPlotStep2fun
#'
#' @description Internal function to help setupPlotStep2, for the setting up of plots for CNVs.
#'
#' @param x List of prepped dataframes for plotting.
#'
#' @noRd
#'
setupPlotStep2fun <- function(x){

  if (missing(x)) stop('x is missing. This is a dataframe within a list being processed by a lapply function in SARC::setupPlotStep2.')

  #reshape the dataframe

  ml <-  reshape2::melt(x, id.vars = "LOC")

  #create new dataframes with needed information

  ml2 <- data.frame(do.call('rbind', strsplit(as.character(ml$LOC),'_',fixed=TRUE)))

  ml3 <- data.frame(do.call('rbind', strsplit(as.character(ml2$X2),':',fixed=TRUE)))

  ml <- cbind(ml$variable, ml$value, ml2$X1, ml3)

  #rename dataframe

  colnames(ml) <- c("SAMPLE", "READ", "START", "GENE", "EXON")

  ml$FEATURE <- paste0(ml$GENE,".",ml$EXON)

  #return object

  return(ml)
}
