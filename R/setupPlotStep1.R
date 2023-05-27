#' setupPlotStep1
#'
#' @description Internal function which will help to set up plots for CNVs.
#'
#' @param namedgranges List of grange objects with exons and genes data added to them. These objects should have exons and genes, and will be in metadata after being created by SARC::addExonsGenes.
#' @noRd
#'
setupPlotStep1 <- function(namedgranges){

  if (missing(namedgranges)) stop('namedgranges is missing. Add a list of named grange objects for plotting. Should be created from SARC::addExonsGenes.')

  #extract samples as data frames

  dflist <- lapply(namedgranges, function(y){

    lapply(y, function(z){

      df <- as.data.frame(z)

    })
  })

  #remove any samples which have 0 values

  dflist <- lapply(dflist, function(x){

    x[sapply(x, nrow)>0]

  })

  #Add columns to dataframes - SYMBOLS and GENE_EXON, from the grange objects

  for (i in seq_len(length(dflist))) {

    y <- dflist[[i]]

    for (j in seq_len(length(y))) {

      dflist[[i]][[j]]$SYMBOL <- as.character(dflist[[i]][[j]]$SYMBOL)

      dflist[[i]][[j]]$GENE_EXON <- paste0(dflist[[i]][[j]]$SYMBOL, ":",

                                           dflist[[i]][[j]]$exon_rank)

      dflist[[i]][[j]] <- dflist[[i]][[j]][!duplicated(dflist[[i]][[j]]$GENE_EXON),]

    }

  }

  dflist <- lapply(dflist, function(y){

    df <- do.call("rbind", y)

  })

  #return list of dataframes

  return(dflist)
}
