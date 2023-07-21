#' regionMean
#'
#' @description Function to calculate MeanScores for each region where a CNV has been detected. Dependent of sequencing quality and WES/WGS input, the estimated MeanScores for CNV types are generally: HOMOZYGOUS DELETION 0-0.3, HETEROZYGOUS DELETION 0.3-0.7, HETEROZYGOUS DUPLICATION 1.3-1.8 HOMOZYGOUS DUPLICATION >1.8.
#'
#' @param RE RaggedExperiment object used to store all information.
#' @param cnv List of CNVs in a dataframe containing CNVs from detection algorithms/ pipelines. It is recomended to use the most recently generated cnv file. Check print(RE) to see if cnv files from other anlyses have been generated.
#' @param splitcov List of small .cov dataframes. Would be stored as metadata after SARC::regionSplit.
#' @param nameofnewdf Name of new dataframe to be saved in metadata(RE)[['CNVlist']]. Default is CNVmeans.
#'
#' @return A new cnv file with an additional column representing the meanScores. A new list of dataframes which contains the meanscores.
#' @export
#' @importFrom utils tail
#'
#' @examples
#' data("test_cnv")
#' data("test_cov")
#' SARC <- regionSet(cnv = test_cnv, cov = test_cov)
#' SARC <- regionSplit(RE = SARC, cnv = metadata(SARC)[['CNVlist']][[1]],
#'                     startlist = metadata(SARC)[[2]],
#'                     endlist = metadata(SARC)[[3]])
#' SARC <- regionMean(RE = SARC, cnv = metadata(SARC)[['CNVlist']][[1]],
#'                    splitcov = metadata(SARC)[[4]])
#'
regionMean <- function(RE, cnv, splitcov, nameofnewdf="CNVmeans"){

  if (missing(RE)) stop('RE is missing. Add a RaggedExperiment object to store data efficiently.')

  if (missing(cnv)) stop('cnv is missing. Add cnv dataframe. Ideally the most recently created cnv file should be used.')

  if (missing(splitcov)) stop('splitcov is missing. Add list of smaller cov dataframes. Will be created from SARC::regionSplit.')

  #create new column in cnv file

  cnv$MeanScore <- ""

  lst <- splitcov

  #calculate mean value for each region and then divide this by the mean in each sample

  split.mean <- list()

  split.mean <- lapply(lst, function(y){

    x <- as.data.frame(y@elementMetadata)

    x$ID <- NULL

    x.m <- colMeans(x)

    x.d <- mean(x.m)

    div <- x.m/x.d

    x.s <- rbind(x, div)

    return(x.s)

  })

  #add the mean scores back to the cnv file

  b <- cnv$SAMPLE

  for (i in seq_len(nrow(cnv))) {

    a <- split.mean[[i]]

    b.i <- b[i]

    c.i <- a[b.i]

    c.t <- round(tail(c.i, 1), 2)

    cnv$MeanScore[i] <- c.t

  }

  cnv$MeanScore <- unlist(cnv$MeanScore)

  #store output as a new dataframe and a new list
  metadata(RE)[["CNVlist"]][[nameofnewdf]] <- cnv

  metadata(RE)[["SPLITMEAN"]] <- split.mean

  #return new RE object

  return(RE)

}
