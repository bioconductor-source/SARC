#' regionMean
#'
#' @description Function to calculate MeanScores for each region where a CNV has been detected. Dependent of sequencing quality and WES/WGS input, the estimated MeanScores for CNV types are generally: HOMOZYGOUS DELETION 0-0.3, HETEROZYGOUS DELETION 0.3-0.7, HETEROZYGOUS DUPLICATION 1.3-1.8 HOMOZYGOUS DUPLICATION >1.8.
#'
#' @param MA MultiAssayExperiment object used to store all information.
#' @param bed CNV bed files stored as dataframes. Chromosomes, Start of CNV, End of CNV, Type of CNV and some other information can be kept here. It is recomended to use the most recently generated bed file. Check print(MA) to see if bed files from other anlyses have been generated.
#' @param splitcov List of small .cov dataframes. Would be stored as metadata after SARC::regionSplit.
#'
#' @return A new bed file with an additional column representing the meanScores. A new list of dataframes which contains the meanscores.
#' @export
#' @importFrom utils tail
#'
#' @examples
#' data("test_bed")
#' data("test_cov")
#' SARC <- regionSet(bed = test_bed, cov = test_cov)
#' SARC <- regionSplit(MA = SARC, bed = test_bed, cov = test_cov,
#'                     startlist = metadata(SARC)[[1]],
#'                     endlist = metadata(SARC)[[2]])
#' SARC <- regionMean(MA = SARC, bed = test_bed, splitcov = metadata(SARC)[[3]])
#'
regionMean <- function(MA, bed, splitcov){

  if (missing(MA)) stop('MA is missing. Add a MultiAssayExperiment object to store data efficiently.')

  if (missing(bed)) stop('bed is missing. Add bed dataframe. Ideally the most recently created bed file should be used.')

  if (missing(splitcov)) stop('splitcov is missing. Add list of smaller cov dataframes. Will be created from SARC::regionSplit.')

  #create new column in bed file

  bed$MeanScore <- ""

  lst <- splitcov

  #calculate mean value for each region and then divide this by the mean in each sample

  split.mean <- list()

  split.mean <- lapply(lst, function(x){

    x.s <- x[5:ncol(x)]

    x.m <- colMeans(x.s)

    x.d <- mean(x.m)

    div <- x.m/x.d

    x.s <- rbind(x.s, div)

    return(x.s)

  })

  #add the mean scores back to the bed file as a new column

  b <- bed$SAMPLE

  for (i in seq_len(nrow(bed))) {

    a <- split.mean[[i]]

    b.i <- b[i]

    c.i <- a[b.i]

    c.t <- round(tail(c.i, 1), 2)

    bed$MeanScore[i] <- c.t

  }

  bed$MeanScore <- unlist(bed$MeanScore)

  #store as a new dataframe and a new list

  MA2 <- suppressWarnings(suppressMessages(MultiAssayExperiment(list(BEDMEAN = bed))))

  metadata(MA)[["SPLITMEAN"]] <- split.mean

  MA <- suppressWarnings(c(MA, MA2))

  #return new MA object

  return(MA)

}
