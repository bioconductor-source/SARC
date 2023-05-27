#' seeDist
#'
#' @description Displays the distribution of mean scores from read depths from all the samples. This is a quick method of checking why some detected CNVs might be false positives. We expect true duplications to have very high read-depths and true deletions to have very low read-depths.
#'
#' @param meanList List of dataframes which show the ranked mean scores for each CNV. Stored in metadata after SARC::setDPlot.
#' @param bed Bed file containing CNVs which the user wishes to generate plots for. It is recommended that the most recently created bed file is used. Check print(MA) to see more bed files created by SARC.
#' @param sample Which CNV/ row from the bed file should be checked. Default is 1.
#' @param plotly Should plotly be used - this could be useful when interested in seeing the samples at each point.
#' @param colourCNV Colour of the sample which had the CNV detected. Default is red.
#' @param size Size if the dots. Default is 2.
#'
#' @return A scatter graph which shows the mean score of read-depths for a particular region of DNA where a CNV was detected. The sample which had the CNV detected is coloured as red/ used selected colour.
#' @export
#'
#' @importFrom plotly ggplotly
#'
#' @examples
#' data("test_bed")
#' data("test_cov")
#' SARC <- regionSet(bed = test_bed, cov = test_cov)
#' SARC <- regionSplit(MA = SARC, bed = test_bed, cov = test_cov,
#'                     startlist = metadata(SARC)[[1]],
#'                     endlist = metadata(SARC)[[2]])
#' SARC <- regionMean(MA = SARC, bed = test_bed, splitcov = metadata(SARC)[[3]])
#' SARC <- regionQuantiles(MA = SARC, bed = experiments(SARC)[[1]],
#'                         meancov = metadata(SARC)[[3]], q1 =.1, q2 = .9)
#' SARC <- setQDplot(MA = SARC, meancov = metadata(SARC)[[4]])
#' p <- seeDist(meanList = metadata(SARC)[[7]], bed = experiments(SARC)[[2]],
#'  plotly=FALSE, sample=1)
seeDist <- function(meanList, bed, sample, plotly=FALSE, colourCNV="red", size=2){

  if (missing(meanList)) stop('meanList is missing. Add list of cov files with region means calculated from SARC::setDPlot.')

  if (missing(bed)) stop('bed is missing. Add bed dataframe. Ideally the most recently created bed file should be used.')


  RANK <- MEANSCORE <- SAMPLE <- NULL

  x <- meanList[[sample]]

  #double check the col names are correct

  colnames(x) <- c("MEANSCORE", "SAMPLE", "RANK")

  #get the sample of interest

  y <- bed$SAMPLE[sample]

  z <- x[which(x$SAMPLE == y),]

  #plot the scatter graph

  if (plotly==FALSE) {
    ggplot(data = x, aes(x = RANK, y=MEANSCORE, group=SAMPLE))+

      geom_point(color="black", size=size)+

      geom_point(data=z, aes(x=RANK, y=MEANSCORE), colour=colourCNV, size=size)+

      labs(title=paste0("Distribution of samples at ", bed$CHROM[sample], ":",
                        bed$START[sample], "-",  bed$END[sample]),

            subtitle = paste0(bed$SAMPLE[sample], "  ", bed$TYPE[sample], "=",
                             bed$VALUE[sample]),

           y = "MEAN SCORE", x= "Samples Ranked by Mean Score")+

      theme_classic()+

      theme(axis.text.x=element_blank(),

            axis.text.y=element_text(size=12, face="bold"),

            axis.ticks.y = element_blank(),

            axis.ticks.x = element_blank(),

            axis.title.x=element_text(size=16, face="bold"),

            axis.title.y=element_text(size=16, face="bold"),

            panel.grid.major = element_blank(),

            panel.grid.minor = element_blank(),

            panel.spacing.y = unit(0, 'lines'),

            plot.title = element_text(size=18, face="bold",
                                      vjust = 1, hjust = 0.1),

            plot.subtitle = element_text(size=15, face="bold",
                                         vjust = -1, hjust = 0.5))
  }

  #plot with plotly

  else if (plotly == TRUE) {

    plotly::ggplotly(

      ggplot(data = x, aes(x = RANK, y=MEANSCORE, group=SAMPLE))+

        geom_point(color="black", size=size)+

        geom_point(data=z, aes(x=RANK, y=MEANSCORE), colour=colourCNV, size=size)+

        labs(title=paste0("Distribution of samples at ", bed$CHROM[sample], ":",
                          bed$START[sample], "-",  bed$END[sample]),

             subtitle = paste0(bed$SAMPLE[sample], "  ", bed$TYPE[sample], "=",
                               bed$VALUE[sample]),

             y = "MEAN SCORE", x= "Samples Ranked by Mean Score")+

        theme_classic()+

        theme(axis.text.x=element_blank(),

              axis.text.y=element_text(size=12, face="bold"),

              axis.ticks.y = element_blank(),

              axis.ticks.x = element_blank(),

              axis.title.x=element_text(size=16, face="bold"),

              axis.title.y=element_text(size=16, face="bold"),

              panel.grid.major = element_blank(),

              panel.grid.minor = element_blank(),

              panel.spacing.y = unit(0, 'lines'),

              plot.title = element_text(size=16, face="bold",
                                        vjust = 1),

              plot.subtitle = element_text(size=15, face="bold",
                                           vjust = -1)))
  }

}
