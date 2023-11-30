#' scan_plot
#'
#' Plot spatial hot-spot detected by scan statistics (Kulldorff, 1997).
#'
#' @param data N (space) X M (time) matrix
#' @param shp shape file with identical space order with data.
#' @param id list of cluster detected by scan.seq function
#' @param wd working directory for saving plots

scan_plot <- function(data, shp, id, wd) {

  library(RColorBrewer)

  if (is.null(id)) {
    stop("id is null, so implement scan statistics using scan.seq function")
  }

  n.space <- nrow(shp)

  setwd(wd)
  col1 <- c(); temp <-0
  for (i in 1:length(id)) {

    temp <- temp + 1
    index <- which(colnames(data) == unique(id[[i]]$week)) # in case of no cluster week
    counts <- data[,index]
    cluster2 <- id[[i]]$shp.order

    jpeg(paste("cluster_",unique(id[[i]]$week),".jpeg",sep = ""), res=600, width=200, height=150, pointsize=9,units='mm')
    brks <- cut(counts, breaks=c(0,210,490,1010,1600,max(counts)+1), include.lowest = TRUE, right=FALSE)
    cols <- brewer.pal(5, "Reds")
    for (k in 1:n.space) {col1[k] <- cols[brks[k]]}
    plot(shp$geometry, col=col1,axes=F, border="gray")
    plot(shp$geometry[cluster2], border="black", add=TRUE)
    title(paste("Cluster map for", unique(id[[i]]$week), sep=" "))
    legend('bottomright',legend=c("[0, 210)", "[210, 490)", "[490, 1010)", "[1010, 1600)", "≥ 1600"),fill=cols,
           cex=0.7)
    dev.off()

    print(temp)
  }

}
