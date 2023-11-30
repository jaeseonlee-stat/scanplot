
#' scan.seq
#'
#' Implement scan statistics (Kulldorff, 1997) for spatio-temporal data.
#'
#'
#' @param data N (space) X M (time) matrix
#' @param shp shape file with identical space order with data.
#' @param shp.name N X 2 matrix with (state name, city name). Order must be identical to shape file!
#' @param centroid N X 2 matrix with (latitude, longitude)
#' @param pop.upper.bound default is 0.2
#' @param n.simulations default is 999
#' @param alpha.level default is 0.05
#' @param wd working directory for saving plots

scan.seq <- function(data, shp, shp.name, centroid, pop.upper.bound = .2, n.simulations = 999, alpha.level = .05, wd) {

  n.space <- nrow(shp)

  library(SpatialEpi)
  library(scanstatistics)

  cent <- centroid
  geo <- SpatialEpi::latlong2grid(cent)
  knn_mat <- scanstatistics::coords_to_knn(geo, 15)
  zones <- scanstatistics::knn_zones(knn_mat)

  setwd(wd)
  temp <-0; id <- list(); scan <- 0
  for (i in 1:(ncol(data))) {

    col1 <- c()
    temp <- temp +1
    cluster2 <- c()
    counts <- as.numeric(data[, i])

    if (length(counts[counts==0])/n.space > 0.3) {
      options(show.error.messages = F)
      try(scan <- scan_eb_zip(counts=counts, zones=zones, n_mcsim=n.simulations,
                              rel_tol=1e-3, population=rep(1,n.space)))
      options(show.error.messages = T)

      if (class(scan)!="scanstatistic") {next}

      cluster <- scan$MLC$locations
      for (j in 1:length(cluster)) {
        if (scan$MC_pvalue<0.05) {cluster2 <- c(cluster2,cluster[j])}
      }
      RR <- round(scan$MLC$relative_risk,2)
      p.val <- scan$MC_pvalue

    } else {
      exp.cases <- expected(rep(1, n.space), counts, 1)
      scan <- kulldorff(geo, counts, rep(1, n.space), expected.cases = exp.cases, pop.upper.bound,
                        n.simulations, alpha.level, FALSE)
      cluster <- scan$most.likely.cluster$location.IDs.included

      for (j in 1:length(cluster)) {
        if (scan$most.likely.cluster$p.value < 0.05) {cluster2 <- c(cluster2,cluster[j])}
      }
      RR <- round(scan$most.likely.cluster$SMR,2)
      p.val <- scan$most.likely.cluster$p.value

    }

    No.cases <- counts[cluster2]
    Exp.cases <- round(length(cluster2)/length(counts)*sum(counts),2)
    total <- sum(counts)
    cluster.cases <- sum(counts[cluster2])
    percent <- round(sum(counts[cluster2])/sum(counts)*100,2)

    if (length(cluster2)!=0) {
      id[[temp]] <- data.frame(week=rep(colnames(data)[i],length(cluster2)),
                               state=shp.name[cluster2, 1], city=shp.name[cluster2, 2],
                               No_of_cases=No.cases, RR=RR, Exp.cases=Exp.cases, pvalue=p.val,
                               total=total, cluster.cases=cluster.cases, percent=percent, shp.order=cluster2)
      id[[temp]] <- id[[temp]][order(id[[temp]][,"state"]),]
    }

    scan <- 0

    print(temp)
  }

  capture.output(id, file="cluster info_covid_korea.txt")

  # for csv file
  tmp_df <- data.frame()
  for (i in 1:length(id)) {
    tmp <- id[[i]]
    if (!is.null(tmp)) {
      tmp_df <- rbind(tmp_df, cbind(id = i, tmp))
    }
  }

  tmp_df[order(tmp_df$week, -tmp_df$No_of_cases), ]
  tmp_table <- as.data.frame(table(tmp_df$id))
  tmp_df <- cbind(frequency = tmp_table$Freq[match(tmp_df$id, tmp_table$Var1)], tmp_df)
  write.csv(tmp_df, 'result.csv', row.names = F, fileEncoding = 'euc-kr')

  return(id)
}

#' scan.plot
#'
#' Plot spatial hot-spot detected by scan statistics (Kulldorff, 1997).
#'
#' @param data N (space) X M (time) matrix
#' @param shp shape file with identical space order with data.
#' @param id list of cluster detected by scan.seq function
#' @param wd working directory for saving plots

scan.plot2 <- function(data, shp, id, wd) {

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
