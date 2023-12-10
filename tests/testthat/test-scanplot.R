test_that("Function works", {

  library(testthat)
  library(SpatialEpi)
  library(scanstatistics)
  library(utils)
  library(sf)
  library(RColorBrewer)
  # library(grDevices)
  library(graphics)
  library(stringr)
  library(stats)
  library(scanplot)

  # Function to implement spatial scan statistics for multiple time points
  scan_seq <- function(data, shp.name, centroid, pop.upper.bound = .2, n.simulations = 999, alpha.level = .05, save = FALSE) {

    n.space <- nrow(data)

    cent <- centroid
    geo <- SpatialEpi::latlong2grid(cent)
    knn_mat <- scanstatistics::coords_to_knn(geo, 15)
    zones <- scanstatistics::knn_zones(knn_mat)

    temp <-0; id <- list(); scan <- 0
    for (i in 1:(ncol(data))) {

      col1 <- c()
      temp <- temp +1
      cluster2 <- c()
      counts <- as.numeric(data[, i])

      if (length(counts[counts==0])/n.space > 0.3) {
        options(show.error.messages = F)
        try(scan <- scanstatistics::scan_eb_zip(counts = matrix(counts, nrow = 1), zones = zones, n_mcsim = n.simulations,
                                                rel_tol = 1e-3, population = matrix(rep(1, n.space), nrow=1)))
        options(show.error.messages = T)

        if (!inherits(scan, "scanstatistic")) {next}

        cluster <- scan$MLC$locations
        for (j in 1:length(cluster)) {
          if (scan$MC_pvalue < 0.05) {cluster2 <- c(cluster2,cluster[j])}
        }
        RR <- round(scan$MLC$relative_risk,2)
        p.val <- scan$MC_pvalue

      } else {
        exp.cases <- SpatialEpi::expected(rep(1, n.space), counts, 1)
        scan <- SpatialEpi::kulldorff(geo, counts, rep(1, n.space), expected.cases = exp.cases, pop.upper.bound,
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

      if (length(cluster2) != 0) {
        id[[temp]] <- data.frame(week = rep(colnames(data)[i], length(cluster2)),
                                 state = shp.name[cluster2, 1], county = shp.name[cluster2, 2],
                                 No_of_cases = No.cases, RR = RR, Exp.cases = Exp.cases, pvalue = p.val,
                                 total = total, cluster.cases = cluster.cases, percent = percent, shp.order = cluster2)
        id[[temp]] <- id[[temp]][order(id[[temp]][,"state"]),]
      }

      scan <- 0

      print(temp)
    }

    if (length(id) == 0) {stop("No detected cluster in given data")}

    if (save == TRUE) {

      utils::capture.output(id, file="cluster info.txt")

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
      utils::write.csv(tmp_df, 'cluster info.csv', row.names = F, fileEncoding = 'euc-kr')

    }

    return(id)
  }

  # Function to visualize the results of scan statistics
  scan_plot <- function(data, shp, id, allmap = FALSE) {

    if (is.null(id)) {
      stop("id is null: implement scan statistics using scan_seq() first")
    }

    n.space <- nrow(shp)

    temp <-0
    for (i in 1:ncol(data)) {

      temp <- temp + 1
      counts <- data[ , i]

      if (tryCatch(is.null(id[[i]]), error = function(e) TRUE) & allmap == TRUE) {
        # No detected cluster

        # grDevices::jpeg(paste("cluster_", colnames(data)[i], ".jpeg", sep = ""),
        #                  res = 600, width = 200, height = 150, pointsize = 9,units = 'mm')

        uniq <- sort(unique(data[,i])) # divide the counts into a few groups to color the map

        if (length(table(data[,i])) < 3) {

          cols <- RColorBrewer::brewer.pal(3, "Reds")[-3]
          brks <- cut(data[,i], breaks=c(0, uniq[2], uniq[2]+1), include.lowest = TRUE, right=FALSE, labels=FALSE)
          leg <- cut(data[,i], breaks=c(0, uniq[2], uniq[2]+1), include.lowest = TRUE, right=FALSE)
          stringr::str_sub(levels(leg)[2],-1) <- ")"

        } else if (length(table(data[,i])) == 3) {

          cols <- RColorBrewer::brewer.pal(3, "Reds")
          brks <- cut(data[,i], breaks=c(0, uniq[2], uniq[3], uniq[3]+1), include.lowest = TRUE, right=FALSE, labels=FALSE)
          leg <- cut(data[,i], breaks=c(0, uniq[2], uniq[3], uniq[3]+1), include.lowest = TRUE, right=FALSE)
          stringr::str_sub(levels(leg)[3],-1) <- ")"

        } else {

          cols <- RColorBrewer::brewer.pal(4, "Reds")
          len <- length(unique(data[,i]))
          brks <- cut(data[,i], breaks=unique(stats::quantile(data[,i])), include.lowest = TRUE, right=FALSE, labels=FALSE)
          leg <- cut(data[,i], breaks=unique(stats::quantile(data[,i])), include.lowest = TRUE, right=FALSE)
          # brks <- cut(data[,i], breaks=c(0, uniq[2], uniq[4], uniq[len]+1), include.lowest = TRUE, right=FALSE, labels=FALSE)
          # leg <- cut(data[,i], breaks=c(0, uniq[2], uniq[4], uniq[len]+1), include.lowest = TRUE, right=FALSE)
          stringr::str_sub(levels(leg)[3], -1) <- ")"

        }

        # plotting the map
        col1 <- c()
        for (k in 1:n.space) {col1[k] <- cols[brks[k]]}
        plot(shp$geometry, col = col1, axes = F, border = "gray")
        graphics::title(paste("Incidence map for", colnames(data)[i], sep = " "))
        graphics::legend('bottomright', legend = c("0", levels(leg)[-1]), fill = cols)
        # grDevices::dev.off()

      } else if (tryCatch(!is.null(id[[i]]), error = function(e) FALSE)) {
        # If there are detected clusters

        shp.order <- id[[i]]$shp.order

        # grDevices::jpeg(paste("cluster_", unique(id[[i]]$week), ".jpeg", sep = ""),
        #                 res = 600, width = 200, height = 150, pointsize = 9, units = 'mm')

        # divide the counts into a few groups to color the map
        uniq <- sort(unique(data[,i]))

        if (length(table(data[,i])) < 3) {

          cols <- RColorBrewer::brewer.pal(3, "Reds")[-3]
          brks <- cut(data[,i], breaks = c(0, uniq[2], uniq[2]+1), include.lowest = TRUE, right=FALSE, labels=FALSE)
          leg <- cut(data[,i], breaks=c(0, uniq[2], uniq[2]+1), include.lowest = TRUE, right=FALSE)
          stringr::str_sub(levels(leg)[2],-1) <- ")"

        } else if (length(table(data[,i])) == 3) {

          cols <- RColorBrewer::brewer.pal(3,"Reds")
          brks <- cut(data[,i], breaks=c(0, uniq[2], uniq[3], uniq[3]+1), include.lowest = TRUE, right=FALSE, labels=FALSE)
          leg <- cut(data[,i], breaks=c(0, uniq[2], uniq[3], uniq[3]+1), include.lowest = TRUE, right=FALSE)
          stringr::str_sub(levels(leg)[3],-1) <- ")"

        } else {

          cols <- RColorBrewer::brewer.pal(4, "Reds")
          len <- length(unique(data[,i]))
          brks <- cut(data[,i], breaks=unique(stats::quantile(data[,i])), include.lowest = TRUE, right=FALSE, labels=FALSE)
          leg <- cut(data[,i], breaks=unique(stats::quantile(data[,i])), include.lowest = TRUE, right=FALSE)
          # brks <- cut(data[,i], breaks=c(0, uniq[2], uniq[4], uniq[len]+1), include.lowest = TRUE, right=FALSE, labels=FALSE)
          # leg <- cut(data[,i], breaks=c(0, uniq[2], uniq[4], uniq[len]+1), include.lowest = TRUE, right=FALSE)
          stringr::str_sub(levels(leg)[3], -1) <- ")"

        }

        # plotting the map
        col1 <- c()
        for (k in 1:n.space) {col1[k] <- cols[brks[k]]}
        plot(shp$geometry, col=col1,axes=F, border="gray")
        plot(shp$geometry[shp.order], border="black", add=TRUE)
        graphics::title(paste("Cluster map with incidence for", unique(id[[i]]$week), sep = " "))
        graphics::legend('bottomright', legend = c("0",levels(leg)[-1]), fill = cols)
        # grDevices::dev.off()

      }

      print(temp)

    }

  }

  # scan_seq
  set.seed(1234)
  simul <- matrix(sample(c(rpois(338 * 2, lambda = 10), rpois(338, lambda = 50))), nrow = 338, ncol = 3)
  colnames(simul) <- c('202301', '202302', '202303')
  data(shp_name); data(centroid)

  id <- scan_seq(data = simul, shp.name = shp_name, centroid = centroid)

  # scan_plot
  shp_US <- sf::st_read(system.file("extdata", "shp_US.shp", package = 'scanplot'))
  scan_plot(data = simul, shp = shp_US, id = id, allmap = F)

})

