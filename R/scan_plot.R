#' scan_plot
#'
#' Plot spatial hot-spot detected by scan statistics (Kulldorff, 1997).
#'
#' @import RColorBrewer
#' @import grDevices
#' @import graphics
#' @import stringr
#' @import sf
#' @import stats
#'
#' @param data N (space) X M (time) matrix
#' @param shp shape file with identical space order with data
#' @param id list of cluster. output of scan.seq function
#' @param allmap logical: if FALSE, plot the map with detected cluster time point only
#'
#' @export
#'
#' @seealso [scan_seq()]
#'
#' @examples
#' # Load the required package
#' library(sf) # st_read
#'
#' # Load the data
#' data(covid_NY)
#' data(shp_name)
#' data(centroid)
#' shp_NY <- sf::st_read(system.file("extdata", "shp_NY.shp", package = 'scanplot'))
#'
#' # Implement spatial scan statistics
#' # Assign appropriate working directory before "save = TRUE"
#' id <- scan_seq(data = covid_NY, shp.name = shp_name, centroid = centroid,
#'                pop.upper.bound = .2, n.simulations = 999, alpha.level = .05, save = FALSE)
#'
#' # Assign appropriate working directory before scan_plot()
#' scan_plot(data = covid_NY, shp = shp_NY, id = id, allmap = FALSE)

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

        grDevices::jpeg(paste("cluster_", colnames(data)[i], ".jpeg", sep = ""),
                        res = 600, width = 200, height = 150, pointsize = 9,units = 'mm')

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
        grDevices::dev.off()

      } else if (tryCatch(!is.null(id[[i]]), error = function(e) FALSE)) {
        # If there are detected clusters

        shp.order <- id[[i]]$shp.order

        grDevices::jpeg(paste("cluster_", unique(id[[i]]$week), ".jpeg", sep = ""),
                        res = 600, width = 200, height = 150, pointsize = 9, units = 'mm')

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
        grDevices::dev.off()

      }

    print(temp)

  }

}
