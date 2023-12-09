#' scan_seq
#'
#' Apply Kulldorff's (1997) scan statistics to detect spatial and temporal patterns in the data at each time point, facilitating cluster detection.
#' We assumed a Zero-inflated Poisson distribution if the zero cases in data exceed 30% of the total incidences, Poisson distribution, otherwise.
#'
#' @import SpatialEpi
#' @import scanstatistics
#' @import utils
#'
#' @param data n (space) X m (time) matrix of spatio-temporal data.
#' @param shp.name n X 2 matrix with (state name, county name). The order must be identical to shape file!
#' @param centroid n X 2 matrix with (latitude, longitude).
#' @param pop.upper.bound The upper bound on the proportion of the total population for each zone. Default is 0.2.
#' @param n.simulations The number of simulate scan statistics to generate to calculate a p-value. Default is 999.
#' @param alpha.level A significance level for kulldorff spatial scan statistics. Default is 0.05.
#' @param save Logical: whether to save the result of scan statistics as csv and txt file. Default is FALSE.
#'
#' @export
#'
#' @seealso [scan_plot()], [scanstatistics::scan_eb_zip()], [SpatialEpi::kulldorff()]
#'
#' @examples
#' # Load the data
#' data(covid_US)
#' data(shp_name)
#' data(centroid)
#'
#' # Implement spatial scan statistics
#' # Assign appropriate working directory before "save = TRUE"
#' # setwd('your working directory')
#' id <- scan_seq(data = covid_US, shp.name = shp_name, centroid = centroid,
#'                pop.upper.bound = .2, n.simulations = 999, alpha.level = .05, save = FALSE)
#' id

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
