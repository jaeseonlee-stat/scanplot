library(sf) # st_read

# shp.kor <- readOGR("D:/covid_korea/shp/bnd_sigungu_00_2018_2018_2Q.shp")
# shp.kor <- read_sf("/Users/leejaeseon/Desktop/Rpackage/shp/bnd_sigungu_00_2018_2018_2Q.shp")
shp.kor <- st_read("/Users/leejaeseon/Desktop/Rpackage/shp/bnd_sigungu_00_2018_2018_2Q.shp")


# SCAN STATISTICS
# INPUT :
# data - domestic.ST.wk (example)
# shp - shape file
# shp.name - (shp file)
# centroid - N X 2 matrix with (latitude, longitude)
# pop.upper.bound -
# n.simulations -
# alpha.level -
# wd - working directory
# OUTPUT :
#

scan.plot <- function(data, shp, shp.name, centroid, pop.upper.bound = .2, n.simulations = 999, alpha.level = .05, wd) {

  n.space <- nrow(shp)
  
  library(SpatialEpi)
  library(scanstatistics)
  library(RColorBrewer)
  
  cent <- centroid
  geo <- SpatialEpi::latlong2grid(cent)
  knn_mat <- scanstatistics::coords_to_knn(geo, 15)
  zones <- scanstatistics::knn_zones(knn_mat)
  
  setwd(wd)  
  temp <-0; id <- list(); scan <- 0
  for (i in 149:(ncol(data)-3)) {

    col1 <- c()
    temp <- temp +1
    cluster2 <- c()
    counts <- as.numeric(data[,(2+i)])
    
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
      id[[temp]] <- data.frame(week=rep(colnames(data)[i+2],length(cluster2)),
                               sido=shp.name$sido2[cluster2], sigungu=shp.name$sigungu[cluster2],
                               No_of_cases=No.cases, RR=RR, Exp.cases=Exp.cases, pvalue=p.val,
                               total=total, cluster.cases=cluster.cases, percent=percent)
      id[[temp]] <- id[[temp]][order(id[[temp]][,"sido"]),]
    }
    
    jpeg(paste("cluster_",colnames(data)[i+2],".jpeg",sep = ""), res=600, width=200, height=150, pointsize=9,units='mm')
    brks <- cut(counts, breaks=c(0,210,490,1010,1600,max(counts)+1), include.lowest = TRUE, right=FALSE)
    cols <- brewer.pal(5, "Reds")
    for (k in 1:250) {col1[k] <- cols[brks[k]]}
    plot(shp$geometry, col=col1,axes=F, border="gray")
    plot(shp$geometry[cluster2], border="black", add=TRUE)
    title(paste("Cluster map for", colnames(data)[i+2], sep=" "))
    legend('bottomright',legend=c("[0, 210)", "[210, 490)", "[490, 1010)", "[1010, 1600)", "≥ 1600"),fill=cols,
           cex=0.7)
    dev.off()
    
    scan <- 0
    
    print(temp)
  }
  capture.output(id, file="cluster info_covid_korea_by20221228.txt")
  
  # for excel file
  # tmp_df <- data.frame()
  # for (i in 1:length(id)) {
  #   tmp <- id[[i]]
  #   if (!is.null(tmp)) {
  #     tmp_df <- rbind(tmp_df, cbind(id = i, tmp)) 
  #   }
  # }
  # 
  # tmp_df[order(tmp_df$week, -tmp_df$No_of_cases), ]
  # tmp_table <- as.data.frame(table(tmp_df$id))
  # tmp_df <- cbind(frequency = tmp_table$Freq[match(tmp_df$id, tmp_table$Var1)], tmp_df)
  # 
  # library(xlsx)
  # write.xlsx(tmp_df, '../result.xlsx', sheetName = 'entire', row.names = F, append = T)
}

domestic.ST.wk <- read.csv("/Users/leejaeseon/Desktop/Rpackage/covid_korea/data/domestic.ST_weekly_FINAL.csv",
                           fileEncoding = 'euc-kr')[,-1]
wd <- "/Users/leejaeseon/Desktop/Rpackage/covid_korea/results_tmp"
shp.name <- read.csv("/Users/leejaeseon/Desktop/Rpackage/shp/sigungu name_newest.csv",
                     fileEncoding = 'euc-kr')
shp.capit <- shp.kor[c(1:25, 50:59, 76:117),]
capit <- domestic.ST.wk[c(1:25, 50:59, 76:117),]
pop <- read.csv("/Users/leejaeseon/Desktop/Rpackage/covid_korea/centroid_pop_2018.csv")
cent <- pop[,2:3]

scan.plot(data = domestic.ST.wk, shp = shp.kor, shp.name = shp.name, centroid = cent, wd = wd)
