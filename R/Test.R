library(sf) # st_read

shp.kor <- st_read("/Users/leejaeseon/Desktop/Rpackage/shp/bnd_sigungu_00_2018_2018_2Q.shp")
shp.capit <- shp.kor[c(1:25, 50:59, 76:117),]

shp.name <- read.csv("/Users/leejaeseon/Desktop/Rpackage/shp/sigungu name_newest.csv",
                     fileEncoding = 'euc-kr')[,-1]

domestic.ST.wk <- read.csv("/Users/leejaeseon/Desktop/Rpackage/covid_korea/data/domestic.ST_weekly_FINAL.csv",
                           fileEncoding = 'euc-kr')[,-c(1:3)]
domestic.ST.wk <- domestic.ST.wk[,101:110]
capit <- domestic.ST.wk[c(1:25, 50:59, 76:117),]

pop <- read.csv("/Users/leejaeseon/Desktop/Rpackage/covid_korea/centroid_pop_2018.csv")
cent <- pop[,2:3]

wd <- "/Users/leejaeseon/Desktop/Rpackage/covid_korea/results_tmp"

scan.plot(data = domestic.ST.wk, shp = shp.kor, shp.name = shp.name, centroid = cent, wd = wd)
