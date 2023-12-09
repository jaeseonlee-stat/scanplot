#' covid_US
#'
#' The covid-19 weekly death counts from 9 US East coast states from March 17 to September 28, 2020.
#' The nine states are Connecticut, Delaware, Maryland, Massachusetts, New Jersey, New York, Pennsylvania, Rhode Island, and Virginia.
#' This dataset was obtained from The New York Times.
#'
#' @format
#' covid_US is a date frame with 338 counties (rows) and 28 weeks (columns).
#' The names of each column indicate the starting day of a given week.
#' For example, the first column '2020.03.17' represents 2020.03.17 - 2020.03.23.
#'
#' \describe{
#'  \item{week}{Each columns (e.g. 2020.03.17) indicate the starting day for given week (e.g. 2020.03.17 - 2020.03.23)}
#' }
#'
"covid_US"

#' centroid
#'
#' The dataset represents 338 counties (rows) with their coordinates, that is, latitude and longitude (columns).
#' Incorporated 338 counties are identical to covid_US data.
#'
#' @format
#' centroid is a data frame with 338 counties (rows) and  2 spatial coordinates (columns).
#'
#' \describe{
#'  \item{X coordinate}{The first column represents 'longitude'}
#'  \item{Y coordinate}{The second column represents 'latitude'}
#' }
#'
"centroid"

#' shp_name
#'
#' shp_name consists of 338 counties (rows) and their state and county names (columns).
#'
#' @format
#' shp_name is a data frame with 338 counties (rows) and 2 spatial names (columns).
#'
#' \describe{
#'  \item{state}{The first column represents the name of state for a given row}
#'  \item{county}{The second column represents the name of county for a given row}
#' }
#'
"shp_name"
