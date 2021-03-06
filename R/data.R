#' Chicago Crime data from 2001 to 2019
#'
#' A randomly shuffled subset of the dataset of reported crimes in the city of Chicago in the period 2001 to 2019.
#' @usage data(crime)
#' @format A data frame with 48,130 rows and 25 variables:
#' \describe{
#'   \item{Date}{date of which the crime was reported}
#'   \item{Primary Type}{type of crime that was reported, one of a finite amount of options}
#'   \item{Description}{description of the reported crime}
#'   \item{Location Description}{type of place where the crime was committed}
#'   \item{Arrest}{logical; if \code{TRUE}, then the crime resulted in an arrest}
#'   \item{Domestic}{logical; if \code{TRUE}, then the incident was domestic-related}
#'   \item{Beat}{the beat where the crime occurred (smallest police geographical area)}
#'   \item{District}{police district where the crime occurred}
#'   \item{Ward}{the ward (City Council district) in Chicago where the crime occurred}
#'   \item{Community Area}{the community area where the incident occurred}
#'   \item{Year}{year when the crime happened}
#'   \item{Latitude}{latitude co-ordinate (in degrees) of the block where the crime occurred}
#'   \item{Longitude}{longitude co-ordinate (in degrees) of the block where the crime occurred}
#'   \item{Hour}{estimated hour of the day when the crime happened}
#'   \item{dist_from_station}{distance in km from the nearest police station (from where the crime happened)}
#'   \item{COMMUNITY AREA NAME}{name of the community area}
#'   \item{PERCENT OF HOUSING CROWDED}{Percentage of occupied housing units, in the community area, with more than one person per room}
#'   \item{PERCENT HOUSEHOLDS BELOW POVERTY}{Percentage of households, in the community area, that are living below the federal poverty level}
#'   \item{PERCENT AGED 16+ UNEMPLOYED}{Percentage of persons, in the community area, that are over the age of 16 years that are unemployed}
#'   \item{PERCENT AGED 25+ WITHOUT HIGH SCHOOL DIPLOMA}{Percentage of persons in the community area that are over the age of 25 years without a high school education}
#'   \item{PERCENT AGED UNDER 18 OR OVER 64}{Percent of the population of the community area under 18 or over 64 years of age}
#'   \item{PER CAPITA INCOME}{Community Area Per capita income; estimated as the sum of tract-level aggregate incomes divided by the total population}
#'   \item{HARDSHIP INDEX}{Score that incorporates each of the six selected socioeconomic indicators}
#'   \item{Crime Density}{crime density as fit by a kernel density estimator}
#'   \item{Population Density}{population density as fit by a kernel density estimator}
#' }
#' @source \url{https://data.cityofchicago.org/Public-Safety/Crimes-2001-to-present/ijzp-q8t2}
"crime"


#' Observed Precipitation Maxima
#'
#' Observed yearly precipitation maxima from 1988 to 2017 in South West UK
#' @usage data(rain_max_obs)
#' @format A data frame with 4731 rows and 8 variables:
#' \describe{
#'   \item{date}{date of observation}
#'   \item{lon}{longitude of observation}
#'   \item{lat}{latitude of observation}
#'   \item{elev}{elevation of observation}
#'   \item{rain}{precipitation value (mm)}
#'   \item{year}{year of observation}
#'   \item{day}{day of month of observation}
#'   \item{month}{month of observation}
#' }
"rain_max_obs"

#' Model Precipitation Maxima
#'
#' Gridded model predicted yearly precipitation maxima from the ERA5 reanalysis from the period 1979 to 2018 in South West UK
#' @usage data(rain_max_obs)
#' @format A data frame with 4200 rows and 7 variables:
#' \describe{
#'   \item{lon}{longitude of centre of grid cell}
#'   \item{year}{year of estimate}
#'   \item{lat}{latitude of centre of grid cell}
#'   \item{rain}{estimated precipitation value (mm)}
#'   \item{id}{id number of the model, corresponding to grid cell}
#'   \item{elev}{average elevation across grid cell}
#'   \item{dist2coast}{distance in km to the coastline from the grid point centre}
#' }
"rain_max_model"


