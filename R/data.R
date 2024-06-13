#' Simulated longitudinal test data
#'
#'
#' @format
#' A data frame with 176 rows and 7 columns:
#' \describe{
#'   \item{id}{Subject id, factor}
#'   \item{arm}{Trial arm, factor}
#'   \item{sex}{Sex, factor}
#'   \item{time}{Measurement time, numeric}
#'   \item{weight}{Weight, numeric}
#'   \item{y}{Response variable, numeric}
#'   \item{f_true}{True signal}
#' }
#' @family built-in datasets
"testdata"

#' Simulated joint longitudinal and time-to-event test data
#'
#'
#' @format
#' A list of two data frames \code{lon} (longitudinal data) and
#' \code{tte} (time-to-event data). Generated using \code{simjm}.
#' @family built-in datasets
"example_data_jm"
