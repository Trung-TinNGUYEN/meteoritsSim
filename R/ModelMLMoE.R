#' A Reference Class which represents a fitted NMoE model.
#'
#' ModelNMoE represents an estimated NMoE model.
#'
#' @field param A [ParamNMoE][ParamNMoE] object. It contains the estimated
#'   values of the parameters.
#' @field stat A [StatNMoE][StatNMoE] object. It contains all the statistics
#'   associated to the NMoE model.
#' @seealso [ParamNMoE], [StatNMoE]
#' @export
#'
#' @examples
#' data(tempanomalies)
#' x <- tempanomalies$Year
#' y <- tempanomalies$AnnualAnomaly
#'
#' nmoe <- emNMoE(X = x, Y = y, K = 2, p = 1, verbose = TRUE)
#'
#' # nmoe is a ModelNMoE object. It contains some methods such as 'summary' and 'plot'
#' nmoe$summary()
#' nmoe$plot()
#'
#' # nmoe has also two fields, stat and param which are reference classes as well
#'
#' # Log-likelihood:
#' nmoe$stat$loglik
#'
#' # Parameters of the polynomial regressions:
#' nmoe$param$beta
ModelMLMoE <- setRefClass(
  "ModelMLMoE",
  fields = list(
    param = "ParamMLMoE"
  )
)
