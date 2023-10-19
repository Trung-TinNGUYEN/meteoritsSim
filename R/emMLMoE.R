#' emMLMoE implements the EM algorithm to fit a Normal Mixture of Experts (NMoE).
#'
#' emMLMoE implements the maximum-likelihood parameter estimation of a Normal
#' Mixture of Experts (NMoE) model by the Expectation-Maximization (EM)
#' algorithm.
#'
#' @details emMLMoE function implements the EM algorithm for the NMoE model. This
#'   function starts with an initialization of the parameters done by the method
#'   `initParam` of the class [ParamNMoE][ParamNMoE], then it alternates between
#'   the E-Step (method of the class [StatNMoE][StatNMoE]) and the M-Step
#'   (method of the class [ParamNMoE][ParamNMoE]) until convergence (until the
#'   relative variation of log-likelihood between two steps of the EM algorithm
#'   is less than the `threshold` parameter).
#'
#' @param X Numeric vector of length \emph{n} representing the covariates/inputs
#'   \eqn{x_{1},\dots,x_{n}}.
#' @param Y Numeric vector of length \emph{n} representing the observed
#'   response/output \eqn{y_{1},\dots,y_{n}}.
#' @param K The number of experts.
#' @param p Optional. The order of the polynomial regression for the experts.
#' @param q Optional. The order of the logistic regression for the gating
#'   network.
#' @param n_tries Optional. Number of runs of the EM algorithm. The solution
#'   providing the highest log-likelihood will be returned.
#' @param max_iter Optional. The maximum number of iterations for the EM
#'   algorithm.
#' @param threshold Optional. A numeric value specifying the threshold for the
#'   relative difference of log-likelihood between two steps of the EM as
#'   stopping criteria.
#' @param verbose Optional. A logical value indicating whether or not values of
#'   the log-likelihood should be printed during EM iterations.
#' @param verbose_IRLS Optional. A logical value indicating whether or not
#'   values of the criterion optimized by IRLS should be printed at each step of
#'   the EM algorithm.
#' @return EM returns an object of class [ModelNMoE][ModelNMoE].
#' @seealso [ModelNMoE], [ParamNMoE], [StatNMoE]
#' @export
#'
#' @examples
#' data(tempanomalies)
#' x <- tempanomalies$Year
#' y <- tempanomalies$AnnualAnomaly
#'
#' nmoe <- emMLMoE(X = x, Y = y, K = 2, p = 1, verbose = TRUE)
#'
#' nmoe$summary()
#'
#' nmoe$plot()
#' To do list:
#' 1. Change alphak_0, betak_0, sigmak_0.
#' 2. param$MStep.
emMLMoE <- function(X, Y, K, R,  p = 1, q = 1, n_tries = 1, max_iter = 1500,
                   threshold = 1e-6, verbose = FALSE, verbose_IRLS = FALSE, update_IRLS = TRUE,
                   nearTrue = FALSE, wk_star = NA, etak_star = NA) {

  top <- 0
  try_EM <- 0
  best_loglik <- -Inf

  while (try_EM < n_tries) {
    try_EM <- try_EM + 1

    if (n_tries > 1 && verbose) {
      message("EM try number: ", try_EM, "\n")
    }

    # Initializations
    param <- ParamMLMoE(X = X, Y = Y, K = K, R = R, p = p, q = q)
    param$initParam(segmental = FALSE, nearTrue, wk_star, etak_star)

    # print(param$alpha)

    iter <- 0
    converge <- FALSE
    prev_loglik <- -Inf

    stat <- StatMLMoE(paramMLMoE = param)

    # ### Test.
    # reg_irls <- param$MStep(stat, verbose_IRLS)
    #
    # stat$computeLikelihood(reg_irls)

    # print("stat$loglik_inital = ")
    # print(stat$loglik)

    # print(param$sigma2)

    while (!converge && (iter <= max_iter)) {

      stat$EStep(param)

      # print(param$alpha)

      reg_irls <- param$MStep(stat, verbose_IRLS, update_IRLS)

      stat$computeLikelihood(reg_irls)

      # print(param$alpha)

      iter <- iter + 1
      if (verbose) {
        message("EM MLMoE: Iteration: ", iter, " | log-likelihood: "  , stat$loglik)
      }

      if (prev_loglik - stat$loglik > 1e-5) {
        if (verbose) {
          warning("EM log-likelihood is decreasing from ", prev_loglik, "to ", stat$loglik, "!")
        }
        top <- top + 1
        if (top > 20)
          break
      }

      # Test of convergence
      converge <- abs((stat$loglik - prev_loglik) / prev_loglik) <= threshold

      if (is.na(converge)) {
        converge <- FALSE
      } # Basically for the first iteration when prev_loglik is Inf

      prev_loglik <- stat$loglik
      stat$stored_loglik <- c(stat$stored_loglik, stat$loglik)
    }# End of an EM loop

    # if (nearTrue){
    #   statSolution <- stat$copy()
    #   paramSolution <- param$copy()
    #   return(ModelNMoE(param = paramSolution, stat = statSolution))
    # }

    if (stat$loglik > best_loglik) {
      # statSolution <- stat$copy()
      paramSolution <- param$copy()

      best_loglik <- stat$loglik
    }

    if (n_tries > 1 && verbose) {
      message("Max value of the log-likelihood: ", stat$loglik, "\n\n")
    }
  }

  # # Computation of c_ig the hard partition of the curves and klas
  # statSolution$MAP()
  #
  # if (n_tries > 1 && verbose) {
  #   message("Max value of the log-likelihood: ", statSolution$loglik, "\n")
  # }
  #
  # statSolution$computeStats(paramSolution)
  #
  # return(ModelNMoE(param = paramSolution, stat = statSolution))

  return(ModelMLMoE(param = paramSolution))
}
