#' Draw a sample from a normal mixture of linear experts model.
#'
#' @param alphak The parameters of the gating network. `alphak` is a matrix of
#'   size \emph{(q + 1, K - 1)}, with \emph{K - 1}, the number of regressors
#'   (experts) and \emph{q} the order of the logistic regression
#' @param betak Matrix of size \emph{(p + 1, K)} representing the regression
#'   coefficients of the experts network.
#' @param sigmak Vector of length \emph{K} giving the standard deviations of
#'   the experts network.
#' @param x A vector og length \emph{n} representing the inputs (predictors).
#'
#' @return A list with the output variable `y` and statistics.
#' \itemize{
#'   \item `y` Vector of length \emph{n} giving the output variable.
#'   \item `zi` A vector of size \emph{n} giving the hidden label of the
#'     expert component generating the i-th observation. Its elements are
#'     \eqn{zi[i] = k}, if the i-th observation has been generated by the
#'     k-th expert.
#'   \item `z` A matrix of size \emph{(n, K)} giving the values of the binary
#'     latent component indicators \eqn{Z_{ik}}{Zik} such that
#'     \eqn{Z_{ik} = 1}{Zik = 1} iff \eqn{Z_{i} = k}{Zi = k}.
#'  \item `stats` A list whose elements are:
#'    \itemize{
#'      \item `Ey_k` Matrix of size \emph{(n, K)} giving the conditional
#'        expectation of Yi the output variable given the value of the
#'        hidden label of the expert component generating the ith observation
#'        \emph{zi = k}, and the value of predictor \emph{X = xi}.
#'      \item `Ey` Vector of length \emph{n} giving the conditional expectation
#'        of Yi given the value of predictor \emph{X = xi}.
#'      \item `Vary_k` Vector of length \emph{k} representing the conditional
#'        variance of Yi given \emph{zi = k}, and \emph{X = xi}.
#'      \item `Vary` Vector of length \emph{n} giving the conditional expectation
#'        of Yi given \emph{X = xi}.
#'    }
#' }
#' @export
#'
#' @examples
#' n <- 500 # Size of the sample
#' alphak <- matrix(c(0, 8), ncol = 1) # Parameters of the gating network
#' betak <- matrix(c(0, -2.5, 0, 2.5), ncol = 2) # Regression coefficients of the experts
#' sigmak <- c(1, 1) # Standard deviations of the experts
#' x <- seq.int(from = -1, to = 1, length.out = n) # Inputs (predictors)
#'
#' # Generate sample of size n
#' sample <- sampleUnivNMoE(alphak = alphak, betak = betak, sigmak = sigmak, x = x)
#'
#' # Plot points and estimated means
#' plot(x, sample$y, pch = 4)
#' lines(x, sample$stats$Ey_k[, 1], col = "blue", lty = "dotted", lwd = 1.5)
#' lines(x, sample$stats$Ey_k[, 2], col = "blue", lty = "dotted", lwd = 1.5)
#' lines(x, sample$stats$Ey, col = "red", lwd = 1.5)
sampleUnivNMoE <- function(alphak, betak, sigmak, x) {

  n <- length(x)

  p <- nrow(betak) - 1
  q <- nrow(alphak) - 1
  K = ncol(betak)

  # Build the regression design matrices
  XBeta <- designmatrix(x, p, q)$XBeta # For the polynomial regression
  XAlpha <- designmatrix(x, p, q)$Xw # For the logistic regression

  y <- rep.int(x = 0, times = n)
  z <- zeros(n, K)
  zi <- rep.int(x = 0, times = n)

  # Calculate the mixing proportions piik:
  piik <- multinomialLogit(alphak, XAlpha, zeros(n, K), ones(n, 1))$piik

  for (i in 1:n) {
    zik <- stats::rmultinom(n = 1, size = 1, piik[i,])

    mu <- as.numeric(XBeta[i,] %*% betak[, zik == 1])
    sigma <- sigmak[zik == 1]

    y[i] <- stats::rnorm(n = 1, mean = mu, sd = sigma)
    z[i, ] <- t(zik)
    zi[i] <- which.max(zik)

  }

  # Statistics (means, variances)

  # E[yi|xi,zi=k]
  Ey_k <- XBeta %*% betak

  # E[yi|xi]
  Ey <- rowSums(piik * Ey_k)

  # Var[yi|xi,zi=k]
  Vary_k <- sigmak ^ 2

  # Var[yi|xi]
  Vary <- rowSums(piik * (Ey_k ^ 2 + ones(n, 1) %*% Vary_k)) - Ey ^ 2

  stats <- list()
  stats$Ey_k <- Ey_k
  stats$Ey <- Ey
  stats$Vary_k <- Vary_k
  stats$Vary <- Vary
  # TestTin
  stats$piik <- piik
  stats$gate <- colSums(piik)
  return(list(y = y, zi = zi, z = z, stats = stats))
}
