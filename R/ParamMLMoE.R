#' A Reference Class which contains parameters of a NMoE model.
#'
#' ParamMLMoE contains all the parameters of a NMoE model.
#'
#' @field X Numeric vector of length \emph{n} representing the covariates/inputs
#'   \eqn{x_{1},\dots,x_{n}}.
#' @field Y Numeric vector of length \emph{n} representing the observed
#'   response/output \eqn{y_{1},\dots,y_{n}}.
#' @field n Numeric. Length of the response/output vector `Y`.
#' @field K The number of experts.
#' @field p The order of the polynomial regression for the experts.
#' @field q The order of the logistic regression for the gating network.
#' @field wk Parameters of the gating network. \eqn{\boldsymbol{\wk} =
#'   (\boldsymbol{\wk}_{1},\dots,\boldsymbol{\wk}_{K-1})}{\wk =
#'   (\wk_{1},\dots,\wk_{K-1})} is a matrix of dimension \eqn{(q + 1, K -
#'   1)}, with `q` the order of the logistic regression for the gating network.
#'   `q` is fixed to 1 by default.
#' @field eta Polynomial regressions coefficients for each expert.
#'   \eqn{\boldsymbol{\eta} =
#'   (\boldsymbol{\eta}_{1},\dots,\boldsymbol{\eta}_{K})}{\eta =
#'   (\beta_{1},\dots,\beta_{K})} is a matrix of dimension \eqn{(p + 1, K)},
#'   with `p` the order of the polynomial regression. `p` is fixed to 1 by
#'   default.
#' @field df The degree of freedom of the NMoE model representing the complexity
#'   of the model.
#' @export
ParamMLMoE <- setRefClass(
  "ParamMLMoE",
  fields = list(

    X = "numeric",
    Y = "numeric",
    n = "numeric",
    phiEta = "list",
    phiWk = "list",

    K = "numeric", # Number of regimes
    R = "numeric", # Number of multinomial response models
    p = "numeric", # Dimension of eta (order of polynomial regression)
    q = "numeric", # Dimension of w (order of logistic regression)
    df = "numeric", # Degree of freedom

    etak = "array",
    wk = "matrix"
    ),

  methods = list(
    initialize = function(X = numeric(), Y = numeric(1), K = 2, R = 2, p = 1, q = 1) {

      X <<- X
      Y <<- Y
      n <<- length(Y)
      phiEta <<- designmatrixML(x = X, p = p)
      phiWk <<- designmatrixML(x = X, p = q)

      K <<- K
      R <<- R
      p <<- p
      q <<- q

      df <<- (q + 1) * (K - 1) + (p + 1) * K + K

      wk <<- matrix(0, q + 1, K - 1)
      # etak <<- matrix(NA, p + 1, K)
      etak <<- array(0, dim = c(p + 1, K, R-1))
      # sigma2 <<- matrix(NA, 1, K)
    },

    initParam = function(segmental = FALSE, nearTrue, wk_star, eta_star) {
      "Method to initialize parameters \\code{wk}, \\code{eta} and
      \\code{sigma2}.

      If \\code{segmental = TRUE} then \\code{wk}, \\code{eta} and
      \\code{sigma2} are initialized by clustering the response \\code{Y}
      uniformly into \\code{K} contiguous segments. Otherwise, \\code{wk},
      \\code{eta} and \\code{sigma2} are initialized by clustering randomly
      the response \\code{Y} into \\code{K} segments."

      ### We initialized the EM algorithm favourably by making
      ### a partition of starting values near the true components.
      if (nearTrue == TRUE){
        ## For testing the code.
        # K <- 3
        # q <- 1
        # p <- 1
        # wk <<- matrix(0, q + 1, K - 1)
        # eta <<- array(0, dim = c(p + 1, K, R-1))
        ###
        Ks <- dim(eta_star)[2]
        qt <- dim(wk_star)[1]
        inds <- seq(from = 1, to = Ks, by = 1)
        while (TRUE){
          s_inds <- sample(1:Ks, K, replace = TRUE)
          s_inds[K] <- Ks
          s_inds_unique <- unique(s_inds)
          if (length(s_inds_unique) == Ks){
            break
          }
        }
        # Original version
        for (i in 1:K){
          if (i < K){
            if (s_inds[i] < Ks){
              wk[, i] <<- wk_star[, s_inds[i]] +
                stats::rnorm(n = 1, mean = 0, sd = 0.007*n^(-0.083))
            } else {
              wk[, i] <<- matrix(0, qt, 1) +
                stats::rnorm(n = 1, mean = 0, sd = 0.007*n^(-0.083))
            }

          }
          etak[, i, ] <<- eta_star[, s_inds[i], ] +
            stats::rnorm(n = 1, mean = 0, sd = 0.007*n^(-0.083))
        }
      # We do not initialize the EM algorithm favourably by making
      ## a partition of starting values near the true components.
      }
      # } else {
      #   if (!segmental) {
      #
      #     klas <- sample(1:K, n, replace = TRUE)
      #?
      #     for (k in 1:K) {
      #
      #       Xk <- phiEta$XEta[klas == k,]
      #       yk <- Y[klas == k]
      #
      #       eta[, k] <<- solve(t(Xk) %*% Xk) %*% t(Xk) %*% yk
      #
      #       sigma2[k] <<- sum((yk - Xk %*% eta[, k]) ^ 2) / length(yk)
      #     }
      #   } else {# Segmental : segment uniformly the data and estimate the parameters
      #
      #     nk <- round(n / K) - 1
      #
      #     klas <- rep.int(0, n)
      #
      #     for (k in 1:K) {
      #       i <- (k - 1) * nk + 1
      #       j <- (k * nk)
      #       yk <- matrix(Y[i:j])
      #       Xk <- phiEta$XEta[i:j, ]
      #
      #       eta[, k] <<- solve(t(Xk) %*% Xk) %*% (t(Xk) %*% yk)
      #
      #       muk <- Xk %*% eta[, k, drop = FALSE]
      #
      #       sigma2[k] <<- t(yk - muk) %*% (yk - muk) / length(yk)
      #
      #       klas[i:j] <- k
      #     }
      #   }
      #
      #   # Intialize the softmax parameters
      #   Z <- matrix(0, nrow = n, ncol = K)
      #   Z[klas %*% ones(1, K) == ones(n, 1) %*% seq(K)] <- 1
      #   tau <- Z
      #   res <- IRLS(phiWk$XEta, tau, ones(nrow(tau), 1), wk)
      #   wk <<- res$W
      # }
    },

    MStep = function(paramMLMoE, verbose_IRLS, update_IRLS) {

      res_irls_Wk <- IRLS(phiWk$XEta, paramMLMoE$tik, ones(nrow(paramMLMoE$tik), 1),
                          wk, verbose_IRLS, update_IRLS)
      paramMLMoE$piikW <- res_irls_Wk$piik
      reg_irls <- res_irls_Wk$reg_irls

      wk <<- res_irls_Wk$W

      etak <<- etak + array(0, dim = c(p + 1, K, R-1))

      ## Note that eta <<- array(0, dim = c(p + 1, K, R-1))
      ## wk <<- matrix(0, q + 1, K - 1)
      ## matrix(paramMLMoE$phiEta$XEta[i, ], nrow = 1)
      ## IRLS(X, Tau, Gamma, Winit, verbose_IRLS, update_IRLS)
      ## ?Maximum of sum_{i=1}^n is not equivalent to sum of maximum.
      ## https://math.stackexchange.com/questions/3803010/when-the-maximum-of-sum-is-equal-to-the-sum-of-maxima
      #
      # for (k in (1:paramMLMoE$K)){
      #   for (i in (1:paramMLMoE$n)) {
      #     res_irls_Eta <<- IRLS(matrix(phiEta$XEta[i, ], nrow = 1),
      #                           statNMoE$tik,
      #                           ones(nrow(statNMoE$tik), 1), Eta[, k ,],
      #                           verbose_IRLS, update_IRLS)
      #     Eta[, k ,] <<- res_irls_Eta$W
      #   }
      # }

            # for (k in (1:paramMLMoE$K)){
      #   for (i in (1:paramMLMoE$n)) {
      #     res_irls_Eta <<- IRLS(matrix(phiEta$XEta[i, ], nrow = 1),
      #                           statNMoE$tik,
      #                           ones(nrow(statNMoE$tik), 1), Eta[, k ,],
      #                           verbose_IRLS, update_IRLS)
      #     Eta[, k ,] <<- res_irls_Eta$W
      #   }
      # }


      # ## Analytically solving a weighted least-squares problem.
      # for (k in 1:K) {
      #
      #   # Update the regression coefficients
      #   Xbeta <- phiEta$XEta * sqrt(statNMoE$tik[, k] %*% ones(1, p + 1))
      #   yk <- Y * sqrt(statNMoE$tik[, k])
      #
      #   eta[, k] <<- solve((t(Xbeta) %*% Xbeta)) %*% (t(Xbeta) %*% yk)
      #
      #   # Update the variances sigma2k
      #   sigma2[k] <<- sum(statNMoE$tik[, k] * ((Y - phiEta$XEta %*% eta[, k]) ^ 2)) / sum(statNMoE$tik[, k])
      #
      # }

      return(reg_irls)
    }
  )
)
