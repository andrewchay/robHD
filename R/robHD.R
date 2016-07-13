#' Robust joint estimator.
#' @rdname robHD
#' @description This is a function to compute the regression coefficients for
#' robust penalized regression. The exponential squared loss function is used to
#' provide the robustness, and penalization is employed to enforce the sparsity
#' of estimators.
#' @details The sequence of models indexed by the regularization parameter
#' \code{lambda} is fit using a coordinate descent algorithm. The loss function
#' is \deqn{\omega exp(-(y - x\beta) ^ 2/(2\theta)) - \lambda|\beta|_1,} where
#' the first term is the exponential squared loss, and the second term is the
#' Lasso penalty. \eqn{\omega} is the Kaplan-Meier weight. To find the maximizer
#' of the loss function, an two-layer MM coordinate descent (MMCD) algorithm is
#' used. The outer layer of the MMCD algorithm is an MM algorithm, in which the
#' exponential term is approximated by a weighted least square. The main
#' contribution of this layer is that it uses the minorization technique to
#' solve the maximization problem using relatively simple approximations in an
#' iterative manner. The inner layer of the MMCD algorithm contains a coordinate
#' descent algorithm. The goal is to find a good enough solution to the minorized
#' penalized regression problem. Because the original loss function is
#' approximated iteratively, it might require a huge amount of time to compute
#' the solution if each inner layer achieves convergence. As a result, the
#' default maximum number of inner iterations is three.
#'
#' The sequence of tuning parameters \eqn{\lambda} is
#' generated based on the largest value of \code{lambda0} and the number of
#' \eqn{\lambda} \code{n_lambda}. The minimum of \eqn{\lambda} equals
#' \code{lambda0} * 0.001, if \eqn{p<n}, and equals \code{lambda0}* 0.05, if
#' \eqn{p>=n}.
#'
#' @param x explanotory variable matrix of size \eqn{n * p}.
#' @param y response variable vector of length \eqn{n}. \code{y} is the logarithm
#' of the observed event/censoring time.
#' @param delta censoring indicator vector of length \eqn{n}. \code{delta} equals
#' 1 meaning event and 0 meaning failure.
#' @param n_lambda number of tuning parameters \eqn{\lambda} to use in the
#' penalized regression. Default is 20.
#' @param lambda0 The user-specified largest value for tuning parameter
#' \eqn{\lambda}. The default will be the \eqn{\lambda} such that all regression
#' coefficients equal zero.
#' @param kappa Tuning parameter in the MCP penalty. \code{kappa} ranges from 0
#' to 1, where 0 corresponds to Lasso penalty, and 1 corresponds to the hard
#' threshold penalty.
#' @param theta Vector of robust tuning parameters.
#' \code{theta} takes positive values.
#' Small \code{theta} downweighs the influence of outliers. Hence, the smaller
#' \code{theta} is, the more robust the estimators are.
#' @param eps Convergence threshold. The algorithm iterates until the
#' change in any coefficient is less than eps. Default is 1e-5.
#' @param max.iter The maximum number of outer layer iterations. Default is 100.
#' @param max.cd The maximum number of inner layer iterations for the coordinate
#' descent algorithm. Default is 5.
#' @param penalty Penalty function. Values could be "LASSO" or "MCP". Default is
#' "MCP".
#' @param init User specified initial values for regression coefficients. If not
#' specified, default value is a length \eqn{p} vector of zeros.
#' @return An object with S3 class "robint" containing:
#' \item{betalasso}{The matrix of regression coefficients using Lasso penalty.
#' The size of the matrix is \eqn{p} * (\code{n_lambda} * n_theta), where
#' n_theta is the number of thetas. The first \code{n_lambda} columns correspond
#' to estimators from \code{theta}[1]. The next \code{n_lambda} columns correspond
#' to estimators from \code{theta}[2], so on and so forth.}
#' \item{betamcp}{The matrix of regression coefficients using MCP penalty.
#' The size of the matrix is \eqn{p} * (\code{n_lambda} * n_theta), where
#' n_theta is the number of thetas. The first \code{n_lambda} columns correspond
#' to estimators from \code{theta}[1]. The next \code{n_lambda} columns correspond
#' to estimators from \code{theta}[2], so on and so forth.}
#' \item{iter}{Vector of length \code{n_lambda} * n_theta for number of
#' iterations in the MM algorithm for the corresponding tuning parameters. The
#' order of tuning parameters is the same as in \code{betalasso}.}
#' \item{lambda}{The vector of \eqn{\lambda}.}
#' \item{kappa}{Tuning parameter in the MCP penalty.}
#' \item{theta}{The vector of robust tuning paramether.}
#' \item{omega}{The vector of Kaplan-Meier weights for each sample.}
#' @examples
#' x = matrix(rnorm(20000), nrow = 100)
#' b = c(0.5, 1, -1, 2, -0.5, rep(0, 195))
#' y = crossprod(t(x), b) + rnorm(100, 0, 1)
#' delta = rbinom(100, 1, .7)
#' robHD(x, y, delta, theta = c(10, 20))
#' @export robHD
#' @useDynLib robHD pathSearch

robHD = function(x, y, delta, n_lambda = 20, lambda0 = -1, kappa = 1 / 3, theta,
                    eps = 1e-5, max.iter = 100, max.cd = 5, penalty = "MCP",
                    init = NaN)
{
  n_theta = length(theta)
  if (length(delta) != length(y)) stop("Length of y and delta do not match.")
  if (length(y) != nrow(as.matrix(x)))
    stop("Dimensionality of y and x do not match.")

  orders = sort(y, index.return = T)$ix
  delta = delta[orders]
  y = y[orders]
  x = as.matrix(x)[orders,]
  p = ncol(x)
  n = length(y)
  tmp = ((n - 1) / n) ^ delta[1]
  for (i in 2:(n - 1))
    tmp[i] = tmp[i - 1] * ((n - i) / (n - i + 1)) ^ delta[i]
  omega = delta[1] / n
  for (i in 2:n)
    omega[i] = delta[i] / (n - i + 1) * tmp[i - 1]
  rm(tmp)
  v = sqrt(omega)
  meanX = as.numeric(crossprod(omega, x) / sum(omega))
  meany = as.numeric(crossprod(y, omega) / sum(omega))
  Y = y - meany
  Y = v * Y
  X = sweep(x, 2, meanX)
  X = sweep(X, 1, v, "*")
  r = sqrt(colSums(X * X))
  X = t(t(X) / r)
  if (lambda0 == -1) lambda0 = max(abs(as.vector(Y %*% X)))
  if (p >= n)
    lambda = (seq(sqrt(lambda0), sqrt(lambda0 * 0.001), len = n_lambda)) ^ 2 else
      lambda = (seq(sqrt(lambda0), sqrt(lambda0 * 0.05), len = n_lambda)) ^ 2

  if (is.nan(init)) init = rep(0, p)
  if (length(kappa) != 1) stop("kappa has to be a single value.")
  res = .Call("pathSearch", x, y, omega, lambda, kappa, theta, n_lambda,
              n_theta, eps, as.integer(max.iter), as.integer(max.cd),
              penalty, init)
  r1 = res[[1]]
  name1 = rep("", p)
  r2 = res[[2]]
  name2 = rep("", p)
  for (i in 1:p)
  {
    name1[i] = paste("las b", i, sep = "")
    name2[i] = paste("mcp b", i, sep = "")
  }
  rownames(r1) = name1
  rownames(r2) = name2
  return(structure(list(betalasso = r1, betamcp = r2, iter = as.vector(res[[3]]),
                        lambda = res[[4]], kappa = res[[5]], theta = res[[6]],
                        omega = omega), class = "rob"))
}
