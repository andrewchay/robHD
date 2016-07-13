### R code from vignette source 'robHD.Rnw'

###################################################
### code chunk number 1: robHD.Rnw:37-44
###################################################
library(robHD)
x = matrix(rnorm(20000), nrow = 100)
b = c(0.5, 1, -1, 2, -0.5, rep(0, 195))
y = crossprod(t(x), b) + rnorm(100, 0, 1)
delta = rbinom(100, 1, .7)
res = robHD(x, y, delta, theta = c(1, 10, 100))
str(res)


###################################################
### code chunk number 2: robHD.Rnw:50-78
###################################################
n.beta0 = 5
P = 200
N = 200
z = matrix(rnorm(N * (P + n.beta0)), nrow = N)
x = z[, 1:n.beta0]
z = z[, (n.beta0 + 1):(n.beta0 + P)]
xz = matrix(apply(z, 2, "*", x), nrow = N)
beta0 = rep(1, n.beta0)
theta0 = c(rep(1, 5), rep(0, P - 5))
tau0 = rep(0, n.beta0 * P)
tau0[c(sort(sample(1:(n.beta0 * 5))[1:20]))] =
  round(runif(20, 0.9, 1.1), 1)
mu = (as.vector(crossprod(t(x), beta0) + crossprod(t(z), theta0) +
                  crossprod(t(xz), tau0)))
e1 = rnorm(N, 0, 1)
e2 = rcauchy(N)
e2 = ifelse(as.logical(rbinom(N, 1, 0.2)), e2, 0)
t = e1 + mu
c = log(rweibull(N, shape = 0.9, scale = exp(mu) * 6))
delta = (t < c)
y = (ifelse(t < c, t, c)) + e2
res = robHD(cbind(x, z, xz), y, delta, theta = c(1, 10, 50))
for (j in 1:P)
{
  res = robHD(cbind(x, z[, j], xz[, (n.beta0 * (j - 1) + 1):(n.beta0 * j)]),
               y, delta, theta = c(1, 10, 50))
}
str(res)


###################################################
### code chunk number 3: robHD.Rnw:84-106 (eval = FALSE)
###################################################
## n.beta0 = 5
## P = 200
## N = 200
## Z = matrix(rnorm(N * (P + n.beta0)), nrow = N)
## X = Z[, 1:n.beta0]
## Z = Z[, (n.beta0 + 1):(n.beta0 + P)]
## XZ = matrix(apply(Z, 2, "*", X), nrow = N)
## beta0 = rep(1, n.beta0)
## theta0 = c(rep(1, 5), rep(0, P - 5))
## tau0 = rep(0, n.beta0 * P)
## tau0[c(sort(sample(1:(n.beta0 * 5))[1:20]))] =
##   round(runif(20, 0.9, 1.1), 1)
## mu = (as.vector(crossprod(t(X), beta0) + crossprod(t(Z), theta0) +
##                   crossprod(t(XZ), tau0)))
## e1 = rnorm(N, 0, 1)
## e2 = rcauchy(N)
## e2 = ifelse(as.logical(rbinom(N, 1, 0.2)), e2, 0)
## t = e1 + mu
## c = log(rweibull(N, shape = 0.9, scale = exp(mu) * 6))
## delta = (t < c)
## y = (ifelse(t < c, t, c)) + e2
## res = robHD(cbind(X, Z, XZ), y, delta, theta = c(1, 10, 50))


