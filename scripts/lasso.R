lasso <- function(x, y, sigma = 1, alpha = 0.05, c = 1.1, rtol = 1e-6,
                  verb = 5)
{
  n <- nrow(x)
  p <- ncol(x)
  lambda <- c * sigma * 2 * sqrt(n) * qnorm(1 - alpha/(2*p))
  P <- list(sense = "min")
  P$c <- c(rep(lambda, 2*p), rep(0, n), 1, 0, 0)/n
  A <- as.matrix.csr(x)
  A <- cbind(A, -A, as(n, "matrix.diag.csr"), as.matrix.csr(0, n, 3))
  A <- rbind(A,cbind(as.matrix.csr(0, 2, 2*p + n),
                     as.matrix.csr(c(-.5,-.5,1,0,0,1), 2, 3)))
  P$A <- as(A,"CsparseMatrix")
  P$bc <- rbind(c(y, -0.5, 0.5), c(y, -0.5, 0.5))
  P$bx <- rbind(c(rep(0, 2 * p), rep(-Inf, n), rep(0, 3)),
                c(rep(Inf, 2 * p + n + 3)))
  P$cones <- matrix(list("QUAD",
                         c(n + 2 * p + 3, (2 * p + 1):(2 * p + n), n + 2 * p + 2)), 2, 1)
  rownames(P$cones) <- c("type", "sub")
  P$dparam$intpnt_nl_tol_rel_gap <- rtol
  z <- mosek(P, opts = list(verbose = verb))
  status <- z$sol$itr$solsta
  f <- z$sol$itr$xx
  coef <- f[1:p] - f[(p + 1):(2 * p)]
  resid <- f[(2 * p + 1):(2 * p + n)]
  list(coef = coef, resid = resid, status = status)
}