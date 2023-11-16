set.seed(116)
n <- 13
p <- 101
k <- 2
q <- 3
is_null       <- rep(FALSE, length = p)
is_null[1:57] <- TRUE

X <- matrix(stats::rnorm(n * q), nrow = n)
B <- matrix(stats::rnorm(q * p), nrow = q)
B[2, is_null] <- 0
Z <- X %*% matrix(stats::rnorm(q * k), nrow = q) +
     matrix(rnorm(n * k), nrow = n)
A <- matrix(stats::rnorm(k * p), nrow = k)
E <- matrix(stats::rnorm(n * p, sd = 1 / 2), nrow = n)
Y <- X %*% B + Z %*% A + E
mout <- mouthwash(Y = Y, X = X, k = k, cov_of_interest = 2,
                  include_intercept = FALSE)

data = rnaseq_data$head
mout <- mouthwash(Y = t(data$l2c), X = data$mod1, k = 1, cov_of_interest = 2,
                  include_intercept = FALSE)
names(mout)
mout$Zhat
residuals = getBatchResiduals(x = data$l2c, tissue = data$tissue, label = "l2c-mouthwash",
                              design = data$design, mod0 = data$mod0, sva = mout$Zhat, col = pull(data$covariates, treatment))
