# Generate example data
     set.seed(985115)
     n <- 350
     p <- 16
     sigmasq_error <- 0.5
     zeroes <- rbinom(p, 1, 0.6)
     beta.true <- rnorm(p, 1, sd = 4)
     beta.true[zeroes] <- 0
     
     X <- cbind(matrix(rnorm(n * p), nrow = n))
     X <- scale(X, center = TRUE, scale = FALSE)
     y <- X %*% matrix(beta.true, ncol = 1) + rnorm(n, 0, sqrt(sigmasq_error))
     y <- scale(y, center = TRUE, scale = FALSE)
     
     # Calculate sufficient statistics
     XtX <- t(X) %*% X
     Xty <- t(X) %*% y
     yty <- t(y) %*% y
     
     # Set the prior
     K <- 9
     sigma0 <- c(0.001, .1, .5, 1, 5, 10, 20, 30, .005)
     omega0 <- rep(1/K, K)
     
     # Calculate summary statistics
     b.hat <- sapply(1:p, function(j) { summary(lm(y ~ X[, j]))$coefficients[-1, 1] })
     s.hat <- sapply(1:p, function(j) { summary(lm(y ~ X[, j]))$coefficients[-1, 2] })
     R.hat <- cor(X)
     var_y <- var(y)
     sigmasq_init <- 1.5
     
     # Run PRS CS
     sumstats = list(BETA=b.hat, MAF=rep(0.5, length(b.hat)))
     LD <- list(blk1 = R.hat)
     write.table(data.frame(sumstats), "sumstats.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
	write.table(LD$blk1, "LD.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
