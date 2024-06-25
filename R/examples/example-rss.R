# Generate training data with mixture error
set.seed(0)
n <- 100
p <- 10
p0 <- 5
ncontam <- 10
beta <- c(rep(1, p0), rep(0, p - p0))
x <- matrix(rnorm(n * p), n, p)
e <- rnorm(n, c(rep(10, ncontam), rep(0, n - ncontam)))
y <- x %*% beta + e

# Robust subset selection
ns_fit = heuristics_ns(x, y, k = p0, h = n - ncontam)

rss_ls(x, y, k = p0, h = n - ncontam, k.mio = p0, h.mio = n - ncontam, l_b = 2, l_a = 2)
fit <- rss(x, y, k.mio = p0, h.mio = n - ncontam)

# Extract model coefficients, generate predictions, and plot cross-validation results
coef(fit, k = p0, h = n - ncontam)
predict(fit, x[1:3, ], k = p0, h = n - ncontam)
plot(fit)


