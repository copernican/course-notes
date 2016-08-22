# SPF example: lecture 4, slide 22
spf <- read.table("data/spf.txt", header = T, sep = "\t")
boxplot(spf, ylab = "Tolerance (min)", cex.lab = 0.75, cex.axis = 0.75)

y <- log(spf$SUNSCREEN / spf$CONTROL)
ybar <- mean(y)
n <- length(y)

# prior hyperparameters
m0 <- 0
k0 <- 0.1
v0 <- 10
s0 <- 4

# posterior hyperparameters
kn <- k0 + n
mn <- (k0 * m0 + n * ybar) / kn
vn <- v0 + n
sn <- (v0 * s0 + (n - 1) * var(y) + (k0 * n / kn) * (ybar - m0)^2) / vn

# draw samples
nsamp <- 1e5
phi <- rgamma(n = nsamp, shape = vn / 2, rate = sn * vn / 2)
sigma2 <- 1 / phi
mu <- rnorm(nsamp, mn, sqrt(sigma2 / kn))
mu.t <- rt(nsamp, vn) * sqrt(sn / kn) + mn
qqplot(mu, mu.t)
abline(a = 0, b = 1, col = "red")
