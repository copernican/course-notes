# why steepest descent gives bad directions
f <- function(x, y) 10^5 * x^2 + y^2
x <- seq(-1e-2, 1e-2, 1e-4)
y <- seq(-1, 1, 1e-2)
z <- outer(x, y, f)
contour(x = x, y = y, z = z, levels = c(0.01, 0.1, 0.25, 0.5, 1, 2.5, 5, 10), 
        drawlabels = F, cex.lab = 0.75, cex.axis = 0.75)

grad.f <- function(x, y) c(2 * 10^5 * x, 2 * y)
x0 <- c(1e-2, -1)
points(x0[1], x0[2])
grad.x0 <- grad.f(x0[1], x0[2])
x1 <- x0 + 5e-6 * -grad.x0
arrows(x0 = x0[1], y0 = x0[2], x1 = 0, y1 = 0, col = "blue", lty = 2,
       length = 0.1)
arrows(x0 = x0[1], y0 = x0[2], x1 = x1[1], y1 = x1[2], col = "red",
       length = 0.1)

# machine representation
dec2bin <- function(x) {   
  q <- x   
  r <- ""   
  while(q > 0 && nchar(q) < 52) {     
    r <- paste0(q %% 2, r)     
    q <- q %/% 2   
  }   
  return(r) 
}

# test the function by comparing its output to a known result
dec2bin(1e4) == "10011100010000"

x <- "10011100010000" # 1e4
x <- "10101" # 21
q <- as.integer(unlist(strsplit(x, split = NULL)))
r <- rev(seq_len(length(q)))[as.logical(q)]
sum(2^(r - 1))

bin2dec <- function(x) {
  q <- as.integer(unlist(strsplit(x, split = NULL)))
  r <- rev(seq_len(length(q)))[as.logical(q)]
  return(sum(2^(r - 1)))
}

bin2dec(x)
bin2dec(dec2bin(1023))
dec2bin(1023)

# change of basis via spectral theorem
lambda <- eigen(matrix(c(1,3,3,2), nr = 2))$vectors

plot.new()
plot.window(xlim = c(-1, 1), ylim = c(-1, 1.2), asp = 1)
arrows(x0 = 0, y0 = 0, x1 = 1, y1 = 0, length = 0.1, col = "black")
arrows(x0 = 0, y0 = 0, x1 = 0, y1 = 1, length = 0.1, col = "black")
arrows(x0 = 0, y0 = 0, x1 = lambda[1,1], y1 = lambda[2,1], length = 0.1, 
       lty = 2, col = "blue")
arrows(x0 = 0, y0 = 0, x1 = lambda[1,2], y1 = lambda[2,2], length = 0.1, 
       lty = 2, col = "blue")
text(x = 1.1, y = 0.1, cex = 0.75, labels = expression(e^(1)))
text(x = 0.1, y = 1.1, cex = 0.75, labels = expression(e^(2)))
text(x = lambda[1,1] + 0.1, y = lambda[2,1] + 0.1, cex = 0.75, 
     labels = expression(q^(1)), col = "blue")
text(x = lambda[1,2] - 0.1, y = lambda[2,2] + 0.1, cex = 0.75, 
     labels = expression(q^(2)), col = "blue")

# example of a convex function
curve(x^2, xlim = c(-2, 2), ylab = "f(x)", cex.lab = 0.75, 
      xaxt = "n", yaxt = "n")
x1 <- c(-0.5, 0.25)
x2 <- c(1, 1)
axis(side = 1, at = c(x1[1], x2[1]), cex.axis = 0.75, 
     labels = c(expression(x[1]), expression(x[2])))
points(x = c(-0.5, 1), y = c(0.25, 1))
segments(x0 = -0.5, y0 = 0.25, x1 = 1, y1 = 1, lty = 2)

# beginning of GS orthogonalization
plot.new()
plot.window(xlim = c(-0.2, 1), ylim = c(0, 1.2), asp = 1)
v1 <- c(0.9, 0.2)
v2 <- c(0.3, 0.8)
q1 <- v1 / norm(v1, type = "2")
arrows(x0 = 0, y0 = 0, x1 = q1[1], y1 = q1[2], length = 0.1, col = "black")
arrows(x0 = 0, y0 = 0, x1 = v2[1], y1 = v2[2], length = 0.1, col = "black")
u <- crossprod(v2, q1) * q1
arrows(x0 = 0, y0 = 0, x1 = u[1], y1 = u[2], length = 0.1, col = "black")
segments(x0 = v2[1], y0 = v2[2], x1 = u[1], y1 = u[2], lty = 2)
w <- v2 - u
arrows(x0 = 0, y0 = 0, x1 = w[1], y1 = w[2], length = 0.1, col = "black")
q2 <- w / norm(w, type = "2")
arrows(x0 = 0, y0 = 0, x1 = q2[1], y1 = q2[2], length = 0.1, col = "black")
arrows(x0 = v2[1], y0 = v2[2], x1 = w[1], y1 = w[2], lty = 2, length = 0.1)
text(x = q1[1] + 0.1, y = q1[2], labels = expression(q^(1)), cex = 0.75)
text(x = u[1], y = u[2] - 0.1, labels = "u", cex = 0.75)
text(x = v2[1] + 0.1, y = v2[2], labels = expression(v^(2)), cex = 0.75)
text(x = w[1] + 0.1, y = w[2] + 0.1, labels = "-u", cex = 0.75)
text(x = q2[1], y = q2[2] + 0.1, labels = expression(q^(2)), cex = 0.75)
text(x = 0.1, y = 0.1, labels = expression(theta), cex = 0.75)
