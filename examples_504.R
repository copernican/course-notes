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
