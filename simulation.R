

## s1 simulation data generating
set.seed(1)
a <- rnorm(200, 2, 0.1)
b <- rnorm(200, 5, 0.1)
c <- rnorm(200, 8, 0.1)
d <- rnorm(200, 10, 0.1)
e <- rnorm(200, 20, 0.1)
s1 <- data.frame(a=a, b=b, c=c, d=d, e=e)
a <- rnorm(500, 8, 0.1)
b <- rnorm(500, 2, 0.1)
c <- rnorm(500, 5, 0.1)
d <- rnorm(500, 10, 0.1)
e <- rnorm(500, 20, 0.1)
s1 <- rbind(s1, data.frame(a=a, b=b, c=c, d=d, e=e))
a <- rnorm(1000, 20, 0.1)
b <- rnorm(1000, 10, 0.1)
c <- rnorm(1000, 8, 0.1)
d <- rnorm(1000, 5, 0.1)
e <- rnorm(1000, 2, 0.1)
s1 <- rbind(s1, data.frame(a=a, b=b, c=c, d=d, e=e))
cluster <- c(rep(1, 200), rep(2, 500), rep(3, 1000))





## s2 simulation data generating
set.seed(2)
a <- rnorm(200, 2, 2)
b <- rnorm(200, 5, 2)
c <- rnorm(200, 8, 2)
d <- rnorm(200, 10, 2)
e <- rnorm(200, 20, 2)
s2 <- data.frame(a=a, b=b, c=c, d=d, e=e)
a <- rnorm(500, 8, 0.5)
b <- rnorm(500, 2, 0.5)
c <- rnorm(500, 5, 0.5)
d <- rnorm(500, 10, 0.5)
e <- rnorm(500, 20, 0.5)
s2 <- rbind(s2, data.frame(a=a, b=b, c=c, d=d, e=e))
a <- rnorm(1000, 20, 0.1)
b <- rnorm(1000, 10, 0.1)
c <- rnorm(1000, 8, 0.1)
d <- rnorm(1000, 5, 0.1)
e <- rnorm(1000, 2, 0.1)
s2 <- rbind(s2, data.frame(a=a, b=b, c=c, d=d, e=e))
cluster <- c(rep(1, 200), rep(2, 500), rep(3, 1000))
















# bandwidth <- seq(1,10, by = 0.1)
# kn <- c()
# ks <- c()
# for(b in bandwidth){
#     agg_cluster <- densityClustX(agg[ ,c(1,2)], dc = b, alpha = 0.01)
#     kn <- c(kn, length(unique(agg_cluster$cluster)))
#     fitrho <- fitdist(agg_cluster$rho, "norm", method="mme")
#     ks <- c(ks, gofstat(fitrho)$ks)
# }