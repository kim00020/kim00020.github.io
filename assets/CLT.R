
n <- 10 

p <- 0.5 

r <- 10000

x <- rbinom(r,n,p)

mu <- n*p

sigma2 <- n*p*(1-p)

Npdf <- dnorm(seq(0,n,0.1), mu, sqrt(sigma2))

hist(x, breaks = seq(0,n,1), prob=T,
     right=F, main ="",
        xlim = c(0,n), ylim=c(0,max(Npdf) ) , 
     mgp = c(2, 0.5, 0), 
       cex = 1.0)

 mtext(side = 3, line=1, outer = F,  
  paste("Histogram for n= ", n, "  p= ", p) ) 
      
lines(seq(0, n, 0.1) + 0.5, Npdf, lty = 2, col="red")



