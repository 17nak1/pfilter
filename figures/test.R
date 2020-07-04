
setwd("~/Git/pfilter/figures/compare20-20000")
js_loglik <- read.csv( "loglik-js.csv")
js_predmean <- read.csv("predmean-js.csv")
loglik_pf <- as.array(js_loglik[,1])
loglik_Js_est <- logmeanexp(loglik_pf)
loglik_js_span <- abs(max(loglik_pf)-min(loglik_pf))


r_loglik <- read.csv( "loglik-r.csv")
r_predmean <- read.csv("predmean-r.csv")
loglik_pf <- as.array(r_loglik[,1])
loglik_r_est <- logmeanexp(loglik_pf)
loglik_r_span <- abs(max(loglik_pf)-min(loglik_pf))

difference <- abs((loglik_Js_est - loglik_r_est)/loglik_Js_est)
js_predmean <- as.data.frame(js_predmean)
plot(r_predmean$S1, type ="l", ylab = "prediction mean", xlab="time step")
s <- seq(2,100,5)
for (i in s){
points(r_predmean[,i],  type ="l")
}

col_js <- seq(1,100,5)
for (i in col_js){
  points(js_predmean[,i], col="red", type ="l")
}
legend("topleft",legend=c("R", "JS"), col=c("black", "red"), lty=1:2)

# On my laptop 20 Times 20,000
Rlog20 <- c(-5527.559, -5565.793, -5381.379, -5598.128, -5506.104, -5501.564,
            -5697.990, -5558.378, -5432.255, -5526.587, -5440.010, -5627.149,
            -5466.653, -5726.683, -5500.616, -5467.925, -5667.815, -5509.453,
            -5757.359 ,-5579.408)
