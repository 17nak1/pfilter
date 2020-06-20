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
plot(r_predmean$H1[-1], type ="l",main= np,col = "red", ylab = "JS")
points(js_predmean$H14[-1],  type ="l")

