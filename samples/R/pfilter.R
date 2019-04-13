cureentfolder <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))

start_time <- Sys.time()
######################################################  Model Snippet
rproc <- Csnippet("
                  double seas, beta, foi;
                  double births, va, tt;
                  double rate[6], trans[6];



                  //Vacination uptake
                  if (t < 1968)
                  va = 0;
                  else if (t>=1968 && t<=1969)
                  va = 0.33;
                  else if (t>=1969 && t<=1970)
                  va = 0.46;
                  else if (t>=1970 && t<=1971)
                  va = 0.51;
                  else if (t>=1971 && t<=1972)
                  va = 0.53;
                  else if (t>=1972 && t<=1973)
                  va = 0.52;
                  else if (t>=1973 && t<=1974)
                  va = 0.46;
                  else if (t>=1974 && t<=1975)
                  va = 0.46;
                  else if (t>=1975 && t<=1976)
                  va = 0.48;
                  else if (t>=1976 && t<=1977)
                  va = 0.48;
                  else if (t>=1977 && t<=1978)
                  va = 0.51;
                  else if (t>=1978 && t<=1979)
                  va = 0.53;
                  else if (t>=1979 && t<=1980)
                  va = 0.55;
                  else if (t>=1980 && t<=1981)
                  va = 0.58;
                  else if (t>=1981 && t<=1982)
                  va = 0.60;
                  else if (t>=1982 && t<=1983)
                  va = 0.63;
                  else if (t>=1983 && t<=1984)
                  va = 0.68;
                  else if (t>=1984 && t<=1985)
                  va = 0.71;
                  else if (t>=1985 && t<=1988)
                  va = 0.76;
                  else if (t>=1988 && t<=1989)
                  va = 0.814;
                  else if (t>=1989 && t<=1990)
                  va = 0.9488;
                  else if (t>=1990 && t<=1991)
                  va = 0.9818;
                  else if (t>=1991 && t<=1992)
                  va = 0.90;
                  else if (t>=1992 && t<=1993)
                  va = 0.92;
                  else if (t>=1993 && t<=1994)
                  va = 0.91;
                  else if (t>=1994 && t<=1995)
                  va = 0.91;
                  else if (t>=1995 && t<=1996)
                  va = 0.92;
                  else if (t>=1996 && t<=1997)
                  va = 0.92;
                  else if (t>=1997 && t<=1998)
                  va = 0.91;
                  else if (t>=1998 && t<=1999)
                  va = 0.88;
                  else if (t>=1999 && t<=2000)
                  va = 0.88;
                  else if (t>=2000 && t<=2001)
                  va = 0.87;
                  else if (t>=2001 && t<=2002)
                  va = 0.84;
                  else if (t>=2002 && t<=2003)
                  va = 0.82;
                  else if (t>=2003 && t<=2004)
                  va = 0.80;
                  else if (t>=2004 && t<=2005)
                  va = 0.81;
                  else if (t>=2005 && t<=2006)
                  va = 0.84;
                  else if (t>=2006 && t<=2007)
                  va = 0.85;
                  else if (t>=2007 && t<=2008)
                  va = 0.85;
                  else if (t>=2008 && t<=2009)
                  va = 0.85;
                  else if (t>=2009 && t<=2010)
                  va = 0.88;
                  else
                  va = 0.89;


                  // term-time seasonality
                  tt = (t-floor(t))*365.25;
                  if ((tt>=7&&tt<=100) || (tt>=115&&tt<=199) || (tt>=252&&tt<=300) || (tt>=308&&tt<=356))
                  seas = 1.0+amplitude*0.2411/0.7589;
                  else
                  seas = 1.0-amplitude;

                  // transmission rate
                  beta = R0*(gamma+mu)*(sigma+mu)*seas/sigma;  //seasonal transmission rate
                  // expected force of infection
                  foi = beta*I/pop;

                  rate[0] = foi;  //         force of infection
                  rate[1] = mu;             // natural S death
                  rate[2] = sigma;        // rate of ending of latent stage
                  rate[3] = mu;             // natural E death
                  rate[4] = gamma;        // recovery
                  rate[5] = mu;             // natural I death

                  //if( t>= 1964.9)
                  //printf(\"%f birthrate %f pop %f\\n\", t, birthrate, pop);
                  // Poisson births
                  births = rpois(birthrate*(1-va)*dt);
                  // transitions between classes
                  reulermultinom(2,S,&rate[0],dt,&trans[0]);
                  reulermultinom(2,E,&rate[2],dt,&trans[2]);
                  reulermultinom(2,I,&rate[4],dt,&trans[4]);
                  S += births - trans[0] - trans[1];
                  E += trans[0] - trans[2] - trans[3];
                  I += trans[2] - trans[4] - trans[5];
                  R = pop - S - E - I;
                  H += trans[4];           // true incidence
                  "
                  )

# The above uses true incidence for reporting

initz <- Csnippet("
                  double m = pop/(S_0+E_0+I_0+R_0);
                  S = nearbyint(m*S_0);
                  E = nearbyint(m*E_0);
                  I = nearbyint(m*I_0);
                  R = nearbyint(m*R_0);
                  H = 0;
                  ")
# Sampling from the normal approximation of the binomial distribution
dmeas <- Csnippet("
                  double m = rho*H;
                  double v = m*(1.0-rho+psi*psi*m);
                  double tol = 1.0e-18;
                  if (R_FINITE(cases)) {
                  if (cases > 0.0) {
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
                  } else {
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
                  }
                  if (give_log) lik = log(lik);
                  } else {
                  lik = (give_log) ? 0 : 1;
                  }
                  ")
rmeas <- Csnippet("

                  double m = rho*H;
                  double v = m*(1.0-rho+psi*psi*m);
                  double tol = 1.0e-18;
                  cases = rnorm(m,sqrt(v)+tol);
                  if (cases > 0.0) {
                  cases = nearbyint(cases);
                  } else {
                  cases = 0.0;
                  }
                  ")




toEst <- Csnippet("
                  Tmu = log(mu);
                  Tpsi = log(psi);
                  Tsigma = log(sigma);
                  Tgamma = log(gamma);
                  TR0 = log(R0);
                  Trho = logit(rho);
                  Tamplitude = logit(amplitude);
                  to_log_barycentric (&TS_0, &S_0, 4);
                  ")

fromEst <- Csnippet("
                    Tmu = exp(mu);
                    Tpsi = exp(psi);
                    Tsigma = exp(sigma);
                    Tgamma = exp(gamma);
                    TR0 = exp(R0);
                    Trho = expit(rho);
                    Tamplitude = expit(amplitude);
                    from_log_barycentric (&TS_0, &S_0, 4);
                    ")

params_mod <- c("R0","amplitude","gamma","mu","sigma","rho","psi")

params_ic <- c("S_0","E_0","R_0","I_0")

statenames <- c("S","E","I","R","H")

zeronames <- c("H")
#########################################            data
cureentfolder <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
London_BiData <- read.csv(file.path(cureentfolder, "London_BiData.csv"))
London_covar <- read.csv(file.path(cureentfolder, "London_covar.csv"))
#########################################       make pomp
pomp(
  data = London_BiData,
  times="time",
  t0=1940,#with(get(paste0(name,"_BiData")),2*time[1]-time[2]),
  rprocess = euler.sim(rproc,delta.t=1/365.25),
  rmeasure=rmeas,
  covar=London_covar,
  tcovar="time",
  dmeasure=dmeas,
  zeronames=zeronames,
  initializer=initz,
  toEstimationScale=toEst,
  fromEstimationScale=fromEst,
  statenames=statenames,
  paramnames=c(params_mod,params_ic)
) -> m1


current_params= c(R0=3.132490e+01 , amplitude=3.883620e-01 , gamma=7.305000e+01 , mu=6.469830e-04 , sigma=4.566000e+01 ,rho= 4.598709e-01 ,psi= 1.462546e-01 ,S_0= 3.399189e-02 ,E_0=2.336327e-04 ,R_0=9.657741e-01,I_0=4.221789e-07)
# pfilter(m1,params=current_params,Np=500) -> ss
# ss ye clasee ke toosh matrix hate mokhtalef dare. Mitooni spred.mean ro baraye moghayese estefade koni ke ye matrix 5* (toole baze zaman)
# hast. satre aval "S" ro mitooni ba in dastoor bebini.

setwd(cureentfolder)
tstart = 1
tend = 548
datasetj <- as.data.frame(read.csv('predmean.csv'))
pfilter(m1,params=current_params,Np=10,filter.mean = T,pred.mean=T, max.fail=3000) -> ss
datapredict <- as.data.frame(ss@pred.mean)
pm=c();for(i in tstart:tend){pm[i]=datapredict[1,i]}
plot(datasetj$S[c(tstart:tend)], type ="l",main="Np",col = "red", ylab = "JS")
points(pm[c(tstart:tend)], type ="l")

end_time <- Sys.time()
end_time - start_time
