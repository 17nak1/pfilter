let pomp = require('./src/pomp/index.js');
let pfilter = require('./src/pfilter/pfilter.js');
fs = require('fs')

pomp.pompclass.t0=1940;
pomp.pompclass.params.names = ['R0','amplitude','gamma','mu','sigma','rho','psi'];
pomp.pompclass.states.names = ['S','E','I','R','H'];
pomp.pompclass.states.zeronames = ['H'];
pomp.pompclass.params.current = { R0:3.132490e+01 , amplitude:3.883620e-01 , gamma:7.305000e+01 , mu:6.469830e-04 , sigma:4.566000e+01 ,rho: 4.598709e-01 ,psi: 1.462546e-01 ,S: 3.399189e-02 ,E:2.336327e-04 ,R:9.657741e-01,I:4.221789e-07};

pomp.pompclass.readData(pomp.pompclass.data,'./src/London_BiData.csv','time');
pomp.pompclass.readData(pomp.pompclass.covar,'./src/London_covar.csv','time');

pomp.pompclass.initializer = obj => {
    res = {};
    m = obj.pop/(obj.S+obj.E+obj.I+obj.R);
    res.S = Math.round(m*obj.S);
    res.E = Math.round(m*obj.E);
    res.I = Math.round(m*obj.I);
    res.R = Math.round(m*obj.R);
    res.H = 0;
    return res;
}

pomp.pompclass.rprocess = obj => {
                  var seas, beta, foi;
                  var births, va, tt;
                  var rate = [6], trans =[6];



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
                  //if( t<=1944.03832991103 && t>=1944.03832991102)
                  //printf(\"%f and %f \\n\", t, foi);
                  // Poisson births
                  births = rpois(birthrate*(1-va)*dt);
                  // transitions between classes
                  reulermultinom(2,S,rate[0],dt,trans[0]);
                  reulermultinom(2,E,rate[2],dt,trans[2]);
                  reulermultinom(2,I,rate[4],dt,trans[4]);
                  S += births - trans[0] - trans[1];
                  E += trans[0] - trans[2] - trans[3];
                  I += trans[2] - trans[4] - trans[5];
                  R = pop - S - E - I;
                  H += trans[4];           // true incidence
}


pfilter(pomp,10);

console.log('Finished');