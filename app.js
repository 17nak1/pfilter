let pomp = require('./similar-to-R');

fs = require('fs')
rootDir = '.'

rproc = function (params, t, del_t, [S,E,I,R,H], pop, birthrate) {
  let seas, beta, beta0, foi, R0, tt, va
  let trans = new Array(6).fill(0)
  let rate = new Array(6) 
  
  R0 = params[0], amplitude = params[1], gamma = params[2], mu = params[3], sigma = params[4] 
  beta0 = R0 * (gamma + mu) * (sigma + mu) / sigma
  
  va = 0;
  tt = (t - Math.floor(t)) * 365.25
  if ((tt >= 7 && tt <= 100) || (tt >= 115 && tt <= 199) || (tt >= 252 && tt <= 300) || (tt >= 308 && tt <= 356)) {
    seas = 1 + amplitude * 0.2411 / 0.7589
  } else {
    seas = 1 - amplitude
  }                 
  beta = R0 * (gamma + mu) * (sigma + mu) * seas / sigma  //seasonal transmission rate
  foi = beta * I / pop
  rate[0] = foi            //force of infection
  rate[1] = mu             // natural S death
  rate[2] = sigma          // rate of ending of latent stage
  rate[3] = mu             // natural E death
  rate[4] = gamma          // recovery
  rate[5] = mu             // natural I death 
   
  let births = rpois.rpoisOne(birthrate * (1 - va) * del_t )// Poisson births
  mathLib.reulermultinom(2, Math.round(S), 0, del_t, 0, rate, trans)
  mathLib.reulermultinom(2, Math.round(E), 2, del_t, 2, rate, trans)
  mathLib.reulermultinom(2, Math.round(I), 4, del_t, 4, rate, trans)
  S += (births - trans[0] - trans[1])
  E += (trans[0] - trans[2] - trans[3]) 
  I += (trans[2] - trans[4] - trans[5]) 
  R = pop - S - E - I
  H += trans[4] 
  return [S, E, I, R, H]
};
initz = function(input) {
  let m = input.pop / (input.S_0 + input.E_0 + input.R_0 + input.I_0),
    S = Math.round(m * input.S_0),
    E = Math.round(m * input.E_0),
    I = Math.round(m * input.I_0),
    R = Math.round(m * input.R_0),
    H = 0
  return {S, E, I, R, H}
};
dmeas = function (rho, psi, H, dCases, giveLog = 1) {
  let lik
  let mn = rho * H
  let v = mn * (1.0 - rho + psi * psi * mn)
  let tol = 1.0e-18
  let modelCases = Number(dCases)
  if(!isNaN(modelCases)){
    if (modelCases > 0.0) {
      lik = mathLib.pnorm(modelCases + 0.5, mn, Math.sqrt(v) + tol, 1, 0) - mathLib.pnorm(modelCases - 0.5, mn, Math.sqrt(v) + tol, 1, 0) + tol
    } else {
      lik = mathLib.pnorm((modelCases + 0.5, mn, Math.sqrt(v) + tol)) + tol
    }
    if (giveLog) lik = Math.log(lik)
  } else {
    lik = (giveLog) ? 0 : 1;
  }
  return lik
}
;
rmeas = function (H, rho, psi) {
  let mn = rho * H
  let v = mn * (1.0 - rho + psi * psi * mn)
  let tol = 1.0e-18
  let cases = mathLib.rnorm(mn, Math.sqrt(v) + tol)
  if (cases > 0) {
    cases = Math.round(cases)
  } else {
    cases = 0
  }
  return cases
};
toEst = {};
fromEst = {};

params_mod = ["R0","amplitude","gamma","mu","sigma","rho","psi"];
params_ic = ["S_0","E_0","R_0","I_0"]
statenames =  ["S","E","I","R","H"]
zeronames = ["H"]

let London_BiData = []
let London_covar = []
let London_covar_file = fs.readFileSync(rootDir+'/samples/London_covar.csv').toString()
London_covar_file = London_covar_file.replace(/['"]+/g, '');
let lines = London_covar_file.split(/\r\n|\n/)
for (let i = 0; i < lines.length; i++) {
  London_covar.push(lines[i].split(','))
}
if (London_covar[London_covar.length - 1].length === 1 ) {
  London_covar.pop()
}

let London_BiData_file = fs.readFileSync(rootDir+'/samples/London_BiDataMain.csv').toString()
London_BiData_file = London_BiData_file.replace(/['"]+/g, '');
lines = London_BiData_file.split(/\r\n|\n/)
for (let i = 0; i < lines.length ; i++) {
  London_BiData.push(lines[i].split(','))
}
if (London_BiData[London_BiData.length - 1].length === 1 ) {
  London_BiData.pop()
}

let m1 = new pomp({
  data: London_BiData,
  times:"time",
  t0: 1940,
  rprocess: rproc,//euler.sim(rproc,delta.t=1/365.25),
  rmeasure: rmeas,
  covar: London_covar,
  tcovar: "time",
  dmeasure: dmeas,
  zeronames: zeronames,
  initializer: initz,
  toEstimationScale: toEst,
  fromEstimationScale: fromEst,
  statenames: statenames,
  paramnames: [...params_mod, ...params_ic]
});

current_params= {R0: 3.132490e+01 , amplitude: 3.883620e-01 , gamma: 7.305000e+01 , mu: 6.469830e-04 , sigma: 4.566000e+01 ,rho:  4.598709e-01 ,psi:  1.462546e-01 ,S_0:  3.399189e-02 ,E_0: 2.336327e-04 ,R_0: 9.657741e-01,I_0: 4.221789e-07}
np= 20;

ss = pomp.pfilter({
  object: m1,
  params: current_params,
  Np: np,
  filter_mean: true,
  pred_mean: true,
  max_fail: 3000
});
console.log('finish');