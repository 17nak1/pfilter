let pomp = require('./similar-to-R');
let rpois = require('./similar-to-R/rpois')
let mathLib = require('./similar-to-R/mathLib')

fs = require('fs')
rootDir = '.'

rproc = function (input) {
  let seas, beta, beta0, foi,  tt, va
  let trans = new Array(6).fill(0)
  let rate = new Array(6) 
  
  beta0 = input.R_0 * (input.gamma + input.mu) * (input.sigma + input.mu) / input.sigma;
  
  va = 0;
  tt = (input.t - Math.floor(input.t)) * 365.25
  if ((tt >= 7 && tt <= 100) || (tt >= 115 && tt <= 199) || (tt >= 252 && tt <= 300) || (tt >= 308 && tt <= 356)) {
    seas = 1 + input.amplitude * 0.2411 / 0.7589
  } else {
    seas = 1 - input.amplitude
  }                 
  beta = input.R_0 * (input.gamma + input.mu) * (input.sigma + input.mu) * seas / input.sigma  //seasonal transmission rate
  foi = beta * input.I / input.pop
  rate[0] = foi            //force of infection
  rate[1] = input.mu             // natural S death
  rate[2] = input.sigma          // rate of ending of latent stage
  rate[3] = input.mu             // natural E death
  rate[4] = input.gamma          // recovery
  rate[5] = input.mu             // natural I death 
   
  let births = rpois.rpoisOne(input.birthrate * (1 - va) * input.dt )// Poisson births
  mathLib.reulermultinom(2, Math.round(input.S), 0, input.dt, 0, rate, trans)
  mathLib.reulermultinom(2, Math.round(input.E), 2, input.dt, 2, rate, trans)
  mathLib.reulermultinom(2, Math.round(input.I), 4, input.dt, 4, rate, trans)
  input.S += (births - trans[0] - trans[1])
  input.E += (trans[0] - trans[2] - trans[3]) 
  input.I += (trans[2] - trans[4] - trans[5]) 
  input.R = input.pop - input.S - input.E - input.I
  input.H += trans[4] 
  return input;
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
dmeas = function (input) {
  let lik
  let mn = input.rho *input. H
  let v = mn * (1.0 - input.rho + input.psi * input.psi * mn)
  let tol = 1.0e-18
  let modelCases = Number(input.modelCases)
  if(!isNaN(modelCases)){
    if (modelCases > 0.0) {
      lik = mathLib.pnorm(modelCases + 0.5, mn, Math.sqrt(v) + tol, 1, 0) - mathLib.pnorm(modelCases - 0.5, mn, Math.sqrt(v) + tol, 1, 0) + tol
    } else {
      lik = mathLib.pnorm((modelCases + 0.5, mn, Math.sqrt(v) + tol)) + tol
    }
    if (input.giveLog) lik = Math.log(lik)
  } else {
    lik = (input.giveLog) ? 0 : 1;
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
np= 3;

ss = pomp.pfilter({
  object: m1,
  params: current_params,
  Np: np,
  filter_mean: true,
  pred_mean: true,
  max_fail: 3000,
  pred_mean: true
});
console.log('finish');