let pomp = require('./similar-to-R');
let rpois = require('./src/rpois');
let mathLib = require('./src/mathLib');

fs = require('fs')
rootDir = '.'

rproc = function (args) {
  let seas, beta, beta0, foi,  tt, va
  let trans = new Array(6).fill(0)
  let rate = new Array(6) 
  
  beta0 = args.R0 * (args.gamma + args.mu) * (args.sigma + args.mu) / args.sigma;
  
  va = 0;
  tt = (args.t - Math.floor(args.t)) * 365.25
  if ((tt >= 7 && tt <= 100) || (tt >= 115 && tt <= 199) || (tt >= 252 && tt <= 300) || (tt >= 308 && tt <= 356)) {
    seas = 1 + args.amplitude * 0.2411 / 0.7589
  } else {
    seas = 1 - args.amplitude
  }                 
  beta = args.R0 * (args.gamma + args.mu) * (args.sigma + args.mu) * seas / args.sigma  //seasonal transmission rate
  foi = beta * args.I / args.pop
  rate[0] = foi            //force of infection
  rate[1] = args.mu             // natural S death
  rate[2] = args.sigma          // rate of ending of latent stage
  rate[3] = args.mu             // natural E death
  rate[4] = args.gamma          // recovery
  rate[5] = args.mu             // natural I death 
   
  let births = rpois.rpoisOne(args.birthrate * (1 - va) * args.dt )// Poisson births
  mathLib.reulermultinom(2, Math.round(args.S), 0, args.dt, 0, rate, trans)
  mathLib.reulermultinom(2, Math.round(args.E), 2, args.dt, 2, rate, trans)
  mathLib.reulermultinom(2, Math.round(args.I), 4, args.dt, 4, rate, trans)
  args.S += (births - trans[0] - trans[1])
  args.E += (trans[0] - trans[2] - trans[3]) 
  args.I += (trans[2] - trans[4] - trans[5]) 
  args.R = args.pop - args.S - args.E - args.I
  args.H += trans[4] 
  return args;
};
initz = function(args) {
  let m = args.pop / (args.S_0 + args.E_0 + args.R_0 + args.I_0),
    S = Math.round(m * args.S_0),
    E = Math.round(m * args.E_0),
    I = Math.round(m * args.I_0),
    R = Math.round(m * args.R_0),
    H = 0
  return {S, E, I, R, H}
};
dmeas = function (args) {
  let lik
  let mn = args.rho *args. H
  let v = mn * (1.0 - args.rho + args.psi * args.psi * mn)
  let tol = 1.0e-18
  let modelCases = Number(args.modelCases)
  if(!isNaN(modelCases)){
    if (modelCases > 0.0) {
      lik = mathLib.pnorm(modelCases + 0.5, mn, Math.sqrt(v) + tol, 1, 0) - mathLib.pnorm(modelCases - 0.5, mn, Math.sqrt(v) + tol, 1, 0) + tol
    } else {
      lik = mathLib.pnorm((modelCases + 0.5, mn, Math.sqrt(v) + tol)) + tol
    }
    if (args.giveLog) lik = Math.log(lik)
  } else {
    lik = (args.giveLog) ? 0 : 1;
  }
  return lik
};
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


myPomp = new pomp({
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
np= 20000;

ss = myPomp.pfilter({
  params: current_params,
  Np: np,
  filter_mean: true,
  pred_mean: true,
  max_fail: 3000,
  pred_mean: true
});

console.log(ss);