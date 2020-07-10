/**
 *  @file        modeSnippet.js
 *               Makes the model and its dependencies.
 *                 
 *  @autor       Nazila Akhavan, nazila@kingsds.network
 *  @date        Feb 2019
 */

snippet = {}
let mathLib = require('./mathLib');
let rpois = require('./rpois');

// Order of index in parameters
let pIndex = {
  "R0": 0,
  "amplitude": 1,
  "gamma": 2,
  "mu": 3,
  "sigma": 4,
  "rho": 5,
  "psi": 6,
  "S": 7,
  "E": 8,
  "I": 9,
  "R":10,
  "H": 11
};

// Order of index in states
let sIndex = {
  "S_0": 0,
  "E_0": 1,
  "I_0": 2,
  "R_0":3,
  "H": 4
}

let dataIndex = {
  "cases": 0
}


snippet.rprocess = function (pomp, states, params, t, deltaT, interpolatorObj) {
  
  let S = states.S;
  let E = states.E;
  let I = states.I;
  let R = states.R;
  let H = states.H;

  let R0 = params.R0;
  let amplitude = params.amplitude;
  let gamma = params.gamma;
  let mu = params.mu;
  let sigma = params.sigma ;

  let pop =interpolatorObj.pop;
  let birthrate = interpolatorObj.birthrate;
  let seas, beta, beta0, foi, tt, va;
  let length = 3;
  let trans = new Array(length * 2).fill(0);
  let rate = new Array(length * 2); 

  
  
  beta0 = R0 * (gamma + mu) * (sigma + mu) / sigma;
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
   
  let births = rpois.rpoisOne(birthrate * (1 - va) * deltaT )// Poisson births
  mathLib.reulermultinom(2, Math.round(S), 0, deltaT, 0, rate, trans)
  mathLib.reulermultinom(2, Math.round(E), 2, deltaT, 2, rate, trans)
  mathLib.reulermultinom(2, Math.round(I), 4, deltaT, 4, rate, trans)
  S += (births - trans[0] - trans[1])
  E += (trans[0] - trans[2] - trans[3]) 
  I += (trans[2] - trans[4] - trans[5]) 
  R = pop - S - E - I
  H += trans[4] 
  return {S, E, I, R, H}
}

snippet.initz = function(args, covar) {
  
  let m = covar.pop / (args.S_0 + args.E_0 + args.R_0 + args.I_0);
  let S_0 = Math.round(m * args.S_0);
  let E_0 = Math.round(m * args.E_0);
  let I_0 = Math.round(m * args.I_0);
  let R_0 = Math.round(m * args.R_0);
  let H_0 = 0;
  return {S:S_0, E:E_0, I:I_0, R:R_0, H:H_0};
}

snippet.dmeasure = function (pomp, data , hidden_state, params, giveLog = 1) {
  
  let lik
  let rho = params.rho;
  let psi = params.psi;
  let H = hidden_state.H;
  let cases = data[0];//[dataIndex['cases

  let tol = 1.0e-18
  let mn = rho * H;
  let v = mn * (1.0 - rho + psi * psi * mn);
  
  let modelCases = Number(cases);
  if(!isNaN(modelCases)){
    if (modelCases > 0.0) {
      lik = mathLib.pnorm(modelCases + 0.5, mn, Math.sqrt(v) + tol, 1, 0) - mathLib.pnorm(modelCases - 0.5, mn, Math.sqrt(v) + tol, 1, 0) + tol
    } else {
      lik = mathLib.pnorm((modelCases + 0.5, mn, Math.sqrt(v) + tol)) + tol;
    }
    if (giveLog) lik = Math.log(lik);
  } else {
    lik = (giveLog) ? 0 : 1;
  }
  return lik
}

snippet.rmeasure = function (H, rho, psi) {
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
}
/**
 * paramnames is an array of all parameters in the model and all code is based on the order
 * of params in this array. The first index is zero i.e R0 = params[0]
 */
snippet.paramnames = ["R0","amplitude","gamma","mu","sigma","rho","psi", "S", "E", "I", "R"];
snippet.zeronames = ["H"];
snippet.statenames = ["S","E","I","R","H"];

snippet.toEst = function(params) {
  let mu = Math.log(params[3]);
  let psi = Math.log(params[6]);
  let sigma = Math.log(params[4]);
  let gamma = Math.log(params[2]);
  let R0 = Math.log(params[0]);
  let rho = mathLib.logit(params[5]);
  let amplitude = mathLib.logit(params[1]);
  let states = mathLib.toLogBarycentric([params[7], params[8], params[9], params[10]]);
  //Parameters order should be the same as paramnames.
  return [R0, amplitude, gamma, mu, sigma, rho, psi, ...states];
}

snippet.fromEst = function(params) {
  let mu = Math.exp(params[3]);
  let psi = Math.exp(params[6]);
  let sigma = Math.exp(params[4]);
  let gamma = Math.exp(params[2]);
  let R0 = Math.exp(params[0]);
  let rho = mathLib.expit(params[5]);
  let amplitude = mathLib.expit(params[1]);
  let states = mathLib.fromLogBarycentric([params[7], params[8], params[9], params[10]]);
  //Parameters order should be the same as paramnames.
  return [R0, amplitude, gamma, mu, sigma, rho, psi, ...states];
}

module.exports = snippet
