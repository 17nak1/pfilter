/*
 * Particle filter
 * A plain vanilla sequential Monte Carlo (particle filter) algorithm.
 * Resampling is performed at each observation.
 * @param Np the number of particles to use.
 * @param tol positive numeric scalar.
 * @param max.fail integer; the maximum number of filtering failures allowed.
 * @param pred.mean logical; if \code{TRUE}, the prediction means are calculated for the state variables and parameters.
 * @param pred.var logical; if \code{TRUE}, the prediction variances are calculated for the state variables and parameters.
 * @param filter.mean logical; if \code{TRUE}, the filtering means are calculated for the state variables and parameters.
 * @param filter.traj logical; if \code{TRUE}, a filtered trajectory is returned for the state variables and parameters.
 * @param save.states logical.
 *  If \code{save.states=TRUE}, the state-vector for each particle at each time is saved.
 * @references
    pfilter.R by Aaron A. King. https://github.com/kingaa/pomp/blob/abc552b335319bd36b21d3313305c89e97f48f68/R/pfilter.R
    M. S. Arulampalam, S. Maskell, N. Gordon, & T. Clapp.
    A Tutorial on Particle Filters for Online Nonlinear, Non-Gaussian Bayesian Tracking.
    IEEE Trans. Sig. Proc. 50:174--188, 2002.
 */
//conditional log liklihood(time) =log(sum(w_i, i in 0:Np) / Np)
var START = new Date()
fs = require('fs')
let fmin = require ('fmin')
let mathLib = require('./mathLib')
let snippet = require('./modelSnippet.js')
let simulator = require('./simulator.js')
var linearInterpolator = require('linear-interpolator/node_main')

//////////////////////////////////////////data///////////////////////////////////////
// var LondonBidata = [], LondonCovar = []
// var rate = new Array(6) 
// // var params = [3.165652e+01 , 3.887624e-01 , 7.305000e+01 , 1.698730e-02  ,4.566000e+01,  4.813669e-01  ,1.963092e-01 , 2.831066e-03 ,3.476483e-04 ,  2.109135e-08,9.968213e-01]
// var params = [3.132490e+01 , 3.883620e-01 , 7.305000e+01 , 6.469830e-04 , 4.566000e+01 , 4.598709e-01 , 1.462546e-01 , 3.399189e-02 ,2.336327e-04 ,4.221789e-07 ,9.657741e-01 ]
// var times =[1940, 1944], maxFail = Infinity
// var Np = 100
// var nvars 
// var toler = 1e-17

// const anonymous = function(){ return -1 }

// //* 1st data set
// var London_covar = fs.readFileSync('./London_covar.csv').toString()
// var lines = London_covar.split('\n')
// for (let i = 1; i < lines.length; i++) {
//   LondonCovar.push(lines[i].split(','))
// }
// var dataCovar = [LondonCovar][0]
// //* 2nd data set
// var London_BiData = fs.readFileSync('./London_BiData1.csv').toString()
// var lines = London_BiData.split('\n')
// for (let i = 1; i < lines.length; i++) {
//   LondonBidata.push(lines[i].split(','))
// }
// var dataCases = [LondonBidata][0]

// var d1 = []// read time and population from 1st data and make interpolation function
// var d2 = []// read time and birthrate from 1st data and make interpolation function
// for (let i = 0; i < dataCovar.length - 1; i++) {
//   d1.push([Number(dataCovar[i][0]), Number(dataCovar[i][1])])
//   d2.push([Number(dataCovar[i][0]), Number(dataCovar[i][2])])
// }
// var interpolPop = linearInterpolator(d1)
// var interpolBirth = linearInterpolator(d2)
// var dt = 1 / 365.25
//////////////////////////////////////////////////////////////////////////////////////* main function//////////////////////////////////////////////////////////////////////
function pfilterCalculation (input) {//filter.traj , save.params
  // {params:inputArr, Np:100,times:times, dt:1 / 365.25,runPredMean:1,  dataCases:dataCases, interpolPop:interpolPopulation, interpolBirth:interpolBirth}
  let START =new Date()
  let defaults = {params:-1, Np:-1, tol:1e-17, maxFail:Infinity, runPredMean:0, runPredVar:0, runFilterMean:0, runSaveStates:0, times:-1, dt:-1,
                   dataCases:0}
  for(let prop in defaults) {
    if(typeof input[prop] == 'undefined') {
      input[prop] = defaults[prop]
    }
  }
  if (input.params === -1 || input.Np === -1 || input.times === -1 || input.dt === -1) {
    throw 'Some required arguments are missed'
  }
  var toler = 1e-17
  var params = input.params
  var Np = input.Np
  var tole = input.tol
  var maxFail = input.maxFail
  var dt = input.dt
  var dataCases = input.dataCases
  var interpolPop = input.interpolPop
  var interpolBirth = input.interpolBirth
  
  let [R0, amplitude, gamma, mu, sigma, rho, psi, S_0, E_0, I_0, R_0] = params
let [t0, tdata] = [1940, 1944]
let nvars = 5
let deltaT = 14 / 365.25
dt = 1 / 365.25
let doPredictionVariance = 0, doPredictionMean = 1, doFilterMean = 0 , allFail = 0

let timeLen = dataCases.length 
let nlost = 0

let rate = [], trans = []
let particles = [], state =[]
let sampleNum = Array.from(Array(Np).keys())
let condLoglik = []
let stateSaved =[]

let timeCountData = 0, ws ,w , vsq, sumsq, ess, loglik = 0, lik 

let predictionMean, predictionVariance, filterMean
let states = Array(Np).fill(null).map(() => Array(nvars))
let weights, normalWeights, S, E, I, R, del_t, ST, simulateValue
let modelCases, likvalue
if (doPredictionMean) {
  predictionMean = Array(timeLen).fill(null).map(() => Array(nvars))
}
if (doPredictionVariance) {
  predictionVariance = Array(timeLen).fill(null).map(() => Array(nvars))
}
if (doFilterMean) {
  filterMean = Array(timeLen).fill(null).map(() => Array(nvars))
}

state = snippet.initz(interpolPop(t0), S_0, E_0, I_0, R_0)
  // time loop
  for (k = t0; k <= Number(dataCases[timeLen - 2][0]) + deltaT / 3 ; k += deltaT){//Number(dataCases[timeLen - 2][0]) + deltaT / 3
  if ( k > tdata - deltaT && k <= tdata) {
    k = tdata
  }
  weights = []; normalWeights = []
  aa = [Np];
  if ( k >= t0) {
    for (np = 0; np < Np; np++) { // copy the particles
      aa[np] = [].concat(particles[sampleNum[np]])
      aa[np][nvars - 1] = 0
    }
  } 
  
  //**PARTICLE LOOP
  for (np = 0; np < Np; np++){ //calc for each particle
    trans = []

    particles[np] = simulator.simulate(aa[np], k, tdata, deltaT, dt, timeCountData, interpolPop, interpolBirth,params, Number(dataCases[dataCases.length - 1]), Number(dataCases[timeCountData + 1][0]))
    
    S = particles[np][0] 
    E = particles[np][1]
    I = particles[np][2] 
    R = particles[np][3] 
    H = particles[np][4] 

    // console.log(k)
     
    //***********weight*************
    if (k >= Number(dataCases[0][0])){
      if (stateSaved) {
        stateSaved.push(particles[np]) //[S,E,I,R,H])
      }
      modelCases = Number(dataCases[timeCountData][1])
      likvalue = snippet.dmeasure(rho, psi, H, modelCases, giveLog = 0)
      weights.push(likvalue)
      
    }
  }//  end particle loop
  
  //normalize
  if (k >= Number(dataCases[0][0])){
    let sumOfWeights = 0
    for (let i = 0; i < Np; i++) {
      sumOfWeights += weights[i]
    }
    for (let i = 0; i < Np; i++) {
      normalWeights[i] = weights[i] / sumOfWeights
    }
    // check the weights and compute sum and sum of squares
    w = 0, ws = 0, nlost = 0
    for (let i = 0; i < Np; i++) {
      if (weights[i] > toler) {
        w += weights[i]
        ws += weights[i] ** 2
      } else { // this particle is lost
        weights[i] = 0;
        nlost++
      }
    }
    if (nlost > maxFail) {
      throw 'execution terminated. The number of filtering failures exceeds the maximum number of filtering failures allowed. '
    }
    if (nlost >= Np) { 
      allFail = 1 // all particles are lost
    } else {
      allFail = 0
    }
    
    if (allFail) {
      lik = Math.log(toler) // minimum log-likelihood
      ess = 0  // zero effective sample size
    } else {
      ess = w * w / ws  // effective sample size
      lik = Math.log(w / Np)// mean of weights is likelihood
    }
    condLoglik[timeCountData] = [timeCountData + 1, lik]
    // the total conditional logliklihood in the time process is loglik
    loglik += lik
    mathLib.nosortResamp(Np, normalWeights, Np, sampleNum, 0)
    
      // Compute outputs
      for (let j = 0; j< nvars; j++) {
      // compute prediction mean
      if (doPredictionMean || doPredictionVariance) {
        let sum = 0, nlost = 0
        for (let nrow =0; nrow < Np; nrow++){
          if (particles[nrow][j]) {
            sum += particles[nrow][j]
          } else {
            nlost++
          }
        }
        sum /= Np
        predictionMean[timeCountData][j] = sum
      }  
      // compute prediction variance
      if (doPredictionVariance) {
        sumsq = 0
        for (let nrow = 0; nrow < Np; nrow++){
          if (particles[nrow][j]) {
            vsq = particles[nrow][j] - sum
            sumsq += Math.pow(vsq, 2)
          }
        }
        predictionVariance[timeCountData][j] = sumsq / (Np - 1) 
      }
      //  compute filter mean
      if (doFilterMean) {
        if (allFail) {   // unweighted average
          ws = 0
          for (let nrow =0; nrow < Np; nrow++){
            if (particles[nrow][j]) {
              ws += particles[nrow][j]
            }
          } 
          filterMean[timeCountData][j] = ws / Np//;console.log(ws / Np)
        } else {      // weighted average
          ws = 0
          for (let nrow =0; nrow < Np; nrow++){
            if (particles[nrow][j]) {
              ws += particles[nrow][j] * weights[nrow]
            }
          }
          filterMean[timeCountData][j] = ws / w
        }
      }
    }
    timeCountData++ 
  }
  }//endTime
  console.log('loglike=',loglik)
  console.log('runing time=', new Date() - START)
  activateDownload ()
   if (input.runPredMean) {
    return predictionMean
  }

  if (input.runPredVar) {
    return predictionVariance
  }

  if (input.runFilterMean) {
    return filterMean
  }
  
  if (input.runSaveStates){
    return stateSaved
  }

}

module.exports = {
  pfilterCalculation
}


    
      


