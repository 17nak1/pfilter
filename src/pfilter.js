
/**
 *  @file        pfilter.js
 *               A plain vanilla sequential Monte Carlo (particle filter) algorithm.
 *               Resampling is performed at each observation.
 *               
 *  @autor       Nazila Akhavan, nazila@kingsds.network
 *  @date        March 2019
 */

let fmin = require ('fmin')
let mathLib = require('./mathLib')
let snippet = require('./modelSnippet.js')
let simulator = require ('./simulator.js')


let pfilter = {}

pfilter.predictionMean = []
pfilter.predictionVariance = []

pfilter.run = function(input){

  let defaults = {params:-1, Np:-1, toler:1e-17, maxFail:Infinity, runPredMean:0, runPredVar:0, runFilterMean:0, runSaveStates:0, timeZero:-1, dt:-1,
    dataCases:0}
  for(let prop in defaults) {
    if(typeof input[prop] == 'undefined') {
      input[prop] = defaults[prop]
    }
  }

  if (input.params === -1 || input.Np === -1 || input.timeZero === -1 || input.dt === -1) {
    throw 'Some required arguments are missed'
  }
  let dataCases = input.dataCases
  let dataCovar = input.dataCovar
  let params = input.params
  let maxFail = input.maxFail
  let Np = input.Np
  let toler = input.toler
  let dt = input.dt
  let t0 = input.timeZero
  let deltaT = 14 / 365.25
  let tdata = dataCases[0][0]
  console.log("Np", Np)

  let d1 = [] // read time and population from 1st data and make interpolation function
  let d2 = [] // read time and birthrate from 1st data and make interpolation function
  for (let i = 0; i < dataCovar.length; i++) {
    d1.push([Number(dataCovar[i][0]), Number(dataCovar[i][1])])
    d2.push([Number(dataCovar[i][0]), Number(dataCovar[i][2])])
  }

  // Define variables
  let interpolPop = mathLib.interpolator (d1)
  let interpolBirth = mathLib.interpolator (d2)
  let START = new Date()
  let [R0, amplitude, gamma, mu, sigma, rho, psi, S_0, E_0, I_0, R_0] = params
  let nvars = 5 // number of states SEIRH
  let doPredictionVariance = 1, doPredictionMean = 1, runFilterMean = 0 , allFail = 0
  let timeLen = dataCases.length 
  let nlost = 0
  let particles = new Array(Np).fill(null).map(() => Array(5))
  let state =[]
  let sampleNum = Array.from(Array(Np).keys())
  let condLoglik = []
  let stateSaved = []
  let temp 
  let ws ,w , vsq, sum, sumsq, ess, loglik = 0, lik 
  let predictionMean, predictionVariance, filterMean
  let weights
  let modelCases, likvalue
  
  // Initiate Matrix for chosen outputs
  if (doPredictionMean) {
    predictionMean = Array(timeLen).fill(null).map(() => Array(nvars))
  }
  if (doPredictionVariance) {
    predictionVariance = Array(timeLen).fill(null).map(() => Array(nvars))
  }
  if (runFilterMean) {
    filterMean = Array(timeLen).fill(null).map(() => Array(nvars))
  }
  
  // Initial states from modelSnippet
  state = snippet.initz(interpolPop(t0), S_0, E_0, I_0, R_0)

  // First matrix of states at time t0; includes Np rows of repeated state 
  temp = new Array(Np).fill(null).map(() => [].concat(state))
  
  // Define t0 as the first time value 
  let k = t0

  /**
   *  Time loops
   *  The first loop simulating from t0 to first time value in reported data.
   *  The second one is based on times in reported data and calculates weights.  
   *  Starting at time t_n we have results at time t_{n+1}. But in the second loop k2 use the 
   *  results from the first loop which is in tdata and start calculing up to timeLen.
   */
  for (k = t0; k < tdata; k += deltaT){
    k2 = k + deltaT
    if ( k2 > tdata) {
      deltaT = tdata - k
      k2 = k + deltaT
    }
    particles = simulator.simulate (Np, temp, dt, interpolPop, interpolBirth, params, k, k2)
    for (np = 0; np < Np; np++) { // copy the particles
      temp[np] = [].concat(particles[sampleNum[np]])
    }
  }
  
  for (let timeCountData = 0; timeCountData < timeLen; timeCountData++){
    weights = []; normalWeights = []
    let k2 =  Number(dataCases[timeCountData][0])
 
    particles = simulator.simulate (Np, temp, dt, interpolPop, interpolBirth, params, k, k2)
    
    // Weights are calculated based on liklihood at each time for each point.
    for (np = 0; np < Np; np++){ 
      if (defaults.runSaveStates) {
        stateSaved.push([k2, ...particles[np]]) //[S,E,I,R,H])
      }
      modelCases = Number(dataCases[timeCountData][1])
      likvalue = snippet.dmeasure(rho, psi, particles[np][4], modelCases, giveLog = 0)
      weights.push(likvalue)           
    }
    
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
      // if (nlost > maxFail) {
      //   throw 'execution terminated. The number of filtering failures exceeds the maximum number of filtering failures allowed. '
      // }
      if (nlost >= Np) { 
        allFail = 1 // all particles are lost
      } else {
        allFail = 0
      }
      
      if (allFail) {
        lik = Math.log(toler)           // minimum log-likelihood
        ess = 0                         // zero effective sample size
      } else {
        ess = w * w / ws               // effective sample size
        lik = Math.log(w / Np)         // mean of weights is likelihood
      }
      condLoglik[timeCountData] = [timeCountData + 1, lik]
      // the total conditional logliklihood in the time process is loglik
      loglik += lik
      
      // Compute outputs
      for (let j = 0; j< nvars; j++) {
        // compute prediction mean
        if (doPredictionMean || doPredictionVariance) {
          sum = 0
          nlost = 0
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
        if (runFilterMean) {
          if (allFail) { //  unweighted average
            ws = 0
            for (let nrow =0; nrow < Np; nrow++){
              if (particles[nrow][j]) {
                ws += particles[nrow][j]
              }
            } 
            filterMean[timeCountData][j] = ws / Np
          } else { //  weighted average
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


      if (!allFail) {
        mathLib.nosortResamp(Np, weights, Np, sampleNum, 0)
        for (np = 0; np < Np; np++) { // copy the particles
          temp[np] = [].concat(particles[sampleNum[np]])
          temp[np][nvars - 1] = 0
        }
      } else {
        for (np = 0; np < Np; np++) { // copy the particles
          temp[np] = [].concat(particles[np])
          temp[np][nvars - 1] = 0
        }
      }

      k = k2
  }//endTime

pfilter.predictionMean = predictionMean
console.log(input.tempi)
console.log(loglik)
  
  let createCsvWriter = require('csv-writer').createArrayCsvWriter;
  let csvWriter = createCsvWriter({
    header: ['i', 'loglik', 'running time'],
    path: rootDir + '/samples/loglik.csv',
    append : true
  })
  csvWriter.writeRecords( [[input.tempi,  loglik, (new Date() - START) /1000]])
  .then(() => {
  console.log('...loglik '+input.tempi)
  })

  csvWriter = createCsvWriter({
    header: ['S', 'E', 'I', 'R', 'H'],
    path: rootDir + '/samples/Result/predmean'+input.tempi+'.csv'
  })  
  csvWriter.writeRecords(predictionMean)
    .then(() => {
    console.log('...predictionMean '+input.tempi)
  })

  csvWriter = createCsvWriter({
    header: ['S', 'E', 'I', 'R', 'H'],
    path: rootDir + '/samples/Result/predvar'+input.tempi+'.csv'
  })  
  csvWriter.writeRecords(predictionVariance)
    .then(() => {
    console.log('...predictionvar '+input.tempi)
  })

  csvWriter = createCsvWriter({
    header: ['lik'],
    path: rootDir + '/samples/Result/condLoglik'+input.tempi+'.csv'
  })  
  csvWriter.writeRecords(condLoglik)
    .then(() => {
    console.log('...condLoglik '+input.tempi)
  })
      
  console.log('running time:',(new Date() - START) /1000)

}

module.exports = pfilter
