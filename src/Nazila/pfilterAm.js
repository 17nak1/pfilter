

fs = require('fs')
let fmin = require ('fmin')
let mathLib = require('./mathLib')
let snippet = require('./modelSnippet.js')

function pfilterCalculation(){

//////////////////////////////////////////data///////////////////////////////////////
var dataCases = [], dataCovar = []
var rate = new Array(6) 
var params = [3.132490e+01, 3.883620e-01, 7.305000e+01, 6.469830e-04, 4.566000e+01, 4.598709e-01, 1.462546e-01, 3.399189e-02, 2.336327e-04, 4.221789e-07, 9.657741e-01 ]
var times =[1940, 1944], maxFail = Infinity
var Np =10
var nvars 
var toler = 1e-17


//* 1st data set
var London_covar = fs.readFileSync('~/../samples/London_covar.csv').toString()
var lines = London_covar.split('\n')
for (let i = 0; i < lines.length - 1; i++) {
  dataCovar.push(lines[i].split(','))
}
dataCovar.shift()
//* 2nd data set
var London_BiData = fs.readFileSync('~/../samples/London_BiData.csv').toString()
var lines = London_BiData.split('\n')
for (let i = 0; i < lines.length - 1; i++) {
  dataCases.push(lines[i].split(','))
}
dataCases.shift()

var d1 = []// read time and population from 1st data and make interpolation function
var d2 = []// read time and birthrate from 1st data and make interpolation function
for (let i = 0; i < dataCovar.length - 1; i++) {
  d1.push([Number(dataCovar[i][0]), Number(dataCovar[i][1])])
  d2.push([Number(dataCovar[i][0]), Number(dataCovar[i][2])])
}
var interpolPop = mathLib.interpolator(d1)
var interpolBirth = mathLib.interpolator(d2)
var START = new Date()


var [R0, amplitude, gamma, mu, sigma, rho, psi, S_0, E_0, I_0, R_0] = params
var paramsIC = [S_0, E_0, I_0, R_0, H = 0]
var [t0, tdata] = times
var nvars = paramsIC.length
var deltaT = 14 / 365.25
var dt = 1 / 365.25
var pop , birthrate
var timeLen = dataCases.length 
var va = 0, seas
var particles = [], state =[]
var nlost = 0
var sampleNum = new Array(Np)
var doPredictionVariance = 0, doPredictionMean = 1, doFilterMean = 0 , allFail = 0
var predictionMean = Array(timeLen).fill(null).map(() => Array(nvars))
var predictionVariance = Array(timeLen).fill(null).map(() => Array(nvars))
var condLoglik = Array(timeLen).fill(null).map(() => Array(1))
var filterMean = Array(timeLen).fill(null).map(() => Array(nvars))
var timeCountData = 0, ws = 0, vsq, sumsq, ess, loglik = 0, lik //condLoglik = []
var stateSaved =[]
var states = Array(Np).fill(null).map(() => Array(nvars))
var weights, normalWeights, trans, S, E, I, R, del_t, ST, simulateValue
var modelCases, likvalue
state = snippet.initz(interpolPop(t0), S_0, E_0, I_0, R_0)
// First Np sets
particles = new Array(Np).fill(null).map(() => [].concat(state))


// Time loop
for (k = t0; k < Number(dataCases[timeLen - 2][0]) + deltaT / 3; k += deltaT){
  if ( k > tdata - deltaT && k <= tdata) {
    k = tdata
  }
  lik = new Array(Np)
  weights = []; normalWeights = []

  //**PARTICLE LOOP
  for (np = 0; np < Np; np++){ //calc for each particle
    trans = new Array(6).fill(0)
    S = particles[np][0]; E = particles[np][1]; I = particles[np][2]; R = particles[np][3]; H = particles[np][4]
      
    // transitions between classes
    if (k <= tdata || k > 1965 - deltaT ) {
      steps = mathLib.numMapSteps(k, k + deltaT, dt)
    } else {
      steps = mathLib.numEulerSteps(k, Number(dataCases[timeCountData + 1][0]), dt)
    }
    del_t = (1 / steps )* deltaT 
    for (let stp = 0; stp < steps; stp++) { // steps in each time interval
      st = k + stp * del_t
      st = st.toFixed(10)
      simulateValue = snippet.rprocess(params, st, del_t, [S,E,I,R,H], interpolPop(st), interpolBirth(st))
      S = simulateValue[0]; E = simulateValue[1], I = simulateValue[2], R = simulateValue[3], H = simulateValue[4]
    }
    particles[np][0] = S
    particles[np][1] = E
    particles[np][2] = I
    particles[np][3] = R
    particles[np][4] = H
   
    states[np][0] = S || 0
    states[np][1] = E || 0
    states[np][2] = I || 0
    states[np][3] = R || 0
    states[np][4] = H || 0
     
    //***********RESAMPLE*************
    if (k >= Number(dataCases[0][0])){
      stateSaved.push([S,E,I,R,H])
      modelCases = Number(dataCases[timeCountData][1])
      likvalue = snippet.dmeasure(rho, psi, H, modelCases, giveLog = 0)
      weights.push(likvalue)
      particles[np][4] = 0
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
    var  w = 0, ws = 0, nlost = 0
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
    mathLib.nosortResamp(Np, normalWeights, Np, sampleNum, 0)//;console.log(k, sampleNum)
    for (np = 0; np < Np; np++) { // copy the particles
      particles[np] = [].concat(particles[sampleNum[np]])
      particles[np][nvars - 1] = 0
    }
    // Compute outputs
    for (let j = 0; j< nvars; j++) {
      // compute prediction mean
      if (doPredictionMean || doPredictionVariance) {
        var sum = 0, nlost = 0
        for (let nrow =0; nrow < Np; nrow++){
          if (states[nrow][j]) {
            sum += states[nrow][j]
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
          if (states[nrow][j]) {
            vsq = states[nrow][j] - sum
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
            if (states[nrow][j]) {
              ws += states[nrow][j]
            }
          } 
          filterMean[timeCountData][j] = ws / Np//;console.log(ws / Np)
        } else {      // weighted average
          ws = 0
          for (let nrow =0; nrow < Np; nrow++){
            if (states[nrow][j]) {
              ws += states[nrow][j] * weights[nrow]
            }
          }
          filterMean[timeCountData][j] = ws / w
        }
      }
    }
    timeCountData++ 
  }
}//endTime
console.log(loglik)
const createCsvWriter = require('csv-writer').createArrayCsvWriter;
const csvWriter = createCsvWriter({
  header: ['S', 'E', 'I', 'R', 'H'],
  path: '~/../samples/predmean.csv'
})
 
csvWriter.writeRecords(predictionMean)
  .then(() => {
  console.log('...predictionMean')
})


  
console.log('running time:',new Date() - START,'ms')
}

module.exports = {
  pfilterCalculation
}