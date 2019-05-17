/**
 *  @file       traj.js        This function attempts to match trajectories of a model's deterministic skeleton to data.
 *                             Trajectory matching is equivalent to maximum likelihood estimatedation under the assumption 
 *                             that process noise is entirely absent, i.e., that all stochasticity is measurement error.
 *                             Accordingly, this method uses only the skeleton and dmeasure components of a POMP model.
 *
 *  @author     Nazila Akhavan
 *  @date       Jan 2019
 */

let snippet = require('./modelSnippet.js')
let mathLib = require('./mathLib')
let fmin    = require('fmin')
let fs = require('fs')

var LondonBidata, LondonCovar
var param = fs.readFileSync('./ParamSet_DeterministicSEIR_run28.csv').toString()
var params = []
var lines = param.split('\n')
for (let i = 1; i < lines.length - 1; i++) {
  params.push(lines[i].split(','))
}

var index =  new Array(12)
index[1] = 1; 
index[3] = 1; index[5] = 1; index[6] = 1

//* 1st data set
var London_covar = fs.readFileSync('./London_covar.csv').toString()
var dataCovar = []
var lines = London_covar.split('\n')
for (let i = 1; i < lines.length; i++) {
  dataCovar.push(lines[i].split(','))
}
//* 2nd data set
dataCases = []
var London_BiData = fs.readFileSync('./London_BiData.csv').toString()
var lines = London_BiData.split('\n')
for (let i = 1; i < lines.length; i++) {
  dataCases.push(lines[i].split(','))
}
let times = [1940, Number(dataCases[0][0]), Number(dataCases[dataCases.length - 2][0])]
var d1 = []// read time and population from 1st data and make interpolation function
var d2 = []// read time and birthrate from 1st data and make interpolation function
for (let i = 0; i < dataCovar.length - 1; i++) {
  d1.push([Number(dataCovar[i][0]), Number(dataCovar[i][1])])
  d2.push([Number(dataCovar[i][0]), Number(dataCovar[i][2])])
}
var interpolPopulation = mathLib.interpolator(d1)
var interpolBirth = mathLib.interpolator(d2)
params = params[0]


function traj_match (interpolPopulation, interpolBirth, dataCases, params, times, index) {
  let deltaT = 0.03832991102
  var tempIndex = 0
  var estimated = []
  var place = []
  var states = []
  var data1 = []
  var data2 = []
  var solution
  
   //*Change the initial values of estimating parameters(with index one) to the log or logit scale.
  // From those amplitude and rho (parameters number 1 and 5) are in logit scale and the rest are in log scale
  for (let i = 0; i < params.length; i++) {
    params[i] = Number(params[i])
    if (index[i] === 1) {
      place.push(i)
      if ((i === 1) || (i === 5)) {
        estimated.push(Math.log(params[i] / (1 - params[i]))) //logit scale
      } else {
        estimated.push(Math.log(params[i])) //log scale
      }
    }
  }


  //* Optimizer function using Nelder Mead method
  solution = fmin.nelderMead(logLik, estimated)
  for (let j = 0;j < params.length; j++) {
    if (index[j] === 1) { // Using exp and expit to get back to the regular scale.
      if ((j === 1) || (j === 5)){
        params[j] = 1/ (1 + Math.exp(-solution.x[tempIndex]))
      } else {
        params[j] = Math.exp(solution.x[tempIndex])
      }
      tempIndex++
    }
  }
  

  //* calculate log likelihood
  function logLik (estimated) {
    var likvalue = 0
    var loglik = 0
    var rho 
    var psi 
    for (let i = 0; i < estimated.length; i++) {
      if ((place[i] === 1) || (place[i] === 5)) { //Change to the exp scale and let optimizer to search all real numbers.
        params[place[i]] = 1 / (1 + Math.exp(-estimated[i]))
      } else {
        params[place[i]] = Math.exp(estimated[i])
      }
    }
    rho = params[5]
    psi = params[6]

    var simH = integrate(interpolPopulation, interpolBirth, params, times, deltaT)
    for (let i = 0; i < simH.length; i++) {
      likvalue = snippet.dmeasure(rho, psi, simH[i], dataCases[i][1])
      loglik = loglik + Math.log(likvalue)
    }
    ;console.log(params, loglik)
    return [-(loglik).toFixed(6)]
  }
  console.log(params, -solution.fx)
  return[params, -solution.fx]
}

//* ODE solver
function integrate (interpolPopulation, interpolBirth, params, times, deltaT) {
  var steps = 1000 // Total number of steps in the Euler method
  var t0 = times[0]
  var dataStartTime = times[1]
  var dataEndTime = times[2]
  var rho = params[5]
  var psi = params[6]
  // var states = []
  var arr = []
  var pop 
  var birthrate
  var timetemp
  var Npre

  var N = snippet.initz(interpolPopulation(t0), params[7], params[8], params[9], params[10])
  // N = [1.030473e+05, 1.223305e+02, 1.294772e+06, 8.176444e+01, 2.049497e+05 ]
  var  k = t0 , count
  var flag = 0
  dt = deltaT
  while ( flag === 0 ) {
    Npre = N
    for (let stp = 0; stp < steps; stp++) { 
      pop = interpolPopulation(k + stp / steps * dt)
      birthrate = interpolBirth(k + stp / steps * dt)
      N45 = mathLib.odeMethod('rkf45', snippet.skeleton, N, k + stp / steps * dt, 1 / steps * dt, params, pop, birthrate)
      N = N45[0]
    }//console.log(k, N)

    timetemp = k
    k += dt
    if (k > dataStartTime) {//console.log(k)
      k = timetemp;// console.log(k)
      dt = dataStartTime - timetemp ;//console.log(k + dt)
      N = Npre
    }
    if (k >= dataStartTime) {  
      k = timetemp + dt
      flag = 1
      arr.push(N[4])
      // states.push([k , N[0], N[1], N[3], N[2], N[4]])
    }
    //console.log(k, N)
  }

  count = 0
  while (k < dataEndTime) {
    
    if (Number(dataCases[count + 1][0]) !== "undefined") {
      dt = Number(dataCases[count + 1][0]) - Number(dataCases[count][0])
    } else {
      dt = deltaT
    }
    
    N[4] = 0
    for (let stp = 0; stp < steps; stp++) { 
      pop = interpolPopulation(k + stp / steps * dt)
      birthrate = interpolBirth(k + stp / steps * dt)
      N45 = mathLib.odeMethod('rkf45', snippet.skeleton, N, k + stp / steps * dt, 1 / steps * dt, params, pop, birthrate)
      N = N45[0]
      H = N[4]
    }//console.log(N)
    k += dt
    count++
    arr.push(H)
    // states.push([k , N[0], N[1], N[3], N[2], H])
  }
  // const createCsvWriter = require('csv-writer').createArrayCsvWriter;
  // const csvWriter = createCsvWriter({
  //   header: ['time', 'S', 'E', 'I', 'R', 'H'],
  //   path: './traj.csv'
  // })   
  // csvWriter.writeRecords(states)
  //   .then(() => {
  //   console.log('...Done')
  // })
  return arr
}

console.log(traj_match (interpolPopulation, interpolBirth, dataCases, params, times, index) )
