
let mathLib = require('./mathLib.js')
let snippet = require('./modelSnippet.js')

let simulator = {}
simulator.simulate = function (Np, temp2, dt, interpolPop, interpolBirth, params, t1,t2 ) {
  let st, steps, del_t, pop, birthrate
  // transitions between classes
      steps = mathLib.numEulerSteps(t1, t2, dt)
      temp = [].concat(temp2)
      del_t = (t2 - t1) / steps
    st = t1
    for (let stp = 0; stp < steps; stp++) { // steps in each time interval
      
      pop = interpolPop(st)
      birthrate = interpolBirth(st)
      for (let np = 0; np < Np; np++){ //calc for each particle
        temp[np] = simulateValue = snippet.rprocess(params, st, del_t, temp[np], pop, birthrate)
      }
      st += del_t
      if (stp == steps - 2) { // penultimate step
        del_t = t2 - st;
        st =  t2 - del_t;
      }
    }
  return temp
}

module.exports = simulator