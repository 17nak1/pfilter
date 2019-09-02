let mathLib = require('./mathLib.js')
let snippet = require('./modelSnippet.js')

simulator = {}
simulator.simulate = function (Np, temp, dt, interpolPop, interpolBirth, params, t1,t2 ) {
  var st, S, E, I, R, H, steps, simulateValue 
  // transitions between classes
      steps = mathLib.numEulerSteps(t1, t2, dt)

      del_t = (t2 - t1) / steps
    st = t1
    for (let stp = 0; stp < steps; stp++) { // steps in each time interval
      
      pop = interpolPop(st)
      birthrate = interpolBirth(st)
      for (let np = 0; np < Np; np++){ //calc for each particle
        trans = []
        S = temp[np][0]; E = temp[np][1]; I = temp[np][2]; R = temp[np][3]; H = temp[np][4]
        simulateValue = snippet.rprocess(params, st, del_t, [S,E,I,R,H], pop, birthrate)
        temp[np][0] = simulateValue[0]; temp[np][1] = simulateValue[1]; temp[np][2] = simulateValue[2]; temp[np][3] = simulateValue[3]; temp[np][4] = simulateValue[4]
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