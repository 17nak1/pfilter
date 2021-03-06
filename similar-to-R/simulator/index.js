
/**
 *  @file        simulator.js
 *               Transitions between classes in the time interval (t1, t2).
 *               The function returns the new value of states at time t2 based on rprocess in modelSnippet.
 *
 *  Inputs       Np : Nomber of points
 *               temp1 : State's values at time t1
 *               temp : State's values at time t2
 *
 *  @autor       Nazila Akhavan, nazila@kingsds.network
 *  @date        March 2019
 */

let mathLib = require('./mathLib.js')

let simulator = function (Np, dt, params, t1, t2 ) {
  let currentTime, steps, del_t
  
  steps = mathLib.numEulerSteps(t1, t2, dt) // Total number of steps in the interval (t1, t2)
  del_t = (t2 - t1) / steps
  currentTime = t1
  for (let i = 0; i < steps; i++) { // steps in each time interval
    args = this.interpolator(currentTime)   
    args.t = currentTime;
    args.dt = del_t;
    for (let np = 0; np < Np; np++){ //calc for each particle
      Object.assign(this.temp[np],args) ;
      this.rprocess(this.temp[np])
    }
    currentTime += del_t
    if (i == steps - 2) { // penultimate step
      del_t = t2 - currentTime;
      currentTime =  t2 - del_t;
    }
  }
}


module.exports = simulator;

