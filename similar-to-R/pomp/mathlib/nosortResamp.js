
let libUnif = require('lib-r-math.js');
const {
    R: { numberPrecision },
    rng: { MersenneTwister, timeseed }
} = libUnif
let U = new MersenneTwister(0) 

module.exports = function (nw, w, np, p, offset) {
    for (j = 1; j < nw; j++) {
     w[j] += w[j-1]
   }
    if (w[nw - 1] <= 0) {
      throw "in 'systematic_resampling': non-positive sum of weight"
    }
    let du = w[nw - 1] / np
    let u = -du * U.unif_rand()//Math.random()
  
    for (j = 0, i = 0; j < np; j++) {
      u += du
      while ((u > w[i]) && (i < nw - 1)) i++;//looking for the low weight
      p[j] = i
    }
    if (offset){// add offset if needed
      for (j = 0; j < np; j++) p[j] += offset
    }
  }