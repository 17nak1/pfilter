let libUnif = require('lib-r-math.js');

const {
    rng: { MersenneTwister },
    rng: { normal: { Inversion } }
} = libUnif
let Uunif = new MersenneTwister(0)

let rngUnif = {}
   
rngUnif.unif_rand = function(){
  return Math.random()
}

 module.exports = rngUnif