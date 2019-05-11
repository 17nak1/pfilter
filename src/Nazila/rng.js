var libUnif = require('lib-r-math.js');

const {
    rng: { MersenneTwister },
    rng: { normal: { Inversion } }
} = libUnif
var Uunif = new MersenneTwister(0)

var rngUnif = {}
   
rngUnif.unif_rand = function(){
  return Math.random()
}

 module.exports = rngUnif