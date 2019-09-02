const libR = require('lib-r-math.js')
const {
 Binomial,
 rng: { KnuthTAOCP2002 }
} = libR
const kn = new KnuthTAOCP2002(1234)
const { rnorm, rbinom } = Binomial(kn)

let reulermultinom = function (m = 1, size, rateAdd, dt, transAdd, rate, trans) {
    var p = 0
    var j, k
    if ((size < 0) || (dt < 0) || (Math.floor(size + 0.5) !== size)) {
      for (k = 0; k < m; k++) trans[k + transAdd] = NaN
      return 0
    }
    for (k = 0; k < m; k++) {
      if (rate[k + rateAdd] < 0.0) {
        for (j = 0; j < m; j++) trans[j + transAdd] = NaN
        return 0
      }
      p += rate[k + rateAdd]// total event rate
    }
    if (p > 0) {
       size = rbinom(1, size, 1 - Math.exp(-p * dt)) // total number of events
      if (!(isFinite(size)))
        throw 'result of binomial draw is not finite.'
      m -= 1
      for (k = 0; k < m; k++) {
        if (rate[k + rateAdd] > p) p = rate[k + rateAdd]
        trans[k + transAdd] = ((size > 0) && (p > 0)) ? rbinom(1, size, rate[k + rateAdd] / p) : 0
        if (!(isFinite(size) && isFinite(p) && isFinite(rate[k + rateAdd]) && isFinite(trans[k + transAdd]))) {
          throw 'result of binomial draw is not finite.'
        }
        size -= trans[k + transAdd]
        p -= rate[k + rateAdd]
      }
      trans[m + transAdd] = size
    } else {
      for (k = 0; k < m; k++) trans[k + transAdd] = 0
    }
   }    

   module.exports = reulermultinom;