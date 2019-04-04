// Note: GetRNGstate(), PutRNGstate()
resample = {}
resample.nosort_resamp  = function (nw, w, np, p, offset) {
  for (j = 1; j < nw; j++) {
   w[j] += w[j-1]
 }
  if (w[nw - 1] <= 0) {
    throw "in 'systematic_resampling': non-positive sum of weight"
  }
  var du = w[nw - 1] / np
  var u = -du * Math.random()// unif_rand()

  for (j = 0, i = 0; j < np; j++) {
    u += du
    while ((u > w[i]) && (i < nw - 1)) i++;//looking for the low weight
    p[j] = i
  }
  if (offset){// add offset if needed
    for (j = 0; j < np; j++) {
      p[j] += offset
    }
  }
}


resample.systematic_resampling = function (weights) {
  let n, perm
  n = weights.length
  perm = n
  // weights = AS_NUMERIC(weights)
  // GetRNGstate();
  nosort_resamp(n, weights, n, perm, 1)
  // PutRNGstate();
  return(perm)
}
module.exports = resample
