


pfilterComputation = function (x, params, p, runPredMean, runPredVar, runFilterMean, trackancestry, doparRS, weights, tol) {
  let allFail = 0, sum, sumsq, wss, loglik
  let st = [], retval = new Array(9)
  let predictionMean = [], predictionVariance = [], filterMean = []
  let sampleNum = Array.from(Array(Np).keys())
  // Check the weights and compute sum and sum of squares
  w = 0, ws = 0, nlost = 0
  for (let i = 0; i < Np; i++) {
    if (weights[i] > toler) {
      w += weights[i]
      ws += weights[i] ** 2
    } else {          // this particle is lost
      weights[i] = 0;
      nlost++
    }
  }
  if (nlost >= Np) { 
    allFail = 1 // all particles are lost
  }
  if (allFail) {
    loglik = Math.log(toler) // minimum log-likelihood
    ess = 0  // zero effective sample size
  } else {
    loglik = Math.log(w / Np)// mean of weights is likelihood
    ess = w * w / ws  // effective sample size
  }
  fail = allFail
  // Compute outputs
  for (let j = 0; j< nvars; j++) {
    // compute prediction mean
    if (doPredictionMean || doPredictionVariance) {
      sum = 0
      for (let nrow =0; nrow < Np; nrow++){
        sum += x[nrow][j]
      }
      sum /= Np
      predictionMean[j] = sum
    }  
    // compute prediction variance
    if (doPredictionVariance) {
      sumsq = 0
      for (let nrow = 0; nrow < Np; nrow++) {
        vsq = x[nrow][j] - sum
        sumsq += Math.pow(vsq, 2)
      }
      predictionVariance[j] = sumsq / (Np - 1) 
    }
    //  compute filter mean
    if (doFilterMean) {
      if (allFail) {   // unweighted average
        wss = 0
        for (let nrow =0; nrow < Np; nrow++){
          wss += x[nrow][j]
        } 
        filterMean[j] = wss / Np//;console.log(ws / Np)
      } else {      // weighted average
        wss = 0
        for (let nrow =0; nrow < Np; nrow++){
          wss += x[nrow][j] * weights[nrow]
        }
        filterMean[j] = wss / w
      }
    }
  }
  if (!allFail) {   // resample the particles unless we have filtering failure
    // resample
    mathLib.nosortResamp(Np, weights, Np, sampleNum, 0)
    for (np = 0; np < Np; np++) { // copy the particles
      st[np] = [].concat(x[sampleNum[np]])
    //   temp[np][nvars - 1] = 0
    }
  } else {  // don't resample: just drop 3rd dimension in x prior to return
    // if (do_ta) {
    //   for (k = 0; k < np; k++) {
    //     xanc[k] = k+1;
    //   }
    // }
  }
  retval[0] = fail
  retval[1] = loglik
  retval[2] = 0 //ess
  if (allFail) {
    retval[3] = x
  } else {
    retval[3] = st
  }
  retval[4] = 0 // newparams
  retval[5] = predictionMean
  retval[6] = predictionVariance
  retval[7] = filterMean
  retval[8] = 0 //anc
  
  return retval
}