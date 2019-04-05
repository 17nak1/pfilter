// Nazila: needs GetRNGstate(), PutRNGstate()

// #include "pomp_internal.h"
// #include <Rdefines.h>

// examines weights for filtering failure.
// computes conditional log likelihood and effective sample size.
// computes (if desired) prediction mean, prediction variance, filtering mean.
// it is assumed that ncol(x) == ncol(params).
// weights are used in filtering mean computation.
// if length(weights) == 1, an unweighted average is computed.
// tracks ancestry of particles if desired.
// returns all of the above in a named list.
resample = require('./resample.js')
mypomp_defines = require('./mypomp_defines')

// weights from dmeasure with logScale = 0
pfilterComputations = function (x, params, Np, predmean = 0, predvar = 0, filtmean = 0,  trackancestry = 0,  doparRS = 0,
  weights, tol) {

  // You use null to handle NULL in Rcpp.
  var pm = null, pv = null, fm = null, anc = null
  var ess, fail, loglik
  var newstates = null, newparams = null
  var retval, retvalnames
  const char = ['variable','rep']
  var xpm = 0, xpv = 0, xfm = 0, xw = 0, xx = 0, xp = 0
  var xanc = 0
  var dimX, dimP, newdim, Xnames, Pnames
  var dim, np // int *dim, np;
  var nvars, npars = 0, nreps, nlost
  var do_pm, do_pv, do_fm, do_ta, do_pr, all_fail = 0
  var sum = 0, sumsq = 0, vsq, ws, w, toler
  var returnValues
  var j, k

  dimX = mypomp_defines.GET_DIM (x)
  dim = dimX
  nvars = dim[0]; nreps = dim[1];console.log('nvars',nvars,'nreps',nreps)
  xx = x // xx = REAL(x);
  // PROTECT(Xnames = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;

  params = mypomp_defines.as_matrix(params);console.log('params',params)
  dimP =  mypomp_defines.GET_DIM (params);console.log('paramsdim',dimP)
  dim = dimP
  npars = dim[0]
  if (nreps % dim[1] !== 0)
    throw "ncol('states') should be a multiple of ncol('params')"
  // PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;

  np = Np// np = *(INTEGER(AS_INTEGER(Np))); // number of particles to resample
  
  // do_pm = *(LOGICAL(AS_LOGICAL(predmean))); // calculate prediction means?
  // do_pv = *(LOGICAL(AS_LOGICAL(predvar)));  // calculate prediction variances?
  // do_fm = *(LOGICAL(AS_LOGICAL(filtmean))); // calculate filtering means?
  // do_ta = *(LOGICAL(AS_LOGICAL(trackancestry))); // track ancestry?
  // // Do we need to do parameter resampling?
  // do_pr = *(LOGICAL(AS_LOGICAL(doparRS)));
  do_pm = predmean
  do_pv = predvar
  do_fm = filtmean
  do_ta = trackancestry
  do_pr = doparRS

  if (do_pr) {
    if (dim[1] !== nreps)
     throw "ncol('states') should be equal to ncol('params')"// # nocov
  }

  // PROTECT(ess = NEW_NUMERIC(1)); nprotect++; // effective sample size
  // PROTECT(loglik = NEW_NUMERIC(1)); nprotect++; // log likelihood
  // PROTECT(fail = NEW_LOGICAL(1)); nprotect++; // particle failure?
  ess = []
  loglik = []
  fail = []

  xw = weights
  toler = tol  // failure tolerance
  // check the weights and compute sum and sum of squares
  var  w = 0, ws = 0, nlost = 0
  for ( i = 0; i < nreps; i++) {
    if (xw[k] > toler) {
      w += xw[k]
      ws += xw[k] * xw[k]
    } else {      // this particle is lost
      xw[k] = 0
      nlost++
    }
  }

  if (nlost >= nreps) { 
    all_fail = 1 // all particles are lost
  }

  if (all_fail) {
    loglik = Math.log(toler) // minimum log-likelihood
    ess = 0  // zero effective sample size
  } else {
    loglik = Math.log(w / nreps)// mean of weights is likelihood
    ess = w * w / ws  // effective sample size
  }
  fail = all_fail
  
  
  if (do_pm || do_pv) {
    pm = new Array(nvars)
    xpm = pm
  }
  
  if (do_pv) {
    pv = new Array(nvars)
    xpv = pv
  }

  if (do_fm) {
    fm = new Array(nvars)
    xfm = fm
  }

  if (do_ta) {
    anc = new Array(np)
    xanc = anc
  }
  
  for (j = 0; j < nvars; j++) { // state variables

    // compute prediction mean
    if (do_pm || do_pv) {
      for (k = 0, sum = 0; k < nreps; k++) {
        sum += xx[j][k]
      } 
      sum /= nreps
      xpm[j] = sum
    }

    // compute prediction variance
    if (do_pv) {
      for (k = 0, sumsq = 0; k < nreps; k++) {
        vsq = xx[j][k] - sum
        sumsq += vsq * vsq
      }
      xpv[j] = sumsq / (nreps - 1)
    }

    //  compute filter mean
    if (do_fm) {
      if (all_fail) {   // unweighted average
        for (k = 0, ws = 0; k < nreps; k++) {
          ws += xx[j][k]
        }
        xfm[j] = ws/ nreps
      } else {      // weighted average
        for (k = 0, ws = 0; k < nreps; k++) {
          ws += xx[j][k] * xw[k]
        }
        xfm[j] = ws / w
      }
    }

  }
  // GetRNGstate();

// all_fail = 0
  // resample the particles unless we have filtering failure
  if (!all_fail) {
    var xdim = new Array(2)
    var sample = new Array(np)
    let ss = 0, st = 0, ps = 0, pt = 0
    // create storage for new states
    xdim[0] = nvars; xdim[1] = np
    newstates = mypomp_defines.makearray(2,xdim)
    // setrownames(newstates,Xnames,2);
    // fixdimnames(newstates,dimnm,2);
    ss = x
    st = newstates

    // create storage for new parameters
    if (do_pr) {
      xdim[0] = npars; xdim[1] = np;
      newparams = mypomp_defines.makearray(2,xdim)
      // setrownames(newparams,Pnames,2);
      // fixdimnames(newparams,dimnm,2);
      ps = params
      pt = newparams
    }
    // resample
    sample = resample.nosort_resamp(nreps,weights,np,sample,0)
    for (k = 0; k < np; k++) { // copy the particles
      for (j = 0; j < nvars; j++) {
        for (g = 0; g < nvars; g++) {
          st[sample[k]][g] = xx[sample[k]] [g]
        }
      }
      if (do_pr) {
        for (j = 0; j < npars; j++) {
          for (g = 0; g < npars; g++) {
            pt[sample[k]][g] = xp[sample[k]][g]
          }
        }
      }
      if (do_ta) {
        xanc[k] = sample[k] + 1
      } 
    }

  } else { // don't resample: just drop 3rd dimension in x prior to return

    newdim = new Array(2)
    dim = newdim
    dim[0] = nvars; dim[1] = nreps;
    x = mypomp_defines.SET_DIM(x,newdim);
    // console.log(dim, newdim)
    // setrownames(x,Xnames,2);
    // fixdimnames(x,dimnm,2);

    if (do_ta){
      for (k = 0; k < np; k++) {
       xanc[k] = k + 1
     }
    }
  }

  // PutRNGstate()
  let returnX = 0; returnPm = 0; returnPv = 0; returnFm = 0; returnAnc = 0
  if (all_fail) {
    returnX = x//;console.log(x)
  } else {
    returnX = newstates//;console.log("!",newstates)
  }

  if (do_pm) {
    returnPm = pm
  }
  if (do_pv) {
    returnPv = xpv
  }
  if (do_fm) {
    returnFm = fm
  }
  if (do_ta) {
    returnAnc = anc
  }

  returnValues = {fail:fail, loglik: loglik, ess: ess, states: returnX, params:0, pm: returnPm, pv: returnPv,fm: returnFm, ancestry:returnAnc }  
  console.log(returnValues)
  return returnValues
}

//Example
pfilterComputations([[1,1,1],[2,2,2]], [[1],[.1],[.1],[.2]], 3, 1,  1,  1,   1,  0,[.01,.9,.01], 1e-17)

