const { initState } = require("./initState.js");
const { rprocessInternal } = require("./rprocessInternal.js");
const { dmeasureInternal } = require("./dmeasureInternal.js");
const { pfilter_computations } = require("./pfilterComputations.js");

exports.pfilter = function (object) {
  
  let params = object.params;
  let Np = object.Np;
  let tol = object.tol ? object.tol : 1e-17;
  let maxFail = object.maxFail ? object.maxFail :  Infinity;
  let predMean = object.predMean ? object.predMean :  false;
  let predVar = object.predVar ? object.predVar :  false;
  let filterMean = object.filterMean ? object.filterMean :  false;
  let filterTraj = object.filterTraj ? object.filterTraj :  false;
  let cooling = object.cooling;
  let coolingm = object.coolingm;
  let verbose = object.verbose ? object.verbose : false;
  let saveStates = object.saveStates ? object.saveStates : false;
  let saveParams = object.saveParams ? object.saveParams : false;
  object = object.object;

  if (params.length === 0)
    throw new Error("In pfilterInternal: params must be specified");
  let onePar = false;
  
  let times = [object.t0, ...object.times];
  let ntimes = times.length - 1;
  
  if (typeof Np === "function" || Np === undefined || Np <= 0) {
    throw new Error(`Number of particles should be a positive number. ${Np} is not translated`)
  }
  if (params.every(x => !Array.isArray(x))) { //there is only one parameter vector
    onePar = true;
    object.coef = params;
    params = [params]; //as.matrix(params)
  }

  let initx = initState(object, params, Np);
  let nvars = initx[0].length;
  let x = initx;
  let xparticles, pparticles, pedigree
    // set up storage for saving samples from filtering distributions
    if (saveStates || filterTraj) {
      xparticles =  new Array(ntimes);
    }
    if (saveParams) {
      pparticles = new Array(ntimes);;
    } else {
      pparticles = [];
    }
    if (filterTraj) {
      pedigree = new Array(ntimes + 1);
    }

  let loglik = new Array(ntimes);
  let effSampleSize = new Array(ntimes).fill(0);
  let nfail = 0;

  // set up storage for prediction means, variances, etc.
  let predm
  if (predMean) {
    predm = new Array(ntimes).fill(null).map(a => Array(nvars).fill(0));
  } else {
    predm = [];
  }

  let predv
  if (predVar) {
    predv = new Array(ntimes).fill(null).map(a => Array(nvars).fill(0));
  } else {
    predv = [];
  }
  
  let filtm
  if (filterMean) {
    filtm  = new Array(ntimes).fill(null).map(a => Array(nvars).fill(0));
  } else {
    filtm = [];
  }

  if (filterTraj) {
    throw new Error ('filterTraj is not implemented')
  }
  // main loop
  let X
  for (nt = 0; nt < ntimes; nt++) {
    try {
      X = rprocessInternal(object, x, [times[nt],times[nt + 1]], params, 1)
    } catch (error) {
      throw new Error(`In pfilterInternal: Process simulation error: ${error}`)
    }
  

    if (predVar) { // check for nonfinite state variables and parameters
      allFinite = X.every(a => a.every(x =>isFinite(x)))
      if (!allFinite) {  // state variables
        throw new Error("In pfilter.js: non-finite state variable(s): ");
      }
    }

    let weights = [];
    try {
      weights = dmeasureInternal(
        object,
        y = object.data[nt],
        X,
        times[nt + 1],
        params,
        log = false
      ); 
    } catch (error) {
      console.error(`In mif2.js: error in calculation of weights: ${error}`);
    }

    let allFinite = weights.map(w => isFinite(w)).reduce((a, b) => a & b, 1);
    if (!allFinite) {
      throw new Error("In dmeasure: weights returns non-finite value");
    }
    
    /** compute prediction mean, prediction variance, filtering mean,
      * effective sample size, log-likelihood
      * also do resampling if filtering has not failed
      */ 
    let xx;
    try {
      xx  = pfilter_computations(
        X,
        params,
        Np = Np,
        0, //rw_sd
        predMean,
        predVar ,
        filterMean ,
        trackancestry = filterTraj,
        onePar,
        weights = weights,
        tol = tol
      );
    } catch (error) {
      console.error(`particle-filter error: ${error}`) 
    } 

    let allFail = xx.fail;
    loglik[nt] = xx.loglik;
    effSampleSize[nt] = xx.ess;
    x = xx.states;
    params = xx.params;

    if (predMean)
      predm[nt] = xx.pm;

    if (predVar)
      predv[nt] = xx.pv;

    if (filterMean)
      filtm[nt] = xx.fm;

    if (filterTraj)
      pedigree[[nt]] = xx.ancestry;

    if (allFail) { // all particles are lost
      nfail = nfail + 1;
      if (nfail > maxFail)
        throw new Error("In mif2Pfilter: too many filtering failures")
    }

    if (saveStates || filterTraj) {
      xparticles[nt] = x;
    }

    if (saveParams) {
      pparticles[nt] = params;
    } 
  } //end of main loop 

  if (filterTraj) { // select a single trajectory
    throw new Error("filterTraj is not translated")
    // b = sample.int(n=length(weights),size=1L,replace=TRUE)
    // filt.t[,1L,ntimes+1] = xparticles[[ntimes]][,b]
    // for (nt in seq.int(from=ntimes-1,to=1L,by=-1L)) {
    //   b = pedigree[[nt+1]][b]
    //   filt.t[,1L,nt+1] = xparticles[[nt]][,b]
    // }
    // if (times[2L] > times[1L]) {
    //   b = pedigree[[1L]][b]
    //   filt.t[,1L,1L] = init.x[,b]
    // } else {
    //   filt.t = filt.t[,,-1L,drop=FALSE]
    // }
  }

  if (!saveStates) xparticles = [];

  if (nfail > 0) {
    console.log("warning! filtering failure occurred.");
  }

  if(predMean) predm.unshift(object.statenames);
  if(predVar) predv.unshift(object.statenames);
  if(filterMean) filtm.unshift(object.statenames);
  
  return {
    object,
    predMean: predm,
    predVar: predv,
    filterMean: filtm,
    filterTraj: false,//filt.t
    paramMatrix: [[]],
    effSamplesize: effSampleSize,
    condLoglik: loglik,
    savedStates: xparticles,
    savedParams: pparticles,
    Np: Np,
    tol: tol,
    nfail: nfail,
    loglik: loglik.reduce((a,b) => a + b, 0)
  }
  
}




