const { initState } = require("./initState.js");
const { rprocessInternal } = require("./rprocessInternal.js");
const { dmeasureInternal } = require("./dmeasureInternal.js");
const { pfilter_computations } = require("./pfilterComputations.js");
const { statenames } = require("./modelSnippet.js");

exports.pfilter = function (args) {
  
  let params = args.params;
  let Np = args.Np;
  let tol = args.tol ? args.tol : 1e-17;
  let maxFail = args.maxFail ? args.maxFail :  Infinity;
  let predMean = args.predMean ? args.predMean :  false;
  let predVar = args.predVar ? args.predVar :  false;
  let filterMean = args.filterMean ? args.filterMean :  false;
  let filterTraj = args.filterTraj ? args.filterTraj :  false;
  let cooling = args.cooling;
  let coolingm = args.coolingm;
  let verbose = args.verbose ? args.verbose : false;
  let saveStates = args.saveStates ? args.saveStates : false;
  let saveParams = args.saveParams ? args.saveParams : false;
  object = args.object;

  if (Object.keys(params).length === 0)
    throw new Error("In pfilterInternal: params must be specified");
  let onePar = false;
  
  let times = [object.t0, ...object.times];
  let ntimes = times.length - 1;
  
  if (typeof Np === "function" || Np === undefined || Np <= 0) {
    throw new Error(`Number of particles should be a positive number. ${Np} is not translated`)
  }
  if (!Array.isArray(params) || params.every(x => !Array.isArray(x))) { //there is only one parameter vector
    onePar = true;
    object.coef = params;
    params = [params]; //as.matrix(params)
  }

  let initx = initState(object, params, Np);
  let nvars = Object.keys(initx[0]).length;
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
    predm = new Array(ntimes).fill(null).map( a => []);
  } else {
    predm = [];
  }

  let predv
  if (predVar) {
    predv = new Array(ntimes).fill({});
    predv = [];
  }
  
  let filtm
  if (filterMean) {
    filtm  = new Array(ntimes).fill({});
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
    } catch (e) {
      console.error(`In pfilterInternal: Process simulation error: ${e}`)
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

    if (predMean) {
      for(let i = 0; i <xx.pm.length; i++) {
        for (let j = 0; j < object.statenames.length; j++) {
          if(Object.keys(xx.pm)[i] === object.statenames[j]) {
            predm[nt][j] =xx.pm[object.statenames[j]];
          }
        }
      }
    }

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




