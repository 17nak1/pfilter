const pfilter = require("./pfilter")


rproc = function(){}
rmeas = function(){}
dmeas = function(){}
initz = function(pomp, initialState){
  let State0 = {}
  let m = pomp.population(pomp.t0) / (initialState.S + initialState.E + initialState.R + initialState.I);
  State0.S = Math.round(m * initialState.S);
  State0.E = Math.round(m * initialState.E);
  State0.I = Math.round(m * initialState.I);
  State0.R = Math.round(m * initialState.R);
  State0.H = 0;
  return State0
}
statenames = ["S","E","I","R","H"]
paramnames = ["R0","amplitude","gamma","mu","sigma","rho","psi", "S", "E", "I", "R"]
euler.sim = function(){}

let pomp = {
data :  London_BiData,
times:  London_BiData.times,
t0: 1940,
rprocess :  euler.sim(snippet.rproc,{deltat: 1/365.25}),
rmeasure: snippet.rmeas,
covar: London_covar,
tcovar: London_covar.time,
dmeasure: snippet.dmeas,
zeronames: zeronames,
initializer: snippet.initz,
toEstimationScale: toEst,//I am not sure how this works
fromEstimationScale: fromEst,//I am not sure how this works
statenames: statenames,
paramnames: paramnames
}

for (let i = 0; i < pomp.covar.length; i++) {
  d1.push([Number(pomp.covar[i][0]), Number(pomp.covar[i][1])])
  d2.push([Number(pomp.covar[i][0]), Number(pomp.covar[i][2])])
}

// Define variables
pomp.population = mathLib.interpolator(d1)
pomp.birth = mathLib.interpolator(d2)

Object.fromEntries = arr => Object.assign({}, ...Array.from(arr, (k) => ({[k]: 0}) ));
let params = Object.fromEntries([...pomp.paramnames, ...pomp.statenames]);
//Now instead of pfilter.run we have
pfilterInternal = function (pomp, params, Np,
  tol, maxFail,
  predMean = false,
  predVar = false,
  filterMean = false,
  filterTraj = false,
  cooling, coolingm,
  verbose = false,
  saveStates = false,
  saveParams = false) 
  {
  let onePar = false;
  if (Object.keys(params).length === 0) throw new Error('Parameters are not specified')
  
  let times = [...pomp.time , pomp.t0];
  let ntimes = times.length - 1;
  
  if (typeof Np === "function" || Np === undefined || Np <= 0) {
    throw new Error(`Np Should be a positive number rather than${Np}`)
  }
  if (Object.values(params).every(function(x){ return !Array.isArray(x)})) {
    onePar = true;
    pomp.coef = params;
  }
  let initx = initState(pomp, params, Np);
  let nvars = initx[0].length;
  let x = initx;
  // set up storage for saving samples from filtering distributions
  if (saveStates || filterTraj) {
    // xparticles =  setNames(vector(mode="list",length=ntimes),time(object))
  }
  if (saveParams) {
    // pparticles = setNames(vector(mode="list",length=ntimes),time(object))
  } else {
    pparticles = {};
  }
  if (filterTraj) {
    pedigree = {};//vector(mode="list",length=ntimes+1)
  }

  let loglik = new Array(ntimes);
  let effSampleSize = new Array(ntimes).fill(0);
  let nfail = 0;

  // set up storage for prediction means, variances, etc.
  let predm
  if (predMean) {
    predm = new Array(ntimes).fill(Array(nvars).fill(0));
  } else {
    predm = [];
  }

  let predv
  if (predVar) {
    predv = new Array(ntimes).fill(Array(nvars).fill(0));
  } else {
    predv = [];
  }
  
  let filtm
  if (filterMean) {
    filtm  = new Array(ntimes).fill(Array(nvars).fill(0));
  } else {
    filtm = [];
  }

  if (filterTraj) {
    throw new Error ('filterTraj is not implemented')
  }
  // main loop
  let X
  for (nt = 0; nt < ntimes.length; nt++) {
    try {
      X = rprocess(pomp, x, [times[nt],times[nt + 1]], params, 1)//How to call rprocess???
    } catch (error) {
      throw new Error(`Process simulation error with${error}`)
    }
  }


  
}

/** Scales initial states using pomp.initializer()
  *  @return {Array} Each element is an array of states ordered based on pomp.statenames.
*/
initState = function (pomp, params, Np) {
  let initialState = {};
  let temp =  pomp.statenames;
  for(let i = 0; i < temp.length; i++) {
    if (params.hasOwnProperty(temp[i])) {
      initialState[temp[i]] = params[temp[i]];
    }
  }
  
  let initVector =  pomp.initializer(pomp, initialState);
  for(let i = 0; i < temp.length; i++) {
    if (initVector.hasOwnProperty(temp[i])) {
      temp[i] = initVector[temp[i]];
    }
  }
  return new Array(Np).fill(temp);
}


