let pfilter = function(input){
    if(input == undefined){
        input = {};
    }
    object = input.object || {};
    params = input.params || {};
    Np = input.Np || {};
    tol = input.tol || 1e-17;
    max_fail = input.max_fail || Infinity;
    pred_mean = input.pred_mean || false;
    pred_var = input.pred_var || false;
    filter_mean = input.filter_mean || false;
    filter_traj = input.filter_traj || false;
    save_states = input.save_states || false;
    save_params = input.save_params || false;
    verbose = input.verbose || false;
    
    toler= 1e-17;
    //times = [...[object.t0],...object.times]
    times = [...object.times]
    //ntimes = times.length - 1;
    ntimes = times.length;
    let START = new Date();
    // Initial states from modelSnippet
    init_x = object.initializer({...object.interpolator(object.t0),...params});
    nvars = Object.keys(init_x).length;

    //loglik = new Array(Np);
    nfail = 0;
    let pred_m, pred_v, filt_m;
    
    // Initiate Matrix for chosen outputs
    if (pred_mean) {
        pred_m = Array(ntimes).fill(null).map(() => Array(nvars))
    }
    if (pred_var) {
        pred_v = Array(ntimes).fill(null).map(() => Array(nvars))
    }
    if (filter_mean) {
        filt_m = Array(ntimes).fill(null).map(() => Array(nvars))
    }
      // First matrix of states at time t0; includes Np rows of repeated state 
    x = new Array(Np).fill(null).map(() => ({...init_x}))
   
    let deltaT = 14 / 365.25;
    let tdata = object.times[0];
    
    let particles = new Array(Np)
    let state =[]
    let sampleNum = Array.from(Array(Np).keys())
    let condLoglik = []
    let stateSaved = []
    let ws ,w , vsq, sum, sumsq, ess, loglik = 0, lik 
    let weights
    let modelCases, likvalue
  /**
   *  Time loops
   *  The first loop simulating from t0 to first time value in reported data.
   *  The second one is based on times in reported data and calculates weights.  
   *  Starting at time t_n we have results at time t_{n+1}. But in the second loop k2 use the 
   *  results from the first loop which is in tdata and start calculing up to timeLen.
   */
    dt= 1 / 365.25;
  for (k = object.t0; k < tdata; k += deltaT){
    k2 = k + deltaT
    if ( k2 > tdata) {
      deltaT = tdata - k
      k2 = k + deltaT
    }
    particles = object.simulator (object, Np, x, dt, params, k, k2)
    for (np = 0; np < Np; np++) { // copy the particles
        x[np] = {...particles[sampleNum[np]]}
    }
  }
  
  for (let timeCountData = 1; timeCountData < ntimes; timeCountData++){
    weights = []; normalWeights = []
    let k2 =  Number(times[timeCountData])
 
    particles = object.simulator (object, Np, x, dt, params, k, k2)
    
    // Weights are calculated based on liklihood at each time for each point.
    for (np = 0; np < Np; np++){ 
      if (save_states) {
        stateSaved.push([k2, ...particles[np]]) //[S,E,I,R,H])
      }
      modelCases = Number(object.data[timeCountData])
      likvalue = object.dmeasure({...params, ...particles[np], ...{modelCases:modelCases, giveLog: 0}})
      weights.push(likvalue)           
    }
    
      w = 0, ws = 0, nlost = 0
      for (let i = 0; i < Np; i++) {
        if (weights[i] > toler) {
          w += weights[i]
          ws += weights[i] ** 2
        } else { // this particle is lost
          weights[i] = 0;
          nfail++
        }
      }
      // if (nlost > maxFail) {
      //   throw 'execution terminated. The number of filtering failures exceeds the maximum number of filtering failures allowed. '
      // }
      if (nfail >= Np) { 
        allFail = 1 // all particles are lost
      } else {
        allFail = 0
      }
      
      if (allFail) {
        lik = Math.log(toler)           // minimum log-likelihood
        ess = 0                         // zero effective sample size
      } else {
        ess = w * w / ws               // effective sample size
        lik = Math.log(w / Np)         // mean of weights is likelihood
      }
      condLoglik[timeCountData] = [timeCountData + 1, lik]
      // the total conditional logliklihood in the time process is loglik
      loglik += lik
      
      // Compute outputs
    //   for (let j = 0; j< nvars; j++) {
    //     // compute prediction mean
    //     if (doPredictionMean || doPredictionVariance) {
    //       sum = 0
    //       nlost = 0
    //       for (let nrow =0; nrow < Np; nrow++){
    //         if (particles[nrow][j]) {
    //           sum += particles[nrow][j]
    //         } else {
    //           nlost++
    //         }
    //       }
    //       sum /= Np
    //       predictionMean[timeCountData][j] = sum
    //     }  
    //     // compute prediction variance
    //     if (doPredictionVariance) {
    //       sumsq = 0
    //       for (let nrow = 0; nrow < Np; nrow++){

    //         if (particles[nrow][j]) {
    //           vsq = particles[nrow][j] - sum
    //           sumsq += Math.pow(vsq, 2)
    //         }
    //       }
    //       predictionVariance[timeCountData][j] = sumsq / (Np - 1) 
    //     }
    //     //  compute filter mean
    //     if (runFilterMean) {
    //       if (allFail) { //  unweighted average
    //         ws = 0
    //         for (let nrow =0; nrow < Np; nrow++){
    //           if (particles[nrow][j]) {
    //             ws += particles[nrow][j]
    //           }
    //         } 
    //         filterMean[timeCountData][j] = ws / Np
    //       } else { //  weighted average
    //         ws = 0
    //         for (let nrow =0; nrow < Np; nrow++){
    //           if (particles[nrow][j]) {
    //             ws += particles[nrow][j] * weights[nrow]
    //           }
    //         }
    //         filterMean[timeCountData][j] = ws / w
    //       }
    //     }
    //   }


      if (!allFail) {
        object.nosortResamp(Np, weights, Np, sampleNum, 0)
        for (np = 0; np < Np; np++) { // copy the particles
          temp[np] = {...particles[sampleNum[np]]}
          temp[np].H = 0
        }
      } else {
        for (np = 0; np < Np; np++) { // copy the particles
          temp[np] = [].concat(particles[np])
          temp[np][nvars - 1] = 0
        }
      }

      k = k2
  }//endTime

    result = {};
    
    //result.predictionMean = predictionMean
    console.log(loglik)
    result.predictionMean = []
    result.predictionVariance = []

    console.log('running time:',(new Date() - START) /1000)
    
    return result;
}

module.exports = pfilter;