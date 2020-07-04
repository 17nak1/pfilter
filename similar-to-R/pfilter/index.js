let mathLib = require('./../simulator/mathLib')
let pfilter = function(args){
    if(args == undefined){
        args = {};
    }
    let params = args.params || {};
    let Np = args.Np || {};
    tol = args.tol || 1e-17;
    max_fail = args.max_fail || Infinity;
    pred_mean = args.pred_mean || false;
    pred_var = args.pred_var || false;
    filter_mean = args.filter_mean || false;
    filter_traj = args.filter_traj || false;
    save_states = args.save_states || false;
    save_params = args.save_params || false;
    verbose = args.verbose || false;

    let START = new Date();
    

    let maxFail = this.maxFail
    let toler = 1e-17
    let dt =  1 / 365.25
    let t0 = this.t0
    let deltaT = 14 / 365.25
    let tdata = this.times[0]
    console.log("Np", Np)
    let nvars = this.statenames.length // number of states SEIRH
    let doPredictionVariance = 1, doPredictionMean = 1, runFilterMean = 0 , allFail = 0
    let timeLen = this.data.length 
    let nlost = 0
    let particles = new Array(Np).fill(null).map(() => Array(5))
    let state =[]
    let sampleNum = Array.from(Array(Np).keys())
    let condLoglik = []
    let stateSaved = []
    let ws ,w , vsq, sum, sumsq, ess, loglik = 0, lik 
    let predictionMean, predictionVariance, filterMean
    let weights
    let modelCases, likvalue

    // Initiate Matrix for chosen outputs
    if (doPredictionMean) {
        predictionMean = Array(timeLen).fill(null).map(() => Array(nvars))
    }
    if (doPredictionVariance) {
        predictionVariance = Array(timeLen).fill(null).map(() => Array(nvars))
    }
    if (runFilterMean) {
        filterMean = Array(timeLen).fill(null).map(() => Array(nvars))
    }

    Object.assign(params, this.interpolator(t0))
    state.t = t0;
    state.dt = dt;
    // Initial states from modelSnippet
    state = this.initializer(params)

    Object.assign(state, params)
    // First matrix of states at time t0; includes Np rows of repeated state 
    this.temp = new Array(Np).fill(null).map(() => Object.assign({}, state))
    
    // Define t0 as the first time value 
    let k = t0

      /**
   *  Time loops
   *  The first loop simulating from t0 to first time value in reported data.
   *  The second one is based on times in reported data and calculates weights.  
   *  Starting at time t_n we have results at time t_{n+1}. But in the second loop k2 use the 
   *  results from the first loop which is in tdata and start calculing up to timeLen.
   */
    for (k = t0; k < tdata; k += deltaT){
        k2 = k + deltaT
        if ( k2 > tdata) {
        deltaT = tdata - k
        k2 = k + deltaT
        }
        this.simulator(Np, dt, params, k, k2)
        for (np = 0; np < Np; np++) { // copy the particles
        this.temp[np] = Object.assign({}, this.temp[sampleNum[np]])
        }
    }

    for (let timeCountData = 0; timeCountData < timeLen; timeCountData++){
        weights = []; normalWeights = []
        let k2 =  Number(this.times[timeCountData])
     
        this.simulator(Np, dt, params, k, k2)
        
        // Weights are calculated based on liklihood at each time for each point.
        for (np = 0; np < Np; np++){ 
        //   if (defaults.runSaveStates) {
        //     stateSaved.push([k2, ...particles[np]]) //[S,E,I,R,H])
        //   }
        this.temp[np].modelCases = Number(this.data[timeCountData])
        this.temp[np].giveLog = 0
          likvalue = this.dmeasure(this.temp[np])
          weights.push(likvalue)           
        }
        
          w = 0, ws = 0, nlost = 0
          for (let i = 0; i < Np; i++) {
            if (weights[i] > toler) {
              w += weights[i]
              ws += weights[i] ** 2
            } else { // this particle is lost
              weights[i] = 0;
              nlost++
            }
          }
          // if (nlost > maxFail) {
          //   throw 'execution terminated. The number of filtering failures exceeds the maximum number of filtering failures allowed. '
          // }
          if (nlost >= Np) { 
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
            mathLib.nosortResamp(Np, weights, Np, sampleNum, 0)
            for (np = 0; np < Np; np++) { // copy the particles
              if(np != sampleNum[np]){
                this.temp[np] = Object.assign({}, this.temp[sampleNum[np]])                
              }
              this.temp[np].H = 0
            }
          } else {
            for (np = 0; np < Np; np++) { // copy the particles
                if(np != sampleNum[np]){
                  this.temp[np] = Object.assign({},this.temp[sampleNum[np]])
                }
                this.temp[np].H = 0
            }
          }
    
          k = k2
      }//endTime
      
      console.log(loglik)
        
    console.log('running time:',(new Date() - START) /1000)
    
    // return result;
}

module.exports = pfilter;