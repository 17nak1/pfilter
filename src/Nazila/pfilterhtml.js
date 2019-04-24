/*
 * Particle filter
 * A plain vanilla sequential Monte Carlo (particle filter) algorithm.
 * Resampling is performed at each observation.
 * @param Np the number of particles to use.
 * @param tol positive numeric scalar.
 * @param max.fail integer; the maximum number of filtering failures allowed.
 * @param pred.mean logical; if \code{TRUE}, the prediction means are calculated for the state variables and parameters.
 * @param pred.var logical; if \code{TRUE}, the prediction variances are calculated for the state variables and parameters.
 * @param filter.mean logical; if \code{TRUE}, the filtering means are calculated for the state variables and parameters.
 * @param filter.traj logical; if \code{TRUE}, a filtered trajectory is returned for the state variables and parameters.
 * @param save.states logical.
 *  If \code{save.states=TRUE}, the state-vector for each particle at each time is saved.
 * @references
    pfilter.R by Aaron A. King. https://github.com/kingaa/pomp/blob/abc552b335319bd36b21d3313305c89e97f48f68/R/pfilter.R
    M. S. Arulampalam, S. Maskell, N. Gordon, & T. Clapp.
    A Tutorial on Particle Filters for Online Nonlinear, Non-Gaussian Bayesian Tracking.
    IEEE Trans. Sig. Proc. 50:174--188, 2002.
 */
//conditional log liklihood(time) =log(sum(w_i, i in 0:Np) / Np)

let pfilter = require ('./pfilter.js')

// main function
function pfilterCalculation (input) {//filter.traj , save.params
    
  pfilter.run(input)

  activateDownload ()
  console.log(pfilter.prediction)
  return pfilter.predictionMean

//    if (input.runPredMean) {
//     return predictionMean
//   }

//   if (input.runPredVar) {
//     return predictionVariance
//   }

//   if (input.runFilterMean) {
//     return filterMean
//   }
  
//   if (input.runSaveStates){
//     return stateSaved
//   }

}

module.exports = {
  pfilterCalculation
}


    
      


