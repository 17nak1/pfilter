




pfilterInternal = function (object, Np, tol, maxFail,
  predMean = FALSE, predVar = FALSE, filterMean = FALSE,
  filterTraj = FALSE, cooling, cooling_m, saveStates = FALSE,
  gnsi = TRUE, verbose = FALSE) {

  // verbose = as.logical(verbose)

  // object = pomp(object,...,verbose=verbose)

  if (typeof object.rprocess === 'undefined' || typeof object.dmeasure === 'undefined')
    throw 'rprocess / dmeasure are needed basic components.'

  // tol = Number(tol)
  // gnsi = as.logical(.gnsi)
  // predMean = as.logical(predMean)
  // predVar = as.logical(predVar)
  // filterMean = as.logical(filterMean)
  // filterTraj = as.logical(filterTraj)
  // saveStates = as.logical(saveStates)

  params = object.coef // coef(object)
  times = object.time  // time(object,t0=TRUE)
  ntimes = times.length-1

  if (Np === undefined) {
    throw 'Np must be specified.'
  // } else if (is.function(Np)) {
  //   Np = tryCatch(
  //     vapply(seq.int(from=0,to=ntimes,by=1),Np,numeric(1)),
  //     error = function (e) {
  //       pStop_("if ",sQuote("Np")," is a function, it must return ",
  //         "a single positive integer.")
  //     }
  //   )
  } else if (!Number.isInteger(Np) && !Np.length) {
    throw 'Np must be a number, a vector of numbers, or a function.'
  }

  if (Np.length == 1) {
    Np = new Array(ntimes + 1).fill(Np)
  } else if (Np.length !== (ntimes + 1)) {
    throw 'Np must have length 1 or length'+ (ntimes + 1)
  }

  if (!all(is.finite(Np)) || any(Np <= 0))
    if (Np.every(function(a) {return isFinite(a)}) || Np.every(function(a) {return a <=0}))

    throw 'Number of particles' + Np + ', must be a positive integer.'

  Np = Number(Math.floor(Np))

  if (tol.length != 1 || !isFinite(tol) || tol < 0)
    throw tol + ' should be a small positive number.'

  // pompLoad(object,verbose=verbose)
  // on.exit(pompUnload(object,verbose=verbose))

  initx = rinit(object,params=params,nsim=Np[0],.gnsi=gnsi)
  statenames = rownames(initx)
  nvars = nrow(initx)
  x = initx

  // ## set up storage for saving samples from filtering distributions
  // if (saveStates || filterTraj) {
  //   xparticles = setNames(vector(mode="list",length=ntimes),time(object))
  // }
  // if (filterTraj) {
  //   pedigree = vector(mode="list",length=ntimes+1)
  // }

  loglik = rep(NA,ntimes)
  effSampleSize = numeric(ntimes)
  nfail = 0

  // ## set up storage for prediction means, variances, etc.
  // if (predMean) {
  //   pred.m = array(data=numeric(1),dim=c(nvars,ntimes),
  //     dimnames=list(variable=statenames,time=time(object)))
  // } else {
  //   pred.m = array(data=numeric(0),dim=c(0,0))
  // }

  // if (predVar) {
  //   pred.v = array(data=numeric(1),dim=c(nvars,ntimes),
  //     dimnames=list(variable=statenames,time=time(object)))
  // } else {
  //   pred.v = array(data=numeric(0),dim=c(0,0))
  // }

  // if (filterMean) {
  //   filt.m = array(data=numeric(1),dim=c(nvars,ntimes),
  //     dimnames=list(variable=statenames,time=time(object)))
  // } else {
  //   filt.m = array(data=numeric(0),dim=c(0,0))
  // }

  // if (filterTraj) {
  //   filt.t = array(data=numeric(1),dim=c(nvars,1,ntimes+1),
  //     dimnames=list(variable=statenames,rep=1,time=times))
  // } else {
  //   filt.t = array(data=numeric(0),dim=c(0,0,0))
  // }

  for (nt = 0; nt < ntimes; nt ++) { // main loop

    // advance the state variables according to the process model
    X = rprocess(object,x0=x,t0=times[nt],times=times[nt+1],params=params,
      .gnsi=gnsi)

  //   if (predVar) { ## check for nonfinite state variables and parameters
  //     problem.indices = unique(which(!is.finite(X),arr.ind=TRUE)[,1L])
  //     nprob = length(problem.indices)
  //     if (nprob > 0)
  //       pStop_("non-finite state variable",ngettext(nprob,"","s"),": ",
  //         paste(rownames(X)[problem.indices],collapse=', '))
  //   }

  //   ## determine the weights
  //   weights = dmeasure(object,y=object@data[,nt,drop=FALSE],x=X,
  //     times=times[nt+1],params=params,log=FALSE,.gnsi=gnsi)

  //   if (!all(is.finite(weights))) {
  //     first = which(!is.finite(weights))[1L]
  //     datvals = object@data[,nt]
  //     weight = weights[first]
  //     states = X[,first,1L]
  //     msg = nonfinite_dmeasure_error(time=times[nt+1],lik=weight,
  //       datvals,states,params)
  //     pStop_(msg)
  //   }

  //   gnsi = FALSE

  //   ## compute prediction mean, prediction variance, filtering mean,
  //   ## effective sample size, log-likelihood
  //   ## also do resampling if filtering has not failed
  //   xx = .Call(P_pfilter_computations,x=X,params=params,Np=Np[nt+1],
  //     predmean=predMean,predvar=predVar,filtmean=filterMean,
  //     trackancestry=filterTraj,doparRS=FALSE,weights=weights,tol=tol)

  //   all.fail = xx$fail
  //   loglik[nt] = xx$loglik
  //   effSampleSize[nt] = xx$ess

  //   x = xx$states
  //   params = xx$params[,1L]

  //   if (predMean) pred.m[,nt] = xx$pm
  //   if (predVar) pred.v[,nt] = xx$pv
  //   if (filterMean) filt.m[,nt] = xx$fm
  //   if (filterTraj) pedigree[[nt]] = xx$ancestry

  //   if (all.fail) { ## all particles are lost
  //     nfail = nfail+1
  //     if (verbose) message("filtering failure at time t = ",times[nt+1])
  //     if (nfail>maxFail) pStop_("too many filtering failures")
  //   }

  //   if (saveStates || filterTraj) {
  //     xparticles[[nt]] = x
  //     dimnames(xparticles[[nt]]) = setNames(dimnames(xparticles[[nt]]),
  //       c("variable","rep"))
  //   }

  //   if (verbose && (nt%%5 == 0))
  //     cat("pfilter timestep",nt,"of",ntimes,"finished\n")

  // } ## end of main loop

  // if (filterTraj) { ## select a single trajectory
  //   b = sample.int(n=length(weights),size=1L,replace=TRUE)
  //   filt.t[,1L,ntimes+1] = xparticles[[ntimes]][,b]
  //   for (nt in seq.int(from=ntimes-1,to=1L,by=-1L)) {
  //     b = pedigree[[nt+1]][b]
  //     filt.t[,1L,nt+1] = xparticles[[nt]][,b]
  //   }
  //   if (times[2L] > times[1L]) {
  //     b = pedigree[[1L]][b]
  //     filt.t[,1L,1L] = initx[,b]
  //   } else {
  //     filt.t = filt.t[,,-1L,drop=FALSE]
  //   }
  // }

  // if (!saveStates) xparticles = list()

  // if (nfail>0)
  //   pWarn("pfilter",nfail," filtering failure",ngettext(nfail,"","s")," occurred.")

  // new(
  //   "pfilterd_pomp",
  //   object,
  //   predMean=pred.m,
  //   predVar=pred.v,
  //   filterMean=filt.m,
  //   filterTraj=filt.t,
  //   paramMatrix=array(data=numeric(0),dim=c(0,0)),
  //   effSampleSize=effSampleSize,
  //   cond.loglik=loglik,
  //   saved.states=xparticles,
  //   Np=as.integer(Np),
  //   tol=tol,
  //   nfail=as.integer(nfail),
  //   loglik=sum(loglik)
  // )
}

}

