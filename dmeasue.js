// dear emacs, please treat this as -*- C++ -*-

// #include <R.h>
// #include <Rmath.h>
// #include <Rdefines.h>
// #include <Rinternals.h>
// #include <R_ext/Rdynload.h>
// #include <R_ext/Arith.h>

// #include "pomp_internal.h"
pompfunmode = require('./pomp_internal.js')

dmeas_args = function (args, Onames, Snames, Pnames, Cnames, log) {
  var variable
  var v

  // we construct the call from end to beginning
  // 'log', covariableiates, parameter, states, observables, then time

  // 'log' is a needed argument
  PROTECT(args = LCONS(AS_LOGICAL(log),args)); nprotect++;
  SET_TAG(args,install("log"));

  // Covariableiates
  for (v = LENGTH(Cnames)-1; v >= 0; v--) {
    PROTECT(variable = NEW_NUMERIC(1)); nprotect++;
    PROTECT(args = LCONS(variable,args)); nprotect++;
    SET_TAG(args,install(CHAR(STRING_ELT(Cnames,v))));
  }

  // Parameters
  for (v = LENGTH(Pnames)-1; v >= 0; v--) {
    PROTECT(variable = NEW_NUMERIC(1)); nprotect++;
    PROTECT(args = LCONS(variable,args)); nprotect++;
    SET_TAG(args,install(CHAR(STRING_ELT(Pnames,v))));
  }

  // Latent state variableiables
  for (v = LENGTH(Snames)-1; v >= 0; v--) {
    PROTECT(variable = NEW_NUMERIC(1)); nprotect++;
    PROTECT(args = LCONS(variable,args)); nprotect++;
    SET_TAG(args,install(CHAR(STRING_ELT(Snames,v))));
  }

  // Observables
  for (v = LENGTH(Onames)-1; v >= 0; v--) {
    PROTECT(variable = NEW_NUMERIC(1)); nprotect++;
    PROTECT(args = LCONS(variable,args)); nprotect++;
    SET_TAG(args,install(CHAR(STRING_ELT(Onames,v))));
  }

  // Time
  PROTECT(variable = NEW_NUMERIC(1)); nprotect++;
  PROTECT(args = LCONS(variable,args)); nprotect++;
  SET_TAG(args,install("t"));

  UNPROTECT(nprotect);
  return args;

}

eval_call = function (fn, args,t, y, nobs, x, nvariable, p, npar, c, ncov) {

  var variable = args, ans
  var v

  *(REAL(CAR(variable))) = *t; variable = CDR(variable);
  for (v = 0; v < nobs; v++, y++, variable=CDR(variable)) *(REAL(CAR(variable))) = *y;
  for (v = 0; v < nvariable; v++, x++, variable=CDR(variable)) *(REAL(CAR(variable))) = *x;
  for (v = 0; v < npar; v++, p++, variable=CDR(variable)) *(REAL(CAR(variable))) = *p;
  for (v = 0; v < ncov; v++, c++, variable=CDR(variable)) *(REAL(CAR(variable))) = *c;

  PROTECT(ans = eval(LCONS(fn,args),CLOENV(fn)));

  UNPROTECT(1);
  return ans;

}

ret_array = function (nreps, ntimes) {
  var dim = [nreps, ntimes]
  const dimnm = ["rep","time"]
  var F
  F = mypomp_defines.makearray(2,dim)
  // fixdimnames(F,dimnm,2);
  // UNPROTECT(1);
  return F;
}

do_dmeasure = function (object, y, x, times, params, log, gnsi) {
  // int nprotect = 0;
  var mode = undefined
  var ntimes, nvariables, npars, ncovariables, nreps, nrepsx, nrepsp, nobs
  var Snames, Pnames, Cnames, Onames
  var cvec, pompfun
  var fn, args, as
  var F
  var dim
  lookup_table_t covariableiate_table;
  var cov

  ntimes = times.length
  if (ntimes < 1) {
    throw "length('times') = 0, no work to do."
  } 

  y = mypomp_defines.as_matrix(y)
  dim = mypomp_defines.GET_DIM(y)
  nobs = dim[0]

  if (ntimes != dim[1]) {
    throw "length of 'times' and 2nd dimension of 'y' do not agree."
  }

  x = mypomp_defines.as_state_array(x)
  dim = mypomp_defines.GET_DIM(x)
  nvariables = dim[0]; nrepsx = dim[1]

  if (ntimes !== dim[2])
    throw "length of 'times' and 3rd dimension of 'x' do not agree."

  params = mypomp_defines.as_matrix(params)
  dim = mypomp_defines.GET_DIM(params)
  npars = dim[0]; nrepsp = dim[1];

  nreps = (nrepsp > nrepsx) ? nrepsp : nrepsx;

  if ((nreps % nrepsp !== 0) || (nreps % nrepsx !== 0))
    throw "larger number of replicates is not a multiple of smaller."

  // PROTECT(Onames = GET_ROWNAMES(GET_DIMNAMES(y))); nprotect++;
  // PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
  // PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  // PROTECT(Cnames = get_covariableiate_names(GET_SLOT(object,install("covariable")))); nprotect++;




  // set up the covariableiate table
  covariableiate_table = make_covariableiate_table(GET_SLOT(object,install("covariable")),&ncovariables);
  PROTECT(cvec = NEW_NUMERIC(ncovariables)); nprotect++;
  cov = REAL(cvec);

  // extract the user-defined function
  PROTECT(pompfun = GET_SLOT(object,install("dmeasure"))); nprotect++;
  PROTECT(fn = pomp_fun_handler(pompfun,gnsi,&mode,Snames,Pnames,Onames,Cnames)); nprotect++;

  // extract 'userdata' as pairlist
  PROTECT(args = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;

  // create array to store results
  PROTECT(F = ret_array(nreps,ntimes)); nprotect++;

  switch (mode) {

  case Rfun: {
    double *ys = REAL(y), *xs = REAL(x), *ps = REAL(params), *time = REAL(times);
    double *ft = REAL(F);
    int j, k;

    // build argument list
    PROTECT(args = dmeas_args(args,Onames,Snames,Pnames,Cnames,log)); nprotect++;

    for (k = 0; k < ntimes; k++, time++, ys += nobs) { // loop over times

      R_CheckUserInterrupt(); // check for user interrupt

      table_lookup(&covariableiate_table,*time,cov); // interpolate the covariableiates

      for (j = 0; j < nreps; j++, ft++) { // loop over replicates

        // evaluate the call
        PROTECT(
          ans = eval_call(
            fn,args,
            time,
            ys,nobs,
            xs+nvariables*((j%nrepsx)+nrepsx*k),nvariables,
            ps+npars*(j%nrepsp),npars,
            cov,ncovariables
          )
        );

        if (k == 0 && j == 0 && LENGTH(ans) != 1)
          errorcall(R_NilValue,"user 'dmeasure' returns a vector of length %d when it should return a scalar.",LENGTH(ans));

        *ft = *(REAL(AS_NUMERIC(ans)));

        UNPROTECT(1);

      }
    }
  }

    break;

  case native: case regNative: {
    int *oidx, *sidx, *pidx, *cidx;
    int give_log;
    pomp_measure_model_density *ff = NULL;
    double *yp = REAL(y), *xs = REAL(x), *ps = REAL(params), *time = REAL(times);
    double *ft = REAL(F);
    double *xp, *pp;
    int j, k;

    // extract state, parameter, covariableiate, observable indices
    sidx = INTEGER(GET_SLOT(pompfun,install("stateindex")));
    pidx = INTEGER(GET_SLOT(pompfun,install("paramindex")));
    oidx = INTEGER(GET_SLOT(pompfun,install("obsindex")));
    cidx = INTEGER(GET_SLOT(pompfun,install("covariableindex")));

    give_log = *(INTEGER(AS_INTEGER(log)));

    // address of native routine
    *((void **) (&ff)) = R_ExternalPtrAddr(fn);

    set_pomp_userdata(args);

    for (k = 0; k < ntimes; k++, time++, yp += nobs) { // loop over times

      R_CheckUserInterrupt(); // check for user interrupt

      // interpolate the covariable functions for the covariableiates
      table_lookup(&covariableiate_table,*time,cov);

      for (j = 0; j < nreps; j++, ft++) { // loop over replicates

        xp = &xs[nvariables*((j%nrepsx)+nrepsx*k)];
        pp = &ps[npars*(j%nrepsp)];

        (*ff)(ft,yp,xp,pp,give_log,oidx,sidx,pidx,cidx,cov,*time);

      }
    }

    unset_pomp_userdata();

  }

    break;

  default: {
    double *ft = REAL(F);
    int j, k;

    for (k = 0; k < ntimes; k++) { // loop over times
      for (j = 0; j < nreps; j++, ft++) { // loop over replicates
        *ft = R_NaReal;
      }
    }

    warningcall(R_NilValue,"'dmeasure' unspecified: likelihood undefined.");

  }

  }

  UNPROTECT(nprotect);
  return F;
}

/////////////////////////////////////////
 // set_tag = function(x = NULL, overwrite = TRUE) {
 //      "set a tag for a file, your tag need to be a list or vector"
 //      if (is.null(x)) stop("please provided your tags")
 //      if (is.character(x)) x <- as.list(x)
 //      if (overwrite) {
 //        auth$api(
 //          path = paste0("files/", id, "/tags"),
 //          method = "PUT",
 //          body = x, ...
 //        )
 //        tags <<- x
 //      } else {
 //        .tags <- tag()
 //        .tags <- c(.tags, x)
 //        auth$api(
 //          path = paste0("files/", id, "/tags"),
 //          method = "PUT",
 //          body = .tags, ...
 //        )
 //        tags <<- .tags
 //      }

 //      tags
 //    },