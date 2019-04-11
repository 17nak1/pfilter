pompInternal = function (data, times, t0, timename,
  rinit, rprocess, dprocess, rmeasure, dmeasure, skeleton, rprior, dprior,
  partrans, params, covar, accumvars, obsnames, statenames,
  paramnames, covarnames, PACKAGE, globals, cdir, cfile, shlibArgs,
  compile,/* .userdata, .solibs = list(), verbose = getOption("verbose", FALSE)*/) {

  // check times
  for (i = 1; i < times.length; i++) {
    if (!isFinite(times[i]) || (times[i - 1] > times[i])) {
      throw 'times must be a non-decreasing sequence of numbers.'
    }
  }
  
  // check t0
  if (!isFinite(t0) || t0 > times[0]) {
    throw 't0 must be a single number not greater than times[1].'
  }

  // if (missing(timename) || is.null(timename))
  //   timename = "time"
  // else
  //   timename = as.character(timename)

  // if (missing(.userdata)) .userdata = list()
  // added.userdata = list(...)
  // if (length(added.userdata)>0) {
  //   message("in ",sQuote("pomp"),": the unrecognized ",
  //     ngettext(length(added.userdata),"argument","arguments")," ",
  //     paste(sQuote(names(added.userdata)),collapse=","),
  //     ngettext(length(added.userdata)," is"," are"),
  //     " available for use by the POMP basic components."
  //   )
  //   .userdata[names(added.userdata)] = added.userdata
  // }

  // if (!is(rprocess,"rprocPlugin")) {
  //   pStop_(sQuote("rprocess"),
  //     " must be specified using one of the plugins:\n",
  //     sQuote("onestep"),", ",sQuote("discrete_time"),
  //     ", ",sQuote("euler"),", ",sQuote("gillespie"),
  //     ", or ",sQuote("gillespie_hl"),".")
  // }

  // if (!is(skeleton,"skelPlugin")) {
  //   pStop_(sQuote("skeleton")," must be specified using either ",
  //     sQuote("map")," or ",sQuote("vectorfield"),".")
  // }

  // if (!is(partrans,"partransPlugin")) {
  //   pStop_(sQuote("partrans")," must be specified using ",
  //     sQuote("parameter_trans"),".")
  // }

  
  // store the data as double-precision matrix
  storage.mode(data) = "double"
  if (is.null(obsnames)) obsnames = rownames(data)

  // check the parameters and force them to be double-precision
  // params = setNames(as.double(params),names(params))
  // if (length(params) > 0) {
  //   if (is.null(names(params)) || !is.numeric(params) ||
  //       !all(nzchar(names(params))))
  //     pStop_(sQuote("params")," must be a named numeric vector.")
  // }

  // if (is(rinit,"Csnippet") && is.null(statenames)) {
  //   pStop_("when ",sQuote("rinit")," is provided as a C snippet, ",
  //     "you must also provide ",sQuote("statenames"),".")
  // }

  // if (is(rmeasure,"Csnippet") && is.null(obsnames)) {
  //   pStop_("when ",sQuote("rmeasure")," is provided as a C snippet, ",
  //          "you must also provide ",sQuote("obsnames"),".")
  // }

  // check and arrange covariates
  // if (is.null(covar)) {
  //   covar = covariate_table()
  // } else if (!is(covar,"covartable")) {
  //   pStop_("bad option for ",sQuote("covar"),".")
  // }

  // if (is.null(covarnames)) covarnames = get_covariate_names(covar)

  // hitches = hitch(
  //   rinit=rinit,
  //   step.fn=rprocess@step.fn,
  //   rate.fn=rprocess@rate.fn,
  //   dprocess=dprocess,
  //   rmeasure=rmeasure,
  //   dmeasure=dmeasure,
  //   rprior=rprior,
  //   dprior=dprior,
  //   toEst=partrans@to,
  //   fromEst=partrans@from,
  //   skeleton=skeleton@skel.fn,
  //   templates=workhorse_templates,
  //   obsnames=obsnames,
  //   statenames=statenames,
  //   paramnames=paramnames,
  //   covarnames=covarnames,
  //   PACKAGE=PACKAGE,
  //   globals=globals,
  //   cfile=cfile,
  //   cdir=cdir,
  //   shlibArgs=shlibArgs,
  //   compile=compile,
  //   verbose=verbose
  // )

  // check to see if covariate times embrace the data times
  covar_time_warning(covar,times,t0,"pomp")

  var pomp = {
    data : data,
    times : times,
    t0 : t0,
    timename : timename,
    rinit : hitches$funs$rinit,
    rprocess : rproc_plugin(
      rprocess,
      step.fn=hitches$funs$step.fn,
      rate.fn=hitches$funs$rate.fn
    ),
    dprocess : hitches$funs$dprocess,
    dmeasure : hitches$funs$dmeasure,
    rmeasure : hitches$funs$rmeasure,
    skeleton : skel_plugin(
      skeleton,
      skel.fn=hitches$funs$skeleton
    ),
    dprior : hitches$funs$dprior,
    rprior : hitches$funs$rprior,
    partrans : parameter_trans(
      toEst=hitches$funs$toEst,
      fromEst=hitches$funs$fromEst
    ),
    params : params,
    covar : covar,
    accumvars : accumvars,
    solibs : if (is.null(hitches$lib)) {
      .solibs
    } else {
      c(list(hitches$lib),.solibs)
    },
    userdata : .userdata
    coef : coef // coef(object)
  }


