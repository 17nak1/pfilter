let interpolatorFactory = require('./mathlib/interpolatorFactory')
Pomp = function(params) {
    if(params == undefined){
        params = {};
    }
    this.data = params.data || [];
    this.times = params.times || [];
    this.t0 = params.t0 || 0;
    this.rprocess = params.rprocess || {};
    this.dprocess = params.dprocess || {};
    this.rmeasure = params.rmeasure || {};
    this.dmeasure = params.dmeasure || {};
    this.measurement_model = params.measurement_model || {};
    this.skeleton = params.skeleton || {};
    this.initializer = params.initializer || {};
    this.rprior = params.rprior || {};
    this.dprior = params.dprior || {};
    this.params = params.params || {};
    this.covar = params.covar || {};
    this.tcovar = params.tcovar || {};
    this.obsnames = params.obsnames || {};
    this.statenames = params.statenames || {};
    this.paramnames = params.paramnames || {};
    this.covarnames = params.covarnames || {};
    this.zeronames = params.zeronames || {};
    this.PACKAGE = params.PACKAGE || {};
    this.fromEstimationScale = params.fromEstimationScale || {};
    this.toEstimationScale = params.toEstimationScale || {};
    this.globals = params.globals || {};
    this.cdir = params.cdir || {};
    this.cfile = params.cfile || {};
    this.shlib_args = params.shlib_args || {};

    let ep = 'in \'pomp\': ';
    if (this.data == undefined || this.data.length == 0) 
        throw new Error(ep + '\'data\' is a required argument.');
    
    let tpos = this.data[0].indexOf(this.times);
    this.times = this.data.map(function(value,index) { return value[tpos]; });
    this.data.map( function(value) { return value.splice(tpos, 1); });
    this.obsnames = this.data.shift();
    this.times.shift();

    tpos = this.covar[0].indexOf(this.tcovar);
    this.tcovar = this.covar.map(function(value,index) { return value[tpos]; });
    this.covar.map( function(value) { return value.splice(tpos, 1); });
    this.covarnames = this.covar.shift();
    this.tcovar.shift();

    this.interpolator = interpolatorFactory(this.covarnames, this.tcovar, this.covar);

}

module.exports = Pomp;

