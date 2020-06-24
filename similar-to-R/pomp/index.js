let interpolatorFactory = require('./mathlib/interpolatorFactory')
let simulator = require('./simulator')
Pomp = function(input) {
    if(input == undefined){
        input = {};
    }
    this.data = input.data || [];
    this.times = input.times || [];
    this.t0 = input.t0 || 0;
    this.rprocess = input.rprocess || {};
    this.dprocess = input.dprocess || {};
    this.rmeasure = input.rmeasure || {};
    this.dmeasure = input.dmeasure || {};
    this.measurement_model = input.measurement_model || {};
    this.skeleton = input.skeleton || {};
    this.initializer = input.initializer || {};
    this.rprior = input.rprior || {};
    this.dprior = input.dprior || {};
    this.params = input.input || {};
    this.covar = input.covar || {};
    this.tcovar = input.tcovar || {};
    this.obsnames = input.obsnames || {};
    this.statenames = input.statenames || {};
    this.paramnames = input.paramnames || {};
    this.covarnames = input.covarnames || {};
    this.zeronames = input.zeronames || {};
    this.PACKAGE = input.PACKAGE || {};
    this.fromEstimationScale = input.fromEstimationScale || {};
    this.toEstimationScale = input.toEstimationScale || {};
    this.globals = input.globals || {};
    this.cdir = input.cdir || {};
    this.cfile = input.cfile || {};
    this.shlib_args = input.shlib_args || {};

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
    this.simulator = simulator;

}

module.exports = Pomp;

