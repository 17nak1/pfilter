let pfilter = require('./pfilter');

Pomp = function(params) {
    if(params == undefined){
        params = {};
    }
    Pomp._parameters.data = params.data || [];
    Pomp._parameters.times = params.times || [];
    Pomp._parameters.t0 = params.t0 || 0;
    Pomp._parameters.rprocess = params.rprocess || {};
    Pomp._parameters.dprocess = params.dprocess || {};
    Pomp._parameters.rmeasure = params.rmeasure || {};
    Pomp._parameters.dmeasure = params.dmeasure || {};
    Pomp._parameters.measurement_model = params.measurement_model || {};
    Pomp._parameters.skeleton = params.skeleton || {};
    Pomp._parameters.initializer = params.initializer || {};
    Pomp._parameters.rprior = params.rprior || {};
    Pomp._parameters.dprior = params.dprior || {};
    Pomp._parameters.params = params.params || {};
    Pomp._parameters.covar = params.covar || {};
    Pomp._parameters.tcovar = params.tcovar || {};
    Pomp._parameters.obsnames = params.obsnames || {};
    Pomp._parameters.statenames = params.statenames || {};
    Pomp._parameters.paramnames = params.paramnames || {};
    Pomp._parameters.covarnames = params.covarnames || {};
    Pomp._parameters.zeronames = params.zeronames || {};
    Pomp._parameters.PACKAGE = params.PACKAGE || {};
    Pomp._parameters.fromEstimationScale = params.fromEstimationScale || {};
    Pomp._parameters.toEstimationScale = params.toEstimationScale || {};
    Pomp._parameters.globals = params.globals || {};
    Pomp._parameters.cdir = params.cdir || {};
    Pomp._parameters.cfile = params.cfile || {};
    Pomp._parameters.shlib_args = params.shlib_args || {};

    let ep = 'in \'pomp\': ';
    if (Pomp._parameters.data == undefined || Pomp._parameters.data.length == 0) 
        throw new Error(ep + '\'data\' is a required argument.');
    
    let tpos = 0;//Pomp._parameters.data[0].indexOf(Pomp._parameters.times);
    Pomp._parameters.times = Pomp._parameters.data.map(function(value,index) { return value[tpos]; });
    Pomp._parameters.data.map( function(value) { return value.splice(tpos, 1); });
    Pomp._parameters.obsnames = Pomp._parameters.data.shift();
    Pomp._parameters.times.shift();

    tpos = 0;//Pomp._parameters.covar[0].indexOf(Pomp._parameters.tcovar);
    Pomp._parameters.tcovar = Pomp._parameters.covar.map(function(value,index) { return value[tpos]; });
    Pomp._parameters.covar.map( function(value) { return value.splice(tpos, 1); });
    Pomp._parameters.covarnames = Pomp._parameters.covar.shift();
    Pomp._parameters.tcovar.shift();

    
}

Pomp._parameters = {};
Pomp.pfilter = pfilter;

module.exports = Pomp;

