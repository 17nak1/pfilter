let pomp = class Pomp {
    constructor(args) {
        if(args == undefined){
            args = {};
        }
        this.data = args.data || [];
        this.times = args.times || [];
        this.t0 = args.t0 || 0;
        this.rprocess = args.rprocess || {};
        this.dprocess = args.dprocess || {};
        this.rmeasure = args.rmeasure || {};
        this.dmeasure = args.dmeasure || {};
        this.measurement_model = args.measurement_model || {};
        this.skeleton = args.skeleton || {};
        this.initializer = args.initializer || {};
        this.rprior = args.rprior || {};
        this.dprior = args.dprior || {};
        this.params = args.args || {};
        this.covar = args.covar || {};
        this.tcovar = args.tcovar || {};
        this.obsnames = args.obsnames || {};
        this.statenames = args.statenames || {};
        this.paramnames = args.paramnames || {};
        this.covarnames = args.covarnames || {};
        this.zeronames = args.zeronames || {};
        this.PACKAGE = args.PACKAGE || {};
        this.fromEstimationScale = args.fromEstimationScale || {};
        this.toEstimationScale = args.toEstimationScale || {};
        this.globals = args.globals || {};
        this.cdir = args.cdir || {};
        this.cfile = args.cfile || {};
        this.shlib_args = args.shlib_args || {};

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

    }
  }

  module.exports = pomp;