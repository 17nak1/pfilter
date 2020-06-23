let pfilter = require('./pfilter');
let fmin = require ('fmin');
let mathLib = require('./mathLib');
let simulator = require ('./simulator.js');

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
    
    let tpos = 0;//this.data[0].indexOf(this.times);
    this.times = this.data.map(function(value,index) { return value[tpos]; });
    this.data.map( function(value) { return value.splice(tpos, 1); });
    this.obsnames = this.data.shift();
    this.times.shift();

    tpos = 0;//this.covar[0].indexOf(this.tcovar);
    this.tcovar = this.covar.map(function(value,index) { return value[tpos]; });
    this.covar.map( function(value) { return value.splice(tpos, 1); });
    this.covarnames = this.covar.shift();
    this.tcovar.shift();

    this.interpolator = function (t) {
        let first, 
            n = this.tcovar.length,
            interpolated,
            leftExtrapolated,
            rightExtrapolated,
            point = {};
      
        if (n === 0) {
            for(let i = 0; i < this.covarnames.length; i++)
                point[this.covarnames[i]] = 0;
        }
      
        if (n === 1) {
            for(let i = 0; i < this.covarnames.length; i++)
                point[this.covarnames[i]] = this.covar[0][i];
        }
      
        first = this.covar[0];
      
        leftExtrapolated = function (t) {
            let t_0 = this.tcovar[0], t_1 = this.tcovar[1];
            let x_0 = this.covar[0], x_1 = this.covar[1];
            for(let i = 0; i < this.covarnames.length; i++)
                point[this.covarnames[i]] = x_0[i] + (t - t_0) * (x_1[i] - x_0[i]) / (t_1 - t_0);

        }
      
        interpolated = function (t, x_0, x_1, t_0, t_1) {            
            for(let i = 0; i < this.covarnames.length; i++)
                point[this.covarnames[i]] = x_0[i] + (t - t_0) * (x_1[i] - x_0[i]) / (t_1 - t_0);
        }
      
        rightExtrapolated = function (t) {            
            let t_0 = this.tcovar[n - 1], t_1 = this.tcovar[n];
            let x_0 = this.covar[n - 1], x_1 = this.covar[n];                     
            for(let i = 0; i < this.covarnames.length; i++)
                point[this.covarnames[i]] = x_1[i] + (t - t_1) * (x_1[i] - x_o[i]) / (t_1 - t_0)
        }

        let i
        if (t <= tcovar[0]) {
            leftExtrapolated(t)
        }
        for (i = 0; i < n; i += 1) {
            if (t > this.tcovar[i] && t <= this.tcovar[i+1]) {
                interpolated(t, this.covar[i], this.covar[i + 1])
            }
        }
        rightExtrapolated(t);
        
      }
    
}


Pomp.pfilter = pfilter;
Pomp.fmin = fmin;
Pomp.mathLib = mathLib;
Pomp.simulator = simulator;

module.exports = Pomp;

