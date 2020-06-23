let pfilter = function(params){
    if(params == undefined){
        params = {};
    }
    object = params.object || {};
    params = params.params || {};
    Np = params.Np || {};
    tol = params.tol || 1e-17;
    max_fail = params.max_fail || Infinity;
    pred_mean = params.pred_mean || false;
    pred_var = params.pred_var || false;
    filter_mean = params.filter_mean || false;
    filter_traj = params.filter_traj || false;
    save_states = params.save_states || false;
    save_params = params.save_params || false;
    verbose = params.verbose || false;
    
    times = [...[object.t0],...object.times]
    ntimes = times.length - 1;

    init_x = object.initializer({covar: object.interpolator(object.t0),params:params});

    result = {};
    result.predictionMean = []
    result.predictionVariance = []

    return result;
}

module.exports = pfilter;