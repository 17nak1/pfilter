let pfilter = function(params){
    if(params == undefined){
        params = {};
    }
    Pomp.pfilter._parameters.params = params.params || {};
    Pomp.pfilter._parameters.Np = params.Np || {};
    Pomp.pfilter._parameters.tol = params.tol || 1e-17;
    Pomp.pfilter._parameters.max_fail = params.max_fail || Number.MAX_VALUE;
    Pomp.pfilter._parameters.pred_mean = params.pred_mean || false;
    Pomp.pfilter._parameters.pred_var = params.pred_var || false;
    Pomp.pfilter._parameters.filter_mean = params.filter_mean || false;
    Pomp.pfilter._parameters.filter_traj = params.filter_traj || false;
    Pomp.pfilter._parameters.save_states = params.save_states || false;
    Pomp.pfilter._parameters.save_params = params.save_params || false;
    Pomp.pfilter._parameters.verbose = params.verbose || false;

    return Pomp.pfilter._result;
}

pfilter._parameters = {};
pfilter._result = {};

module.exports = pfilter;