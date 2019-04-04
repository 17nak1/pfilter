let pomp = require('./src/pomp.js');
let pfilter = require('./src/pfilter.js');

let gompertz = new pomp({  
    data : [
        [0,1,1.07952841598157],
        [1,0.950636618475229,1.11861594539405]],
    times:'time',
    t0:0,
    params:{K:1,r:0.1,sigma:0.1,tau:0.1,X_0:1},
    rprocess:'discrete.time.sim('+
        'step.fun=_gompertz_simulator'+
        ')',
    rmeasure:'_gompertz_normal_rmeasure',
    dmeasure:'_gompertz_normal_dmeasure',
    skeleton:'map(_gompertz_skeleton,delta_t=1)',
    paramnames:['r','K','sigma','tau'],
    statenames:['X'],
    userdata : { 
        fromEstimationScale: function(params){
            exp(params)
        },
        toEstimationScale:function(params){
            log(params)
        }}
});
let pf = new pfilter({
    data : gompertz,
    np : 3
});

console.log(gompertz.constructor.name);
console.log(pf.constructor.name);
console.log('Finished');