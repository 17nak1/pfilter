let pomp = require('./src/pomp.js');
fs = require('fs')

let pf = new pomp();

pf.t0=1940;
pf.params.names = ['R0','amplitude','gamma','mu','sigma','rho','psi'];
pf.states.names = ['S','E','I','R','H'];
pf.states.zeronames = ['H'];
pf.params.current = { R0:3.132490e+01 , amplitude:3.883620e-01 , gamma:7.305000e+01 , mu:6.469830e-04 , sigma:4.566000e+01 ,rho: 4.598709e-01 ,psi: 1.462546e-01 ,S: 3.399189e-02 ,E:2.336327e-04 ,R:9.657741e-01,I:4.221789e-07};

pf.readData(pf.data,'./src/London_BiData.csv','time');
pf.readData(pf.covar,'./src/London_covar.csv','time');

pf.initializer = obj => {
    res = {};
    m = obj.pop/(obj.S+obj.E+obj.I+obj.R);
    res.S = Math.round(m*obj.S);
    res.E = Math.round(m*obj.E);
    res.I = Math.round(m*obj.I);
    res.R = Math.round(m*obj.R);
    res.H = 0;
    return res;
}

pf.pfilter();

console.log('Finished');