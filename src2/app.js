/**
 * 
 *
 */
// let mif2 = require('./pfilter.js');
let snippet = require('./modelSnippet.js');
let fs = require('fs');
let mathLib = require('./mathLib.js');
let { pfilter } = require('./pfilter.js');


rootDir = '.'

let dataCases = [];
let dataCasesTimes = [];
let dataCovar = [];
let dataCovarTimes = [];
let currentParams = {}; 

// 1st data set; read all rows and delete last one if it is ['']
let temp, file;
file = fs.readFileSync(rootDir+'/samples/London_covar.csv').toString();
let lines = file.split('\n');
let dataCovar_name = lines[0].split(',');
for (let i = 1; i < lines.length; i++) {
  temp = lines[i].split(',');
  if(temp.length > 1) {
    temp = temp.map(x => Number(x));
    dataCovarTimes.push(temp[0]);
    dataCovar.push(temp.slice(1));
  }
}

//* 2nd data set
file = fs.readFileSync(rootDir+'/samples/London_BiDataMain.csv').toString()
lines = file.split('\n');
let dataCases_name = lines[0].split(',');
for (let i = 1; i < lines.length ; i++) {
  temp = lines[i].split(',');
  if(temp.length > 1) {
    temp = temp.map(x => Number(x));
    dataCasesTimes.push(temp[0]);
    dataCases.push(temp.slice(1));
  }
}

//* 3nd data set and names
file = fs.readFileSync(rootDir+'/samples/initial_parameters.csv').toString()
lines = file.split(/\r\n|\n/);
let currentParams_name = lines[0].split(',');
for (let i = 1; i < lines.length ; i++) {
  temp = lines[i].split(',');
  if(temp.length > 1) {
    temp = temp.map(function (x) {return Number(x)});
    for(let j = 0; j < temp.length; j++){
      currentParams[currentParams_name[j]] = temp[j];
    }
    
  }
}

let params_ic_fit = [];
let params_mod_fit = ["R0", "amplitude", "mu", "rho", "psi"];


///////////////////////////////////////////////////

/////////////////////////////////////
const pomp = {
  data :  dataCases,
  times:  dataCasesTimes,
  t0: 1940,
  rprocess :  { type:"euler_sim", stepFunction: snippet.rprocess, deltaT: 1/365.25 },
  rmeasure: snippet.rmeas,
  covar: dataCovar,
  tcovar: dataCovarTimes,
  dmeasure: snippet.dmeasure,
  zeronames: snippet.zeronames,
  initializer: snippet.initz,
  toEstimationScale: snippet.toEst,
  fromEstimationScale: snippet.fromEst,
  statenames: snippet.statenames,
  paramnames: snippet.paramnames,
  coef: currentParams,
}
let d1 = [], d2 = [];
for (let i = 0; i < pomp.covar.length; i++) {
  d1.push([Number(pomp.tcovar[i]), Number(pomp.covar[i][0])])
  d2.push([Number(pomp.tcovar[i]), Number(pomp.covar[i][1])])
}

pomp.population = mathLib.interpolator(d1);
pomp.birthrate = mathLib.interpolator(d2);
let t = new Date()
let pf = pfilter({object: pomp, params: currentParams, Np: 200, filterMean: true, predMean: true, maxFail: 3000})

console.log(new Date() - t, pf.loglik)

let createCsvWriter = require('csv-writer').createArrayCsvWriter;
  let csvWriter = createCsvWriter({
    header: [],
    path: './samples/predmean.csv',
    append : true
  })
  csvWriter.writeRecords( pf.predMean)
  .then(() => {
  console.log('...loglik ')
  })