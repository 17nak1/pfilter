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
let currentParams = []; 

// 1st data set;
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
lines = file.split('\n');
let currentParams_name = lines[0].split(',');
for (let i = 1; i < lines.length ; i++) {
  temp = lines[i].split(',');
  if(temp.length > 1) {
    temp = temp.map(function (x) {return Number(x)});
    currentParams.push(temp);
  }
}

//remove "\r" from the last variable's name.
let st = currentParams_name[currentParams_name.length - 1];
currentParams_name[currentParams_name.length - 1] = st.substring(0, st.length - 1);
st = dataCovar_name[dataCovar_name.length - 1];
dataCovar_name[dataCovar_name.length - 1] = st.substring(0, st.length - 1);
st = dataCases_name[dataCases_name.length - 1];
dataCases_name[dataCases_name.length - 1] = st.substring(0, st.length - 1);

let sortedCurrentParams = new Array(currentParams.length).fill(Array(currentParams[0].length));
// sortedCurrentParams is sorted currentParams based on snippet.paramnames.
for(let i = 0; i < snippet.paramnames.length; i++) {
  for( let j = 0; j < currentParams_name.length; j++) {
    if(snippet.paramnames[i] === currentParams_name[j] || snippet.paramnames[i] +"_0" === currentParams_name[j]) {
      for (let k = 0; k < currentParams.length; k++) {
        sortedCurrentParams[k][i] = currentParams[k][j];
      }
    }
  }       
}

let params_ic_fit = [];
let params_mod_fit = ["R0", "amplitude", "mu", "rho", "psi"];
let current_params = sortedCurrentParams[0];



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
  coef: current_params,
}
let d1 = [], d2 = [];
for (let i = 0; i < pomp.covar.length; i++) {
  d1.push([Number(pomp.tcovar[i]), Number(pomp.covar[i][0])])
  d2.push([Number(pomp.tcovar[i]), Number(pomp.covar[i][1])])
}

pomp.population = mathLib.interpolator(d1);
pomp.birthrate = mathLib.interpolator(d2);
let t = new Date()
let pf = pfilter({object: pomp, params: current_params, Np: 200, filterMean: true, predMean: true, maxFail: 3000})

console.log(new Date() - t, pf.loglik)

let createCsvWriter = require('csv-writer').createArrayCsvWriter;
  let csvWriter = createCsvWriter({
    header: [],
    path: rootDir+'/samples/predmean.csv',
    append : true
  })
  csvWriter.writeRecords( pf.predMean)
  .then(() => {
  console.log('...loglik ')
  })
