/**
 * 
 *
 */

fs = require('fs')
rootDir = '..'

let pfilter = require('./pfilter.js')

let dataCases = []
let dataCovar = []

//* 1st data set; read all rows and delete last one if it is ['']
let London_covar = fs.readFileSync(rootDir+'/../samples/London_covar.csv').toString()
let lines = London_covar.split('\n')
for (let i = 1; i < lines.length; i++) {
  dataCovar.push(lines[i].split(','))
}
if (dataCovar[dataCovar.length - 1].length === 1 ) {
  dataCovar.pop()
}

//* 2nd data set
let London_BiData = fs.readFileSync(rootDir+'/../samples/London_BiDataMain.csv').toString()
lines = London_BiData.split('\n')
for (let i = 1; i < lines.length ; i++) {
  dataCases.push(lines[i].split(','))
}
if (dataCases[dataCases.length - 1].length === 1 ) {
  dataCases.pop()
}

pfilter.run({
    maxFail : 500,
    dataCases : dataCases,
    dataCovar : dataCovar,
    runSaveStates : 1,
    params : [3.132490e+01, 3.883620e-01, 7.305000e+01, 6.469830e-04, 4.566000e+01, 4.598709e-01, 1.462546e-01, 3.399189e-02, 2.336327e-04, 4.221789e-07, 9.657741e-01 ],
    Np : 10000 ,
    dt : 1 / 365.25,// Input from pomp model
    timeZero : 1940})

