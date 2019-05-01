fs = require('fs')
rootDir = '..'

let pfilter = require('./pfilter.js')

let dataCases = []
let dataCovar = []


  //* 1st data set
  let London_covar = fs.readFileSync(rootDir+'/../samples/London_covar.csv').toString()
  var lines = London_covar.split('\n')
  for (let i = 1; i < lines.length - 1; i++) {
    dataCovar.push(lines[i].split(','))
  }

  //* 2nd data set
  let London_BiData = fs.readFileSync(rootDir+'/../samples/London_BiData.csv').toString()
  var lines = London_BiData.split('\n')
  for (let i = 1; i < lines.length - 1; i++) {
    dataCases.push(lines[i].split(','))
  }

pfilter.run({
    dataCases : dataCases,
    dataCovar : dataCovar,
    params : [3.132490e+01, 3.883620e-01, 7.305000e+01, 6.469830e-04, 4.566000e+01, 4.598709e-01, 1.462546e-01, 3.399189e-02, 2.336327e-04, 4.221789e-07, 9.657741e-01 ],
    Np : 100,
    dt : 1 / 365.25,
    times : [1940, 1944]})

// const createCsvWriter = require('csv-writer').createArrayCsvWriter;
//   const csvWriter = createCsvWriter({
//     header: ['S', 'E', 'I', 'R', 'H'],
//     path: rootDir +'/../samples/predmean.csv'
//   })
  
//   csvWriter.writeRecords(pfilter.predictionMean)
//     .then(() => {
//     console.log('...predictionMean')
//   })