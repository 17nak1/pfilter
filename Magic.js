/**
 *  @file       Magic.js  
 */

let pfilterCalculation = require('pfilter').pfilterCalculation
var dataCovar = [], dataCases = [], inputArr = [], init = [], res = []
var indx = new Array(12)

function start () {
  let req = new XMLHttpRequest()
  req.open('GET', './demo-compiled.js')
  req.onload = function () {
    code = req.responseText
  }
  var fileChooser = document.getElementById('file1-upload')
  fileChooser.onclick = function () {
    this.value = ''
  }
  document.getElementById('file1-upload').onchange = function () {
    document.getElementById('label-file1').innerHTML = 'Uploaded'
    document.getElementById('label-file1').style.backgroundColor = '#ffbf00'
    var file = this.files[0]
    dataCovar = []
    var reader = new FileReader()
    reader.onload = function () {
      var lines = this.result.split('\n')
      for (var line = 1; line < lines.length; line++) {
        dataCovar.push(lines[line].split(','))
      }
      console.log(dataCovar)
    }
    reader.readAsText(file)
  }
  var fileChooser = document.getElementById('file2-upload')
  fileChooser.onclick = function () {
    this.value = ''
  }
  document.getElementById('file2-upload').onchange = function () {
    document.getElementById('label-file2').innerHTML = 'Uploaded'
    document.getElementById('label-file2').style.backgroundColor = '#ffbf00'
    var file = this.files[0]
    dataCases = []
    var reader = new FileReader ()
    reader.onload = function () {
      var lines = this.result.split('\n')
      for (var line = 1; line < lines.length; line++) {
        dataCases.push(lines[line].split(','))
      }
      console.log(dataCases)
    }
    reader.readAsText(file)
  }

  let computeButton = document.querySelector('button#calc')
  let downloadButton = document.querySelector('button#download')
  downloadButton.style.display = 'none'
  computeButton.onclick = function () {
    inputArr = [], res = []
    computeButton.style.display = 'none'
    downloadButton.style.display = ''
    downloadButton.style.backgroundColor = '#D3D3D3'
    let setup1 = document.querySelector('#setup')
    let table = setup1.querySelector('table#table')
    let rows = table.querySelectorAll('tr')
    
    for (let i = 1; i < rows.length; i++) {
      let row = rows[i]
      let cols = row.querySelectorAll('td')
      let cell = cols[cols.length - 1]
      var input = cell.querySelector('input#valInp').value;
      inputArr.push(Number(input))// Read parameters from the table
      // if (init.length) {
      //   cell.querySelector('input#valInp').disabled = 'true'
      //   if (i !== 12) {
      //     cell.querySelector('input#valInp').value = ''
      //   } 
      // }
    }
    var data1 = []
    var data2 = []
    for (let i = 0; i < dataCovar.length - 1; i++) {
    data1.push([Number(dataCovar[i][0]), Number(dataCovar[i][1])])
    data2.push([Number(dataCovar[i][0]), Number(dataCovar[i][2])])
    }
    var interpolPopulation = interpolator(data1)
    var interpolBirth = interpolator(data2)

    let times = [inputArr[11], Number(dataCases[0][0])]
    inputArr.pop()
    setTimeout(function () {
      // var tem = inputArr[9]
      // inputArr[9] = inputArr[10]
      // inputArr[10] = tem
      res.push(pfilterCalculation({params:inputArr, Np:50,times:times, dt:1 / 365.25,runPredMean:1,  dataCases:dataCases, interpolPop:interpolPopulation, interpolBirth:interpolBirth}))
        
      console.log(res)
      res.splice(0, 0, ['S', 'E', 'I', 'R', 'H'])
    }, 0)
  }
  downloadButton.onclick = function () {
    Csv()
  }
}

function Csv () {
  var csv = ''
  res.forEach(function (row) {
    csv += row.join(',')
    csv += '\n'
  })
  var hiddenElement = document.createElement('a')

  hiddenElement.href = 'data:text/csv;charset=utf-8,' + encodeURI(csv)
  hiddenElement.setAttribute('download', 'result.csv')
  hiddenElement.click()
}
function activateDownload () {
  document.querySelector('button#download').disabled = false
  document.querySelector('button#download').style.backgroundColor = '#2ed573'
}

function interpolator(points) {
  var first, n = points.length - 1,
    interpolated,
    leftExtrapolated,
    rightExtrapolated;

  if (points.length === 0) {
    return function () {
      return 0
    }
  }

  if (points.length === 1) {
    return function () {
      return points[0][1]
    }
  }

  points = points.sort(function (a, b) {
    return a[0] - b[0]
  })
  first = points[0]

  leftExtrapolated = function (x) {
    var a = points[0], b = points[1];
    return a[1] + (x - a[0]) * (b[1] - a[1]) / (b[0] - a[0])
  }

  interpolated = function (x, a, b) {
    return a[1] + (x - a[0]) * (b[1] - a[1]) / (b[0] - a[0])
  }

  rightExtrapolated = function (x) {
    var a = points[n - 1], b = points[n];
    return b[1] + (x - b[0]) * (b[1] - a[1]) / (b[0] - a[0])
  }

  return function (x) {
    var i
    if (x <= first[0]) {
      return leftExtrapolated(x)
    }
    for (i = 0; i < n; i += 1) {
      if (x > points[i][0] && x <= points[i + 1][0]) {
        return interpolated(x, points[i], points[i + 1])
      }
    }
    return rightExtrapolated(x);
  }
}