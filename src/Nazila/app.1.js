fs = require('fs')
rootDir = '~'

let rpois = require('./rpois.js')


var res = new Array()
for(let i = 0; i< 500000; i++){
  // console.log(i)
  var tres = rpois.rpoisOne(0.5)
    res.push(tres);
}


var file = fs.createWriteStream(rootDir +'/../samples/rpoisjs.txt');
file.on('error', function(err) { /* error handling */ });
res.forEach(function(v) { file.write(Math.round (v) + '\n'); });
file.end();
