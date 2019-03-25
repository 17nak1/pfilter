/* 
 *                                    User needs to provide values for 
 *                                    dataset : DeterministicSEIR_all.csv
 *                                    indices : LogLikIndex, [R0Index, AMPLITUDE, MU, RHO, PSI]
 *                                              the index of the  LogLik column and also indices of those parameters that we want to generate sets for them.
 *                                    tolerance : Determine how far from the best liklihood is acceptable.
 *                                    number of profile : Number of points to be consider in the interval of paramLimits.
 *                                    number of points : Number of  points to be generated.
 *                                    s : The value that is used to add noice in the best set and generate more points.
 *                                    indexMult : Defines how to generate new values ('divide & multiply' or 'subtract & add')
 *                                    indexInc : Defines how to generate new values (divide or multiply or both)(subtract or add or both)
 * 
 *                                    
 *                                    Also user needs to define the following for each estimating parameter through the function 'determineRunProperties';
 *                                    paramLimits : Lower and upper bound.
 *                                    logScale : If we consider calculation in the log scale, this value equals one.
 *                                    flagBound : If the generated values should be in the interval (0,1), this value equals one.
*/

let generateSets = require('./generateSets.js')
fs = require('fs')

/** dataset :
 *  Read the dataset and order the inputs based on ['R0', 'amplitude', 'gamma', 'mu', 'sigma', 'rho', 'psi', 'S_0', 'E_0', 'R_0', 'I_0']
 */
var dataset = []
var file = fs.readFileSync('./DeterministicSEIR_all.csv').toString()
var lines = file.split('\n')
for (let i = 1; i < lines.length; i++) {
  dataset.push(lines[i].split(','))
}
dataset.pop()
if(JSON.parse(lines[0].split(',')[9]) !== "R_0"){
  for (let i = 0; i < dataset.length; i++) {
    var tem = dataset[i][9]
    dataset[i][9] = dataset[i][10]
    dataset[i][10] = tem
  }
}
var tolerance = 50
var numberOfProfile = 50
var numberOfPoints = 2000

var indexInc = 0
var indexMult = 1
var s = 0.01

// Indicies for parameters. 
const paramObject = {
R0Index : 0,
AMPLITUDE : 1,
GAMMA : 2,
MU : 3,
SIGMA : 4,
RHO : 5,
PSI : 6,
S_0 : 7,
E_0 : 8,
R_0 : 9,
I_0 : 10,
LogLikIndex : 11
}

// This function should be define by the user for logScale, paramLimits and flagBound for each parameters.
determineRunProperties = function  (run, paramObject) {
  if (run == paramObject.R0Index) {
    logScale = 0
    paramLimits = [0,80]
    flagBound = 0
  } else if (run == paramObject.AMPLITUDE) {
    logScale = 0
    paramLimits = [0,1]
    flagBound = 1
  } else if (run == paramObject.MU) {
    logScale = 0
    paramLimits = [0,1]
    flagBound = 1
  } else if (run == paramObject.RHO) {
    logScale = 0
    paramLimits = [0,1]
    flagBound = 1
  } else if (run == paramObject.PSI) {
    logScale = 0   
    paramLimits = [0,1]
    flagBound = 0
  } 
  return [paramLimits, logScale, flagBound ]  
}
/** Main function
 *   For each parameters it call function 'generate' which provides csv file includes initial values for calculation in trajMatch.
 */
generateSet = function (dataset, paramObject, paramIndexArray, tolerance, numberOfProfile, numberOfPoints, indexInc, indexMult, s) {
  var Maxloglik
  // Reorder dataset descending based on LogLik column and find the maximum LogLik
  dataset.sort(sortFunction)
  Maxloglik = dataset[0][paramObject.LogLikIndex]

  for ( index = 0; index < paramIndexArray.length; index++) {
    var paramIndex = paramIndexArray[index]
    var paramLimits, logScale, flagBound
    [paramLimits, logScale, flagBound] = determineRunProperties (paramIndex, paramObject)
    generateSets.generate (dataset,Maxloglik, paramObject, paramIndex, tolerance, numberOfProfile, numberOfPoints, paramLimits, logScale, flagBound, indexInc, indexMult, s)
  }
}
// Helper function
function sortFunction(a, b) {
  if (Number(a[a.length - 1]) === Number(b[a.length - 1])) {
    return 0
  }
  else {
    return (Number(a[a.length - 1]) < Number(b[a.length - 1])) ? 1 : -1;
  }
}

generateSet (dataset, paramObject, [paramObject.R0Index, paramObject.AMPLITUDE, paramObject.MU, paramObject.RHO, paramObject.PSI], tolerance, numberOfProfile, numberOfPoints, indexInc, indexMult, s)