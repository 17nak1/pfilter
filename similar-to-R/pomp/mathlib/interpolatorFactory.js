let interpolatorFactory = function (covarnames, tcovar, covar) {
    let n = tcovar.length,
        nCovarnames = covarnames.length;                
    
    return function(t){
        point = {}; 
        if (n === 0) {
            for(let k = 0; k < nCovarnames; k++)
                point[covarnames[k]] = 0;
        }
      
        if (n === 1) {
            for(let k = 0; k < nCovarnames; k++)
                point[covarnames[k]] = +covar[0][k];
        }
    
        if (t <= tcovar[0]) {
            for(let k = 0; k < nCovarnames; k++)
                point[covarnames[k]] = +covar[0][k];
        } else if(t >= tcovar[n - 1]){
            for(let k = 0; k < nCovarnames; k++)
                point[covarnames[k]] = +covar[n-1][k];
        } else
            for (let i = 0; i < n; i++) {
                if (t > tcovar[i] && t <= tcovar[i+1]) {
                    for(let k = 0; k < nCovarnames; k++)
                        point[covarnames[k]] = +covar[i][k] + (t - tcovar[i]) * (covar[i+1][k] - covar[i][k]) / (tcovar[i+1] -tcovar[i])
                }
            }
        
        return point;
    }
  }

  
module.exports = interpolatorFactory;