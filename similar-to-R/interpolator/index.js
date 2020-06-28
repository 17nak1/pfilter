let interpolator = function(t){
    let n = this.covar.length,
        nCovarnames = this.covarnames.length;  
    point = {}; 
    if (n === 0) {
        for(let k = 0; k < nCovarnames; k++)
            point[this.covarnames[k]] = 0;
    }
  
    if (n === 1) {
        for(let k = 0; k < nCovarnames; k++)
            point[this.covarnames[k]] = +this.covar[0][k];
    }

    if (t <= this.tcovar[0]) {
        for(let k = 0; k < nCovarnames; k++)
            point[this.covarnames[k]] = +this.covar[0][k];
    } else if(t >= this.tcovar[n - 1]){
        for(let k = 0; k < nCovarnames; k++)
            point[this.covarnames[k]] = +this.covar[n-1][k];
    } else
        for (let i = 0; i < n; i++) {
            if (t > this.tcovar[i] && t <= this.tcovar[i+1]) {
                for(let k = 0; k < nCovarnames; k++)
                    point[this.covarnames[k]] = +this.covar[i][k] + (t - this.tcovar[i]) * (this.covar[i+1][k] - this.covar[i][k]) / (this.tcovar[i+1] -this.tcovar[i])
            }
        }
    
    return point;
}

module.exports = interpolator;