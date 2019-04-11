let pfilter = function(np){ 
    this.pompclass.covar.curent = this.interpolator(this.pompclass.covar.rows)(this.pompclass.t0);
    let res = this.pompclass.initializer( this.pompclass.states.current, this.pompclass.params.current, this.pompclass.covar.curent);
    for(let i=0; i<np; i++){
        this.pompclass.xstart.push(res.clone());
    }
    this.pompclass.loglik = [this.pompclass.data.rows.length]
    
    let ntimes = this.pompclass.data.rows.length;
    let rproctimes = [this.pompclass.t0,this.pompclass.t0];
    for(let i = 0; i < ntimes; i++){
        rproctimes[0] = rproctimes[1];
        rproctimes[1] = this.pompclass.data.rows[i]['time'];
        let tempx = [this.pompclass.xstart.clone(),this.pompclass.xstart.clone()];
        nsteps = Math.ceil((rproctimes[1] - rproctimes[0]) / (this.pompclass.dt/(1+1.490116e-08)));
        ndt = (rproctimes[1] - rproctimes[0]) / nsteps;
        t = rproctimes[0];
        // zero col shoule be zero
        let cov = this.interpolator(this.pompclass.covar.rows)(t);
        for( k = 0; k < nsteps; k++){
            for( j = 0; j < np; j++){
                this.pompclass.rprocess( t, ndt, tempx[1][j], this.pompclass.params.current, cov);
            }
            t += ndt;
            if( k == nsteps - 2){
                ndt = rproctimes[1]-t;
                t = rproctimes[1]-ndt;
            }
        }
    }
}

module.exports = pfilter;