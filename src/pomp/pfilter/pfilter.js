let pfilter = function(np){ 
    let covartimeindex = this.pompclass.covar.cname.indexOf(this.pompclass.covar.timesname);
    if(this.pompclass.t0 == this.pompclass.covar.rows[0][covartimeindex]){
        this.pompclass.covar.curent = this.pompclass.covar.rows[0];
    } else{
        // function tu create t0 value
        this.pompclass.covar.curent.pop = 1326748;
        this.pompclass.covar.curent.birthrate = 2;
    }

    let res = this.pompclass.initializer(this.pompclass.states.current, this.pompclass.params.current, this.pompclass.covar.curent);
    for(let i=0; i<np; i++){
        this.pompclass.xstart.push(res);
    }
    this.pompclass.loglik = [this.pompclass.data.rows.length]
    
    this.pompclass.x = this.pompclass.xstart;
    for(let i = 0; i < np; i++){
        this.pompclass.rprocess( 0, this.pompclass.dt, this.pompclass.x[i], this.pompclass.params.current, this.pompclass.covar.curent);
    }
}

module.exports = pfilter;