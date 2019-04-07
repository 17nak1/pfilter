var pfilter = function(obj,np){ 
    var pomp = obj.pompclass;       
    let covartimeindex = pomp.covar.cname.indexOf(pomp.covar.timesname);
    if(pomp.t0 == pomp.covar.rows[0][covartimeindex]){
        //foreach read and make object covar curent
    } else{
        // function tu create t0 value
        pomp.covar.curent.pop = 1326748;
    }
    var obj = Object.assign({}, pomp.params.current, pomp.covar.curent);

    var res = pomp.initializer(obj);
    for(var i=0; i<np; i++){
        pomp.xstart.push(res);
    }
    pomp.loglik = [pomp.data.rows.length]
}

module.exports = pfilter;