pomp = class pomp {
    constructor(inputs){
        this.data = [];
        this.times = [];
        this.t0 = 0;
        this.timename = 'time';
        this.rinit= { def :'pomp_fun(slotname=rinit)'};
        this.rprocess = { def :'rproc_plugin()'};
        this.dprocess = { def :'pomp_fun(slotname=dprocess)'};
        this.rmeasure = { def :'pomp_fun(slotname=rmeasure)'};
        this.dmeasure = { def :'pomp_fun(slotname=dmeasure)'};
        this.rprior = { def :'pomp_fun(slotname=rprior)'};
        this.dprior = { def :'pomp_fun(slotname=dprior)'};
        this.skeleton = { def :'skel_plugin()'};
        this.partrans = { def :'parameter_trans()'};
        this.params = 0;
        this.states = [];
        this.covar = { def :'covariate_table()'};
        this.accumvars = '';
        this.solibs = {};
        this.userdata = {}
    }
    construct_pomp(inputs){

    }
}

module.exports = pomp;