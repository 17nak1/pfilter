pomp = class pomp {
    constructor(){

        this.t0 = 0;
        this.params = {};
        this.states = {};
        this.data = {};
        this.covar = {};
        this.covar.curent = {};
        this.initializer = {};
        this.x = [];
    }
    readData(obj,fileaddress,timesname){
        obj.timesname = timesname;
        obj.cname = []
        obj.rows = []
        var temp = fs.readFileSync(fileaddress).toString()
        temp = temp.replace(/\r/g,'');
        temp = temp.replace(/\"/g,'');
        var lines = temp.split('\n')
        obj.cname = lines[0].split(',')
        for (let i = 1; i < lines.length; i++) {
            obj.rows.push(lines[i].split(','))
        }
    }

    pfilter(np){        
        let covartimeindex = this.covar.cname.indexOf(this.covar.timesname);
        if(this.t0 == this.covar.rows[0][covartimeindex]){
            //foreach read and make object covar curent
        } else{
            // function tu create t0 value
            this.covar.curent.pop = 1326748;
        }
        var obj = Object.assign({}, this.params.current, this.covar.curent);
        var res = this.initializer(obj);
    }
}

module.exports = pomp;