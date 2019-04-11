    
let readData = function(obj,fileaddress,timesname){
        obj.timesname = timesname;
        obj.cname = []
        obj.rows = []
        let temp = fs.readFileSync(fileaddress).toString()
        temp = temp.replace(/\r/g,'');
        temp = temp.replace(/\"/g,'');
        let lines = temp.split('\n')
        obj.cname = lines[0].split(',')
        for (let i = 1; i < lines.length; i++) {
            let temprow = {};
            let temprowarr = lines[i].split(',');
            for(let j =0; j < obj.cname.length; j++){
                temprow[obj.cname[j]] = +temprowarr[j]
            }
            obj.rows.push(temprow);
        }
    }

module.exports = {
    t0 : 0,
    dt : 0,
    params : { names:[],current:{}},
    states : { names:[],zeronames:[]},
    data : { curent : {},timesname:'',cname:'',rows:''},
    covar : { curent : {},timesname:'',cname:'',rows:'' },
    initializer : {},
    x : [],
    xstart : [],
    loglik : [],

    readData : readData,
    initializer : {},
    rprocess : {}
};