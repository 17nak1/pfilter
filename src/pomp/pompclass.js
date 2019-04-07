    
var readData = function(obj,fileaddress,timesname){
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

module.exports = {
    t0 : 0,
    params : { names:[],current:{}},
    states : { names:[],zeronames:[]},
    data : { curent : {},timesname:'',cname:'',rows:''},
    covar : { curent : {},timesname:'',cname:'',rows:'' },
    initializer : {},
    x : [],
    xstart : [],
    loglik : [],

    readData : readData
};