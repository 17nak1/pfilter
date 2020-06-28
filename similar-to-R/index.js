let pomp = require('./pomp')
let pfilter = require('./pfilter')
let interpolator = require('./interpolator')
let simulator = require('./simulator')

pomp.prototype.pfilter = pfilter;
pomp.prototype.interpolator = interpolator;
pomp.prototype.simulator = simulator;

        // this.simulator = simulator;
        // this.nosortResamp = nosortResamp;
module.exports = pomp;