

let pompclass = require('./pompclass.js');
let interpolator = require('./interpolator.js');
let reulermultinom = require('./reulermultinom.js');
let pfilter = require('./pfilter/pfilter.js');
Object.prototype.clone = function() {
    var newObj = (this instanceof Array) ? [] : {};
    for (i in this) {
      if (i == 'clone') continue;
      if (this[i] && typeof this[i] == "object") {
        newObj[i] = this[i].clone();
      } else newObj[i] = this[i]
    } return newObj;
  };

module.exports = {
    pompclass : pompclass,
    interpolator : interpolator,
    pfilter : pfilter,
    reulermultinom : reulermultinom
};
