//OMLY FOR PFILTER
/** Scales initial states using object.initializer()
  *  @return {Array} Each element is an array of states ordered based on object.statenames.
*/
exports.initState = function (object, params, nsim) {
  
  if(params === undefined) params = object.coef;
  if(isNaN(nsim)  || nsim === null) nsim = params.length;
  
  return do_init_state(object, params, nsim)
}

const do_init_state = function (object, params, ns) {
  
  params = params[0];
  // if (ns % nrep != 0) 
  //   throw new Error("in 'init.state': number of desired state-vectors 'nsim' is not a multiple of ncol('params')");

  let sNames =  object.statenames;
  let pNames = object.paramnames;
  let x = new Array(sNames.length);
  
  for(let i = 0; i < sNames.length; i++) {
    for(let j = 0; j < params.length; j++) {
      if (pNames[j] === sNames[i]) {
        x[i] = params[j];
      }
    }
  }
  let initVector =  object.initializer(object, x);  
  
  return  new Array(ns).fill(null).map(a => initVector);
}
