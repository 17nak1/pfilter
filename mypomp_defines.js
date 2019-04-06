 mypomp_defines = {}

 mypomp_defines.as_matrix = function(x) {
  let dim, names
  let xdim, nrow, ncol
  dim = mypomp_defines.GET_DIM(x)
  if (typeof dim[0] === 'undefined' ) {
    // PROTECT(names = GET_NAMES(x)); nprotect++;
    dim = new Array(2)
    xdim = dim
    xdim[0] = x.length; xdim[1] = 1;
    // SET_DIM(x,dim);
    // SET_NAMES(x,R_NilValue);
    // setrownames(x,names,2);
  } else if (dim.length ===  1) {
    // PROTECT(names = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
    dim = new Array(2)
    xdim = dim
    xdim[0] = x.length; xdim[1] = 1;
    // SET_DIM(x,dim);
    // SET_NAMES(x,R_NilValue);
    // setrownames(x,names,2);
  } else if (dim.length > 2) {
    // PROTECT(names = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
    nrow = dim[0]
    ncol = x.length / nrow
    dim = new Array(2)
    xdim = dim
    xdim[0] = nrow; xdim[1] = ncol;
    // SET_DIM(x,dim);
    // SET_NAMES(x,R_NilValue);
    // setrownames(x,names,2);
  }
  // UNPROTECT(nprotect);
  return x;
}

mypomp_defines.GET_DIM = function (x) {
  
  if (x === null || typeof x === 'undefined' || x.length === 0) {
    return 0
  } else if (x[0].length === undefined) {
    return [1,x.length]
  } else {
    return [x.length, x[0].length]
  }
}

mypomp_defines.SET_DIM = function (x,newdim) {
  var rows = newdim[0]
  var cols = newdim[1]
  var matrix = []
  if (typeof x[0] === 'object') {
    matrix = x
  } else if (rows * cols === x.length) {
    for(var i = 0; i < x.length; i+= cols) {
      matrix.push(x.slice(i, cols + i));
    }
  } else {
    throw 'dimantions do not match the length of object'
  }
  return matrix;

}
mypomp_defines.makearray = function(rank, dim) {
  if ( dim[0] === 1) {
    return Array(dim[1])
  } else {
    return Array(dim[0]).fill(null).map(() => Array(dim[1]))
  }
}


module.exports = mypomp_defines
 // console.log(mypomp_defines.GET_DIM([1,1]))
 // console.log(mypomp_defines.as_matrix ([1,1]))
 // console.log(mypomp_defines.SET_DIM([[1,2],[3,4]],[2,2]))