 as_matrix = function(x) {
  let dim, names
  let *xdim, nrow, ncol;
  PROTECT(dim = GET_DIM(x)); nprotect++;
  if (isNull(dim)) {
    PROTECT(x = duplicate(x)); nprotect++;
    PROTECT(names = GET_NAMES(x)); nprotect++;
    dim = NEW_INTEGER(2);
    xdim = INTEGER(dim);
    xdim[0] = LENGTH(x); xdim[1] = 1;
    SET_DIM(x,dim);
    SET_NAMES(x,R_NilValue);
    setrownames(x,names,2);
  } else if (LENGTH(dim) == 1) {
    PROTECT(x = duplicate(x)); nprotect++;
    PROTECT(names = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
    dim = NEW_INTEGER(2);
    xdim = INTEGER(dim);
    xdim[0] = LENGTH(x); xdim[1] = 1;
    SET_DIM(x,dim);
    SET_NAMES(x,R_NilValue);
    setrownames(x,names,2);
  } else if (LENGTH(dim) > 2) {
    PROTECT(x = duplicate(x)); nprotect++;
    PROTECT(names = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
    nrow = INTEGER(dim)[0];
    ncol = LENGTH(x)/nrow;
    dim = NEW_INTEGER(2);
    xdim = INTEGER(dim);
    xdim[0] = nrow; xdim[1] = ncol;
    SET_DIM(x,dim);
    SET_NAMES(x,R_NilValue);
    setrownames(x,names,2);
  }
  UNPROTECT(nprotect);
  return x;
}