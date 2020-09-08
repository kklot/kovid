#include "utils.hpp"

SEXP get_value(SEXP list, const char *str) {
  SEXP out = R_NilValue, names = GET_NAMES(list);
  int l_len = GET_LENGTH(list);
  for (int i = 0; i < l_len; i++ ) {
    if ( strcmp(CHAR(STRING_ELT(names, i)), str) == 0 ) {
      out = VECTOR_ELT(list, i);
        break;
    }
  }
  if ( out == R_NilValue )
    Rf_warning("%s is NULL", str);
  return out;
}

bool has_value(SEXP list, const char *str) {
  SEXP names = GET_NAMES(list);
  int l_len = GET_LENGTH(list);
  for (int i = 0; i < l_len; i++ )
    if ( strcmp(CHAR(STRING_ELT(names, i)), str) == 0 ) {
      if (VECTOR_ELT(list, i) == R_NilValue) {
        Rf_warning("%s is NULL", str);
        return false;
      }
      else
        return true;
    }
  Rf_warning("%s does not exist.", str);
  return false;
}

boost1I array_dim(SEXP array) {
  SEXP r_dims = Rf_protect(Rf_getAttrib(array, R_DimSymbol));
  boost1I dims(extents[Rf_length(r_dims)]);
  int *rdims = INTEGER(r_dims);
  for (int i = 0; i < dims.num_elements(); ++i)
    dims[i] = rdims[i];
  UNPROTECT(1);
  return dims;
}

boost::array<boost1D_ptr::index, 1> get_dim_1D(SEXP r_in, const char *str) {
  int len = GET_LENGTH(get_value(r_in, str));
  boost::array<boost1D_ptr::index, 1> out = {{ len }};
  return out;
}

boost::array<boost2D_ptr::index, 2> get_dim_2D(SEXP r_in, const char *str) {
  SEXP r_dims = Rf_protect(Rf_getAttrib(get_value(r_in, str), R_DimSymbol));
  if (r_dims == R_NilValue)
    Rf_error("%s dimension is wrong! expecting a 2-D", str);
  int * rdims = INTEGER(r_dims);
  UNPROTECT(1);
  boost::array<boost2D_ptr::index, 2> out = {{ rdims[1], rdims[0] }};
  return out;
}

boost::array<boost3D_ptr::index, 3> get_dim_3D(SEXP r_in, const char *str) {
  SEXP r_dims = Rf_protect(Rf_getAttrib(get_value(r_in, str), R_DimSymbol));
  if (r_dims == R_NilValue || Rf_length(r_dims) != 3)
    Rf_error("%s dimension is wrong! expecting a 3-D", str);
  int *rdims = INTEGER(r_dims);
  UNPROTECT(1);
  boost::array<boost3D_ptr::index, 3> out = {{ rdims[2], rdims[1], rdims[0] }};
  return out;
}

boost::array<boost4D_ptr::index, 4> get_dim_4D(SEXP r_in, const char *str) {
  SEXP r_dims = Rf_protect(Rf_getAttrib(get_value(r_in, str), R_DimSymbol));
  if (r_dims == R_NilValue || Rf_length(r_dims) != 4)
    Rf_error("%s dimension is wrong! expecting a 4-D", str);
  int *rdims = INTEGER(r_dims);
  UNPROTECT(1);
  boost::array<boost4D_ptr::index, 4> 
    out = {{ rdims[3], rdims[2], rdims[1], rdims[0] }};
  return out;
}