// Copyright (C) 2019  Kinh Nguyen

// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.

// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.

// You should have received a copy of the GNU General Public License along
// with this program.  If not, see <http://www.gnu.org/licenses/>.

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <armadillo>
using namespace arma;

SEXP get_value(SEXP list, const char *str) {
  SEXP out = R_NilValue, names = Rf_getAttrib(list, R_NamesSymbol);
  int i;
  for ( i = 0; i < Rf_length(list); i++ ) {
    if ( strcmp(CHAR(STRING_ELT(names, i)), str) == 0 ) {
      out = VECTOR_ELT(list, i);
        break;
    }
  }
  if ( out == R_NilValue )
    Rf_error("%s missing from list, check ?prepare_fp_for_Cpp", str);
  return out;
}

// NOTES:
// uvec, double, int, and bool are copied, these are cheap
// matrix and array are not copied with arma aux_mem=false
// uvec is used for indexing, default shifting offset to zero-based, vec isn't

uvec arma_uvec(SEXP sexp_vector, bool shift_index=true) {
  // check if arma update has advanced constructor for uvec to simplify this
  if (!Rf_isInteger(sexp_vector))
    Rf_error("Index need to be integer");
	int *values = INTEGER(sexp_vector);
  std::vector<int> ans(values, values + Rf_length(sexp_vector));
  if (shift_index)
    for (int i = 0; i < ans.size(); ++i) 
      --ans[i];
  uvec out = conv_to<uvec>::from(ans);
  return out;
}

vec arma_vec(SEXP sexp_vector) {
  vec out(REAL(sexp_vector), Rf_length(sexp_vector), false);
  return out;
}

mat arma_mat(SEXP sexp_mat) {
  if (!Rf_isMatrix(sexp_mat))
    Rf_error("armat_mat needs input to be a matrix");
  mat out(REAL(sexp_mat), Rf_nrows(sexp_mat), Rf_ncols(sexp_mat), false);
  return out;
}

uvec array_dim(SEXP array) {
  SEXP r_dims = Rf_protect(Rf_getAttrib(array, R_DimSymbol));
  uvec dims(Rf_length(r_dims));
  int *rdims = INTEGER(r_dims);
  for (int i = 0; i < dims.n_elem; ++i)
    dims(i) = rdims[i];
  UNPROTECT(1);
  return dims;
}

cube arma_cube(SEXP sexp_cube) {
  uvec dimcube = array_dim(sexp_cube);
  cube out(REAL(sexp_cube), dimcube(0), dimcube(1), dimcube(2), false);
  return out;
}

field<cube> arma_4D(SEXP sexp_field) {
  uvec dimfield = array_dim(sexp_field);
  int N = dimfield(2) * dimfield(3);
  cube f3D(REAL(sexp_field), dimfield(0), dimfield(1), N, false);
  field<cube> out(dimfield(3));
  int start = 0, length_cube = dimfield(2);
  for (int i = 0; i < dimfield(3); ++i) {
    out(i) = f3D.slices(start, start + length_cube - 1);
    start += length_cube;
  }
  return out;
}

field<cube> arma_5D(SEXP sexp_field) {
  uvec dimfield = array_dim(sexp_field);
  int N = dimfield(2) * dimfield(3) * dimfield(4);
  cube f3D(REAL(sexp_field), dimfield(0), dimfield(1), N, false);
  field<cube> out(dimfield(4), dimfield(3)); // year in row, sex in col
  int start = 0;
  int length_cube = dimfield(2);
  for (int i = 0; i < dimfield(4); ++i) {
    for (int j = 0; j < dimfield(3); ++j) {
      out(i,j) = f3D.slices(start, start + length_cube - 1);
      start += length_cube;
    }
  }
  return out;
}

// sexp to arma pointer

vec vec2ptr (SEXP& sexp_obj) {
  vec V = vec(REAL(sexp_obj), GET_LENGTH(sexp_obj), false);
	V.zeros();
	return V;
}

cube cube2ptr (SEXP& sexp_obj) {
  uvec dimSEXP = array_dim(sexp_obj);
  if (dimSEXP.n_elem != 3)
  	Rf_error("cube2ptr requires 3D array"); // outputs 3D
	cube C = cube(REAL(sexp_obj), dimSEXP(0), dimSEXP(1), dimSEXP(2), false);
	C.zeros();
  return C;
}

field<cube> fourD2ptr (SEXP& sexp_obj) {
  uvec dimSEXP = array_dim(sexp_obj);
  if (dimSEXP.n_elem != 4) 
  	Rf_error("fourD2ptr requires 4D array"); // e.g. pop, hivpop
	mat M(REAL(sexp_obj), dimSEXP(0) * dimSEXP(1) * dimSEXP(2), dimSEXP(3), false);
	field<cube> F(dimSEXP(3)); // years
	for (int i = 0; i < dimSEXP(3); ++i)
		F(i) = cube(M.colptr(i), dimSEXP(0), dimSEXP(1), dimSEXP(2), false);
	F.for_each( [&](cube& X) { X.zeros(); } );
	return F;
}

field<cube> fiveD2ptr (SEXP& sexp_obj) {
  uvec dimSEXP = array_dim(sexp_obj);
  if (dimSEXP.n_elem != 5) 
  	Rf_error("fiveD2ptr requires 5D array"); // art pop
	mat M(REAL(sexp_obj), dimSEXP(0) * dimSEXP(1) * dimSEXP(2), 
												dimSEXP(3) * dimSEXP(4), false);
	field<cube> F(dimSEXP(4), dimSEXP(3)); // years, sexes
	uword col_count = 0;
	for (int year = 0; year < dimSEXP(4); ++year)
		for (uword sex = 0; sex < 2; ++sex) {
			F(year, sex) = cube(M.colptr(col_count), 
													dimSEXP(0), dimSEXP(1), dimSEXP(2), false);
			++col_count;
		}
	F.for_each( [&](cube& X) { X.zeros(); } );
	return F;
}