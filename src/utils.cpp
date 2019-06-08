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
#include "utils.hpp"

SEXP get_value(SEXP list, const char *str) {
  SEXP out = R_NilValue, names = GET_NAMES(list);
  for (int i = 0; i < GET_LENGTH(list); i++ ) {
    if ( strcmp(CHAR(STRING_ELT(names, i)), str) == 0 ) {
      out = VECTOR_ELT(list, i);
        break;
    }
  }
  if ( out == R_NilValue )
    Rf_warning("%s missing from list, check ?prepare_fp_for_Cpp", str);
  return out;
}

bool has_value(SEXP list, const char *str) {
  SEXP names = GET_NAMES(list);
  for (int i = 0; i < GET_LENGTH(list); i++ )
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

boost::array<boost2D_ptr::index, 2> get_extents_2D(SEXP array) {
  SEXP r_dims = Rf_protect(Rf_getAttrib(array, R_DimSymbol));
  int *rdims = INTEGER(r_dims);
  UNPROTECT(1);
  boost::array<boost2D_ptr::index, 2> out = {{ rdims[1], rdims[0] }};
  return out;
}

boost::array<boost3D_ptr::index, 3> get_extents_3D(SEXP array) {
  SEXP r_dims = Rf_protect(Rf_getAttrib(array, R_DimSymbol));
  int *rdims = INTEGER(r_dims);
  UNPROTECT(1);
  boost::array<boost3D_ptr::index, 3> out = {{ rdims[2], rdims[1], rdims[0] }};
  return out;
}

boost::array<boost4D_ptr::index, 4> get_extents_4D(SEXP array) {
  SEXP r_dims = Rf_protect(Rf_getAttrib(array, R_DimSymbol));
  int *rdims = INTEGER(r_dims);
  UNPROTECT(1);
  boost::array<boost4D_ptr::index, 4> 
    out = {{ rdims[3], rdims[2], rdims[1], rdims[0] }};
  return out;
}

boost2D_ptr sexp_2D_to_boost(SEXP sexp_mat) {
  boost1I dimcube = array_dim(sexp_mat);
  boost2D_ptr out(REAL(sexp_mat), boost::extents[dimcube[1]][dimcube[0]]);
  return out;
}

// age_of_interest example is ag_idx, not shifted
boost2D sumByAG (const boost2D& B, const boost1I& age_of_interest, int new_size) 
{
  int current_age_group = age_of_interest[0],
      row = B.shape()[0], col = B.shape()[1]; // first age group
  boost2D A(extents[ row ][ new_size ]);
  for (int j = 0; j < row; ++j) {
    for (int i = 0; i < col; ++i) {
      if ( age_of_interest[i] == current_age_group)
        A[j][current_age_group - 1] += B[j][i];
      else {
        ++current_age_group;
        A[j][current_age_group - 1] += B[j][i];
      }
    } // end age-groups
    current_age_group = age_of_interest[0]; // reset
  } // end columns
  return A;
}

boost1D sumByAG (const boost1D& B, const boost1I& age_of_interest, int new_size) 
{
  boost1D A(extents[ new_size ]);
  int current_age_group = age_of_interest[0]; // first age group
  for (int i = 0; i < B.num_elements(); ++i) {
    if ( age_of_interest[i] == current_age_group)
      A[current_age_group - 1] += B[i];
    else {
      ++current_age_group;
      A[current_age_group - 1] += B[i];
    }
  } // end age-groups
  return A;
}

dvec sumByAG (const dvec& B, const ivec& age_of_interest, int new_size) {
  dvec A(new_size);
  int current_age_group = age_of_interest[0]; // first age group
  for (int i = 0; i < B.size(); ++i) {
    if ( age_of_interest[i] == current_age_group)
      A[current_age_group - 1] += B[i];
    else {
      ++current_age_group;
      A[current_age_group - 1] += B[i];
    }
  } // end age-groups
  return A;
}