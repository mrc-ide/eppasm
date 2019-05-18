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
#include <armadillo>
using namespace arma;

SEXP get_value(SEXP list, const char *str);

// NOTES:
// uvec, double, int, and bool are copied, these are cheap
// matrix and array are not copied with arma aux_mem=false
// uvec is used for indexing, default shifting offset to zero-based, vec isn't

uvec arma_uvec(SEXP sexp_vector, bool shift_index=true);
vec arma_vec(SEXP sexp_vector);
mat arma_mat(SEXP sexp_mat);
uvec array_dim(SEXP array);
cube arma_cube(SEXP sexp_cube);
field<cube> arma_4D(SEXP sexp_field);
field<cube> arma_5D(SEXP sexp_field);
vec vec2ptr (SEXP& sexp_obj);
cube cube2ptr (SEXP& sexp_obj);
field<cube> fourD2ptr (SEXP& sexp_obj);
field<cube> fiveD2ptr (SEXP& sexp_obj);