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

vec sumByAG (vec B, uvec age_of_interest);
mat sumByAG (mat B, uvec age_of_interest);
cube sweepX23(cube A, mat B);
// R sweep(array, 3:4, mat, "*")
field<cube> sweepX34 (field<cube> A, mat B, uword max_slice=0);
mat rep_each_row_by_span(mat M0, vec v_span);
vec col_sum_cube_to_vec(cube C);

// Pmin
template <class K>
K p_min (K A, double B, bool na_rm = true) {
  A.elem( find(A > B)).fill(B);
  if (na_rm)
    A.replace(datum::nan, B);
  return A;
}
// Pmax
template <class K>
K p_max (K A, double B, bool na_rm = true) {
  A.elem( find(A < B)).fill(B);
  if (na_rm)
    A.replace(datum::nan, B);
  return A;
}