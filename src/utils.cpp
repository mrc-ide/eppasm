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

vec sumByAG (vec B, uvec age_of_interest) {
  uvec C = unique(age_of_interest);
  vec A = zeros(C.n_elem);
  uword current_age_group = age_of_interest(0); // first age group
  for (uword i = 0; i < age_of_interest.n_elem; ++i) {
    if (age_of_interest(i)==current_age_group) 
      A(current_age_group) += B(i);
    else {
      ++current_age_group;
      A(current_age_group) += B(i);
    }
  }
  return A;
}

mat sumByAG (mat B, uvec age_of_interest) {
  uvec C = unique(age_of_interest);
  mat A = zeros(C.n_elem, B.n_cols);
  uword current_age_group = age_of_interest(0); // first age group
  for (uword j = 0; j < B.n_cols; ++j) {
    for (uword i = 0; i < age_of_interest.n_elem; ++i) {
      if (age_of_interest(i)==current_age_group) 
        A(current_age_group, j) += B(i, j);
      else {
        ++current_age_group;
        A(current_age_group, j) += B(i, j);
      }
    } // end age-groups
    current_age_group = age_of_interest(0); // reset
  } // end columns
  return A;
}

cube sweepX23(cube A, mat B) {
  mat current(A.n_cols, A.n_slices);
  if (size(current) != size(B))
    Rf_error("Mat to sweep is of wrong size.");
  cube out(size(A));
  for (uword i = 0; i < A.n_rows; ++i) {
    current = A.row(i);
    out.row(i) = current % B;
  }
  return out;
}

// R sweep(array, 3:4, mat, "*")
field<cube> sweepX34 (field<cube> A, mat B, uword max_slice) {
  uword n_slice = (max_slice > 0) ? max_slice : A(0).n_slices;
  field<cube> FC(size(A)); //copied
  FC.for_each([&] (cube& X) { 
    X.zeros(A(0).n_rows, A(0).n_cols, n_slice); // resize the cube if needed
  });
  for (uword sex = 0; sex < A.n_cols; ++sex)
    for (uword agr = 0; agr < n_slice; ++agr)
      FC(sex).slice(agr) = A(sex).slice(agr) * B(agr, sex);
  return FC;
}

// rep each mat row by span
mat rep_each_row_by_span(mat M0, vec v_span) {
  mat M1( accu(v_span), M0.n_cols);
  vec my_leng = cumsum(v_span) - 1;
  uword j = 0;
  for (int i = 0; i < M0.n_rows; ++i) {
    M1.rows(j, my_leng(i)) = repmat(M0.row(i), v_span(i), 1);
    j = my_leng(i) + 1;
  }
  return M1;
}

vec col_sum_cube_to_vec(cube C) {
  mat A = sum(C);
  return trans(sum(A));
}