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
#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
using namespace boost;
#pragma once

SEXP get_value(SEXP list, const char *str);
bool has_value(SEXP list, const char *str);

typedef multi_array<double, 1> boost1D;
typedef multi_array<int,    1> boost1I;
typedef multi_array<double, 2> boost2D;
typedef multi_array<double, 3> boost3D;
typedef multi_array<double, 4> boost4D;
typedef multi_array<double, 5> boost5D;

typedef std::vector<double>    dvec;
typedef std::vector<int>       ivec;

typedef multi_array_ref<double, 1> boost1D_ptr;
typedef multi_array_ref<int,    1> boost1I_ptr;
typedef multi_array_ref<double, 2> boost2D_ptr;
typedef multi_array_ref<double, 3> boost3D_ptr;
typedef multi_array_ref<double, 4> boost4D_ptr;
typedef multi_array_ref<double, 5> boost5D_ptr;

typedef multi_array_types::index_range in;

array<boost2D_ptr::index, 2> get_dim_2D(SEXP array, const char *str);
array<boost3D_ptr::index, 3> get_dim_3D(SEXP array, const char *str);
array<boost4D_ptr::index, 4> get_dim_4D(SEXP array, const char *str);

boost2D sumByAG (const boost2D& B, const boost1I& age_of_interest, int new_size);
boost1D sumByAG (const boost1D& B, const boost1I& age_of_interest, int new_size);
dvec sumByAG (const dvec& B, const ivec& age_of_interest, int new_size);

// Boost array NA/INF to zero: num/0.0 or 0.0/0.0
template <class K>
void replace_na_with (K& A, double B = 0) {
  for (auto i = A.data(); i < (A.data() + A.num_elements()); ++i)
    if (std::isnan(*i) || std::isinf(*i))
      *i = B;
}

// Boost array all to zero
template <class K>
void zeroing (K& A) {
  for (auto i = A.data(); i < (A.data() + A.num_elements()); ++i)
      *i = 0.0;
}

// Boost array sum
template <class K>
double sumArray (K& A) {
  double sum = 0;
  for (auto i = A.data(); i < (A.data() + A.num_elements()); ++i)
    sum += *i;
  return sum;
}

// stand vector sum
template <class K>
double sum_vector (K& A) {
  double sum = 0;
  for (auto i = A.data(); i < (A.data() + A.size()); ++i)
    sum += *i;
  return sum;
}

// Boost array multiply each cell
template <class K>
K multiply_with (K A, double X) { // not by reference
  for (auto i = A.data(); i < (A.data() + A.num_elements()); ++i)
      *i *= X;
  return A;
}

// Boost array multiply each cell
template <class K>
void multiply_with_inplace (K& A, double X) { // by reference
  for (auto i = A.data(); i < (A.data() + A.num_elements()); ++i)
      *i *= X;
}

// Boost array multiply each cell
template <class K>
K add_to_each (K A, double X) { // not by reference
  for (auto i = A.data(); i < (A.data() + A.num_elements()); ++i)
      *i += X;
  return A;
}

// Boost array multiply each cell
template <class K>
void add_to_each_inplace (K& A, double X) { // by reference
  for (auto i = A.data(); i < (A.data() + A.num_elements()); ++i)
      *i += X;
}

// Boost array substract each cell
template <class K>
K substract_from_each (K A, double X) { // not by reference
  for (auto i = A.data(); i < (A.data() + A.num_elements()); ++i)
      *i -= X;
  return A;
}

// Boost array fill
template <class K>
K fill_with (K A, double B) {
  for (auto i = A.data(); i < (A.data() + A.num_elements()); ++i)
      *i = B;
  return A;
}

// Boost array fill inplace
template <class K>
void fill_with_inplace (K& A, double B) {
  for (auto i = A.data(); i < (A.data() + A.num_elements()); ++i)
      *i = B;
}