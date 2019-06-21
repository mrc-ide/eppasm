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
#include "Classes.hpp"

// calculate ART initiation distribution
void popC::art_distribute (const dvec& art_need) {
  if (!p.med_cd4init_input[year]) {
    if (p.art_alloc_method == 4L) { // by lowest CD4
      // Calculate proportion to be initiated in each CD4 category
      dvec init_pr(NG);
      for (int cd4 = hDS - 1; cd4 > 0; --cd4) { //6->0
        dvec elig_hm(NG);
        for (int sex = 0; sex < NG; sex++)
          for (int age = 0; age < hAG; age++)
            elig_hm[sex] += art_elig_[sex][age][cd4];
        if ( elig_hm[m_idx] == 0 & elig_hm[f_idx] == 0 )
          init_pr = elig_hm;
        else {
          double x;
          for (int sex = 0; sex < NG; ++sex) {
            x = art_need[sex] / elig_hm[sex];
            init_pr[sex] = ( (x > 1) | std::isnan(x) | std::isinf(x)) ? 1 : x;
          }
        }
        for (int sex = 0; sex < NG; ++sex)
          for (int agr = 0; agr < hAG; ++agr)
            art_init_[sex][agr][cd4] = art_elig_[sex][agr][cd4] * init_pr[sex];
      }
    } 
    else { // Spectrum Manual p168--p169, 
      int A = h_age15plus_idx[0] - 1;
      dvec artX(NG), artY(NG);
      for (int sex = 0; sex < NG; sex++)
        for (int agr = A; agr < hAG_15plus; agr++)
          for (int cd4 = 0; cd4 < hDS; cd4++) {
            artX[sex] += art_elig_[sex][agr][cd4] * p.cd4_mort[sex][agr][cd4];
            artY[sex] += art_elig_[sex][agr][cd4];
          }
      double xx;
      for (int sex = 0; sex < NG; sex++)
        for (int agr = A; agr < hAG_15plus; agr++)
          for (int cd4 = 0; cd4 < hDS; cd4++) {
            xx = (p.cd4_mort[sex][agr][cd4] / artX[sex] * p.art_alloc_mxweight +
                  ((1 - p.art_alloc_mxweight) / artY[sex]) ) *
                art_elig_[sex][agr][cd4] * art_need[sex];
            art_init_[sex][agr][cd4] =
              (xx > art_elig_[sex][agr][cd4]) ? art_elig_[sex][agr][cd4] : xx;
          }
    }
  }
  else {
    int CD4_LO[] = {500,  350, 250, 200, 100, 50,  0 };
    int CD4_UP[] = {1000, 500, 350, 250, 200, 100, 50};
    int j = p.med_cd4init_cat[year] - 1; // R to C++
    double pr_below = (p.median_cd4init[year] - CD4_LO[j]) / 
                      (CD4_UP[j] - CD4_LO[j]);
    dvec elig_below(NG);
    for (int sex = 0; sex < NG; sex++)
      for (int agr = 0; agr < hAG; agr++)
        elig_below[sex] += art_elig_[sex][agr][j] * pr_below;
    dvec A(NG);
    if (j < (hDS - 1)) {
      for (int sex = 0; sex < NG; sex++) {
        for (int agr = 0; agr < hAG; agr++)
          for (int cd4 = j+1; cd4 < hDS; cd4++)
            A[sex] += art_elig_[sex][agr][cd4];
        elig_below[sex] += A[sex];
      }
    }
    dvec elig_above(NG);
    dvec B(NG);
    for (int sex = 0; sex < NG; sex++) {
      for (int agr = 0; agr < hAG; agr++)
          B[sex] += art_elig_[sex][agr][j] * (1.0 - pr_below);
      elig_above[sex] += B[sex];
    }
    if (j > 1) {
      dvec C(NG);
      for (int sex = 0; sex < NG; sex++) {
        for (int agr = 0; agr < hAG; agr++)
          for (int cd4 = 0; cd4 < j-1; cd4++)
            C[sex] = art_elig_[sex][agr][cd4];
        elig_above[sex] += C[sex];
      }
    }
    dvec initpr_below(NG), initpr_above(NG), initpr_medcat(NG);
    double x, y;
    for (int sex = 0; sex < NG; ++sex) {
      x = art_need[sex] * 0.5 / elig_below[sex];
      y = art_need[sex] * 0.5 / elig_above[sex];
      initpr_below[sex] = (x < 1) ? x : 1;
      initpr_above[sex] = (y < 1) ? y : 1;
      initpr_medcat[sex] = initpr_below[sex] *      pr_below + 
                           initpr_above[sex] * (1 - pr_below);
    }
    if (j < (hDS - 1)) {
      for (int sex = 0; sex < NG; sex++)
        for (int agr = 0; agr < hAG; agr++)
          for (int cd4 = j + 1; cd4 < hDS; cd4++)
            art_init_[sex][agr][cd4] = 
              art_elig_[sex][agr][cd4] * initpr_below[sex];
    }
    for (int sex = 0; sex < NG; sex++)
      for (int agr = 0; agr < hAG; agr++)
        art_init_[sex][agr][j] =
          art_elig_[sex][agr][j] * initpr_medcat[sex];
    if (j > 0) {
      for (int sex = 0; sex < NG; sex++)
        for (int agr = 0; agr < hAG; agr++)
          for (int cd4 = 0; cd4 < j - 1; cd4++)
            art_init_[sex][agr][cd4] = 
              art_elig_[sex][agr][cd4] * initpr_above[sex];
    }
  }
}