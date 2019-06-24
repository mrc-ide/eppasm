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
void popC::art_distribute (const dvec& art_need, const Parameters& p,
                           const StateSpace& s) {
  if (!p.ad.med_cd4init_input[s.year]) {
    if (p.ad.art_alloc_method == 4L) { // by lowest CD4
      // Calculate proportion to be initiated in each CD4 category
      dvec init_pr(s.NG);
      for (int cd4 = s.hDS - 1; cd4 > 0; --cd4) { //6->0
        dvec elig_hm(s.NG);
        for (int sex = 0; sex < s.NG; sex++)
          for (int age = 0; age < s.hAG; age++)
            elig_hm[sex] += art_elig_[sex][age][cd4];
        if ( (elig_hm[s.m_idx] == 0) && (elig_hm[s.f_idx] == 0) )
          init_pr = elig_hm;
        else {
          double x;
          for (int sex = 0; sex < s.NG; ++sex) {
            x = art_need[sex] / elig_hm[sex];
            init_pr[sex] = ( (x > 1) | std::isnan(x) | std::isinf(x)) ? 1 : x;
          }
        }
        for (int sex = 0; sex < s.NG; ++sex)
          for (int agr = 0; agr < s.hAG; ++agr)
            art_init_[sex][agr][cd4] = art_elig_[sex][agr][cd4] * init_pr[sex];
      }
    } 
    else { // Spectrum Manual p168--p169, 
      int A = s.h_age15plus_idx[0] - 1;
      boost1D artX(extents[s.NG]), artY(extents[s.NG]);
      for (int sex = 0; sex < s.NG; sex++)
        for (int agr = A; agr < s.hAG_15plus; agr++)
          for (int cd4 = 0; cd4 < s.hDS; cd4++) {
            artX[sex] += art_elig_[sex][agr][cd4] * p.nh.cd4_mort[sex][agr][cd4];
            artY[sex] += art_elig_[sex][agr][cd4];
          }
      double xx;
      for (int sex = 0; sex < s.NG; sex++)
        for (int agr = A; agr < s.hAG_15plus; agr++)
          for (int cd4 = 0; cd4 < s.hDS; cd4++) {
            xx = (p.nh.cd4_mort[sex][agr][cd4] / artX[sex] *
                  p.ad.art_alloc_mxweight +
                  ((1 - p.ad.art_alloc_mxweight) / artY[sex]) ) *
                art_elig_[sex][agr][cd4] * art_need[sex];
            art_init_[sex][agr][cd4] =
              (xx > art_elig_[sex][agr][cd4]) ? art_elig_[sex][agr][cd4] : xx;
          }
    }
  }
  else {
    int CD4_LO[] = {500,  350, 250, 200, 100, 50,  0 };
    int CD4_UP[] = {1000, 500, 350, 250, 200, 100, 50};
    int j = p.ad.med_cd4init_cat[s.year] - 1; // R to C++
    double pr_below = (p.ad.median_cd4init[s.year] - CD4_LO[j]) / 
                      (CD4_UP[j] - CD4_LO[j]);
    dvec elig_below(s.NG);
    for (int sex = 0; sex < s.NG; sex++)
      for (int agr = 0; agr < s.hAG; agr++)
        elig_below[sex] += art_elig_[sex][agr][j] * pr_below;
    dvec A(s.NG);
    if (j < (s.hDS - 1)) {
      for (int sex = 0; sex < s.NG; sex++) {
        for (int agr = 0; agr < s.hAG; agr++)
          for (int cd4 = j+1; cd4 < s.hDS; cd4++)
            A[sex] += art_elig_[sex][agr][cd4];
        elig_below[sex] += A[sex];
      }
    }
    dvec elig_above(s.NG);
    dvec B(s.NG);
    for (int sex = 0; sex < s.NG; sex++) {
      for (int agr = 0; agr < s.hAG; agr++)
          B[sex] += art_elig_[sex][agr][j] * (1.0 - pr_below);
      elig_above[sex] += B[sex];
    }
    if (j > 1) {
      dvec C(s.NG);
      for (int sex = 0; sex < s.NG; sex++) {
        for (int agr = 0; agr < s.hAG; agr++)
          for (int cd4 = 0; cd4 < j-1; cd4++)
            C[sex] = art_elig_[sex][agr][cd4];
        elig_above[sex] += C[sex];
      }
    }
    dvec initpr_below(s.NG), initpr_above(s.NG), initpr_medcat(s.NG);
    double x, y;
    for (int sex = 0; sex < s.NG; ++sex) {
      x = art_need[sex] * 0.5 / elig_below[sex];
      y = art_need[sex] * 0.5 / elig_above[sex];
      initpr_below[sex] = (x < 1) ? x : 1;
      initpr_above[sex] = (y < 1) ? y : 1;
      initpr_medcat[sex] = initpr_below[sex] *      pr_below + 
                           initpr_above[sex] * (1 - pr_below);
    }
    if (j < (s.hDS - 1)) {
      for (int sex = 0; sex < s.NG; sex++)
        for (int agr = 0; agr < s.hAG; agr++)
          for (int cd4 = j + 1; cd4 < s.hDS; cd4++)
            art_init_[sex][agr][cd4] = 
              art_elig_[sex][agr][cd4] * initpr_below[sex];
    }
    for (int sex = 0; sex < s.NG; sex++)
      for (int agr = 0; agr < s.hAG; agr++)
        art_init_[sex][agr][j] =
          art_elig_[sex][agr][j] * initpr_medcat[sex];
    if (j > 0) {
      for (int sex = 0; sex < s.NG; sex++)
        for (int agr = 0; agr < s.hAG; agr++)
          for (int cd4 = 0; cd4 < j - 1; cd4++)
            art_init_[sex][agr][cd4] = 
              art_elig_[sex][agr][cd4] * initpr_above[sex];
    }
  }
}