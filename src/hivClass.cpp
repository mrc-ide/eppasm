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

void hivC::aging(const boost2D& ag_prob) {
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hAG; agr++)
      for (int cd4 = 0; cd4 < hDS; cd4++)
        data[year][sex][agr][cd4] = data[year-1][sex][agr][cd4];
  double nHup;
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hAG - 1; agr++)
      for (int cd4 = 0; cd4 < hDS; cd4++) {
        nHup = data[year-1][sex][agr][cd4] * ag_prob[sex][agr];
        data[year][sex][agr][cd4]   -= nHup;
        data[year][sex][agr+1][cd4] += nHup;
      }
  if (MODEL == 2) {
    for (int sex = 0; sex < NG; sex++)
      for (int agr = 0; agr < hDB; agr++)
        for (int cd4 = 0; cd4 < hDS; cd4++)
          data_db[year][sex][agr][cd4] = data_db[year-1][sex][agr][cd4];
    for (int sex = 0; sex < NG; sex++)
      for (int agr = 0; agr < hDB - 1; agr++)
        for (int cd4 = 0; cd4 < hDS; cd4++) {
          nHup = data_db[year-1][sex][agr][cd4] * ag_prob[sex][agr];
          data_db[year][sex][agr][cd4]   -= nHup;
          data_db[year][sex][agr+1][cd4] += nHup;
        }
  }
}

void hivC::add_entrants(const dvec& artYesNo) { // see pop.entrant_art
  if (MODEL == 1)
    for (int sex = 0; sex < NG; sex++)
      for (int cd4 = 0; cd4 < hDS; cd4++)
        data[year][sex][0][cd4] +=
          p.paedsurv_cd4dist[year][sex][cd4] * artYesNo[sex+2];
  if (MODEL == 2) // add to virgin then debut
    for (int sex = 0; sex < NG; sex++)
      for (int cd4 = 0; cd4 < hDS; cd4++)
        data_db[year][sex][0][cd4] +=
          p.paedsurv_cd4dist[year][sex][cd4] * artYesNo[sex+2];
}

void hivC::sexual_debut() {
  for (int sex = 0; sex < NG; sex++)
    for (int adb = 0; adb < hDB; adb++)
      for (int cd4 = 0; cd4 < hDS; cd4++) {
        double n_db = data_db[year][sex][adb][cd4] * p.db_pr[sex][adb];
        data[year][sex][adb][cd4]    += n_db;
        data_db[year][sex][adb][cd4] -= n_db;
      }
}

void hivC::deaths (const boost2D& survival_pr) {
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hAG; agr++)
      for (int cd4 = 0; cd4 < hDS; cd4++)
        data[year][sex][agr][cd4] *= survival_pr[sex][agr];
  if (MODEL == 2)
    for (int sex = 0; sex < NG; sex++)
      for (int adb = 0; adb < hDB; adb++)
        for (int cd4 = 0; cd4 < hDS; cd4++)
          data_db[year][sex][adb][cd4] *= survival_pr[sex][adb];
}

void hivC::migration (const boost2D& migration_pr) {
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hAG; agr++)
      for (int cd4 = 0; cd4 < hDS; cd4++)
        data[year][sex][agr][cd4] *= migration_pr[sex][agr];
  if (MODEL == 2)
    for (int sex = 0; sex < NG; sex++)
      for (int adb = 0; adb < hDB; adb++)
        for (int cd4 = 0; cd4 < hDS; cd4++)
          data_db[year][sex][adb][cd4] *= migration_pr[sex][adb];
}

void hivC::update_infection (const boost2D& new_infect) {
  zeroing(grad); // reset every time step
  boost2D infectAG = sumByAG(new_infect, ag_idx, hAG);
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hAG; agr++)
      for (int cd4 = 0; cd4 < hDS; cd4++)
        grad[sex][agr][cd4] += 
          p.cd4_initdist[sex][agr][cd4] * infectAG[sex][agr];
}

void hivC::grad_progress (const boost3D& mortality_rate) { // HIV gradient progress
  if (p.eppmod == 2)
    zeroing(grad); // reset every time step
  // remove cd4 stage progression (untreated)
  double nHup;
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hAG; agr++) {
      for (int cd4 = 0; cd4 < hDS - 1; cd4++) {
        nHup = data[year][sex][agr][cd4] * p.cd4_prog[sex][agr][cd4];
        _death[sex][agr][cd4] = 
          data[year][sex][agr][cd4] * mortality_rate[sex][agr][cd4];
        grad[sex][agr][cd4]   -= (nHup + _death[sex][agr][cd4]);
        grad[sex][agr][cd4+1] += nHup;
      }
      _death[sex][agr][hDS-1] = 
        data[year][sex][agr][hDS-1] * mortality_rate[sex][agr][hDS-1];
      grad[sex][agr][hDS-1] -= _death[sex][agr][hDS-1];
    }
  if (MODEL == 2) {
    zeroing(grad_db); // reset, this's the 1st time grad_db is used
    for (int sex = 0; sex < NG; sex++)
      for (int agr = 0; agr < hDB; agr++) {
        for (int cd4 = 0; cd4 < hDS - 1; cd4++) {
          nHup = data_db[year][sex][agr][cd4] * p.cd4_prog[sex][agr][cd4];
          _death_db[sex][agr][cd4] =
            data_db[year][sex][agr][cd4] * mortality_rate[sex][agr][cd4];
          grad_db[sex][agr][cd4]   -= (nHup + _death_db[sex][agr][cd4]);
          grad_db[sex][agr][cd4+1] += nHup;
        }
        _death_db[sex][agr][hDS-1] =
          data_db[year][sex][agr][hDS-1] * mortality_rate[sex][agr][hDS-1];
        grad_db[sex][agr][hDS-1] -= _death_db[sex][agr][hDS-1];
      }
  }
}

dvec hivC::eligible_for_art () { // this one does not depend on model state
                                 // can just do this in fp and have a new par
  dvec A(hDS);
  for (int i = p.artcd4elig_idx[year] - 1; i < hDS; ++i) 
    A[i] = 1;
  dvec B(hDS);
  for (int i = 2; i < hDS; ++i) 
    B[i] = p.who34percelig;
  boost1D C(extents[hDS]);
  C = add_to_each(C, p.specpop_percelig[year]);
  dvec D(hDS);
  for (int i = 0; i < hDS; ++i) 
    D[i] = 1 - (1 - A[i]) * (1 - B[i]) * (1 - C[i]);
  return D;
}

void hivC::distribute_artinit (boost3D& artinit, artC& artpop) {
    double debut_now, all_hivpop, pr_weight_db, n_artinit_db;
    boost3D artinit_db(extents[NG][hAG][hDS]);
    for (int sex = 0; sex < NG; sex++)
      for (int agr = 0; agr < hAG; agr++)
        for (int cd4 = 0; cd4 < hDS; cd4++) {
          debut_now =
            data_db[year][sex][agr][cd4] + DT * grad_db[sex][agr][cd4];
          all_hivpop =
            (data[year][sex][agr][cd4] + DT * grad[sex][agr][cd4]) + debut_now;
          if (artinit[sex][agr][cd4] > all_hivpop)
            artinit[sex][agr][cd4]   = all_hivpop;
          pr_weight_db               = debut_now / all_hivpop;
          n_artinit_db               = artinit[sex][agr][cd4] * pr_weight_db;
          artinit_db[sex][agr][cd4]  = n_artinit_db;
          artinit[sex][agr][cd4]    -= n_artinit_db;
          grad_db[sex][agr][cd4]    -= n_artinit_db / DT;
        }
    artpop.grad_db_init(artinit_db);
}

void hivC::add_grad_to_pop () { // this is over-extraction
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hAG; agr++)
      for (int cd4 = 0; cd4 < hDS; cd4++)
        data[year][sex][agr][cd4] += DT * grad[sex][agr][cd4];
  if (MODEL == 2)
    for (int sex = 0; sex < NG; sex++)
      for (int adb = 0; adb < hDB; adb++)
        for (int cd4 = 0; cd4 < hDS; cd4++)
          data_db[year][sex][adb][cd4] += DT * grad_db[sex][adb][cd4];
}

void hivC::adjust_pop (const boost2D& adj_prob) {
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hAG; agr++)
      for (int cd4 = 0; cd4 < hDS; cd4++)
        data[year][sex][agr][cd4] *= adj_prob[sex][agr];
  if (MODEL == 2)
    for (int sex = 0; sex < NG; sex++)
      for (int adb = 0; adb < hDB; adb++)
        for (int cd4 = 0; cd4 < hDS; cd4++)
          data_db[year][sex][adb][cd4] *= adj_prob[sex][adb];
}