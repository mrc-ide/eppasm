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

void artC::aging (const boost2D& ag_prob) {
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hAG; agr++)
      for (int cd4 = 0; cd4 < hDS; cd4++)
        for (int dur = 0; dur < hTS; dur++) {
          data[year][sex][agr][cd4][dur] = data[year-1][sex][agr][cd4][dur];
          if (MODEL == 2 && agr < hDB)
            data_db[year][sex][agr][cd4][dur] =
              data_db[year-1][sex][agr][cd4][dur];
        }
  double nARTup;
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hAG-1; agr++)
      for (int cd4 = 0; cd4 < hDS; cd4++)
        for (int dur = 0; dur < hTS; dur++) {
          nARTup = data[year-1][sex][agr][cd4][dur] * ag_prob[sex][agr];
          data[year][sex][agr][cd4][dur]   -= nARTup;
          data[year][sex][agr+1][cd4][dur] += nARTup;
          if (MODEL == 2 && agr < hDB - 1) {
            nARTup = data_db[year-1][sex][agr][cd4][dur] * ag_prob[sex][agr];
            data_db[year][sex][agr][cd4][dur]   -= nARTup;
            data_db[year][sex][agr+1][cd4][dur] += nARTup;
          }
        }
}

void artC::add_entrants (const dvec& artYesNo) {
  for (int sex = 0; sex < NG; sex++)
    for (int cd4 = 0; cd4 < hDS; cd4++)
      for (int dur = 0; dur < hTS; dur++) {
        double add = p.paedsurv_artcd4dist[year][sex][cd4][dur] * artYesNo[sex];
        if (MODEL == 1)
          data[year][sex][0][cd4][dur] += add;
        if (MODEL == 2) // add to virgin then debut
          data_db[year][sex][0][cd4][dur] += add;
      }
}

void artC::sexual_debut () {
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hDB; agr++)
      for (int cd4 = 0; cd4 < hDS; cd4++)
        for (int dur = 0; dur < hTS; dur++) {
          double n_db = data_db[year][sex][agr][cd4][dur] * p.db_pr[sex][agr];
          data[year][sex][agr][cd4][dur]    += n_db;
          data_db[year][sex][agr][cd4][dur] -= n_db;
        }
}

void artC::deaths (const boost2D& survival_pr) {
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hAG; agr++)
      for (int cd4 = 0; cd4 < hDS; cd4++)
        for (int dur = 0; dur < hTS; dur++) {
          data[year][sex][agr][cd4][dur] *= survival_pr[sex][agr];
          if (MODEL == 2 && agr < hDB)
            data_db[year][sex][agr][cd4][dur] *= survival_pr[sex][agr];
        }
}

void artC::migration (const boost2D& migration_pr) {
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hAG; agr++)
      for (int cd4 = 0; cd4 < hDS; cd4++)
        for (int dur = 0; dur < hTS; dur++) {
          data[year][sex][agr][cd4][dur] *= migration_pr[sex][agr];
          if (MODEL == 2 && agr < hDB)
            data_db[year][sex][agr][cd4][dur] *= migration_pr[sex][agr];
        }
}

void artC::grad_progress () {
  zeroing(gradART); // reset gradient
  if (MODEL == 2)
    zeroing(gradART_db); // reset gradient
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hAG; agr++)
      for (int cd4 = 0; cd4 < hDS; cd4++) {
        for (int dur = 0; dur < hTS - 1; dur++) {
          double art_up = 2.0 * data[year][sex][agr][cd4][dur];
          gradART[sex][agr][cd4][dur] -= (art_up + _death[sex][agr][cd4][dur]);
          gradART[sex][agr][cd4][dur+1] += art_up;
          if (MODEL == 2 && agr < hDB) {
            art_up = 2.0 * data_db[year][sex][agr][cd4][dur];
            gradART_db[sex][agr][cd4][dur] -=
              (art_up + _death_db[sex][agr][cd4][dur]);
            gradART_db[sex][agr][cd4][dur+1] += art_up;
          }
        }
        gradART[sex][agr][cd4][hTS-1] -= _death[sex][agr][cd4][hTS-1];
        if (MODEL == 2 && agr < hDB)
          gradART_db[sex][agr][cd4][hTS-1] -= _death_db[sex][agr][cd4][hTS-1];
      }
}

void artC::art_dropout (hivC& hivpop) {
  double n_dropout;
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hAG; agr++)
      for (int cd4 = 0; cd4 < hDS; cd4++)
        for (int dur = 0; dur < hTS; dur++) {
          n_dropout = data[year][sex][agr][cd4][dur] * p.art_dropout[year];
          hivpop.grad[sex][agr][cd4]  += n_dropout;
          gradART[sex][agr][cd4][dur] -= n_dropout;
          if (MODEL == 2 && agr < hDB) {
            n_dropout = data_db[year][sex][agr][cd4][dur] * p.art_dropout[year];
            hivpop.grad_db[sex][agr][cd4]  += n_dropout;
            gradART_db[sex][agr][cd4][dur] -= n_dropout;
          }
        }
}

dvec artC::current_on_art () {
  dvec art_curr(NG);
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hAG; agr++)
      for (int cd4 = 0; cd4 < hDS; cd4++)
        for (int dur = 0; dur < hTS; dur++) {
          art_curr[sex] += (data[year][sex][agr][cd4][dur] + 
                               gradART[sex][agr][cd4][dur] * DT);
          if (MODEL == 2 && agr < hDB)  // add art from virgin pop
            art_curr[sex] += (data_db[year][sex][agr][cd4][dur] + 
                                 gradART_db[sex][agr][cd4][dur] * DT);          
        }
  return art_curr;
}

void artC::grad_init (const boost3D& artinit) { // 7x9x2
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hAG; agr++)
      for (int cd4 = 0; cd4 < hDS; cd4++) {
        gradART[sex][agr][cd4][0] += artinit[sex][agr][cd4] / DT;
        for (int dur = 0; dur < hTS; dur++)
          data[year][sex][agr][cd4][dur] += DT * gradART[sex][agr][cd4][dur];
      }
}

void artC::grad_db_init (const boost3D& artinit_db) {
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hDB; agr++)
      for (int cd4 = 0; cd4 < hDS; cd4++) {
        gradART_db[sex][agr][cd4][0] += artinit_db[sex][agr][cd4] / DT;
        for (int dur = 0; dur < hTS; dur++)
          data_db[year][sex][agr][cd4][dur] += 
            DT * gradART_db[sex][agr][cd4][dur];
      }
}

void artC::adjust_pop (const boost2D& adj_prob) {
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hAG; agr++)
      for (int cd4 = 0; cd4 < hDS; cd4++)
        for (int dur = 0; dur < hTS; dur++) {
          data[year][sex][agr][cd4][dur] *= adj_prob[sex][agr];
          if (MODEL == 2 && agr < hDB)
            data_db[year][sex][agr][cd4][dur] *= adj_prob[sex][agr];
        }
}

void artC::count_death () {
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hAG; agr++)
      for (int cd4 = 0; cd4 < hDS; cd4++)
        for (int dur = 0; dur < hTS; dur++) {
          double x = p.art_mort[sex][agr][cd4][dur] * p.artmx_timerr[year][dur];
          _death[sex][agr][cd4][dur] = data[year][sex][agr][cd4][dur] * x;
          if (MODEL == 2 && agr < hDB)
            _death_db[sex][agr][cd4][dur] = data_db[year][sex][agr][cd4][dur] * x;
        }
}