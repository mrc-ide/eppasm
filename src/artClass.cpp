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
#include "Classes.h"

void artC::aging (mat ag_prob) {
  data.row(year) = data.row(year-1);
  field<cube> nARTup = sweepX34(data.row(year), ag_prob, hAG - 1);
  for (uword sex = 0; sex < 2; ++sex) {
    data.row(year)(sex).slices(0, hAG-2) -= nARTup(sex);
    data.row(year)(sex).slices(1, hAG-1) += nARTup(sex);
  }
  if (MODEL==2) {
    data_db.row(year) = data_db.row(year-1);
    nARTup = sweepX34(data_db.row(year), ag_prob, hAG - 1);
    for (uword sex = 0; sex < 2; ++sex) {
      data_db.row(year)(sex).slices(0, hAG-2) -= nARTup(sex);
      data_db.row(year)(sex).slices(1, hAG-1) += nARTup(sex);
    }
  }
}

void artC::add_entrants (vec artYesNo) {
  vec artYes = artYesNo.head(2);
  cube n_in = p.paedsurv_artcd4dist(year); // 3x7x2
  for (uword sex = 0; sex < 2; ++sex)
    n_in.slice(sex) = p.paedsurv_artcd4dist(year).slice(sex) * artYes(sex);
  if (MODEL==1)
    for (uword sex = 0; sex < 2; ++sex)
      data.row(year)(sex).slice(0) += n_in.slice(sex);
  if (MODEL==2) // add to virgin then debut
    for (uword sex = 0; sex < 2; ++sex)
      data_db.row(year)(sex).slice(0) += n_in.slice(sex);
}

void artC::sexual_debut () {
  field<cube> debut_art = sweepX34(data_db.row(year), p.db_pr, pDB);
  data.row(year)(0).slices(0, pDB-1) += debut_art(0);
  data.row(year)(1).slices(0, pDB-1) += debut_art(1);
  data_db.row(year)(0).slices(0, pDB-1) -= debut_art(0);
  data_db.row(year)(1).slices(0, pDB-1) -= debut_art(1);
}

void artC::deaths (mat survival_pr) {
  data.row(year) = sweepX34(data.row(year), survival_pr);
  if (MODEL==2)
    data_db.row(year) = sweepX34(data_db.row(year), survival_pr);
}

void artC::migration (mat migration_pr) {
  data.row(year) = sweepX34(data.row(year), migration_pr);
  if (MODEL==2)
    data_db.row(year) = sweepX34(data_db.row(year), migration_pr);
}

void artC::grad_progress () {
  gradART.for_each([](cube& X) { X.zeros(); } ); // reset gradient
  if (MODEL==2)
    gradART_db.for_each([](cube& X) { X.zeros(); } ); // reset gradient

  // progression and mortality (HARD CODED 6 months duration)
  for (uword sex = 0; sex < 2; ++sex)
    gradART(sex).rows(0, 1) -= 2.0 * data.row(year)(sex).rows(0, 1);
  for (uword sex = 0; sex < 2; ++sex)
    gradART(sex).rows(1, 2) += 2.0 * data.row(year)(sex).rows(0, 1);
  if (MODEL==2) {
    for (uword sex = 0; sex < 2; ++sex)
      gradART_db(sex).rows(0, 1) -= 2.0 * data_db.row(year)(sex).rows(0, 1);
    for (uword sex = 0; sex < 2; ++sex)
      gradART_db(sex).rows(1, 2) += 2.0 * data_db.row(year)(sex).rows(0, 1);
  }

  // ART mortality
  field<cube> art_mort_now = p.art_mort;
  for (uword sex = 0; sex < 2; ++sex) {
    for (uword agr = 0; agr < hAG; ++agr)
      art_mort_now(sex).slice(agr).each_col() %= p.artmx_timerr.col(year);
    gradART(sex) -= (art_mort_now(sex) % data.row(year)(sex));
    if (MODEL==2)
      gradART_db(sex) -= art_mort_now(sex) % data_db.row(year)(sex);
  }
}

void artC::art_dropout (hivC& hivpop) {
  double dropout_rate = p.art_dropout(year);
  cube cube_now(size(hivpop.grad)); // 7x9x2
  for (uword sex = 0; sex < 2; ++sex) {
    mat mat_now  = sum(data.row(year)(sex)); //7x9
    cube_now.slice(sex) = mat_now;
  }
  hivpop.grad += cube_now * dropout_rate;
  for (uword sex = 0; sex < 2; ++sex)
    gradART(sex) -= data.row(year)(sex) * dropout_rate;

  if (MODEL==2) {
    cube_now.zeros();
    for (uword sex = 0; sex < 2; ++sex) {
      mat mat_now  = sum(data_db.row(year)(sex));
      cube_now.slice(sex) = mat_now;
    }
    hivpop.grad_db += cube_now * dropout_rate;
    for (uword sex = 0; sex < 2; ++sex)
      gradART_db(sex) -= data.row(year)(sex)  * dropout_rate;
  }
}

vec artC::current_on_art () {
  vec art_curr = { accu(data.row(year)(0)), accu(data.row(year)(1)) }; 
  vec grad_curr = { accu(gradART(0)), accu(gradART(1)) };
  art_curr += (grad_curr * DT);
  if (MODEL==2) {
    vec art_curr_2 = {accu(data_db.row(year)(0)), accu(data_db.row(year)(1))};
    vec grad_curr_2 = {accu(gradART_db(0)), accu(gradART_db(1))};
    art_curr += (art_curr_2 + grad_curr_2 * DT);
  }
  return art_curr;
}

void artC::grad_init (cube artinit) {
  for (uword sex = 0; sex < 2; ++sex) {
    gradART(sex).row(0) += artinit.slice(sex) / DT;
    data.row(year)(sex) += DT * gradART(sex);
  }
}

void artC::grad_db_init (cube artinit_db) {
  for (uword sex = 0; sex < 2; ++sex) {
    gradART_db(sex).row(0) += artinit_db.slice(sex) / DT;
    data_db.row(year)(sex) += DT * gradART_db(sex);
  }
}

void artC::adjust_pop (mat adj_prob) {
  data.row(year) = sweepX34(data.row(year), adj_prob);
  if (MODEL==2)
    data_db.row(year) = sweepX34(data_db.row(year), adj_prob);
}

// set_data (FUN="+", x, TS=T, DS=T, AG=T, NG=T, YEAR=NULL) {
//     FUN <- match.fun(FUN)
//     if (is.null(YEAR))
//         data[TS, DS, AG, NG] <<- FUN(data[TS, DS, AG, NG], x)
//     else
//         data[TS, DS, AG, NG, YEAR] <<- FUN(data[TS, DS, AG, NG, YEAR], x)
// } 

// get (YEAR=NULL, TS=T, DS=T, AG=T, NG=T) {
//     'get does not change anything, just return the data, use set to change'
//     if (is.null(YEAR))
//         data[TS, DS, AG, NG, TRUE]
//     else
//         data[TS, DS, AG, NG, YEAR]
// }

// sweep_sex (FUN="*", x, year) {
//     data[,,,,year] <<- sweep(data[,,,,year], 3:4, x, FUN)
// }