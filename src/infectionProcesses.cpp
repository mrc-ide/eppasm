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

void popC::infect_spec (const hivC& hivpop, const artC& artpop, int time_step) {
  int ts = (year-1)/ DT + time_step,
      p_lo = p_age15to49_idx[0] - 1, h_lo = h_age15to49_idx[0] - 1;
  double dt_ii = 1 - DT * time_step, // transition of population in 1 year
         n_neg_mf = 0, n_pos_mf = 0, 
         n_pos_lo = 0, n_pos_up = 0, n_neg_lo = 0, n_neg_up = 0,
         n_hiv_lo = 0, n_art_lo = 0, n_hiv_up = 0, n_art_up = 0, art_ii = 0;
  update_active_pop_to(year);
  for (int sex = 0; sex < NG; sex++) {
    for (int age = p_lo; age < pAG_1549; age++) {
      n_neg_mf += data_active[hivn_idx][sex][age];
      n_pos_mf += data_active[hivp_idx][sex][age];
    }
    n_neg_lo += data_active[hivn_idx][sex][p_lo];
    n_neg_up += data_active[hivn_idx][sex][pAG_1549];
    n_pos_lo += data_active[hivp_idx][sex][p_lo];
    n_pos_up += data_active[hivp_idx][sex][pAG_1549];
    for (int agr = h_lo; agr < hAG_1549; agr++)
      for (int cd4 = 0; cd4 < hDS; cd4++)
        for (int dur = 0; dur < hTS; dur++) {
          art_ii += artpop.data[year][sex][agr][cd4][dur];
          n_art_lo += artpop.data[year][sex][h_lo][cd4][dur];
          n_art_up += artpop.data[year][sex][hAG_1549][cd4][dur];
        }
    for (int cd4 = 0; cd4 < hDS; cd4++) {
      n_hiv_lo += hivpop.data[year][sex][h_lo][cd4];
      n_hiv_up += hivpop.data[year][sex][hAG_1549][cd4];
    }
  }
  double hivn_ii = n_neg_mf - n_neg_lo * dt_ii + n_neg_up * dt_ii;
  double hivp_ii = n_pos_mf - n_pos_lo * dt_ii + n_pos_up * dt_ii;
  if ( n_hiv_lo + n_art_lo > 0) {
    for (int sex = 0; sex < NG; sex++) {
      double art_trans = 0, hiv_trans = 0;
      for (int cd4 = 0; cd4 < hDS; cd4++) {
        for (int dur = 0; dur < hTS; dur++)
          art_trans += artpop.data[year][sex][h_lo][cd4][dur];
        hiv_trans += hivpop.data[year][sex][h_lo][cd4];
      }
      art_ii -= ( data_active[hivp_idx][sex][p_lo] * 
                  art_trans / (hiv_trans + art_trans) ) * dt_ii;
    }
  }
  if ( n_hiv_up + n_art_up > 0) {
    for (int sex = 0; sex < NG; sex++) {
      double art_trans = 0, hiv_trans = 0;
      for (int cd4 = 0; cd4 < hDS; cd4++) {
        for (int dur = 0; dur < hTS; dur++)
          art_trans += artpop.data[year][sex][hAG_1549][cd4][dur];
        hiv_trans += hivpop.data[year][sex][hAG_1549][cd4];
      }
      art_ii += ( data_active[hivp_idx][sex][pAG_1549] * 
                  art_trans / (hiv_trans + art_trans) ) * dt_ii;
    }
  }
  double
  transm_prev = (hivp_ii - art_ii * (1 - p.relinfectART)) / (hivn_ii + hivp_ii);
  double w = (p.proj_steps[ts] == p.tsEpidemicStart) ? p.iota : 0.0;
  double inc_rate = rvec[ts] * transm_prev + w;

  // Incidence: male = negative / female negative * sexratio + male negative; 
  //          female = male * sexratio
  double n_neg_m = 0, n_neg_f = 0;
  for (int age = p_lo; age < pAG_1549; age++) {
    n_neg_m += data_active[hivn_idx][m_idx][age];
    n_neg_f += data_active[hivn_idx][f_idx][age];
  }
  double adj_sex = (n_neg_m + n_neg_f) / (n_neg_m + n_neg_f * p.incrr_sex[year]);
  double sex_inc[2] = {inc_rate * adj_sex, 
                       inc_rate * adj_sex * p.incrr_sex[year]};
  // New infections distributed by age: ratio age_i/ 25-29 age
  for (int sex = 0; sex < NG; sex++) {
    double n_neg = 0, n_neg_rr = 0, adj_age;
    for (int age = p_lo; age < pAG_1549; age++) {
      n_neg += data_active[hivn_idx][sex][age]; 
      n_neg_rr += p.incrr_age[year][sex][age] * data_active[hivn_idx][sex][age];
    }
    adj_age = sex_inc[sex] / ( n_neg_rr / n_neg );
    for (int age = 0; age < pAG; age++) {
      if (sex==m_idx) // age-specific incidence among circumcised men
        adj_age *= (1 - p.circ_incid_rr * p.circ_prop[year][age]);
      infections_[sex][age] = p.incrr_age[year][sex][age] * adj_age * 
        data_active[hivn_idx][sex][age];
    }
  }
  // saving
  incrate15to49_ts[ts] = inc_rate;
  prev_last = hivp_ii / (hivn_ii + hivp_ii);
  prev15to49_ts[ts] = prev_last;
}

void popC::infect_mix (int ii) {
  update_active_pop_to(year);
  int ts = (year-1)/DT + ii;
  boost2D transm_prev(extents[NG][pAG]);
  double N_hivp;
  for (int sex = 0; sex < NG; sex++)
    for (int age = 0; age < pAG; age++) {
      N_hivp = data_active[hivp_idx][sex][age];
      transm_prev[sex][age] = ((N_hivp * (1 - artcov[sex])) + 
                               (N_hivp * artcov[sex] * (1 - p.relinfectART)))/
                               (data_active[hivn_idx][sex][age] + N_hivp);
      }
  //+intervention effects and time epidemic start
  double w = (p.proj_steps[ts] == p.tsEpidemicStart) ? p.iota : 0.0;
  multiply_with_inplace(transm_prev, rvec[ts]);
  add_to_each_inplace(transm_prev, w);

  // sweep over sexual mixing matrices
  zeroing(infections_);
  for (int my_age = 0; my_age < pAG; ++my_age)
    for (int partner_age = 0; partner_age < pAG; ++partner_age) {
      infections_[m_idx][my_age] +=
        p.mat_m[partner_age][my_age] * transm_prev[f_idx][partner_age];
      infections_[f_idx][my_age] +=
        p.mat_f[partner_age][my_age] * transm_prev[m_idx][partner_age];
    }
  // if (exists("f_fun", fp)) // that fun
  //   ir = ir * fp.f_fun
  for (int sex = 0; sex < NG; sex++)
    for (int age = 0; age < pAG; age++)
      infections_[sex][age] *= data_active[hivn_idx][sex][age];
  // incrate15to49_ts_m.slice(ts) = ir_mf;
  // prev15to49_ts_m should use this one! now just store as below
  boost2D n_pos = data_active[ indices[hivp_idx][in(0, NG)][in(0, pAG)] ];
  prev15to49_ts[ts] = sumArray(n_pos) / sumArray(data_active);
  prev_last = prev15to49_ts[ts];
}

void popC::epp_disease_model_direct  (hivC& hivpop, artC& artpop) {
  int a_l, a_r;
  if (p.incidpopage) { // incidence for 15+ population
    a_l = p_age15plus_idx[0] - 1;
    a_r = pAG_15plus;
  }
  else { // incidence for 15 -49 population
    a_l = p_age15to49_idx[0] - 1;
    a_r = pAG_1549;
  }
  update_active_last_year();
  double n_m = 0, n_f = 0;
  for (int age = a_l; age < a_r; ++age) {
    n_m += active_last_year_[hivn_idx][m_idx][age];
    n_f += active_last_year_[hivn_idx][f_idx][age];
  }
  dvec sex_inc(NG); 
  sex_inc[m_idx] = (n_m + n_f) * p.incidinput[year] / 
                   (n_m + n_f  * p.incrr_sex[year]);
  sex_inc[f_idx] = (n_m + n_f) * p.incidinput[year] * p.incrr_sex[year] /
                   (n_m + n_f  * p.incrr_sex[year]);
  dvec ageinc(NG);
  for (int sex = 0; sex < NG; sex++) {
    double neg_sa = 0, inc_sa = 0;
    for (int age = a_l; age < a_r; age++) {
      neg_sa += active_last_year_[hivn_idx][sex][age];
      inc_sa += active_last_year_[hivn_idx][sex][age] * p.incrr_age[year][sex][age];
    }
    ageinc[sex] = inc_sa / neg_sa;
  }
  double new_infect = 0;
  for (int sex = 0; sex < NG; sex++)
    for (int age = 0; age < pAG; age++) {
      new_infect = p.incrr_age[year][sex][age] * ( sex_inc[sex] / ageinc[sex]) *
                   active_last_year_[hivn_idx][sex][age];
      infections[year][sex][age]      = new_infect;
      data[year][hivn_idx][sex][age] -= new_infect;
      data[year][hivp_idx][sex][age] += new_infect;
    }
  boost2D infect_agrp = 
    sumByAG(infections[ indices[year][in(0, NG)][in(0, pAG)]], ag_idx, hAG);
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hAG; agr++)
      for (int cd4 = 0; cd4 < hDS; cd4++)
        hivpop.data[year][sex][agr][cd4] +=
          p.cd4_initdist[sex][agr][cd4] * infect_agrp[sex][agr];
  for (int sex = 0; sex < NG; sex++)
    for (int age = p_age15to49_idx[0] - 1; age < pAG_1549; age++)
      incid15to49[year] += infections[year][sex][age];
}