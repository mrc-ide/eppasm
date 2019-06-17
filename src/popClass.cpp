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

void popC::update_active_pop_to (int when) {
  data_active = data[ indices[when][in(0,pDS)][in(0,NG)][in(0, pAG)] ];
  if (MODEL==2)
    for (int ds = 0; ds < pDS; ds++)
      for (int sex = 0; sex < NG; sex++)
        for (int age = 0; age < pDB; age++)
          data_active[ds][sex][age] -= data_db[when][ds][sex][age];
}

boost3D popC::get_active_pop_in (int when) {
  boost3D out = data[ indices[when][in(0,pDS)][in(0,NG)][in(0, pAG)] ];
  if (MODEL==2)
    for (int ds = 0; ds < pDS; ds++)
      for (int sex = 0; sex < NG; sex++)
        for (int age = 0; age < pDB; age++)
          out[ds][sex][age] -= data_db[when][ds][sex][age];
  return out;
}

void popC::aging () { // open ended
  for (int ds = 0; ds < pDS; ds++)
    for (int sex = 0; sex < NG; sex++) {
      for (int age = 1; age < pAG; age++) 
        data[year][ds][sex][age]  = data[year-1][ds][sex][age-1];
      data[year][ds][sex][pAG-1] += data[year-1][ds][sex][pAG-1];
    }
  if (MODEL==2) // age the debut pop
    for (int ds = 0; ds < pDS; ds++)
      for (int sex = 0; sex < NG; sex++)
        for (int age = 1; age < pDB; age++)
          data_db[year][ds][sex][age] = data_db[year-1][ds][sex][age-1];
}

void popC::add_entrants () { // Add lagged births into youngest age group
  double healthy[2], positiv[2];
  // if (exists("popadjust", where=p) & p.popadjust) {
  if (p.popadjust) {
    if (MODEL==0)
      for (int sex = 0; sex < NG; ++sex)
        healthy[sex] = p.entrantpop[year-1][sex];
    else {
      for (int sex = 0; sex < NG; ++sex) {
        healthy[sex] = p.entrantpop[year-1][sex] * (1-p.entrantprev[year][sex]);
        positiv[sex] = p.entrantpop[year-1][sex] *    p.entrantprev[year][sex];
      }
    }
  }
  else {
    if (MODEL==0)
      for (int sex = 0; sex < NG; ++sex)
        healthy[sex] = birthslag[year-1][sex] * p.cumsurv[year-1][sex] / 
                       p.paedsurv_lag[year-1] + p.cumnetmigr[year-1][sex];
    else {
      for (int sex = 0; sex < NG; ++sex) {
        healthy[sex] = birthslag[year-1][sex] * p.cumsurv[year-1][sex] * 
          (1 - p.entrantprev[year][sex] / p.paedsurv_lag[year-1]) +
          p.cumnetmigr[year-1][sex] * (1 - pregprevlag[year-1] * p.netmig_hivprob);
        positiv[sex] = (birthslag[year-1][sex] * p.cumsurv[year-1][sex] + 
          p.cumnetmigr[year-1][sex]) * p.entrantprev[year][sex];
      }
    }
  }
  // save and update pop
  double sum_p = 0, sum_h = 0;
  for (int sex = 0; sex < NG; ++sex) {
    data[year][hivn_idx][sex][0] = healthy[sex];
    if (MODEL != 0) {
      data[year][hivp_idx][sex][0] = positiv[sex];
      sum_p += positiv[sex]; sum_h += healthy[sex];
      hivp_entrants_out[year][sex] = positiv[sex];
    }
    if (MODEL==2) { // add to virgin to record
      data_db[year][hivn_idx][sex][0] = healthy[sex];
      data_db[year][hivp_idx][sex][0] = positiv[sex];
    }
  }
  
  if (MODEL!=0)
    entrantprev[year] = sum_p / (sum_p + sum_h);
}

void popC::sexual_debut () {
  for (int ds = 0; ds < pDS; ds++)
    for (int sex = 0; sex < NG; sex++)
      for (int age = 0; age < pDB; age++)
        data_db[year][ds][sex][age] *= (1 - p.db_pr[sex][age]);
}

boost2D popC::hiv_aging_prob () {
  boost2D pop_hivp =
    sumByAG(data[indices[year-1][hivp_idx][in(0, NG)][in(0, pAG)]], ag_idx, hAG);
  boost2D out(extents[NG][hAG]);
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hAG; agr++)
      out[sex][agr] = data[year-1][hivp_idx][sex][aglast_idx[agr]-1] / 
                      pop_hivp[sex][agr];
  replace_na_with(out, 0.0); // inplace
  return out;
}

dvec popC::entrant_art () { // return these for updating HIV and ART pop
  dvec out(4);
  for (int sex = 0; sex < NG; ++sex) {
    out[sex]   = data[year][hivp_idx][sex][0] *      p.entrantartcov[year][sex];
    out[sex+2] = data[year][hivp_idx][sex][0] * (1 - p.entrantartcov[year][sex]);
  }
  return out; // 1:2 ART+, 3:4 ART-
}

void popC::deaths () {
  boost2D p_n =
    sumByAG(data[indices[year][hivp_idx][in(0, NG)][in(0, pAG)]], ag_idx, hAG);
  boost3D death_now(extents[pDS][NG][pAG]);
  for (int ds = 0; ds < pDS; ds++)
    for (int sex = 0; sex < NG; sex++)
      for (int age = 0; age < pAG; age++) {
        death_now[ds][sex][age] =
          data[year][ds][sex][age] * (1 - p.Sx[year][sex][age]);
        data[year][ds][sex][age] *= p.Sx[year][sex][age];
      }

  boost2D d_n =
    sumByAG(death_now[indices[hivp_idx][in(0, NG)][in(0, pAG)]], ag_idx, hAG);
  for (int sex = 0; sex < NG; sex++)
    for (int age = 0; age < hAG; age++)
      hiv_sx_prob[sex][age] = 1 - (d_n[sex][age] / p_n[sex][age]);
  replace_na_with(hiv_sx_prob, 0);

  if (MODEL == 2)
    for (int ds = 0; ds < pDS; ds++)
      for (int sex = 0; sex < NG; sex++)
        for (int age = 0; age < pDB; age++)
          data_db[year][ds][sex][age] *= p.Sx[year][sex][age];

  for (int sex = 0; sex < NG; ++sex)
    for (int age = 0; age < pAG; ++age)
      natdeaths[year][sex][age] = 
        death_now[0][sex][age] + death_now[1][sex][age];
}

void popC::migration () {
  boost2D netmigsurv(extents[NG][pAG]);
  for (int sex = 0; sex < NG; sex++)
    for (int age = 0; age < pAG; age++)
      netmigsurv[sex][age] = 
        p.netmigr[year][sex][age] * (1 + p.Sx[year][sex][age]) / 2;
  boost2D mr_prob(extents[NG][pAG]);
  boost2D nH_mr(extents[NG][pAG]);
  for (int sex = 0; sex < NG; ++sex)
    for (int age = 0; age < pAG; ++age) {
      mr_prob[sex][age] = 1 + netmigsurv[sex][age] /
        ( data[year][hivn_idx][sex][age] + data[year][hivp_idx][sex][age] );
      nH_mr[sex][age] = mr_prob[sex][age] * data[year][hivp_idx][sex][age];
    }
  boost2D nH_mr_agr = sumByAG(nH_mr, ag_idx, hAG);
  boost2D n_pop_agr = 
    sumByAG(data[indices[year][hivp_idx][in(0, NG)][in(0, pAG)]], ag_idx, hAG);
  for (int sex = 0; sex < NG; sex++)
    for (int age = 0; age < hAG; age++)
      hiv_mr_prob[sex][age] = nH_mr_agr[sex][age] / n_pop_agr[sex][age];
  replace_na_with(hiv_mr_prob, 0);
  if (MODEL==2)
    for (int ds = 0; ds < pDS; ds++)
      for (int sex = 0; sex < NG; sex++)
        for (int age = 0; age < pDB; age++)
          data_db[year][ds][sex][age] *= mr_prob[sex][age];
  for (int ds = 0; ds < pDS; ds++)
    for (int sex = 0; sex < NG; sex++)
      for (int age = 0; age < pAG; age++)
        data[year][ds][sex][age] *= mr_prob[sex][age];
}

void popC::update_fertile () { // only on active pop
  update_active_pop_to(year);
  boost3D last_year = get_active_pop_in(year-1);
  for (int age = 0; age < pAG_FERT; age++)
    birth_age[age] =
      ((data_active[hivp_idx][f_idx][age] + data_active[hivn_idx][f_idx][age] +
          last_year[hivp_idx][f_idx][age] + last_year[hivn_idx][f_idx][age]) /
      2) * p.asfr[year][age];
  ivec sub_id(ag_idx.begin() + p_fert_idx[0] - 1, ag_idx.begin() + pAG_FERT);
  birth_agrp = sumByAG(birth_age, sub_id, hAG_FERT);
  double n_births = sum_vector(birth_agrp);
  if ( (year + AGE_START) <= (PROJ_YEARS - 1) )
    for (int sex = 0; sex < NG; ++sex)
      birthslag[year + AGE_START - 1][sex] = p.srb[year][sex] * n_births;
}

void popC::adjust_pop () {
  for (int sex = 0; sex < NG; sex++)
    for (int age = 0; age < pAG; age++)
      popadjust[year][sex][age] = p.targetpop[year][sex][age] / 
        ( data[year][hivn_idx][sex][age] + data[year][hivp_idx][sex][age] );
  if (MODEL!=0) {
    boost2D n_adjust(extents[NG][pAG]);
    for (int sex = 0; sex < NG; sex++)
      for (int age = 0; age < pAG; age++)
        n_adjust[sex][age] = 
          popadjust[year][sex][age] * data[year][hivp_idx][sex][age];
    boost2D adj_num = sumByAG(n_adjust, ag_idx, hAG);
    boost2D adj_dem = 
      sumByAG(data[year][indices[hivp_idx][in(0, NG)][in(0,pAG)]], ag_idx, hAG);
    for (int sex = 0; sex < NG; sex++)
      for (int age = 0; age < hAG; age++)
        adj_prob[sex][age] = adj_num[sex][age] / adj_dem[sex][age];
    replace_na_with(adj_prob, 0);
  }
  for (int ds = 0; ds < pDS; ds++)
    for (int sex = 0; sex < NG; sex++) {
      for (int age = 0; age < pAG; age++)
        data[year][ds][sex][age] *= popadjust[year][sex][age];
      if (MODEL==2)
        for (int adb = 0; adb < pDB; adb++)
          data_db[year][ds][sex][adb] *= popadjust[year][sex][adb];
    }
}

void popC::cal_prev_pregant (const hivC& hivpop, const artC& artpop) { // only on active pop
  dvec n_mean(pAG_FERT); // 1 X 35
  update_active_pop_to(year); boost3D last_year = get_active_pop_in(year-1);
  for (int age = 0; age < pAG_FERT; ++age)
    n_mean[age] = (last_year[hivn_idx][f_idx][age] +
                   data_active[hivn_idx][f_idx][age]) / 2;
  ivec sub_id(ag_idx.begin() + p_fert_idx[0] - 1, ag_idx.begin() + pAG_FERT);
  dvec hivn = sumByAG(n_mean, sub_id, hAG_FERT); // 1 x 8
  double frp = 0, fra = 0, frap = 0;
  for (int agr = 0; agr < hAG_FERT; ++agr) {
    for (int cd4 = 0; cd4 < hDS; ++cd4) {
      frp += (hivpop.data[year-1][f_idx][agr][cd4] +
              hivpop.data[year  ][f_idx][agr][cd4] ) / 2 * 
              p.frr_cd4[year][agr][cd4];
      for (int dur = 0; dur < hTS; dur++)
        fra += (artpop.data[year-1][f_idx][agr][cd4][dur] +
                artpop.data[year  ][f_idx][agr][cd4][dur]) / 2 *
                p.frr_art[year][agr][cd4][dur];
    }
    frap += birth_agrp[agr] * (1 - hivn[agr] / (hivn[agr] + frp + fra));
    frp = 0; fra = 0;
  }
  pregprevlag[year + AGE_START - 1] = frap / sum_vector(birth_age);
}

void popC::save_prev_n_inc () {
  boost3D last_year = get_active_pop_in(year-1);
  double n_positive = 0, everyone_now = 0, s_previous = 0;
  for (int sex = 0; sex < NG; sex++)
    for (int age = p_age15to49_idx[0] - 1; age < pAG_1549; age++) {
      n_positive += data[year][hivp_idx][sex][age]; // +virgin
      s_previous += last_year[hivn_idx][sex][age]; // susceptible -virgin
      for (int ds = 0; ds < pDS; ds++)
        everyone_now += data[year][ds][sex][age];
    }
  prev15to49[year] = n_positive / everyone_now;
  incid15to49[year] /= s_previous;
  // prev(year) = accu(data_all.slice(hivp_idx)) / accu(data_all);
  // incid(year) = incid15to49(year) / accu(data(year-1).slice(hivn_idx)); // toBfixed
}

boost2D popC::infect_mix (int ii) {
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
  boost2D infections_ts(extents[NG][pAG]);
  for (int my_age = 0; my_age < pAG; ++my_age)
    for (int partner_age = 0; partner_age < pAG; ++partner_age) {
      infections_ts[m_idx][my_age] +=
        p.mat_m[partner_age][my_age] * transm_prev[f_idx][partner_age];
      infections_ts[f_idx][my_age] +=
        p.mat_f[partner_age][my_age] * transm_prev[m_idx][partner_age];
    }
  // if (exists("f_fun", fp)) // that fun
  //   ir = ir * fp.f_fun
  for (int sex = 0; sex < NG; sex++)
    for (int age = 0; age < pAG; age++)
      infections_ts[sex][age] *= data_active[hivn_idx][sex][age];
  // incrate15to49_ts_m.slice(ts) = ir_mf;
  // prev15to49_ts_m should use this one! now just store as below
  boost2D n_pos = data_active[ indices[hivp_idx][in(0, NG)][in(0, pAG)] ];
  prev15to49_ts[ts] = sumArray(n_pos) / sumArray(data_active);
  prev_last = prev15to49_ts[ts];
  return infections_ts;
}

boost2D popC::infect_spec (const hivC& hivpop, const artC& artpop, int time_step) {
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
  boost2D infections_ts(extents[NG][pAG]);
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
      infections_ts[sex][age] = p.incrr_age[year][sex][age] * adj_age * 
        data_active[hivn_idx][sex][age];
    }
  }
  // saving
  incrate15to49_ts[ts] = inc_rate;
  prev_last = hivp_ii / (hivn_ii + hivp_ii);
  prev15to49_ts[ts] = prev_last;
  return infections_ts;
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
  boost3D last_year = get_active_pop_in(year-1);
  double n_m = 0, n_f = 0;
  for (int age = a_l; age < a_r; ++age) {
    n_m += last_year[hivn_idx][m_idx][age];
    n_f += last_year[hivn_idx][f_idx][age];
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
      neg_sa += last_year[hivn_idx][sex][age];
      inc_sa += last_year[hivn_idx][sex][age] * p.incrr_age[year][sex][age];
    }
    ageinc[sex] = inc_sa / neg_sa;
  }
  double new_infect = 0;
  for (int sex = 0; sex < NG; sex++)
    for (int age = 0; age < pAG; age++) {
      new_infect = p.incrr_age[year][sex][age] * ( sex_inc[sex] / ageinc[sex]) *
                   last_year[hivn_idx][sex][age];
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

double popC::calc_rtrend_rt (int ts, double time_step) {
  double rveclast = rvec[ts-1];
  Rf_error("K not write for r_trend, no fp template to write");
  // double dtii =  1 - DT * (time_step - 1);
  // int a_l = p_age15to49_idx[0] -1, a_r = a_l + p_age15to49_idx.num_elements();
  // boost2D A = data[ indices[year][hivn_idx][in(0, NG)][in(a_l, a_r)] ];
  // boost1D B = data[ indices[year][hivn_idx][in(0, NG)][a_l] ];
  // boost1D C = data[ indices[year][hivn_idx][in(0, NG)][a_r] ];
  // double hivn_ii = sumArray(A) - sumArray(B) * dtii + sumArray(C) * dtii;
  // A = data[ indices[year][hivp_idx][in(0, NG)][in(a_l, a_r)] ];
  // B = data[ indices[year][hivp_idx][in(0, NG)][a_l] ];
  // C = data[ indices[year][hivp_idx][in(0, NG)][a_r] ];
  // double hivp_ii = sumArray(A) - sumArray(B) * dtii + sumArray(C) * dtii;
  // double prevcurr = hivp_ii / (hivn_ii + hivp_ii);
  // double t_ii     = p.proj_steps[ts];
  // if (t_ii > p.tsEpidemicStart) {
    // par = p.rtrend;
    // if (t_ii < par.tStabilize)
    //   gamma_t = 0
    // else 
    //   gamma_t = (prevcurr-prevlast) * (t_ii - par$tStabilize) / (DT * prevlast)
    // logr.diff =  par$beta[2] * (par$beta[1] - rveclast) +
    //              par$beta[3] * prevlast + 
    //              par$beta[4] * gamma_t
    // return(exp(log(rveclast) + logr.diff))
  // }
  // else
    // return(p.rtrend$r0)
  return rveclast; // remove this
}

void popC::update_rvec (double time_step) {
  int ts = (year-1) / DT + time_step;
  // if (p.eppmod %in% c("rtrend", "rtrend_rw")) <<- need to be fixed in FP
  if (p.eppmod == 1) // rtrend see prepare_fp_for_Cpp
    rvec[ts] = calc_rtrend_rt(ts, time_step);
  else
    rvec[ts] = p.rvec[ts];
}

void popC::update_infection (const boost2D& infect) {
  double n_infect;
  for (int sex = 0; sex < NG; sex++)
    for (int age = 0; age < pAG; age++) {
      n_infect = infect[sex][age] * DT;
      data[year][hivn_idx][sex][age] -= n_infect;
      data[year][hivp_idx][sex][age] += n_infect;
      infections[year][sex][age]     += n_infect;
      if ( (age >=  p_age15to49_idx[0] - 1) & (age < pAG_1549) )
        incid15to49[year] += n_infect;
    }
}

void popC::remove_hiv_death (const boost3D& cd4_mx,
                             const hivC& hivpop, const artC& artpop) {
  boost2D dbyAG(extents[NG][hAG]);  // death by age group
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hAG; agr++) {
      double hivD = 0, artD = 0;
      for (int cd4 = 0; cd4 < hDS; cd4++) {
        hivD += hivpop._death[sex][agr][cd4];
        if (MODEL==2) // add hiv deaths from inactive population
          hivD += hivpop._death_db[sex][agr][cd4];
        for (int dur = 0; dur < hTS; dur++) {
          artD += artpop._death[sex][agr][cd4][dur];
          if (MODEL==2) // add art deaths from inactive population
            artD += artpop._death_db[sex][agr][cd4][dur];
        }
      }
      dbyAG[sex][agr] = DT * (hivD + artD); // deaths by single-year
    }
  boost2D n_hiv =
    sumByAG(data[indices[year][hivp_idx][in(0,NG)][in(0,pAG)]], ag_idx, hAG);
  double hiv_mx;
  for (int sex = 0; sex < NG; sex++) {
    int age = 0;
    for (int agr = 0; agr < hAG; agr++) {
      if (n_hiv[sex][agr] != 0) {
        hiv_mx = dbyAG[sex][agr] / n_hiv[sex][agr];
        for (int i = 0; i < h_ag_span[agr]; ++i) {
          hivdeaths[year][sex][age] += data[year][hivp_idx][sex][age] * hiv_mx;
          data[year][hivp_idx][sex][age] *= (1 - hiv_mx);
          if (age < pDB)
            data_db[year][hivp_idx][sex][age] *= (1 - hiv_mx);
          age++;
        }
      } else 
        age += h_ag_span[agr] - 1;
    }
  }
}

void popC::update_preg (boost3D& art_elig,
                        const hivC& hivpop, const artC& artpop) {
  int h_lo = h_fert_idx[0] - 1, // 0 9
      p_lo = p_fert_idx[0] - 1; // 0 35
  update_active_pop_to(year);
  boost1I sub_id = ag_idx[indices[in(p_lo, pAG_FERT)]];
  boost1D hivn =
    sumByAG(data_active[indices[hivn_idx][f_idx][in(p_lo, pAG_FERT)]],
            sub_id, hAG_FERT); // 1 x 8
  dvec all_art(hAG_FERT);
  for (int agr = h_lo; agr < hAG_FERT; agr++)
    for (int cd4 = 0; cd4 < hDS; cd4++) {
      for (int dur = 0; dur < hTS; dur++)
        all_art[agr] += artpop.data[year][f_idx][agr][cd4][dur] * 
          p.frr_art[year][agr][cd4][dur];
      all_art[agr] +=
        hivpop.data[year][f_idx][agr][cd4] * p.frr_cd4[year][agr][cd4];
    }
  for (int agr = h_lo; agr < hAG_FERT; agr++)
    for (int cd4 = 0; cd4 < (p.artcd4elig_idx[year] - 1); cd4++)
      art_elig[f_idx][agr][cd4] += 
        hivpop.data[year][f_idx][agr][cd4] * p.frr_cd4[year][agr][cd4] *
          (birth_agrp[agr] / (hivn[agr] + all_art[agr]));
}

dvec popC::art_initiate (const dvec& art_curr, const boost3D& art_elig,
                         int time_step) {
  dvec out(NG);
  int year_w = (DT * (time_step + 1) < 0.5) ? 0 : 1;
  double trans = DT * (time_step + 1) + 0.5 - year_w;
  int year_l = year - (2 - year_w), year_r = year - (1 - year_w);
  for (int sex = 0; sex < NG; ++sex) {
    if( (p.art15plus_isperc[year_l][sex] == 0) & 
        (p.art15plus_isperc[year_r][sex] == 0) ) { // both number
      out[sex] = p.art15plus_num[year_l][sex] * (1 - trans) + 
                 p.art15plus_num[year_r][sex] *      trans;
      if (MIX) {
        boost2D art_elig_sex = art_elig[ indices[sex][in(0, hAG)][in(0, hDS)] ];
        double cov = out[sex] / (sumArray(art_elig_sex) + art_curr[sex]);
        artcov[sex] = (cov <= 1) ? cov : 1;
      }
    } 
    else if ( (p.art15plus_isperc[year_l][sex] == 1) &
              (p.art15plus_isperc[year_r][sex] == 1) ) { // both percentage
      double cov = p.art15plus_num[year_l][sex] * (1 - trans) +
                   p.art15plus_num[year_r][sex] *      trans;
      if (MIX)
        artcov[sex] = (cov <= 1) ? cov : 1;
      boost2D art_elig_sex = art_elig[ indices[sex][in(0, hAG)][in(0, hDS)] ];
      out[sex] = cov * (sumArray(art_elig_sex) + art_curr[sex]);
    }
    else if ( (p.art15plus_isperc[year_l][sex] == 0) & 
              (p.art15plus_isperc[year_r][sex] == 1) ) { // transition number to percentage
      boost2D art_elig_sex = art_elig[ indices[sex][in(0, hAG)][in(0, hDS)] ];
      double actual_cov =
        art_curr[sex] / (sumArray(art_elig_sex) + art_curr[sex]);
      double diff_cov = p.art15plus_num[year_r][sex] - actual_cov;
      double cov = actual_cov + diff_cov * DT / (0.5 + year_w - DT * time_step);
      if (MIX)
        artcov[sex] = (cov < 1) ? cov : 1;
      out[sex] = cov * ( sumArray(art_elig_sex) + art_curr[sex] );
    }
  }
  return out;
} 

// calculate ART initiation distribution
boost3D popC::art_distribute (const boost3D& art_elig, const dvec& art_need) {
  if (!p.med_cd4init_input[year]) {
    if (p.art_alloc_method == 4L) { // by lowest CD4
      // Calculate proportion to be initiated in each CD4 category
      dvec init_pr(NG);
      boost3D art_real(extents[NG][hAG][hDS]);
      for (int cd4 = hDS - 1; cd4 > 0; --cd4) { //6->0
        dvec elig_hm(NG);
        for (int sex = 0; sex < NG; sex++)
          for (int age = 0; age < hAG; age++)
            elig_hm[sex] += art_elig[sex][age][cd4];
        if ( elig_hm[m_idx] == 0 & elig_hm[f_idx] == 0 )
          init_pr = elig_hm;
        else {
          double x;
          for (int sex = 0; sex < NG; ++sex) {
            x = art_need[sex] / elig_hm[sex];
            init_pr[sex] = ( (x > 1) | std::isnan(x) | std::isinf(x)) ? 1 : x;
          }
        }
        boost2D current_m = art_elig[ indices[in(0, NG)][in(0, hAG)][cd4] ];
        for (int sex = 0; sex < NG; ++sex)
          for (int agr = 0; agr < hAG; ++agr)
            art_real[sex][agr][cd4] = current_m[sex][agr] * init_pr[sex];
      }
      return art_real;
    } 
    else { // Spectrum Manual p168--p169, 
      int A = h_age15plus_idx[0] - 1;
      dvec artX(NG), artY(NG);
      for (int sex = 0; sex < NG; sex++)
        for (int agr = A; agr < hAG_15plus; agr++)
          for (int cd4 = 0; cd4 < hDS; cd4++) {
            artX[sex] += art_elig[sex][agr][cd4] * p.cd4_mort[sex][agr][cd4];
            artY[sex] += art_elig[sex][agr][cd4];
          }
      boost3D art_real(extents[NG][hAG_15plus][hDS]);
      double xx;
      for (int sex = 0; sex < NG; sex++)
        for (int agr = A; agr < hAG_15plus; agr++)
          for (int cd4 = 0; cd4 < hDS; cd4++) {
            xx = (p.cd4_mort[sex][agr][cd4] / artX[sex] * p.art_alloc_mxweight +
                  ((1 - p.art_alloc_mxweight) / artY[sex]) ) *
                art_elig[sex][agr][cd4] * art_need[sex];
            art_real[sex][agr][cd4] =
              (xx > art_elig[sex][agr][cd4]) ? art_elig[sex][agr][cd4] : xx;
          }
      return art_real;
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
        elig_below[sex] += art_elig[sex][agr][j] * pr_below;
    dvec A(NG);
    if (j < (hDS - 1)) {
      for (int sex = 0; sex < NG; sex++) {
        for (int agr = 0; agr < hAG; agr++)
          for (int cd4 = j+1; cd4 < hDS; cd4++)
            A[sex] += art_elig[sex][agr][cd4];
        elig_below[sex] += A[sex];
      }
    }
    dvec elig_above(NG);
    dvec B(NG);
    for (int sex = 0; sex < NG; sex++) {
      for (int agr = 0; agr < hAG; agr++)
          B[sex] += art_elig[sex][agr][j] * (1.0 - pr_below);
      elig_above[sex] += B[sex];
    }
    if (j > 1) {
      dvec C(NG);
      for (int sex = 0; sex < NG; sex++) {
        for (int agr = 0; agr < hAG; agr++)
          for (int cd4 = 0; cd4 < j-1; cd4++)
            C[sex] = art_elig[sex][agr][cd4];
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
    boost3D art_real(extents[NG][hAG][hDS]);
    if (j < (hDS - 1)) {
      for (int sex = 0; sex < NG; sex++)
        for (int agr = 0; agr < hAG; agr++)
          for (int cd4 = j + 1; cd4 < hDS; cd4++)
            art_real[sex][agr][cd4] = 
              art_elig[sex][agr][cd4] * initpr_below[sex];
    }
    for (int sex = 0; sex < NG; sex++)
      for (int agr = 0; agr < hAG; agr++)
        art_real[sex][agr][j] =
          art_elig[sex][agr][j] * initpr_medcat[sex];
    if (j > 0) {
      for (int sex = 0; sex < NG; sex++)
        for (int agr = 0; agr < hAG; agr++)
          for (int cd4 = 0; cd4 < j - 1; cd4++)
            art_real[sex][agr][cd4] = 
              art_elig[sex][agr][cd4] * initpr_above[sex];
    }
    return art_real;
  }
}

boost3D popC::scale_cd4_mort (hivC& hivpop, artC& artpop) {
  if (p.scale_cd4_mort) {
    boost3D cd4mx(extents[NG][hAG][hDS]);
    double num, den = 0;
    for (int sex = 0; sex < NG; sex++)
      for (int agr = 0; agr < hAG; agr++)
        for (int cd4 = 0; cd4 < hDS; cd4++) {
          num = hivpop.data[year][sex][agr][cd4] + 
                hivpop.data_db[year][sex][agr][cd4];
          for (int dur = 0; dur < hTS; dur++)
            den += artpop.data[year][sex][agr][cd4][dur] +
                   artpop.data_db[year][sex][agr][cd4][dur];
          num = (num + den == 0.0) ? 1 : num / (num + den);
          cd4mx[sex][agr][cd4] = num * p.cd4_mort[sex][agr][cd4];
          den = 0;
        }
    return cd4mx;
  } 
  else
    return p.cd4_mort;
}

void popC::finalize (hivC& hivpop, artC& artpop) {
  int np = 0;

  SEXP pop_sexp_dim = PROTECT(NEW_INTEGER(4)); ++np;
  INTEGER(pop_sexp_dim)[0] = pAG;
  INTEGER(pop_sexp_dim)[1] = NG;
  INTEGER(pop_sexp_dim)[2] = pDS;
  INTEGER(pop_sexp_dim)[3] = PROJ_YEARS;
  SET_DIM(pop_sexp, pop_sexp_dim);

  if (MODEL!=0) {
    SEXP age_sex_year_dim = PROTECT(NEW_INTEGER(3)); ++np;
    INTEGER(age_sex_year_dim)[0] = pAG;
    INTEGER(age_sex_year_dim)[1] = NG;
    INTEGER(age_sex_year_dim)[2] = PROJ_YEARS;
    SET_DIM(infections_sexp, age_sex_year_dim);
    SET_DIM(hivdeaths_sexp, age_sex_year_dim);
    SET_DIM(natdeaths_sexp, age_sex_year_dim);
    SET_DIM(popadjust_sexp, age_sex_year_dim);

    SEXP hiv_sexp_dim = PROTECT(NEW_INTEGER(4)); ++np;
    INTEGER(hiv_sexp_dim)[0] = hDS;
    INTEGER(hiv_sexp_dim)[1] = hAG;
    INTEGER(hiv_sexp_dim)[2] = NG;
    INTEGER(hiv_sexp_dim)[3] = PROJ_YEARS;
    SET_DIM(hivpop.hiv_sexp, hiv_sexp_dim);

    SEXP art_sexp_dim = PROTECT(NEW_INTEGER(5)); ++np;
    INTEGER(art_sexp_dim)[0] = hTS;
    INTEGER(art_sexp_dim)[1] = hDS;
    INTEGER(art_sexp_dim)[2] = hAG;
    INTEGER(art_sexp_dim)[3] = NG;
    INTEGER(art_sexp_dim)[4] = PROJ_YEARS;
    SET_DIM(artpop.art_sexp, art_sexp_dim);

    SET_ATTR(pop_sexp, Rf_install("infections"), infections_sexp);
    SET_ATTR(pop_sexp, Rf_install("hivdeaths"), hivdeaths_sexp);
    SET_ATTR(pop_sexp, Rf_install("natdeaths"), natdeaths_sexp);
    SET_ATTR(pop_sexp, Rf_install("popadjust"), popadjust_sexp);
    SET_ATTR(pop_sexp, Rf_install("pregprevlag"), pregprevlag_sexp);
    SET_ATTR(pop_sexp, Rf_install("rvec_ts"), rvec_sexp);
    SET_ATTR(pop_sexp, Rf_install("prev15to49"), prev15to49_sexp);
    SET_ATTR(pop_sexp, Rf_install("incid15to49"), incid15to49_sexp);
    SET_ATTR(pop_sexp, Rf_install("entrantprev"), entrantprev_sexp);
    SET_ATTR(pop_sexp, Rf_install("incrate15to49_ts"), inci15to49_ts_sexp);
    SET_ATTR(pop_sexp, Rf_install("prev15to49_ts_sexp"), prev15to49_ts_sexp);
    SET_ATTR(pop_sexp, Rf_install("artpop"), artpop.art_sexp);
    SET_ATTR(pop_sexp, Rf_install("hivpop"), hivpop.hiv_sexp);
    SET_CLASS(pop_sexp, Rf_mkString("spec"));
    if (MODEL==2) {
      SEXP pop_db_sexp_dim = PROTECT(NEW_INTEGER(4)); ++np;
      INTEGER(pop_db_sexp_dim)[0] = pDB;
      INTEGER(pop_db_sexp_dim)[1] = NG;
      INTEGER(pop_db_sexp_dim)[2] = pDS;
      INTEGER(pop_db_sexp_dim)[3] = PROJ_YEARS;
      SET_DIM(data_db_sexp, pop_db_sexp_dim);
      SET_ATTR(pop_sexp, Rf_install("debut_pop"), data_db_sexp);
    }
  }
  else 
    SET_CLASS(pop_sexp, Rf_mkString("dempp"));
  UNPROTECT(np + 15); // pop 
}