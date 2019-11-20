#include "Classes.hpp"

void popC::infect_spec (const hivC& hivpop, const artC& artpop, int time_step,
                        Views& v, const Parameters& p, const StateSpace& s) {
  int ts = (s.year-1)/ s.DT + time_step,
      p_lo = s.p_age15to49_[0] - 1, h_lo = s.h_age15to49_[0] - 1;
  double dt_ii = 1 - s.DT * time_step, // transition of population in 1 year
         n_neg_mf = 0, n_pos_mf = 0, 
         n_pos_lo = 0, n_pos_up = 0, n_neg_lo = 0, n_neg_up = 0,
         n_hiv_lo = 0, n_art_lo = 0, n_hiv_up = 0, n_art_up = 0, art_ii = 0;

  update_active_pop_to(s.year, v, s); // substract virgin when needed

  for (int sex = 0; sex < s.NG; sex++) {
    for (int age = p_lo; age < s.pAG_1549; age++) {
      n_neg_mf += data_active[s.N][sex][age];
      n_pos_mf += data_active[s.P][sex][age];
    }
    n_neg_lo += data_active[s.N][sex][p_lo];
    n_neg_up += data_active[s.N][sex][s.pAG_1549];
    n_pos_lo += data_active[s.P][sex][p_lo];
    n_pos_up += data_active[s.P][sex][s.pAG_1549];
    for (int agr = h_lo; agr < s.hAG_1549; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        for (int dur = 0; dur < s.hTS; dur++) {
          art_ii   += v.now_art[sex][agr][cd4][dur];
          n_art_lo += v.now_art[sex][h_lo][cd4][dur];
          n_art_up += v.now_art[sex][s.hAG_1549][cd4][dur];
        }
    for (int cd4 = 0; cd4 < s.hDS; cd4++) {
      n_hiv_lo += v.now_hiv[sex][h_lo][cd4];
      n_hiv_up += v.now_hiv[sex][s.hAG_1549][cd4];
    }
  }
  double hivn_ii = n_neg_mf - n_neg_lo * dt_ii + n_neg_up * dt_ii;
  double hivp_ii = n_pos_mf - n_pos_lo * dt_ii + n_pos_up * dt_ii;
  if ( n_hiv_lo + n_art_lo > 0) {
    for (int sex = 0; sex < s.NG; sex++) {
      double art_trans = 0, hiv_trans = 0;
      for (int cd4 = 0; cd4 < s.hDS; cd4++) {
        for (int dur = 0; dur < s.hTS; dur++)
          art_trans += v.now_art[sex][h_lo][cd4][dur];
        hiv_trans += v.now_hiv[sex][h_lo][cd4];
      }
      art_ii -= ( data_active[s.P][sex][p_lo] * 
                  art_trans / (hiv_trans + art_trans) ) * dt_ii;
    }
  }
  if ( n_hiv_up + n_art_up > 0) {
    for (int sex = 0; sex < s.NG; sex++) {
      double art_trans = 0, hiv_trans = 0;
      for (int cd4 = 0; cd4 < s.hDS; cd4++) {
        for (int dur = 0; dur < s.hTS; dur++)
          art_trans += v.now_art[sex][s.hAG_1549][cd4][dur];
        hiv_trans += v.now_hiv[sex][s.hAG_1549][cd4];
      }
      art_ii += ( data_active[s.P][sex][s.pAG_1549] * 
                  art_trans / (hiv_trans + art_trans) ) * dt_ii;
    }
  }
  double
  transm_prev = (hivp_ii - art_ii * (1 - p.ic.relinfectART)) / (hivn_ii + hivp_ii);
  double w = (p.ic.proj_steps[ts] == p.ic.tsEpidemicStart) ? p.ic.iota : 0.0;
  double inc_rate = rvec[ts] * transm_prev + w;

  // Incidence: male = negative / female negative * sexratio + male negative; 
  //          female = male * sexratio
  double n_neg_m = 0, n_neg_f = 0;
  for (int age = p_lo; age < s.pAG_1549; age++) {
    n_neg_m += data_active[s.N][s.M][age];
    n_neg_f += data_active[s.N][s.F][age];
  }
  double adj_sex = (n_neg_m + n_neg_f) / (n_neg_m + n_neg_f * p.ic.incrr_sex[s.year]);
  double sex_inc[2] = {inc_rate * adj_sex, inc_rate * adj_sex * p.ic.incrr_sex[s.year]};
  // New infections distributed by age: ratio age_i/ 25-29 age
  for (int sex = 0; sex < s.NG; sex++) {
    double n_neg = 0, n_neg_rr = 0, adj_age;
    for (int age = p_lo; age < s.pAG_1549; age++) {
      n_neg += data_active[s.N][sex][age]; 
      n_neg_rr += p.ic.incrr_age[s.year][sex][age] * data_active[s.N][sex][age];
    }
    adj_age = sex_inc[sex] / ( n_neg_rr / n_neg );
    for (int age = 0; age < s.pAG; age++) {
      if (sex == s.M) // age-specific incidence among circumcised men
        adj_age *= (1 - p.ic.circ_incid_rr * p.ic.circ_prop[s.year][age]);
      infections_[sex][age] = p.ic.incrr_age[s.year][sex][age] * adj_age * 
        data_active[s.N][sex][age];
    }
  }
  // saving
  incrate15to49_ts[ts] = inc_rate;
  prev_last = hivp_ii / (hivn_ii + hivp_ii);
  prev15to49_ts[ts] = prev_last;
}

void popC::infect_mix (int ii, Views& v, const Parameters& p, const StateSpace& s) {
  update_active_pop_to(s.year, v, s);
  int ts = (s.year-1)/s.DT + ii;
  boost2D transm_prev(extents[s.NG][s.pAG]);
  double N_hivp;
  for (int sex = 0; sex < s.NG; sex++)
    for (int age = 0; age < s.pAG; age++) {
      N_hivp = data_active[s.P][sex][age];
      transm_prev[sex][age] = ((N_hivp * (1 - artcov[sex])) + 
                               (N_hivp * artcov[sex] * (1 - p.ic.relinfectART)))/
                               (data_active[s.N][sex][age] + N_hivp);
      }
  //+intervention effects and time epidemic start
  double w = (p.ic.proj_steps[ts] == p.ic.tsEpidemicStart) ? p.ic.iota : 0.0;
  multiply_with_inplace(transm_prev, rvec[ts]);
  add_to_each_inplace(transm_prev, w);

  // sweep over sexual mixing matrices
  zeroing(infections_);
  for (int my_age = 0; my_age < s.pAG; ++my_age)
    for (int partner_age = 0; partner_age < s.pAG; ++partner_age) {
      infections_[s.M][my_age] +=
        p.ic.mat_m[partner_age][my_age] * transm_prev[s.F][partner_age];
      infections_[s.F][my_age] +=
        p.ic.mat_f[partner_age][my_age] * transm_prev[s.M][partner_age];
    }
  // if (exists("f_fun", fp)) // that fun
  //   ir = ir * fp.f_fun

  // Match IRR by age // incrr_age was scaled in R
  for (int sex = 0; sex < s.NG; sex++)
    for (int age = 0; age < s.pAG; age++)
      infections_[sex][age] *= p.ic.incrr_age[s.year][sex][age];

  // Match IRR by Sex
  double inc_M = 0, inc_F = 0, S_M = 0, S_F = 0;
  for (int age = 0; age < s.pAG; age++) {
    S_M   += data_active[s.N][s.M][age];
    S_F   += data_active[s.N][s.F][age];
    inc_M += infections_[s.M][age] * data_active[s.N][s.M][age];
    inc_F += infections_[s.F][age] * data_active[s.N][s.F][age];
  }
  double IIR_estimated = (inc_F/S_F) / (inc_M/S_M);
  if (!std::isnan(IIR_estimated)) {
    for (int age = 0; age < s.pAG; age++) {
      infections_[s.F][age] *= p.ic.incrr_sex[s.year];
      infections_[s.M][age] *= IIR_estimated;
    }
  }

  for (int sex = 0; sex < s.NG; sex++)
    for (int age = 0; age < s.pAG; age++)
      infections_[sex][age] *= data_active[s.N][sex][age];
  // incrate15to49_ts_m.slice(ts) = ir_mf;
  // prev15to49_ts_m should use this one! now just store as below
  boost2D n_pos = data_active[ indices[s.P][_all][_all] ];
  prev15to49_ts[ts] = sumArray(n_pos) / sumArray(data_active);
  prev_last = prev15to49_ts[ts];
}

void popC::epp_disease_model_direct(hivC& hivpop, artC& artpop, Views& v,
                                    const Parameters& p, const StateSpace& s) {
  int a_l, a_r;
  if (p.ic.incidpopage) { // incidence for 15+ population
    a_l = s.p_age15plus_[0] - 1;
    a_r = s.pAG_15plus;
  } else { // incidence for 15 -49 population
    a_l = s.p_age15to49_[0] - 1;
    a_r = s.pAG_1549;
  }
  update_active_last_year(v, s);
  double n_m = 0, n_f = 0;
  for (int age = a_l; age < a_r; ++age) {
    n_m += active_last_year_[s.N][s.M][age];
    n_f += active_last_year_[s.N][s.F][age];
  }
  dvec sex_inc(s.NG);
  sex_inc[s.M] = (n_m + n_f) * p.ic.incidinput[s.year] / 
                     (n_m + n_f  * p.ic.incrr_sex[s.year]);
  sex_inc[s.F] = (n_m + n_f) * p.ic.incidinput[s.year] * p.ic.incrr_sex[s.year] /
                     (n_m + n_f  * p.ic.incrr_sex[s.year]);
  dvec ageinc(s.NG);
  for (int sex = 0; sex < s.NG; sex++) {
    double neg_sa = 0, inc_sa = 0;
    for (int age = a_l; age < a_r; age++) {
      neg_sa += active_last_year_[s.N][sex][age];
      inc_sa += active_last_year_[s.N][sex][age] * p.ic.incrr_age[s.year][sex][age];
    }
    ageinc[sex] = inc_sa / neg_sa;
  }
  double new_infect = 0;
  for (int sex = 0; sex < s.NG; sex++)
    for (int age = 0; age < s.pAG; age++) {
      new_infect =
        p.ic.incrr_age[s.year][sex][age] * ( sex_inc[sex] / ageinc[sex]) * 
        active_last_year_[s.N][sex][age];
      infections[s.year][sex][age]        = new_infect;
      v.now_pop[s.N][sex][age] -= new_infect;
      v.now_pop[s.P][sex][age] += new_infect;
    }
  boost2D infect_agrp = 
    sumByAG(infections[ indices[s.year][_all][_all]], s.ag_, s.hAG);
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        v.now_hiv[sex][agr][cd4] +=
          p.nh.cd4_initdist[sex][agr][cd4] * infect_agrp[sex][agr];
  for (int sex = 0; sex < s.NG; sex++)
    for (int age = s.p_age15to49_[0] - 1; age < s.pAG_1549; age++)
      incid15to49[s.year] += infections[s.year][sex][age];
}