#include "Classes.hpp"

void popC::infect_spec (const hivC& hivpop, const artC& artpop, int time_step,
                        Views& v, const Parameters& p, const StateSpace& s) {
  int ts = (s.year-1)/ s.DT + time_step,
      p_lo = s.p_age15to49_[0] - 1, h_lo = s.h_age15to49_[0] - 1;
  double dt_ii = 1 - s.DT * time_step, // transition of population in 1 year
         n_neg_mf = 0, n_pos_mf = 0, n_pos_inactive = 0, n_pos_inactive_lo = 0,
         n_pos_lo = 0, n_pos_up = 0, n_neg_lo = 0, n_neg_up = 0,
         n_hiv_lo = 0, n_art_lo = 0, n_hiv_up = 0, n_art_up = 0, art_ii = 0;

  update_active_pop_to(s.year, v, s); // substract virgin when needed

  for (int sex = 0; sex < s.NG; sex++) {
    for (int age = p_lo; age < s.pAG_1549; age++) {
      n_neg_mf += v.now_pop[s.N][sex][age]; //  this includes debut neg
      if (s.MODEL == 2 && age < s.pDB)
        n_pos_inactive += data_db[s.year][s.P][sex][age]; // "safe" positive
      n_pos_mf += data_active[s.P][sex][age]; // transmissible positive
    }
    n_neg_lo += v.now_pop[s.N][sex][p_lo];
      if (s.MODEL == 2)
        n_pos_inactive_lo += data_db[s.year][s.P][sex][p_lo]; // "safe" positive
    n_neg_up += v.now_pop[s.N][sex][s.pAG_1549];
    n_pos_lo += data_active[s.P][sex][p_lo];
    n_pos_up += data_active[s.P][sex][s.pAG_1549];
    if (s.year >= s.tARTstart-1) {
      for (int agr = h_lo; agr < s.hAG_1549; agr++)
        for (int cd4 = 0; cd4 < s.hDS; cd4++)
          for (int dur = 0; dur < s.hTS; dur++) {
            art_ii   += v.now_art[sex][agr][cd4][dur];
            n_art_lo += v.now_art[sex][h_lo][cd4][dur];
            n_art_up += v.now_art[sex][s.hAG_1549][cd4][dur];
          }
    }
    for (int cd4 = 0; cd4 < s.hDS; cd4++) {
      n_hiv_lo += v.now_hiv[sex][h_lo][cd4];
      n_hiv_up += v.now_hiv[sex][s.hAG_1549][cd4];
    }
  }
  
  double hivp_inactive = n_pos_inactive - n_pos_inactive_lo * dt_ii;
  double hivn_both = n_neg_mf - n_neg_lo * dt_ii + n_neg_up * dt_ii;
  double hivp_active = n_pos_mf - n_pos_lo * dt_ii + n_pos_up * dt_ii;
  
  if (s.year >= s.tARTstart-1) {
    if (n_hiv_lo + n_art_lo > 0) {
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
    if (n_hiv_up + n_art_up > 0) {
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
  }
  
  // Prob of contacting a sexual active, H+ in total pop
  double
  transm_prev = (hivp_active - art_ii * (1 - p.ic.relinfectART)) / 
                (hivn_both + hivp_active + hivp_inactive);

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
  prev_last = (hivp_active + hivp_inactive) / (hivn_both + hivp_active + hivp_inactive);
  prev15to49_ts[ts] = prev_last;
}

void popC::infect_mix (hivC& hivpop, artC& artpop, int ii, Views& v, const Parameters& p, const StateSpace& s) {
  update_active_pop_to(s.year, v, s);
  // for (int ds = 0; ds < s.pDS; ds++)
  //   for (int sex = 0; sex < s.NG; sex++)
  //     for (int age = 0; age < s.pAG; age++)
  //       data_active[ds][sex][age] *= p.ic.est_senesence[sex][age];

  dvec prop_n_m(s.pAG), prop_n_f(s.pAG);

  for (int r = 0; r < s.pAG; ++r) {
    double Ma = data_active[s.N][s.M][r] + data_active[s.P][s.M][r];
    double Fa = data_active[s.N][s.F][r] + data_active[s.P][s.F][r];
    prop_n_m[r] = data_active[s.N][s.M][r] / Ma;
    prop_n_f[r] = data_active[s.N][s.F][r] / Fa;
  }
  
  boost3D actual_active = data_active;

  // number of sex acts
  // for (int ds = 0; ds < s.pDS; ds++)
  //   for (int sex = 0; sex < s.NG; sex++)
  //     for (int age = 0; age < s.pAG; age++)
  //       data_active[ds][sex][age] *= (1 + p.ic.est_pcr[sex][age]);

  // balancing number of sex acts
  boost2D
    nc_m(extents[s.pAG][s.pAG]),
    nc_f(extents[s.pAG][s.pAG]),
    nc_m_adj(extents[s.pAG][s.pAG]),
    nc_f_adj(extents[s.pAG][s.pAG]),
    n_m_active_negative(extents[s.pAG][s.pAG]),
    n_f_active_negative(extents[s.pAG][s.pAG]);

  for (int r = 0; r < s.pAG; ++r) {
    double Ma = data_active[s.N][s.M][r] + data_active[s.P][s.M][r];
    double Fa = data_active[s.N][s.F][r] + data_active[s.P][s.F][r];
    for (int c = 0; c < s.pAG; ++c) {
      nc_m[c][r] = Ma * p.ic.mixmat[s.M][c][r];
      nc_f[c][r] = Fa * p.ic.mixmat[s.F][c][r];
    }
  }

  for (int r = 0; r < s.pAG; ++r) 
    for (int c = 0; c < s.pAG; ++c) {
      double ratio_mf = nc_m[c][r] / nc_f[r][c];
      double rr = ratio_mf - p.ic.balancing * (ratio_mf  - 1);
      nc_m_adj[c][r] = nc_m[c][r] * rr / ratio_mf;
      nc_f_adj[c][r] = nc_f[r][c] * rr;
    }

  // Number of sex acts in HIV negative pop
  for (int r = 0; r < s.pAG; ++r)
    for (int c = 0; c < s.pAG; ++c) {
      n_m_active_negative[c][r] = nc_m_adj[c][r] * prop_n_m[r];
      n_f_active_negative[c][r] = nc_m_adj[r][c] * prop_n_f[r];
    }
  boost2D art_cov(extents[s.NG][s.pAG]);
  if (s.year >= s.tARTstart-1)
    art_cov = age_sex_cov(hivpop, artpop, v, p, s);

  int ts = (s.year-1)/s.DT + ii;

  boost1D transm_prev(extents[s.NG]);
  for (int sex = 0; sex < s.NG; sex++) {
    double all_pop=0, hiv_treated=0, hiv_not_treated=0;
    for (int age = 0; age < s.pAG; age++) {
      hiv_treated     += data_active[s.P][sex][age] * art_cov[sex][age];
      hiv_not_treated += data_active[s.P][sex][age] * (1 - art_cov[sex][age]);
      all_pop += actual_active[s.N][sex][age] + actual_active[s.P][sex][age];
    }
    transm_prev[sex] = (hiv_not_treated + hiv_treated * (1 - p.ic.relinfectART)) / all_pop;
  }

  //+intervention effects and time epidemic start
  double w = (p.ic.proj_steps[ts] == p.ic.tsEpidemicStart) ? p.ic.iota : 0.0;
  multiply_with_inplace(transm_prev, rvec[ts]);
  transm_prev[s.F] *= p.ic.incrr_sex[s.year];
  add_to_each_inplace(transm_prev, w);

  boost2D inc_m(extents[s.pAG][s.pAG]), inc_f(extents[s.pAG][s.pAG]);

  // adjusted to IRRa
  for (int r = 0; r < s.pAG; ++r)
    for (int c = 0; c < s.pAG; ++c) {
      inc_m[c][r] = n_m_active_negative[c][r] * transm_prev[s.F] * p.ic.incrr_age[s.year][s.M][r];
      inc_f[c][r] = n_f_active_negative[c][r] * transm_prev[s.M] * p.ic.incrr_age[s.year][s.F][r];
    }

  boost1D inc_mv = rowSums(inc_m), inc_fv = rowSums(inc_f);

  for (int age = 0; age < s.pAG; ++age) {
    infections_[s.M][age] = inc_mv[age];
    infections_[s.F][age] = inc_fv[age];
  }

  // prev15to49_ts_m should use this one! now just store as below
  boost2D n_pos = v.now_pop[ indices[s.P][_all][_all] ];
  prev15to49_ts[ts] = sumArray(n_pos) / sumArray(v.now_pop);
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