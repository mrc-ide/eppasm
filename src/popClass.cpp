#include "Classes.hpp"

void popC::initiate (const Parameters& p, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int age = 0; age < s.pAG; age++)
      data[0][s.hivn_idx][sex][age] = p.dm.basepop[sex][age];

  birthslag = p.dm.birthslag;

  if (p.ic.eppmod == 0)
    for (int i = 0; i < s.n_steps; ++i)
      rvec[i] = p.ic.rvec[i];
}

void popC::update_active_pop_to (int when, const StateSpace& s) {
  data_active = data[ indices[when][_all][_all][_all] ];
  if (s.MODEL == 2)
    for (int ds = 0; ds < s.pDS; ds++)
      for (int sex = 0; sex < s.NG; sex++)
        for (int age = 0; age < s.pDB; age++)
          data_active[ds][sex][age] -= data_db[when][ds][sex][age];
}

void popC::update_active_last_year (const StateSpace& s) {
  active_last_year_ = data[ indices[s.year-1][_all][_all][_all] ];
  if (s.MODEL == 2)
    for (int ds = 0; ds < s.pDS; ds++)
      for (int sex = 0; sex < s.NG; sex++)
        for (int age = 0; age < s.pDB; age++)
          active_last_year_[ds][sex][age] -= data_db[s.year-1][ds][sex][age];
}

void popC::aging (const StateSpace& s) { // open ended
  for (int ds = 0; ds < s.pDS; ds++)
    for (int sex = 0; sex < s.NG; sex++) {
      for (int age = 1; age < s.pAG; age++) {
        data[s.year][ds][sex][age]  = data[s.year-1][ds][sex][age-1];
        if (s.MODEL == 2 && age < s.pDB)
          data_db[s.year][ds][sex][age] = data_db[s.year-1][ds][sex][age-1];
      }
      data[s.year][ds][sex][s.pAG-1] += data[s.year-1][ds][sex][s.pAG-1];
    }
}

void popC::add_entrants (const Parameters& p, const StateSpace& s) { // Add lagged births into youngest age group
  double healthy[2], positiv[2];
  // if (exists("popadjust", where=p) & p.popadjust) {
  if (p.dm.flag_popadjust) {
    if (s.MODEL==0) {
      for (int sex = 0; sex < s.NG; ++sex)
        healthy[sex] = p.dm.entrantpop[s.year-1][sex];
    } else {
      for (int sex = 0; sex < s.NG; ++sex) {
        healthy[sex] = p.dm.entrantpop[s.year-1][sex] * (1-p.ph.entrantprev[s.year][sex]);
        positiv[sex] = p.dm.entrantpop[s.year-1][sex] *    p.ph.entrantprev[s.year][sex];
      }
    }
  } else {
    if (s.MODEL==0) {
      for (int sex = 0; sex < s.NG; ++sex)
        healthy[sex] = birthslag[s.year-1][sex] * p.dm.cumsurv[s.year-1][sex] / 
                       p.ph.paedsurv_lag[s.year-1] + p.dm.cumnetmigr[s.year-1][sex];
    } else {
      for (int sex = 0; sex < s.NG; ++sex) {
        healthy[sex] = birthslag[s.year-1][sex] * p.dm.cumsurv[s.year-1][sex] * 
          (1 - p.ph.entrantprev[s.year][sex] / p.ph.paedsurv_lag[s.year-1]) +
          p.dm.cumnetmigr[s.year-1][sex] * (1 - pregprevlag[s.year-1] * p.ph.netmig_hivprob);
        positiv[sex] = (birthslag[s.year-1][sex] * p.dm.cumsurv[s.year-1][sex] + 
          p.dm.cumnetmigr[s.year-1][sex]) * p.ph.entrantprev[s.year][sex];
      }
    }
  }
  // save and update pop
  double sum_p = 0, sum_h = 0;
  for (int sex = 0; sex < s.NG; ++sex) {
    data[s.year][s.hivn_idx][sex][0] = healthy[sex];
    if (s.MODEL != 0) {
      data[s.year][s.hivp_idx][sex][0] = positiv[sex];
      sum_p += positiv[sex]; sum_h += healthy[sex];
      hivp_entrants_out[s.year][sex] = positiv[sex];
      // 1:2 ART+, 3:4 ART-
      entrant_art_[sex]   = positiv[sex] *      p.ph.entrantartcov[s.year][sex];
      entrant_art_[sex+2] = positiv[sex] * (1 - p.ph.entrantartcov[s.year][sex]);
    }
    if (s.MODEL==2) { // add to virgin to record
      data_db[s.year][s.hivn_idx][sex][0] = healthy[sex];
      data_db[s.year][s.hivp_idx][sex][0] = positiv[sex];
    }
  }
  
  if (s.MODEL != 0)
    entrantprev[s.year] = sum_p / (sum_p + sum_h);
}

void popC::sexual_debut (const Parameters& p, const StateSpace& s) {
  for (int ds = 0; ds < s.pDS; ds++)
    for (int sex = 0; sex < s.NG; sex++)
      for (int age = 0; age < s.pDB; age++)
        data_db[s.year][ds][sex][age] *= (1 - p.ic.db_pr[sex][age]);
}

void popC::update_hiv_aging_prob (const StateSpace& s) {
  zeroing(hiv_by_agrp_);
  for (int sex = 0; sex < s.NG; ++sex) {
    int current_age_group = s.ag_idx[0]; // first age group
    for (int age = 0; age < s.pAG; ++age) {
      if ( s.ag_idx[age] != current_age_group)
        ++current_age_group;
      hiv_by_agrp_[sex][current_age_group-1] +=
        data[s.year-1][s.hivp_idx][sex][age];
    } // end age-groups
  }
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      hiv_aging_prob_[sex][agr] = (hiv_by_agrp_[sex][agr] == 0) ? 0 :
        data[s.year-1][s.hivp_idx][sex][s.aglast_idx[agr]-1] / hiv_by_agrp_[sex][agr];
}

void popC::deaths (const Parameters& p, const StateSpace& s) {
  zeroing(hiv_by_agrp_);
  for (int sex = 0; sex < s.NG; ++sex) {
    int current_age_group = s.ag_idx[0]; // first age group
    for (int age = 0; age < s.pAG; ++age) {
      if ( s.ag_idx[age] != current_age_group)
        ++current_age_group;
      hiv_by_agrp_[sex][current_age_group-1] += data[s.year][s.hivp_idx][sex][age];
    } // end age-groups
  }
  for (int ds = 0; ds < s.pDS; ds++)
    for (int sex = 0; sex < s.NG; sex++)
      for (int age = 0; age < s.pAG; age++) {
        num_death_[ds][sex][age] =
          data[s.year][ds][sex][age] * (1 - p.dm.Sx[s.year][sex][age]);
        data[s.year][ds][sex][age] *= p.dm.Sx[s.year][sex][age];
        if (s.MODEL == 2 && age < s.pDB)
          data_db[s.year][ds][sex][age] *= p.dm.Sx[s.year][sex][age];
      }
  // calculate survival prob for hivpop and artpop
  zeroing(death_by_agrp_);
  for (int sex = 0; sex < s.NG; ++sex) {
    int current_age_group = s.ag_idx[0]; // first age group
    for (int age = 0; age < s.pAG; ++age) {
      if ( s.ag_idx[age] != current_age_group)
        ++current_age_group;
      death_by_agrp_[sex][current_age_group-1] += num_death_[s.hivp_idx][sex][age];
    } // end age-groups
  }
  for (int sex = 0; sex < s.NG; sex++)
    for (int age = 0; age < s.hAG; age++)
      hiv_sx_prob[sex][age] = (hiv_by_agrp_[sex][age] == 0) ? 0 :
        1 - (death_by_agrp_[sex][age] / hiv_by_agrp_[sex][age]);
  //  save natural death outputs
  for (int sex = 0; sex < s.NG; ++sex)
    for (int age = 0; age < s.pAG; ++age)
      natdeaths[s.year][sex][age] = 
        num_death_[0][sex][age] + num_death_[1][sex][age];
}

void popC::migration (const Parameters& p, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; ++sex)
    for (int age = 0; age < s.pAG; ++age) {
      double netmigsurv =
        p.dm.netmigr[s.year][sex][age] * (1 + p.dm.Sx[s.year][sex][age]) / 2;
      migrate_prob_[sex][age] = 1 + netmigsurv /
        ( data[s.year][s.hivn_idx][sex][age] + data[s.year][s.hivp_idx][sex][age] );
      num_migrate_[sex][age] =
        migrate_prob_[sex][age] * data[s.year][s.hivp_idx][sex][age];
    }
  zeroing(migrant_by_agrp_);
  for (int sex = 0; sex < s.NG; ++sex) {
    int current_age_group = s.ag_idx[0]; // first age group
    for (int age = 0; age < s.pAG; ++age) {
      if ( s.ag_idx[age] != current_age_group)
        ++current_age_group;
      migrant_by_agrp_[sex][current_age_group-1] += num_migrate_[sex][age];
    } // end age-groups
  }
  zeroing(hiv_by_agrp_);
  for (int sex = 0; sex < s.NG; ++sex) {
    int current_age_group = s.ag_idx[0]; // first age group
    for (int age = 0; age < s.pAG; ++age) {
      if ( s.ag_idx[age] != current_age_group)
        ++current_age_group;
      hiv_by_agrp_[sex][current_age_group-1] += data[s.year][s.hivp_idx][sex][age];
    } // end age-groups
  }
  for (int sex = 0; sex < s.NG; sex++)
    for (int age = 0; age < s.hAG; age++) {
      hiv_mr_prob[sex][age] = (hiv_by_agrp_[sex][age] == 0) ? 0 :
        migrant_by_agrp_[sex][age] / hiv_by_agrp_[sex][age];
    }
  for (int ds = 0; ds < s.pDS; ds++)
    for (int sex = 0; sex < s.NG; sex++)
      for (int age = 0; age < s.pAG; age++) {
        data[s.year][ds][sex][age] *= migrate_prob_[sex][age];
        if (s.MODEL == 2 && age < s.pDB)
          data_db[s.year][ds][sex][age] *= migrate_prob_[sex][age];
      }
}

void popC::update_fertile (const Parameters& p, const StateSpace& s) { // only on active pop
  update_active_pop_to(s.year, s);
  update_active_last_year(s);
  for (int age = 0; age < s.pAG_FERT; age++)
    birth_age[age] =
      ((data_active[s.hivp_idx][s.f_idx][age] + data_active[s.hivn_idx][s.f_idx][age] +
        active_last_year_[s.hivp_idx][s.f_idx][age] +
        active_last_year_[s.hivn_idx][s.f_idx][age]) /
      2) * p.dm.asfr[s.year][age];
  ivec sub_id(s.ag_idx.begin() + s.p_fert_idx[0] - 1, s.ag_idx.begin() + s.pAG_FERT);
  birth_agrp = sumByAG(birth_age, sub_id, s.hAG_FERT);
  double n_births = sum_vector(birth_agrp);
  if ( (s.year + s.AGE_START) <= (s.PROJ_YEARS - 1) )
    for (int sex = 0; sex < s.NG; ++sex)
      birthslag[s.year + s.AGE_START - 1][sex] = p.dm.srb[s.year][sex] * n_births;
}

void popC::adjust_pop (const Parameters& p, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int age = 0; age < s.pAG; age++)
      popadjust[s.year][sex][age] = p.dm.targetpop[s.year][sex][age] / 
        ( data[s.year][s.hivn_idx][sex][age] + data[s.year][s.hivp_idx][sex][age] );
  if (s.MODEL!=0) {  // calculate asjust prob for hiv and art pops
    for (int sex = 0; sex < s.NG; sex++)
      for (int age = 0; age < s.pAG; age++)
        num_adjust_[sex][age] = 
          popadjust[s.year][sex][age] * data[s.year][s.hivp_idx][sex][age];
    zeroing(num_adjust_by_agrp_);
    for (int sex = 0; sex < s.NG; ++sex) {
      int current_age_group = s.ag_idx[0]; // first age group
      for (int age = 0; age < s.pAG; ++age) {
        if ( s.ag_idx[age] != current_age_group)
          ++current_age_group;
        num_adjust_by_agrp_[sex][current_age_group-1] += num_adjust_[sex][age];
      } // end age-groups
    }
    zeroing(hiv_by_agrp_);
    for (int sex = 0; sex < s.NG; ++sex) {
      int current_age_group = s.ag_idx[0]; // first age group
      for (int age = 0; age < s.pAG; ++age) {
        if ( s.ag_idx[age] != current_age_group)
          ++current_age_group;
        hiv_by_agrp_[sex][current_age_group-1] += data[s.year][s.hivp_idx][sex][age];
      } // end age-groups
    }
    for (int sex = 0; sex < s.NG; sex++)
      for (int age = 0; age < s.hAG; age++) {
        adj_prob[sex][age] = (hiv_by_agrp_[sex][age] == 0) ? 0 :
          num_adjust_by_agrp_[sex][age] / hiv_by_agrp_[sex][age];
      }
  }
  for (int ds = 0; ds < s.pDS; ds++)  // only then adjust myself
    for (int sex = 0; sex < s.NG; sex++)
      for (int age = 0; age < s.pAG; age++) {
        data[s.year][ds][sex][age] *= popadjust[s.year][sex][age];
        if (s.MODEL == 2 && age < s.pDB)
          data_db[s.year][ds][sex][age] *= popadjust[s.year][sex][age];
      }
}

void popC::cal_prev_pregant (const hivC& hivpop, const artC& artpop,
                             const Parameters& p, const StateSpace& s) { // only on active pop
  dvec n_mean(s.pAG_FERT); // 1 X 35
  update_active_pop_to(s.year, s); 
  update_active_last_year(s);
  for (int age = 0; age < s.pAG_FERT; ++age)
    n_mean[age] = (active_last_year_[s.hivn_idx][s.f_idx][age] +
                   data_active[s.hivn_idx][s.f_idx][age]) / 2;
  ivec sub_id(s.ag_idx.begin() + s.p_fert_idx[0] - 1, s.ag_idx.begin() + s.pAG_FERT);
  dvec hivn = sumByAG(n_mean, sub_id, s.hAG_FERT); // 1 x 8
  double frap = 0;
  for (int agr = 0; agr < s.hAG_FERT; ++agr) {
    double frp = 0, fra = 0;
    for (int cd4 = 0; cd4 < s.hDS; ++cd4) {
      frp += (hivpop.data[s.year-1][s.f_idx][agr][cd4] +
              hivpop.data[s.year  ][s.f_idx][agr][cd4] ) / 2 * 
              p.nh.frr_cd4[s.year][agr][cd4];
      if (s.year >= s.tARTstart - 1)
        for (int dur = 0; dur < s.hTS; dur++)
          fra += (artpop.data[s.year-1][s.f_idx][agr][cd4][dur] +
                  artpop.data[s.year  ][s.f_idx][agr][cd4][dur]) / 2 *
                  p.nh.frr_art[s.year][agr][cd4][dur];
    }
    frap += birth_agrp[agr] * (1 - hivn[agr] / (hivn[agr] + frp + fra));
  }
  pregprevlag[s.year + s.AGE_START - 1] = frap / sum_vector(birth_age);
}

void popC::save_prev_n_inc (const StateSpace& s) {
  if (s.year + s.AGE_START > s.PROJ_YEARS - 1) // otherwise did in cal_prev_pregant
    update_active_last_year(s);
  double n_positive = 0, everyone_now = 0, s_previous = 0;
  for (int sex = 0; sex < s.NG; sex++)
    for (int age = s.p_age15to49_idx[0] - 1; age < s.pAG_1549; age++) {
      n_positive += data[s.year][s.hivp_idx][sex][age]; // +virgin
      s_previous += active_last_year_[s.hivn_idx][sex][age]; // susceptible -virgin
      for (int ds = 0; ds < s.pDS; ds++)
        everyone_now += data[s.year][ds][sex][age];
    }
  prev15to49[s.year] = n_positive / everyone_now;
  incid15to49[s.year] /= s_previous;
  // prev(s.year) = accu(data_all.slice(s.hivp_idx)) / accu(data_all);
  // incid(s.year) = incid15to49(s.year) / accu(data(s.year-1).slice(s.hivn_idx)); // toBfixed
}

double popC::calc_rtrend_rt (int ts, double time_step, const StateSpace& s) {
  double rveclast = rvec[ts-1];
  Rf_error("K not write for r_trend, no fp template to write");
  // double dtii =  1 - s.DT * (time_step - 1);
  // int a_l = s.p_age15to49_idx[0] -1, a_r = a_l + s.p_age15to49_idx.num_elements();
  // boost2D A = data[ indices[s.year][s.hivn_idx][in(0, s.NG)][in(a_l, a_r)] ];
  // boost1D B = data[ indices[s.year][s.hivn_idx][in(0, s.NG)][a_l] ];
  // boost1D C = data[ indices[s.year][s.hivn_idx][in(0, s.NG)][a_r] ];
  // double hivn_ii = sumArray(A) - sumArray(B) * dtii + sumArray(C) * dtii;
  // A = data[ indices[s.year][s.hivp_idx][in(0, s.NG)][in(a_l, a_r)] ];
  // B = data[ indices[s.year][s.hivp_idx][in(0, s.NG)][a_l] ];
  // C = data[ indices[s.year][s.hivp_idx][in(0, s.NG)][a_r] ];
  // double hivp_ii = sumArray(A) - sumArray(B) * dtii + sumArray(C) * dtii;
  // double prevcurr = hivp_ii / (hivn_ii + hivp_ii);
  // double t_ii     = p.proj_steps[ts];
  // if (t_ii > p.tsEpidemicStart) {
    // par = p.rtrend;
    // if (t_ii < par.tStabilize)
    //   gamma_t = 0
    // else 
    //   gamma_t = (prevcurr-prevlast) * (t_ii - par$tStabilize) / (s.DT * prevlast)
    // logr.diff =  par$beta[2] * (par$beta[1] - rveclast) +
    //              par$beta[3] * prevlast + 
    //              par$beta[4] * gamma_t
    // return(exp(log(rveclast) + logr.diff))
  // }
  // else
    // return(p.rtrend$r0)
  return rveclast; // remove this
}

void popC::update_rvec (double time_step, const Parameters& p, const StateSpace& s) {
  int ts = (s.year-1) / s.DT + time_step;
  // if (p.eppmod %in% c("rtrend", "rtrend_rw")) <<- need to be fixed in FP
  if (p.ic.eppmod == 1) // rtrend see prepare_fp_for_Cpp
    rvec[ts] = calc_rtrend_rt(ts, time_step, s);
}

void popC::update_infection (const StateSpace& s) {
  double n_infect;
  for (int sex = 0; sex < s.NG; sex++)
    for (int age = 0; age < s.pAG; age++) {
      n_infect = infections_[sex][age] * s.DT;
      data[s.year][s.hivn_idx][sex][age] -= n_infect;
      data[s.year][s.hivp_idx][sex][age] += n_infect;
      infections[s.year][sex][age]     += n_infect;
      if ( (age >=  s.p_age15to49_idx[0] - 1) & (age < s.pAG_1549) )
        incid15to49[s.year] += n_infect;
    }
}

void popC::remove_hiv_death (const hivC& hivpop, const artC& artpop,
                             const Parameters& p, const StateSpace& s) {
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++) {
      double hivD = 0, artD = 0;
      for (int cd4 = 0; cd4 < s.hDS; cd4++) {
        hivD += hivpop.death_[sex][agr][cd4];
        if (s.MODEL==2) // add hiv deaths from inactive population
          hivD += hivpop.death_db_[sex][agr][cd4];
        if (s.year >= s.tARTstart - 1) {
          for (int dur = 0; dur < s.hTS; dur++) {
            artD += artpop.death_[sex][agr][cd4][dur];
            if (s.MODEL==2) // add art deaths from inactive population
              artD += artpop.death_db_[sex][agr][cd4][dur];
          }
        }
      }
      death_by_agrp_[sex][agr] = s.DT * (hivD + artD); // deaths by single-year
    }
  zeroing(hiv_by_agrp_);
  for (int sex = 0; sex < s.NG; ++sex) {
    int current_age_group = s.ag_idx[0]; // first age group
    for (int age = 0; age < s.pAG; ++age) {
      if ( s.ag_idx[age] != current_age_group)
        ++current_age_group;
      hiv_by_agrp_[sex][current_age_group-1] += data[s.year][s.hivp_idx][sex][age];
    } // end age-groups
  }
  double hiv_mx;
  for (int sex = 0; sex < s.NG; sex++) {
    int age = 0;
    for (int agr = 0; agr < s.hAG; agr++) {
      if (hiv_by_agrp_[sex][agr] != 0) {
        hiv_mx = death_by_agrp_[sex][agr] / hiv_by_agrp_[sex][agr];
        for (int i = 0; i < s.h_ag_span[agr]; ++i) {
          hivdeaths[s.year][sex][age] += data[s.year][s.hivp_idx][sex][age] * hiv_mx;
          data[s.year][s.hivp_idx][sex][age] *= (1 - hiv_mx);
          if (age < s.pDB)
            data_db[s.year][s.hivp_idx][sex][age] *= (1 - hiv_mx);
          age++;
        }
      } else {
        age += s.h_ag_span[agr] - 1;
      }
    }
  }
}

void popC::update_preg (const hivC& hivpop, const artC& artpop,
                        const Parameters& p, const StateSpace& s) {
  int h_lo = s.h_fert_idx[0] - 1, // 0 9
      p_lo = s.p_fert_idx[0] - 1; // 0 35
  update_active_pop_to(s.year, s);
  dvec hivn(s.hAG_FERT);
  int current_age_group = s.ag_idx[p_lo]; // first age group
  for (int age = p_lo; age < s.pAG_FERT; ++age) {
    if ( s.ag_idx[age] != current_age_group)
      ++current_age_group;
    hivn[current_age_group - 1] += data_active[s.hivn_idx][s.f_idx][age];
  } // end age-groups
  dvec all_art(s.hAG_FERT);
  for (int agr = h_lo; agr < s.hAG_FERT; agr++)
    for (int cd4 = 0; cd4 < s.hDS; cd4++) {
      for (int dur = 0; dur < s.hTS; dur++)
        all_art[agr] += artpop.data[s.year][s.f_idx][agr][cd4][dur] * 
          p.nh.frr_art[s.year][agr][cd4][dur];
      all_art[agr] +=
        hivpop.data[s.year][s.f_idx][agr][cd4] * p.nh.frr_cd4[s.year][agr][cd4];
    }
  for (int agr = h_lo; agr < s.hAG_FERT; agr++)
    for (int cd4 = 0; cd4 < (p.ad.artcd4elig_idx[s.year] - 1); cd4++)
      art_elig_[s.f_idx][agr][cd4] += 
        hivpop.data[s.year][s.f_idx][agr][cd4] * p.nh.frr_cd4[s.year][agr][cd4] *
          (birth_agrp[agr] / (hivn[agr] + all_art[agr]));
}

dvec popC::art_initiate (const dvec& art_curr, int time_step,
                         const Parameters& p, const StateSpace& s) {
  dvec out(s.NG);
  int year_w = (s.DT * (time_step + 1) < 0.5) ? 0 : 1;
  double trans = s.DT * (time_step + 1) + 0.5 - year_w;
  int year_l = s.year - (2 - year_w), year_r = s.year - (1 - year_w);
  for (int sex = 0; sex < s.NG; ++sex) {
    if( (p.ad.art15plus_isperc[year_l][sex] == 0) & 
        (p.ad.art15plus_isperc[year_r][sex] == 0) ) { // both number
      out[sex] = p.ad.art15plus_num[year_l][sex] * (1 - trans) + 
                 p.ad.art15plus_num[year_r][sex] *      trans;
      if (s.MIX) {
        boost2D art_elig_sex = art_elig_[ indices[sex][_all][_all] ];
        double cov = out[sex] / (sumArray(art_elig_sex) + art_curr[sex]);
        artcov[sex] = (cov <= 1) ? cov : 1;
      }
    } 
    else if ( (p.ad.art15plus_isperc[year_l][sex] == 1) &
              (p.ad.art15plus_isperc[year_r][sex] == 1) ) { // both percentage
      double cov = p.ad.art15plus_num[year_l][sex] * (1 - trans) +
                   p.ad.art15plus_num[year_r][sex] *      trans;
      if (s.MIX)
        artcov[sex] = (cov <= 1) ? cov : 1;
      boost2D art_elig_sex = art_elig_[ indices[sex][_all][_all] ];
      out[sex] = cov * (sumArray(art_elig_sex) + art_curr[sex]);
    }
    else if ( (p.ad.art15plus_isperc[year_l][sex] == 0) & 
              (p.ad.art15plus_isperc[year_r][sex] == 1) ) { // transition number to percentage
      boost2D art_elig_sex = art_elig_[ indices[sex][_all][_all] ];
      double actual_cov =
        art_curr[sex] / (sumArray(art_elig_sex) + art_curr[sex]);
      double diff_cov = p.ad.art15plus_num[year_r][sex] - actual_cov;
      double cov = actual_cov + diff_cov * s.DT / (0.5 + year_w - s.DT * time_step);
      if (s.MIX)
        artcov[sex] = (cov < 1) ? cov : 1;
      out[sex] = cov * ( sumArray(art_elig_sex) + art_curr[sex] );
    }
  }
  return out;
}
