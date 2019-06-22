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

void popC::initiate (const Parameters& p) {
  for (int sex = 0; sex < NG; sex++)
    for (int age = 0; age < pAG; age++)
      data[0][hivn_idx][sex][age] = p.basepop[sex][age];

  birthslag = p.birthslag;

  if (p.eppmod == 0)
    for (int i = 0; i < n_steps; ++i)
      rvec[i] = p.rvec[i];
}

void popC::update_active_pop_to (int when) {
  data_active = data[ indices[when][_all][_all][_all] ];
  if (MODEL == 2)
    for (int ds = 0; ds < pDS; ds++)
      for (int sex = 0; sex < NG; sex++)
        for (int age = 0; age < pDB; age++)
          data_active[ds][sex][age] -= data_db[when][ds][sex][age];
}

boost3D popC::get_active_pop_in (int when) {
  boost3D out = data[ indices[when][in(0,pDS)][in(0,NG)][in(0, pAG)] ];
  if (MODEL == 2)
    for (int ds = 0; ds < pDS; ds++)
      for (int sex = 0; sex < NG; sex++)
        for (int age = 0; age < pDB; age++)
          out[ds][sex][age] -= data_db[when][ds][sex][age];
  return out;
}

void popC::update_active_last_year () {
  active_last_year_ = data[ indices[year-1][_all][_all][_all] ];
  if (MODEL == 2)
    for (int ds = 0; ds < pDS; ds++)
      for (int sex = 0; sex < NG; sex++)
        for (int age = 0; age < pDB; age++)
          active_last_year_[ds][sex][age] -= data_db[year-1][ds][sex][age];
}

void popC::aging () { // open ended
  for (int ds = 0; ds < pDS; ds++)
    for (int sex = 0; sex < NG; sex++) {
      for (int age = 1; age < pAG; age++) {
        data[year][ds][sex][age]  = data[year-1][ds][sex][age-1];
        if (MODEL == 2 && age < pDB)
          data_db[year][ds][sex][age] = data_db[year-1][ds][sex][age-1];
      }
      data[year][ds][sex][pAG-1] += data[year-1][ds][sex][pAG-1];
    }
}

void popC::add_entrants (const Parameters& p) { // Add lagged births into youngest age group
  double healthy[2], positiv[2];
  // if (exists("popadjust", where=p) & p.popadjust) {
  if (p.popadjust) {
    if (MODEL==0) {
      for (int sex = 0; sex < NG; ++sex)
        healthy[sex] = p.entrantpop[year-1][sex];
    } else {
      for (int sex = 0; sex < NG; ++sex) {
        healthy[sex] = p.entrantpop[year-1][sex] * (1-p.entrantprev[year][sex]);
        positiv[sex] = p.entrantpop[year-1][sex] *    p.entrantprev[year][sex];
      }
    }
  } else {
    if (MODEL==0) {
      for (int sex = 0; sex < NG; ++sex)
        healthy[sex] = birthslag[year-1][sex] * p.cumsurv[year-1][sex] / 
                       p.paedsurv_lag[year-1] + p.cumnetmigr[year-1][sex];
    } else {
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
      // 1:2 ART+, 3:4 ART-
      entrant_art_[sex]   = positiv[sex] *      p.entrantartcov[year][sex];
      entrant_art_[sex+2] = positiv[sex] * (1 - p.entrantartcov[year][sex]);
    }
    if (MODEL==2) { // add to virgin to record
      data_db[year][hivn_idx][sex][0] = healthy[sex];
      data_db[year][hivp_idx][sex][0] = positiv[sex];
    }
  }
  
  if (MODEL!=0)
    entrantprev[year] = sum_p / (sum_p + sum_h);
}

void popC::sexual_debut (const Parameters& p) {
  for (int ds = 0; ds < pDS; ds++)
    for (int sex = 0; sex < NG; sex++)
      for (int age = 0; age < pDB; age++)
        data_db[year][ds][sex][age] *= (1 - p.db_pr[sex][age]);
}

void popC::update_hiv_aging_prob () {
  zeroing(hiv_by_agrp_);
  for (int sex = 0; sex < NG; ++sex) {
    int current_age_group = ag_idx[0]; // first age group
    for (int age = 0; age < pAG; ++age) {
      if ( ag_idx[age] != current_age_group)
        ++current_age_group;
      hiv_by_agrp_[sex][current_age_group-1] +=
        data[year-1][hivp_idx][sex][age];
    } // end age-groups
  }
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hAG; agr++)
      hiv_aging_prob_[sex][agr] = (hiv_by_agrp_[sex][agr] == 0) ? 0 :
        data[year-1][hivp_idx][sex][aglast_idx[agr]-1] / hiv_by_agrp_[sex][agr];
}

void popC::deaths (const Parameters& p) {
  zeroing(hiv_by_agrp_);
  for (int sex = 0; sex < NG; ++sex) {
    int current_age_group = ag_idx[0]; // first age group
    for (int age = 0; age < pAG; ++age) {
      if ( ag_idx[age] != current_age_group)
        ++current_age_group;
      hiv_by_agrp_[sex][current_age_group-1] += data[year][hivp_idx][sex][age];
    } // end age-groups
  }
  for (int ds = 0; ds < pDS; ds++)
    for (int sex = 0; sex < NG; sex++)
      for (int age = 0; age < pAG; age++) {
        num_death_[ds][sex][age] =
          data[year][ds][sex][age] * (1 - p.Sx[year][sex][age]);
        data[year][ds][sex][age] *= p.Sx[year][sex][age];
        if (MODEL == 2 && age < pDB)
          data_db[year][ds][sex][age] *= p.Sx[year][sex][age];
      }
  // calculate survival prob for hivpop and artpop
  zeroing(death_by_agrp_);
  for (int sex = 0; sex < NG; ++sex) {
    int current_age_group = ag_idx[0]; // first age group
    for (int age = 0; age < pAG; ++age) {
      if ( ag_idx[age] != current_age_group)
        ++current_age_group;
      death_by_agrp_[sex][current_age_group-1] += num_death_[hivp_idx][sex][age];
    } // end age-groups
  }
  for (int sex = 0; sex < NG; sex++)
    for (int age = 0; age < hAG; age++)
      hiv_sx_prob[sex][age] = (hiv_by_agrp_[sex][age] == 0) ? 0 :
        1 - (death_by_agrp_[sex][age] / hiv_by_agrp_[sex][age]);
  //  save natural death outputs
  for (int sex = 0; sex < NG; ++sex)
    for (int age = 0; age < pAG; ++age)
      natdeaths[year][sex][age] = 
        num_death_[0][sex][age] + num_death_[1][sex][age];
}

void popC::migration (const Parameters& p) {
  for (int sex = 0; sex < NG; ++sex)
    for (int age = 0; age < pAG; ++age) {
      double netmigsurv =
        p.netmigr[year][sex][age] * (1 + p.Sx[year][sex][age]) / 2;
      migrate_prob_[sex][age] = 1 + netmigsurv /
        ( data[year][hivn_idx][sex][age] + data[year][hivp_idx][sex][age] );
      num_migrate_[sex][age] =
        migrate_prob_[sex][age] * data[year][hivp_idx][sex][age];
    }
  zeroing(migrant_by_agrp_);
  for (int sex = 0; sex < NG; ++sex) {
    int current_age_group = ag_idx[0]; // first age group
    for (int age = 0; age < pAG; ++age) {
      if ( ag_idx[age] != current_age_group)
        ++current_age_group;
      migrant_by_agrp_[sex][current_age_group-1] += num_migrate_[sex][age];
    } // end age-groups
  }
  zeroing(hiv_by_agrp_);
  for (int sex = 0; sex < NG; ++sex) {
    int current_age_group = ag_idx[0]; // first age group
    for (int age = 0; age < pAG; ++age) {
      if ( ag_idx[age] != current_age_group)
        ++current_age_group;
      hiv_by_agrp_[sex][current_age_group-1] += data[year][hivp_idx][sex][age];
    } // end age-groups
  }
  for (int sex = 0; sex < NG; sex++)
    for (int age = 0; age < hAG; age++) {
      hiv_mr_prob[sex][age] = (hiv_by_agrp_[sex][age] == 0) ? 0 :
        migrant_by_agrp_[sex][age] / hiv_by_agrp_[sex][age];
    }
  for (int ds = 0; ds < pDS; ds++)
    for (int sex = 0; sex < NG; sex++)
      for (int age = 0; age < pAG; age++) {
        data[year][ds][sex][age] *= migrate_prob_[sex][age];
        if (MODEL == 2 && age < pDB)
          data_db[year][ds][sex][age] *= migrate_prob_[sex][age];
      }
}

void popC::update_fertile (const Parameters& p) { // only on active pop
  update_active_pop_to(year);
  update_active_last_year();
  for (int age = 0; age < pAG_FERT; age++)
    birth_age[age] =
      ((data_active[hivp_idx][f_idx][age] + data_active[hivn_idx][f_idx][age] +
        active_last_year_[hivp_idx][f_idx][age] +
        active_last_year_[hivn_idx][f_idx][age]) /
      2) * p.asfr[year][age];
  ivec sub_id(ag_idx.begin() + p_fert_idx[0] - 1, ag_idx.begin() + pAG_FERT);
  birth_agrp = sumByAG(birth_age, sub_id, hAG_FERT);
  double n_births = sum_vector(birth_agrp);
  if ( (year + AGE_START) <= (PROJ_YEARS - 1) )
    for (int sex = 0; sex < NG; ++sex)
      birthslag[year + AGE_START - 1][sex] = p.srb[year][sex] * n_births;
}

void popC::adjust_pop (const Parameters& p) {
  for (int sex = 0; sex < NG; sex++)
    for (int age = 0; age < pAG; age++)
      popadjust[year][sex][age] = p.targetpop[year][sex][age] / 
        ( data[year][hivn_idx][sex][age] + data[year][hivp_idx][sex][age] );
  if (MODEL!=0) {  // calculate asjust prob for hiv and art pops
    for (int sex = 0; sex < NG; sex++)
      for (int age = 0; age < pAG; age++)
        num_adjust_[sex][age] = 
          popadjust[year][sex][age] * data[year][hivp_idx][sex][age];
    zeroing(num_adjust_by_agrp_);
    for (int sex = 0; sex < NG; ++sex) {
      int current_age_group = ag_idx[0]; // first age group
      for (int age = 0; age < pAG; ++age) {
        if ( ag_idx[age] != current_age_group)
          ++current_age_group;
        num_adjust_by_agrp_[sex][current_age_group-1] += num_adjust_[sex][age];
      } // end age-groups
    }
    zeroing(hiv_by_agrp_);
    for (int sex = 0; sex < NG; ++sex) {
      int current_age_group = ag_idx[0]; // first age group
      for (int age = 0; age < pAG; ++age) {
        if ( ag_idx[age] != current_age_group)
          ++current_age_group;
        hiv_by_agrp_[sex][current_age_group-1] += data[year][hivp_idx][sex][age];
      } // end age-groups
    }
    for (int sex = 0; sex < NG; sex++)
      for (int age = 0; age < hAG; age++) {
        adj_prob[sex][age] = (hiv_by_agrp_[sex][age] == 0) ? 0 :
          num_adjust_by_agrp_[sex][age] / hiv_by_agrp_[sex][age];
      }
  }
  for (int ds = 0; ds < pDS; ds++)  // only then adjust myself
    for (int sex = 0; sex < NG; sex++)
      for (int age = 0; age < pAG; age++) {
        data[year][ds][sex][age] *= popadjust[year][sex][age];
        if (MODEL == 2 && age < pDB)
          data_db[year][ds][sex][age] *= popadjust[year][sex][age];
      }
}

void popC::cal_prev_pregant (const hivC& hivpop, const artC& artpop,
                             const Parameters& p) { // only on active pop
  dvec n_mean(pAG_FERT); // 1 X 35
  update_active_pop_to(year); 
  update_active_last_year();
  for (int age = 0; age < pAG_FERT; ++age)
    n_mean[age] = (active_last_year_[hivn_idx][f_idx][age] +
                   data_active[hivn_idx][f_idx][age]) / 2;
  ivec sub_id(ag_idx.begin() + p_fert_idx[0] - 1, ag_idx.begin() + pAG_FERT);
  dvec hivn = sumByAG(n_mean, sub_id, hAG_FERT); // 1 x 8
  double frap = 0;
  for (int agr = 0; agr < hAG_FERT; ++agr) {
    double frp = 0, fra = 0;
    for (int cd4 = 0; cd4 < hDS; ++cd4) {
      frp += (hivpop.data[year-1][f_idx][agr][cd4] +
              hivpop.data[year  ][f_idx][agr][cd4] ) / 2 * 
              p.frr_cd4[year][agr][cd4];
      if (year >= p.tARTstart - 1)
        for (int dur = 0; dur < hTS; dur++)
          fra += (artpop.data[year-1][f_idx][agr][cd4][dur] +
                  artpop.data[year  ][f_idx][agr][cd4][dur]) / 2 *
                  p.frr_art[year][agr][cd4][dur];
    }
    frap += birth_agrp[agr] * (1 - hivn[agr] / (hivn[agr] + frp + fra));
  }
  pregprevlag[year + AGE_START - 1] = frap / sum_vector(birth_age);
}

void popC::save_prev_n_inc () {
  if (year + AGE_START > PROJ_YEARS - 1) // otherwise did in cal_prev_pregant
    update_active_last_year();
  double n_positive = 0, everyone_now = 0, s_previous = 0;
  for (int sex = 0; sex < NG; sex++)
    for (int age = p_age15to49_idx[0] - 1; age < pAG_1549; age++) {
      n_positive += data[year][hivp_idx][sex][age]; // +virgin
      s_previous += active_last_year_[hivn_idx][sex][age]; // susceptible -virgin
      for (int ds = 0; ds < pDS; ds++)
        everyone_now += data[year][ds][sex][age];
    }
  prev15to49[year] = n_positive / everyone_now;
  incid15to49[year] /= s_previous;
  // prev(year) = accu(data_all.slice(hivp_idx)) / accu(data_all);
  // incid(year) = incid15to49(year) / accu(data(year-1).slice(hivn_idx)); // toBfixed
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

void popC::update_rvec (double time_step, const Parameters& p) {
  int ts = (year-1) / DT + time_step;
  // if (p.eppmod %in% c("rtrend", "rtrend_rw")) <<- need to be fixed in FP
  if (p.eppmod == 1) // rtrend see prepare_fp_for_Cpp
    rvec[ts] = calc_rtrend_rt(ts, time_step);
}

void popC::update_infection () {
  double n_infect;
  for (int sex = 0; sex < NG; sex++)
    for (int age = 0; age < pAG; age++) {
      n_infect = infections_[sex][age] * DT;
      data[year][hivn_idx][sex][age] -= n_infect;
      data[year][hivp_idx][sex][age] += n_infect;
      infections[year][sex][age]     += n_infect;
      if ( (age >=  p_age15to49_idx[0] - 1) & (age < pAG_1549) )
        incid15to49[year] += n_infect;
    }
}

void popC::remove_hiv_death (const hivC& hivpop, const artC& artpop,
                             const Parameters& p) {
  for (int sex = 0; sex < NG; sex++)
    for (int agr = 0; agr < hAG; agr++) {
      double hivD = 0, artD = 0;
      for (int cd4 = 0; cd4 < hDS; cd4++) {
        hivD += hivpop.death_[sex][agr][cd4];
        if (MODEL==2) // add hiv deaths from inactive population
          hivD += hivpop.death_db_[sex][agr][cd4];
        if (year >= p.tARTstart - 1) {
          for (int dur = 0; dur < hTS; dur++) {
            artD += artpop.death_[sex][agr][cd4][dur];
            if (MODEL==2) // add art deaths from inactive population
              artD += artpop.death_db_[sex][agr][cd4][dur];
          }
        }
      }
      death_by_agrp_[sex][agr] = DT * (hivD + artD); // deaths by single-year
    }
  zeroing(hiv_by_agrp_);
  for (int sex = 0; sex < NG; ++sex) {
    int current_age_group = ag_idx[0]; // first age group
    for (int age = 0; age < pAG; ++age) {
      if ( ag_idx[age] != current_age_group)
        ++current_age_group;
      hiv_by_agrp_[sex][current_age_group-1] += data[year][hivp_idx][sex][age];
    } // end age-groups
  }
  double hiv_mx;
  for (int sex = 0; sex < NG; sex++) {
    int age = 0;
    for (int agr = 0; agr < hAG; agr++) {
      if (hiv_by_agrp_[sex][agr] != 0) {
        hiv_mx = death_by_agrp_[sex][agr] / hiv_by_agrp_[sex][agr];
        for (int i = 0; i < h_ag_span[agr]; ++i) {
          hivdeaths[year][sex][age] += data[year][hivp_idx][sex][age] * hiv_mx;
          data[year][hivp_idx][sex][age] *= (1 - hiv_mx);
          if (age < pDB)
            data_db[year][hivp_idx][sex][age] *= (1 - hiv_mx);
          age++;
        }
      } else {
        age += h_ag_span[agr] - 1;
      }
    }
  }
}

void popC::update_preg (const hivC& hivpop, const artC& artpop,
                        const Parameters& p) {
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
      art_elig_[f_idx][agr][cd4] += 
        hivpop.data[year][f_idx][agr][cd4] * p.frr_cd4[year][agr][cd4] *
          (birth_agrp[agr] / (hivn[agr] + all_art[agr]));
}

dvec popC::art_initiate (const dvec& art_curr, int time_step,
                         const Parameters& p) {
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
        boost2D art_elig_sex = art_elig_[ indices[sex][_all][_all] ];
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
      boost2D art_elig_sex = art_elig_[ indices[sex][_all][_all] ];
      out[sex] = cov * (sumArray(art_elig_sex) + art_curr[sex]);
    }
    else if ( (p.art15plus_isperc[year_l][sex] == 0) & 
              (p.art15plus_isperc[year_r][sex] == 1) ) { // transition number to percentage
      boost2D art_elig_sex = art_elig_[ indices[sex][_all][_all] ];
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