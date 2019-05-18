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

void popC::my_all (int when) { // get all data when there are extra pops
  data_all = data(when);
  if (MODEL==2)
    data_all.rows(0, pDB-1) += data_db(when);
}

void popC::aging () { // open ended
  data(year).rows(1, pAG-2) = data(year-1).rows(0, pAG-3);
  data(year).row(pAG-1) = data(year-1).row(pAG-1) + data(year-1).row(pAG-2); 
  if (MODEL==2) { // @debut age the debut pop
    data(year).row(0).zeros(); // clean for adding debut entrant
    data_db(year).rows(1, pDB-1) = data_db(year-1).rows(0, pDB-2);
    data_db(year).row(0).zeros(); // clean for adding debut entrant
  }
}

void popC::add_entrants () { // Add lagged births into youngest age group
  vec healthy, positiv;
  // if (exists("popadjust", where=p) & p.popadjust) {
  if (p.popadjust) {
    if (MODEL==0)
      healthy = p.entrantpop.col(year-1);
    else {
      healthy = p.entrantpop.col(year-1) % (1-p.entrantprev.col(year));
      positiv = p.entrantpop.col(year-1) %    p.entrantprev.col(year);
    }
  }
  else {
    if (MODEL==0)
      healthy = birthslag.col(year-1) % p.cumsurv.col(year-1) / 
                    p.paedsurv_lag(year-1) + p.cumnetmigr.col(year-1);
    else {
      vec survivedBirth = birthslag.col(year-1) % p.cumsurv.col(year-1);
      vec prev_now = p.entrantprev.col(year); // avoid repeat accesses
      vec imm_now  = p.cumnetmigr.col(year-1);
      healthy = survivedBirth % (1 - prev_now / p.paedsurv_lag.col(year-1)) + 
                imm_now % (1 - pregprevlag.col(year-1) * p.netmig_hivprob);
      positiv = (survivedBirth + imm_now) % prev_now;
    }
  }
  // save and update pop
  if (MODEL==0)
    data(year).slice(hivn_idx).row(0) = healthy.t();
  if (MODEL==1) {
    data(year).slice(hivn_idx).row(0) = healthy.t();
    data(year).slice(hivp_idx).row(0) = positiv.t();
  }
  if (MODEL==2) { // add to virgin then debut
    data_db(year).slice(hivn_idx).row(0) = healthy.t();
    data_db(year).slice(hivp_idx).row(0) = positiv.t();
  }
  if (MODEL!=0) {
    entrantprev(year)           = accu(positiv)/accu(healthy+positiv);
    hivp_entrants_out.col(year) = positiv;
  }
}

void popC::sexual_debut () {
  cube debut_now = data_db(year).each_slice() % p.db_pr;
  data(year).rows(0, pDB-1) += debut_now;
  data_db(year)             -= debut_now;
}

mat popC::hiv_aging_prob () {
    my_all(year-1); // now data_all refers to total pop: debut+not
    mat hiv_ag_prob = data_all.slice(hivp_idx).rows(aglast_idx) / 
                      sumByAG(data_all.slice(hivp_idx), ag_idx);
    hiv_ag_prob.replace(datum::nan, 0);
    return hiv_ag_prob;
}

vec popC::entrant_art () { // return these for updating HIV and ART pop
  my_all(year);
  vec out(4), hivp_here = data_all.slice(hivp_idx).row(0).t();
  out.rows(0,1) = hivp_here %      p.entrantartcov.col(year);
  out.rows(2,3) = hivp_here % (1 - p.entrantartcov.col(year));
  return out; // 1:2 ART+, 3:4 ART-
}

void popC::deaths () {
  cube death_now = data(year).each_slice() % (1 - p.Sx.slice(year));
  if (MODEL==1) {
    hiv_sx_prob = 1 - sumByAG(death_now.slice(hivp_idx), ag_idx) / 
                      sumByAG(data(year).slice(hivp_idx), ag_idx);
    hiv_sx_prob.replace(datum::nan, 0);
  }
  cube death_all;
  if (MODEL==2) {
    my_all(year);
    cube death_all = data_all.each_slice() % (1 - p.Sx.slice(year));
    hiv_sx_prob = 1 - sumByAG(death_all.slice(hivp_idx), ag_idx) / 
                      sumByAG(data_all.slice(hivp_idx), ag_idx);
    hiv_sx_prob.replace(datum::nan, 0);
    cube death_db = data_db(year).each_slice() % (1 - p.Sx.slice(year).rows(0, pDB));
    data_db(year) -= death_db;
  }
  data(year) -= death_now;
  if (MODEL!=2)
    natdeaths.slice(year) = sum(death_now, 2);
  else
    natdeaths.slice(year) = sum(death_all, 2);
}

void popC::migration () {
  mat netmigsurv = p.netmigr.slice(year) % (1 + p.Sx.slice(year)) / 2;
  mat mr_prob, current_pop;
  if (MODEL!=2) {
    mat current_pop = sum(data(year), 2);
    mr_prob = 1 + netmigsurv / current_pop;
  }
  if (MODEL==1) {
    mat n_mr_hivp = mr_prob % data(year).slice(hivp_idx);
    hiv_mr_prob = sumByAG(n_mr_hivp, ag_idx) / 
                  sumByAG(data(year).slice(hivp_idx), ag_idx);
    hiv_mr_prob.replace(datum::nan, 0);
  }
  if (MODEL==2) {
    my_all(year);
    mat current_pop = sum(data_all, 2);
    mr_prob = 1 + netmigsurv / current_pop;
    mat n_mr_hivp = mr_prob % data(year).slice(hivp_idx);
    hiv_mr_prob = sumByAG(n_mr_hivp, ag_idx) / 
                  sumByAG(data_all.slice(2), ag_idx);
    hiv_mr_prob.replace(datum::nan, 0);
    data_db(year).each_slice() %= mr_prob.rows(0, pDB);
  }
  data(year).each_slice() %= mr_prob;
}

void popC::update_fertile () { // only on active pop
  vec ave_fert_pop = (sum(data(year).tube(span(0, 34), span(f_idx)),2) + 
                    sum(data(year-1).tube(span(0, 34), span(f_idx)),2)) / 2;
  birth_age = ave_fert_pop % p.asfr.col(year);
  birth_agrp = sumByAG(birth_age, ag_idx(p_fert_idx));
  vec births = p.srb.col(year) * accu(birth_agrp);
  if (year + AGE_START <= PROJ_YEARS - 1)
    birthslag.col(year + AGE_START - 1) = births;
}

void popC::adjust_pop () {
  my_all(year);
  mat all_pop_by_sex = sum(data_all,2);
  popadjust.slice(year) = p.targetpop.slice(year) / all_pop_by_sex;
  if (MODEL!=0) {
    mat n_adjust = popadjust.slice(year) % data_all.slice(hivp_idx);
    adj_prob = sumByAG(n_adjust, ag_idx) / sumByAG(data_all.slice(hivp_idx), ag_idx);
    adj_prob.replace(datum::nan, 0);
  }
  data(year).each_slice() %= popadjust.slice(year);
  if (MODEL==2)
    data_db(year).each_slice() %= popadjust.slice(year).rows(0, pDB);
}

void popC::cal_prev_pregant (hivC& hivpop, artC& artpop) { // only on active pop
  vec n_mean = (data(year-1).slice(hivn_idx)(span(0, 34), span(hivn_idx)) + 
                data(year  ).slice(hivn_idx)(span(0, 34), span(hivn_idx))) / 2;
  vec hivn = sumByAG(n_mean, ag_idx(p_fert_idx));
  mat hivp = (hivpop.data(year-1).slice(f_idx).cols(h_fert_idx) + 
              hivpop.data(year  ).slice(f_idx).cols(h_fert_idx)) / 2;
  vec frp  = sum(p.frr_cd4.slice(year) % hivp).t();
  cube art = (artpop.data.row(year-1)(f_idx).slices(h_fert_idx) + 
              artpop.data.row(year-1)(f_idx).slices(h_fert_idx)) / 2;
  vec fra  = sum(sum(p.frr_art(year) % art));
  double pregprev = accu( birth_agrp % (1 - hivn / ( hivn + frp + fra)) ) / 
                    accu(birth_age);
  pregprevlag(year + AGE_START - 1) = pregprev;
}

void popC::save_prev_n_inc () {
  my_all(year);
  prev15to49(year) = accu(data_all.slice(hivp_idx).rows(p_age15to49_idx)) / 
                     accu(data_all.rows(as_scalar(p_age15to49_idx.head(1)),
                                        as_scalar(p_age15to49_idx.tail(1))));
  incid15to49(year) /= accu(data(year-1).slice(hivn_idx).rows(p_age15to49_idx));
  // prev(year) = accu(data_all.slice(hivp_idx)) / accu(data_all);
  // incid(year) = incid15to49(year) / accu(data(year-1).slice(hivn_idx)); // toBfixed
}

mat popC::infect_mix (uword ii) {
  uword ts = (year-1)/DT + ii;
  mat hivp_now = data(year).slice(hivp_idx);
  mat hiv_treated = hivp_now.each_row() % artcov.t();
  mat hiv_not_treated = hivp_now - hiv_treated;
  mat total_pop = sum(data(year), 2);
  mat transm_prev = (hiv_not_treated + hiv_treated * (1 - p.relinfectART)) / 
                    total_pop; // prevalence adjusted for art
  // +intervention effects and time epidemic start
  double w = (p.proj_steps(ts) == p.tsEpidemicStart) ? p.iota : 0;
  mat ir = rvec(ts) * transm_prev + w;
  // sweep over sexual mixing matrices
  vec ir_m = sum(p.mat_m.each_row() % ir.col(m_idx).t(), 1); // male
  vec ir_f = sum(p.mat_f.each_row() % ir.col(f_idx).t(), 1); // male
  mat ir_mf = join_rows(ir_m, ir_f);
  // if (exists("f_fun", fp)) // that fun
  //   ir = ir * fp.f_fun
  mat infections_ts = ir_mf % data(year).slice(hivn_idx);
  incrate15to49_ts_m.slice(ts) = ir_mf;
  // prev15to49_ts_m should use this one! now just store as below
  prev15to49_ts(ts) = accu(data(year).slice(hivp_idx)) / accu(data(year));
  prev_last = prev15to49_ts(ts);
  return infections_ts;
}

mat popC::infect_spec (hivC& hivpop, artC& artpop, uword time_step) {
  uword ts = (year-1)/ DT + time_step;
  double dt_ii = 1 - DT * time_step; // transition of population in 1 year
  uword first_age  = p_age15to49_idx(0),
        first_agrp = h_age15to49_idx(0);
  uword last_age   = as_scalar(p_age15to49_idx.tail(1)) + 1,
        last_agrp  = as_scalar(h_age15to49_idx.tail(1)) + 1;
  cube now = data(year);
  double hivn_ii = accu(now.slice(hivn_idx).rows(p_age15to49_idx)) -
                   accu(now.slice(hivn_idx).row(first_age)) * dt_ii +
                   accu(now.slice(hivn_idx).row(last_age))  * dt_ii;
  double hivp_ii = accu(now.slice(hivp_idx).rows(p_age15to49_idx)) -
                   accu(now.slice(hivp_idx).row(first_age)) * dt_ii +
                   accu(now.slice(hivp_idx).row(last_age)) * dt_ii;
  double art_ii = 0.0;
  for (uword sex = 0; sex < 2; ++sex)
    art_ii += accu(artpop.data.row(year)(sex).slices(h_age15to49_idx));
  if ((accu(hivpop.data(year).row(first_agrp)) + 
       accu(artpop.data.row(year)(0).slice(first_agrp)) + // female
       accu(artpop.data.row(year)(1).slice(first_agrp))) > 0) { // male
    rowvec art_first = {
      sum(sum(artpop.data.row(year)(0).slice(first_agrp))),
      sum(sum(artpop.data.row(year)(1).slice(first_agrp)))
    }, 
    hiv_first = sum(hivpop.data(year).col(first_agrp));
    art_ii -= accu( now.slice(hivp_idx).row(first_age) % // rowvec 2
        art_first / (hiv_first + art_first) ) * dt_ii;
  }
  if ((accu(hivpop.data(year).col(last_agrp)) + 
       accu(artpop.data.row(year)(0).slice(last_agrp)) +
       accu(artpop.data.row(year)(1).slice(last_agrp))) > 0) {
    rowvec art_last = {
      sum(sum(artpop.data.row(year)(0).slice(last_agrp))),
      sum(sum(artpop.data.row(year)(1).slice(last_agrp)))
    }, 
    hiv_last = sum(hivpop.data(year).col(last_agrp));
    art_ii += accu( now.slice(hivp_idx).row(last_age) %
        art_last / (hiv_last + art_last) ) * dt_ii;
  }
  double transm_prev = (hivp_ii - art_ii * (1 - p.relinfectART)) / 
                       (hivn_ii + hivp_ii);
  double w = (p.proj_steps(ts) == p.tsEpidemicStart) ? p.iota : 0.0;
  double incrate15to49_ts_now = rvec(ts) * transm_prev + w;
    
  // Incidence: male = negative / female negative * sexratio + male negative; 
  //          female = male * sexratio
  mat sus_by_age_sex = now.slice(hivn_idx).rows(p_age15to49_idx);
  double adj_sex =  accu(sus_by_age_sex) /
    ( accu(now.slice(hivn_idx)( span(first_age, last_age-1), span(m_idx) )) +
      accu(now.slice(hivn_idx)( span(first_age, last_age-1), span(f_idx) )) * 
      p.incrr_sex(year) );
  rowvec sexinc15to49_ts {1, p.incrr_sex(year)};
  sexinc15to49_ts *= incrate15to49_ts_now * adj_sex;

  // New infections distributed by age: ratio age_i/ 25-29 age
  rowvec adj_age = sexinc15to49_ts /
    ( sum(sus_by_age_sex % p.incrr_age.slice(year).rows(p_age15to49_idx) ) / 
      sum(sus_by_age_sex) );
  mat agesex_inc = p.incrr_age.slice(year).each_row() % adj_age;

  // Adjust age-specific incidence among men for circumcision coverage
  agesex_inc.col(m_idx) %= (1 - p.circ_incid_rr * p.circ_prop.col(year));
  mat infections_ts = agesex_inc % now.slice(hivn_idx);
  // saving
  incrate15to49_ts(ts) = incrate15to49_ts_now;
  prev_last = hivp_ii / (hivn_ii + hivp_ii);
  prev15to49_ts(ts) = prev_last;
  return infections_ts;
}

// epp_disease_model_direct (hivpop, artpop) {
//   if (p.incidpopage == 0L) // incidence for 15 -49 population
//     age_id <- p_age15to49_idx
//   else if (p.incidpopage == 1L) // incidence for 15+ population
//     age_id <- p.age15plus_idx
//   num <- c(1, p.incrr_sex[year]) * sum(data[age_id,, hivn_idx, year-1])
//   den <- sum(data[age_id, m_idx, hivn_idx, year-1]) + 
//          sum(data[age_id, f_idx, hivn_idx, year-1]) * p.incrr_sex[year]
//   sexinc <- p.incidinput[year] * num / den
//   ageinc <- colSums(data[age_id,,hivn_idx,year-1] * p.incrr_age[age_id,,year]) /
//             colSums(data[age_id,,hivn_idx,year-1])

//   agesex.inc <- sweep(p.incrr_age[,,year], 2, sexinc / ageinc, "*")
//   infections[,,year]     <<- agesex.inc * data[,,hivn_idx,year-1]
//   data[,,hivn_idx,year]  <<- data[,,hivn_idx,year] - infections[,,year]
//   data[,,hivp_idx,year]  <<- data[,,hivp_idx,year] + infections[,,year]
//   infect_agrp            <<- sumByAGs(infections[,,year], ag_idx)
//   hivpop.data[,,,year] <- hivpop.data[,,,year] + 
//                              sweep(p.cd4_initdist, 2:3, infect_agrp, "*")
//   incid15to49[year]      <<- sum(infections[p_age15to49_idx,,year])
// }

double popC::calc_rtrend_rt (uword ts, double time_step) {
  double rveclast = rvec(ts-1);
  // double dtii =  1 - DT * (time_step - 1);
  // hivn.ii <- sum(data[p_age15to49_idx,,hivn_idx,year]) - 
  //            sum(data[p_age15to49_idx[1],,hivn_idx,year]) * dtii + 
  //            sum(data[tail(p_age15to49_idx,1)+1,,hivn_idx,year]) * dtii

  // hivp.ii <- sum(data[p_age15to49_idx,,hivp_idx,year]) -
  //            sum(data[p_age15to49_idx[1],,hivp_idx,year]) * dtii +
  //            sum(data[tail(p_age15to49_idx,1)+1,,hivp_idx,year]) * dtii
  // prevcurr <- hivp.ii / (hivn.ii + hivp.ii)
  // t.ii     <- p.proj.steps[ts]
  // if (t.ii > p.tsEpidemicStart){
  //   par <- p.rtrend
  //   if (t.ii < par$tStabilize)
  //     gamma.t <- 0
  //   else 
  //     gamma.t <- (prevcurr-prevlast) * (t.ii - par$tStabilize) / (DT * prevlast)
  //   logr.diff <- par$beta[2] * (par$beta[1] - rveclast) +
  //                par$beta[3] * prevlast + 
  //                par$beta[4] * gamma.t
  //   return(exp(log(rveclast) + logr.diff))
  // }
  // else
  //   return(p.rtrend$r0)
  return rveclast;
}

void popC::update_rvec (double time_step) {
  uword ts = (year-1) / DT + time_step;
  // if (p.eppmod %in% c("rtrend", "rtrend_rw")) <<- need to be fixed in FP
  if (p.eppmod == 2)
    rvec(ts) = calc_rtrend_rt(ts, time_step);
  else
    rvec(ts) = p.rvec(ts);
}

void popC::update_infection (mat infect) {
    mat change_ts = infect * DT;
    data(year).slice(hivn_idx) -= change_ts;
    data(year).slice(hivp_idx) += change_ts;
    infections.slice(year) += change_ts;
    incid15to49(year) += accu(change_ts.rows(p_age15to49_idx));
}

void popC::remove_hiv_death (cube cd4_mx, hivC& hivpop, artC& artpop) {
  // death by age group
  field<cube> artmx_now = p.art_mort; // 3x7x9x2
  for (uword sex = 0; sex < 2; ++sex)
    for (int agr = 0; agr < artmx_now(0).n_slices; ++agr)
      artmx_now(sex).slice(agr).each_col() %= p.artmx_timerr.col(year);
  field<cube> artD = p.art_mort;
  for (uword sex = 0; sex < 2; ++sex)
    artD(sex) = artpop.data.row(year)(sex) % artmx_now(sex);
  cube hivD = cd4_mx % hivpop.data(year);
  if (MODEL==2) { // add deaths from inactive population
    for (uword sex = 0; sex < 2; ++sex)
      artD(sex) += artpop.data_db.row(year)(sex) % artmx_now(sex);
    hivD += cd4_mx % hivpop.data_db(year);
  }
  mat hivDbyAG = sum(hivD), artDbyAG = hivDbyAG;
  for (uword sex = 0; sex < 2; ++sex)
    artDbyAG.col(sex) = col_sum_cube_to_vec (artD(sex));
  mat dbyAG = DT * ( hivDbyAG + artDbyAG );
  // deaths by single-year
  mat pA = data(year).slice(hivp_idx) / 
           rep_each_row_by_span(sumByAG(data(year).slice(hivp_idx), ag_idx),
                                h_ag_span);
  pA.replace(datum::nan, 0);
  mat dbyA = rep_each_row_by_span(dbyAG, h_ag_span) % pA;
  data(year).slice(hivp_idx) -= dbyA;
  hivdeaths.slice(year) += dbyA;
}

cube popC::update_preg (cube art_elig, hivC& hivpop, artC& artpop) {
  mat hivp = hivpop.data(year).slice(f_idx).cols(h_fert_idx) % p.frr_cd4.slice(year);
  int a_lo = as_scalar(p_fert_idx.head(1)), a_up = as_scalar(p_fert_idx.tail(1));
  vec hivn_vec = data(year).slice(hivn_idx)(span(a_lo, a_up), span(f_idx));
  vec hivn = sumByAG(hivn_vec, ag_idx(p_fert_idx));
  vec art_now = col_sum_cube_to_vec(
    artpop.data.row(year)(f_idx).slices(h_fert_idx) % p.frr_art(year));
  vec hiv_bagr = birth_agrp / (hivn + conv_to<vec>::from(sum(hivp)) + art_now);
  mat birthdist = hivp.each_row() % conv_to<rowvec>::from(hiv_bagr);
  uvec elDS = linspace<uvec>(0, p.artcd4elig_idx(year) - 1, p.artcd4elig_idx(year));
  uvec elAG = h_fert_idx - min(h_age15plus_idx);
  art_elig.slice(f_idx)(elDS, elAG) += birthdist.rows(elDS);
  return art_elig;
}

vec popC::artInit (vec art_curr, cube art_elig, int time_step) {
  vec out = zeros(2);
  int year_w = (DT * (time_step + 1) < 0.5) ? 0 : 1;
  double trans = DT * (time_step + 1) + 0.5 - year_w;
  vec transition {1 - trans, trans};
  uword year_l = year - (2 - year_w), year_r = year - (1 - year_w);
  for (uword sex = 0; sex < 2; ++sex) {
    if( !any(p.art15plus_isperc(sex, span(year_l, year_r))) ) { // both number
      out(sex) = accu(p.art15plus_num(sex, span(year_l, year_r)) * transition);
      if (MIX) {
        double cov = out(sex) / (accu(art_elig.slice(sex)) + art_curr(sex));
        artcov(sex) = (cov < 1) ? cov : 1;
      }
    } 
    else if (all(p.art15plus_isperc(sex, span(year_l, year_r)) == 1)) { // both percentage
      double cov = accu(p.art15plus_num(sex, span(year_l, year_r)) * transition);
      if (MIX)
        artcov(sex) = (cov < 1) ? cov : 1;
      out(sex) = cov * (accu(art_elig.slice(sex)) + art_curr(sex));
    }
    else if ((p.art15plus_isperc(sex, year_l)==0) & 
             (p.art15plus_isperc(sex, year_r)==1)) { // transition number to percentage
      double curr_n_cov = art_curr(sex) / (accu(art_elig.slice(sex)) + art_curr(sex));
      double diff_n_cov = p.art15plus_num(sex, year_r) - curr_n_cov;
      double cov = curr_n_cov + diff_n_cov * DT / (0.5 + year_w - DT * time_step);
      if (MIX)
        artcov(sex) = (cov < 1) ? cov : 1;
      out(sex) = cov * ( accu(art_elig.slice(sex)) + art_curr(sex) );
    }
  }
  return out;
} 

// calculate ART initiation distribution
cube popC::artDist (cube art_elig, vec art_need) {
  cube art_real = zeros(size(art_elig));
  if (!p.med_cd4init_input(year)) {
    if (p.art_alloc_method == 4L) { // by lowest CD4
      // Calculate proportion to be initiated in each CD4 category
      vec init_pr(2);
      mat current_m(art_elig.n_cols, art_elig.n_slices);
      for (uword m = hDS - 1; m > 0; --m) {
        vec elig_hm = sum(art_elig.row(m));
        if (all(elig_hm == 0))
          init_pr = elig_hm;
        else {
          init_pr = art_need / elig_hm;
          init_pr.elem( find(init_pr > 1) ).ones();
          init_pr.replace(datum::nan, 1); // pmin(x, 1, na.rm)
        }
        current_m = art_elig.row(m);
        art_real.row(m) = current_m.each_row() % init_pr;
      }
    } 
    else { // Spectrum Manual p168--p169, 
      int A = as_scalar(h_age15plus_idx.head(1)),
          Z = as_scalar(h_age15plus_idx.tail(1)); // add these in fp for easier
      cube expect_mort_w = p.cd4_mort.cols(A, Z);
      vec artmort = col_sum_cube_to_vec( art_elig % expect_mort_w );
      for (uword sex = 0; sex < 2; ++sex) expect_mort_w.slice(sex) /= artmort(sex);
      cube init_w = p.art_alloc_mxweight * expect_mort_w;
      vec mx_w = (1 - p.art_alloc_mxweight) / col_sum_cube_to_vec(art_elig);
      for (uword sex = 0; sex < 2; ++sex) init_w.slice(sex) += mx_w(sex);
      art_real = init_w % art_elig;
      for (uword sex = 0; sex < 2; ++sex) 
        art_real.slice(sex) *= art_need(sex);
      uvec greater = find(art_real > art_elig);
      art_real(greater) = art_elig(greater);
    }
  }
  else {
    vec CD4_LO = {500, 350, 250, 200, 100, 50, 0};
    vec CD4_UP = {1000, 500, 350, 250, 200, 100, 50};
    uword j = p.med_cd4init_cat(year) - 1; // R to C++
    double pr_below = (p.median_cd4init(year) - CD4_LO(j)) / 
                      (CD4_UP(j) - CD4_LO(j));
    vec elig_below = sum(art_elig.row(j)) * pr_below;
    if (j < (hDS - 1))
      elig_below += sum(art_elig.rows(j+1, hDS-1));
    vec elig_above = sum(art_elig.row(j)) * (1.0 - pr_below);
    if (j > 1) 
      elig_above += sum(art_elig.rows(0, j-1));

    vec initpr_below = art_need * 0.5 / elig_below;
    initpr_below = p_min(initpr_below, 1);
    vec initpr_above = art_need * 0.5 / elig_above;
    initpr_above = p_min(initpr_above, 1);
    vec initpr_medcat = initpr_below * pr_below + initpr_above * (1-pr_below);
    if (j < (hDS - 1)) {
      for (uword sex = 0; sex < 2; ++sex)
        art_real.slice(sex).rows(j+1, hDS-1) = 
          art_elig.slice(sex).rows(j+1, hDS-1) * initpr_below(sex);
    }
    for (uword sex = 0; sex < 2; ++sex)
      art_real.slice(sex).row(j) = 
        art_elig.slice(sex).row(j) * initpr_medcat(sex);
    if (j > 0) { // ? the condition might not correct
      for (uword sex = 0; sex < 2; ++sex)
        art_real.slice(sex).rows(0, j-1) = 
          art_elig.slice(sex).rows(0, j-1) * initpr_above(sex);
    }
  }
  return art_real;
}

cube popC::scale_cd4_mort (hivC& hivpop, artC& artpop) {
  cube cd4_mort_ts;
  if (p.scale_cd4_mort) {
    cube num = hivpop.data(year) + hivpop.data_db(year);
    cube den(size(num));
    for (uword sex = 0; sex < 2; ++sex)
      den.slice(sex) = sum(artpop.data.row(year)(sex) + 
                           artpop.data_db.row(year)(sex));
    cube cd4mx = num / (num + den);
    cd4mx.elem( find_nonfinite(cd4mx) ).ones();
    cd4_mort_ts = cd4mx % p.cd4_mort;
  } 
  else
    cd4_mort_ts = p.cd4_mort;
  return cd4_mort_ts;
}

void popC::finalize (hivC& hivpop, artC& artpop) {
  SET_ATTR(data_sexp, Rf_install("infections"), infections_sexp);
  SET_ATTR(data_sexp, Rf_install("hivdeaths"), hivdeaths_sexp);
  SET_ATTR(data_sexp, Rf_install("natdeaths"), natdeaths_sexp);
  SET_ATTR(data_sexp, Rf_install("popadjust"), popadjust_sexp);
  SET_ATTR(data_sexp, Rf_install("pregprevlag"), pregprevlag_sexp);
  SET_ATTR(data_sexp, Rf_install("incrate15to49_ts"), incrate15to49_ts_sexp);
  SET_ATTR(data_sexp, Rf_install("prev15to49_ts"), prev15to49_ts_sexp);
  SET_ATTR(data_sexp, Rf_install("rvec_ts"), rvec_sexp);
  SET_ATTR(data_sexp, Rf_install("prev15to49"), prev15to49_sexp);
  SET_ATTR(data_sexp, Rf_install("incid15to49"), incid15to49_sexp);
  SET_ATTR(data_sexp, Rf_install("entrantprev"), entrantprev_sexp);
  if (MODEL!=0) {
    SET_ATTR(data_sexp, Rf_install("hivpop"), hivpop.data_sexp);
    SET_ATTR(data_sexp, Rf_install("artpop"), artpop.data_sexp);
    SET_CLASS(data_sexp, Rf_mkString("spec"));
  }
  else 
    SET_CLASS(data_sexp, Rf_mkString("dempp"));
}