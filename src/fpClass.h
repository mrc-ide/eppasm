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
#include "arma.h"
#pragma once
// Keep the same order as in R to keep track
// 
// All *_idx parameters are shift to zero index
// 
class fp_rt {
public:
  void init_me(SEXP fp);
  fp_rt(){};
public:
  vec proj_steps;
  double rw_start;
  double rw_trans;
  vec rlogistic_steps;
  vec rw_steps;
  double n_rw;
  double rw_dk;
  vec rw_knots;
  uvec rw_idx;
  double n_param;
  vec rw_transition;
  // dt             : num 0.1
  // eppmod         : chr "rhybrid"
};

class fp_main {
public:
  void init_me(SEXP fp);
  fp_main(){};
public:
  int         SIM_YEARS;
  vec         proj_steps;
  mat         basepop;
  cube        Sx;
  mat         asfr;
  mat         srb;
  cube        netmigr;
  mat         birthslag;
  mat         cumsurv;
  mat         cumnetmigr;
  bool        popadjust;
  mat         entrantpop;
  cube        targetpop;
  vec         births;
  double      relinfectART;
  vec         incrr_sex;
  cube        incrr_age;
  cube        cd4_initdist;
  cube        cd4_prog;
  cube        cd4_mort;
  field<cube> art_mort;
  mat         artmx_timerr;
  cube        frr_cd4;
  field<cube> frr_art;
  mat         art15plus_num;
  mat        art15plus_isperc;
  vec         specpop_percelig;
  uvec        artcd4elig_idx;
  vec         pw_artelig;
  double      who34percelig;
  vec         art_dropout;
  vec         median_cd4init;
  uvec        med_cd4init_input;
  uvec        med_cd4init_cat;
  int         tARTstart;
  int         art_alloc_method;
  double      art_alloc_mxweight;
  int         scale_cd4_mort;
  vec         verttrans_lag;
  vec         paedsurv_lag;
  mat         entrantprev;
  mat         entrantartcov;
  cube        paedsurv_cd4dist;
  field<cube> paedsurv_artcd4dist;
  double      netmig_hivprob;
  double      netmighivsurv;
  double      circ_incid_rr;
  mat         circ_prop;
  double      tsEpidemicStart;
  int         eppmod;
  bool        ancsitedata;
  int         ancrt;
  bool        logitiota;
  double      rw_start;
  int         incidmod;
  vec         rvec;
  double      iota;
  double      ancbias;
  double      v_infl;
  double      log_frr_adjust;
  double      ancrtcens_vinfl;
  double      ancrtsite_beta;
  mat         mat_m;
  mat         mat_f;
  mat         db_pr;
  // rt parameters
  fp_rt       rt;
};

// Master parameters class
class CeppFP {
public:
  CeppFP(SEXP fp);
  // ~CeppFP();
public:
  int      proj_start;
  int      PROJ_YEARS;
  int      AGE_START;
  int      hiv_steps_per_year;
  double   time_epi_start;
  double   NG;
  double   pDS;
  double   m_idx;
  double   f_idx;
  double   hivn_idx;
  double   hivp_idx;
  double   pAG;
  double   ag_rate;
  uvec     p_fert_idx;
  uvec     p_age15to49_idx;
  uvec     p_age15plus_idx;
  vec      h_ag_span;
  int      hAG;
  double   hDS;
  double   hTS;
  uvec     ag_idx;
  uvec     agfirst_idx;
  uvec     aglast_idx;
  uvec     h_fert_idx;
  uvec     h_age15to49_idx;
  uvec     h_age15plus_idx;
  double   DT;
  double   pDB;
  // vec      db_aid;
  fp_main p;
};