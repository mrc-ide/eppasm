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
#include "utils.hpp"
#pragma once
// Keep the same order as in R to keep track
// 
class fp_rt {
public:
  void init_me(SEXP fp);
  fp_rt() {};
public:
  double * proj_steps;
  double rw_start;
  double rw_trans;
  double * rlogistic_steps;
  double * rw_steps;
  double n_rw;
  double rw_dk;
  double * rw_knots;
  int * rw_idx;
  double n_param;
  double * rw_transition;
  // dt             : num 0.1
  // eppmod         : chr "rhybrid"
};

class fp_main {
public:
  fp_main(SEXP fp);
public:
  int         SIM_YEARS;
  double *    proj_steps;
  boost2D_ptr basepop;
  boost2D_ptr asfr;
  boost2D_ptr srb;
  boost2D_ptr birthslag;
  boost2D_ptr cumsurv;
  boost2D_ptr cumnetmigr;
  boost2D_ptr entrantpop;
  boost2D_ptr artmx_timerr;
  boost2D_ptr art15plus_num;
  boost2D_ptr art15plus_isperc;
  boost2D_ptr entrantprev;
  boost2D_ptr entrantartcov;
  boost2D_ptr circ_prop;
  boost2D_ptr mat_m;
  boost2D_ptr mat_f;
  boost2D_ptr db_pr;
  boost3D_ptr Sx;
  boost3D_ptr netmigr;
  boost3D_ptr targetpop;
  boost3D_ptr incrr_age;
  boost3D_ptr cd4_initdist;
  boost3D_ptr cd4_prog;
  boost3D_ptr cd4_mort;
  boost3D_ptr frr_cd4;
  boost3D_ptr paedsurv_cd4dist;
  boost4D_ptr art_mort;
  boost4D_ptr frr_art;
  boost4D_ptr paedsurv_artcd4dist;
  bool        popadjust;
  double *    births;
  double      relinfectART;
  double *    incrr_sex;
  double *    specpop_percelig;
  int *       artcd4elig_idx;
  double *    pw_artelig;
  double      who34percelig;
  double *    art_dropout;
  double *    median_cd4init;
  int *       med_cd4init_input;
  int *       med_cd4init_cat;
  int         tARTstart;
  int         art_alloc_method;
  double      art_alloc_mxweight;
  int         scale_cd4_mort;
  double *    verttrans_lag;
  double *    paedsurv_lag;
  double      netmig_hivprob;
  double      netmighivsurv;
  double      circ_incid_rr;
  double      tsEpidemicStart;
  int         eppmod;
  double *    incidinput;
  int         incidpopage = 0;
  bool        ancsitedata;
  int         ancrt;
  bool        logitiota;
  double      rw_start;
  int         incidmod;
  double *    rvec;
  double      iota;
  double      ancbias;
  double      v_infl;
  double      log_frr_adjust;
  double      ancrtcens_vinfl;
  double      ancrtsite_beta;
  // rt parameters
  fp_rt       rt;
};

// Master parameters class
class CeppFP {
public:
  CeppFP(SEXP fp);
public:
  SEXP         fp_ss;
  int          proj_start;
  int          PROJ_YEARS;
  int          AGE_START;
  int          hiv_steps_per_year;
  double       time_epi_start;
  int          NG;
  int          pDS;
  int          m_idx;
  int          f_idx;
  int          hivn_idx;
  int          hivp_idx;
  int          pAG;
  double       ag_rate;
  boost1I_ptr  p_fert_idx;
  boost1I_ptr  p_age15to49_idx;
  boost1I_ptr  p_age15plus_idx;
  boost1D_ptr  h_ag_span;
  int          hAG;
  int          hDS;
  int          hTS;
  boost1I_ptr  ag_idx;
  boost1I_ptr  agfirst_idx;
  boost1I_ptr  aglast_idx;
  boost1I_ptr  h_fert_idx;
  boost1I_ptr  h_age15to49_idx;
  boost1I_ptr  h_age15plus_idx;
  double       DT;
  int          pDB;
  int          hDB;
  fp_main      p;
  int          n_steps;
  // 
  int pAG_FERT, hAG_FERT, pAG_1549, hAG_1549, pAG_15plus, hAG_15plus;
};