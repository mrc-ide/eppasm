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
#include "fpClass.hpp"

void fp_rt::init_me(SEXP fp) {
  SEXP fp_rt      = get_value(fp, "rt");
  proj_steps      = REAL(get_value(fp_rt, "proj_steps"));
  rw_start        = *REAL(get_value(fp_rt, "rw_start"));
  rw_trans        = *REAL(get_value(fp_rt, "rw_trans"));
  rlogistic_steps = REAL(get_value(fp_rt, "rlogistic_steps"));
  rw_steps        = REAL(get_value(fp_rt, "rw_steps"));
  n_rw            = *REAL(get_value(fp_rt, "n_rw"));
  rw_dk           = *REAL(get_value(fp_rt, "rw_dk"));
  rw_knots        = REAL(get_value(fp_rt, "rw_knots"));
  rw_idx          = INTEGER(get_value(fp_rt, "rw_idx"));
  n_param         = *REAL(get_value(fp_rt, "n_param"));
  rw_transition   = REAL(get_value(fp_rt, "rw_transition"));
}

fp_main::fp_main(SEXP fp) :
// init list for boost array class
  basepop(REAL(get_value(fp, "basepop")), get_dim_2D(fp, "basepop")),
  asfr(REAL(get_value(fp, "asfr")), get_dim_2D(fp, "asfr")),
  srb(REAL(get_value(fp, "srb")), get_dim_2D(fp, "srb")),
  birthslag(REAL(get_value(fp, "birthslag")), get_dim_2D(fp, "birthslag")),
  cumsurv(REAL(get_value(fp, "cumsurv")), get_dim_2D(fp, "cumsurv")),
  cumnetmigr(REAL(get_value(fp, "cumnetmigr")), get_dim_2D(fp, "cumnetmigr")),
  entrantpop(REAL(get_value(fp, "entrantpop")), get_dim_2D(fp, "entrantpop")),
  artmx_timerr(REAL(get_value(fp, "artmx_timerr")),
               get_dim_2D(fp, "artmx_timerr")),
  art15plus_num(REAL(get_value(fp, "art15plus_num")),
                get_dim_2D(fp, "art15plus_num")), 
  art15plus_isperc(REAL(get_value(fp, "art15plus_isperc")),
                   get_dim_2D(fp, "art15plus_isperc")), 
  entrantprev(REAL(get_value(fp, "entrantprev")),
              get_dim_2D(fp, "entrantprev")), 
  entrantartcov(REAL(get_value(fp, "entrantartcov")),
                get_dim_2D(fp, "entrantartcov")), 
  circ_prop(REAL(get_value(fp, "circ_prop")), get_dim_2D(fp, "circ_prop")),
  mat_m(REAL(get_value(fp, "mat_m")), get_dim_2D(fp, "mat_m")),
  mat_f(REAL(get_value(fp, "mat_f")), get_dim_2D(fp, "mat_f")),
  db_pr(REAL(get_value(fp, "db_pr")), get_dim_2D(fp, "db_pr")),
  Sx(REAL(get_value(fp, "Sx")), get_dim_3D(fp, "Sx")),
  netmigr(REAL(get_value(fp, "netmigr")), get_dim_3D(fp, "netmigr")),
  targetpop(REAL(get_value(fp, "targetpop")), get_dim_3D(fp, "targetpop")),
  incrr_age(REAL(get_value(fp, "incrr_age")), get_dim_3D(fp, "incrr_age")),
  cd4_initdist(REAL(get_value(fp, "cd4_initdist")),
               get_dim_3D(fp, "cd4_initdist")),
  cd4_prog(REAL(get_value(fp, "cd4_prog")), get_dim_3D(fp, "cd4_prog")),
  cd4_mort(REAL(get_value(fp, "cd4_mort")), get_dim_3D(fp, "cd4_mort")),
  frr_cd4(REAL(get_value(fp, "frr_cd4")), get_dim_3D(fp, "frr_cd4")),
  paedsurv_cd4dist(REAL(get_value(fp, "paedsurv_cd4dist")),
                   get_dim_3D(fp, "paedsurv_cd4dist")),
  art_mort(REAL(get_value(fp, "art_mort")), get_dim_4D(fp, "art_mort")),
  frr_art(REAL(get_value(fp, "frr_art")), get_dim_4D(fp, "frr_art")),
  paedsurv_artcd4dist(REAL(get_value(fp, "paedsurv_artcd4dist")),
                      get_dim_4D(fp, "paedsurv_artcd4dist"))
  {
// Non class init
    SIM_YEARS          = *INTEGER(get_value(fp, "SIM_YEARS"));
    proj_steps         = REAL(get_value(fp, "proj_steps"));
    popadjust          = *LOGICAL(get_value(fp, "popadjust"));
    births             = REAL(get_value(fp, "births"));
    relinfectART       = *REAL(get_value(fp, "relinfectART"));
    incrr_sex          = REAL(get_value(fp, "incrr_sex"));         
    specpop_percelig   = REAL(get_value(fp, "specpop_percelig"));         
    artcd4elig_idx     = INTEGER(get_value(fp, "artcd4elig_idx"));         
    pw_artelig         = REAL(get_value(fp, "pw_artelig"));         
    who34percelig      = *REAL(get_value(fp, "who34percelig"));      
    art_dropout        = REAL(get_value(fp, "art_dropout"));         
    median_cd4init     = REAL(get_value(fp, "median_cd4init"));         
    med_cd4init_input  = INTEGER(get_value(fp, "med_cd4init_input"));
    med_cd4init_cat    = INTEGER(get_value(fp, "med_cd4init_cat"));
    tARTstart          = *INTEGER(get_value(fp, "tARTstart"));
    art_alloc_method   = *INTEGER(get_value(fp, "art_alloc_method"));
    art_alloc_mxweight = *REAL(get_value(fp, "art_alloc_mxweight"));
    scale_cd4_mort     = *INTEGER(get_value(fp, "scale_cd4_mort"));
    verttrans_lag      = REAL(get_value(fp, "verttrans_lag"));
    paedsurv_lag       = REAL(get_value(fp, "paedsurv_lag"));
    netmig_hivprob     = *REAL(get_value(fp, "netmig_hivprob"));
    netmighivsurv      = *REAL(get_value(fp, "netmighivsurv"));
    circ_incid_rr      = *REAL(get_value(fp, "circ_incid_rr"));
    if (has_value(fp, "tsEpidemicStart"))
      tsEpidemicStart  = *REAL(get_value(fp, "tsEpidemicStart"));
    eppmod             = *INTEGER(get_value(fp, "eppmodInt")); // char to int
    if (eppmod==2) {  // direct incidence input
      if (has_value(fp, "incidpopage")) {
        incidinput     = REAL(get_value(fp, "incidinput"));
        incidpopage    = *INTEGER(get_value(fp, "incidpopage"));
      }
    }
    if (has_value(fp, "ancsitedata")) {
      ancsitedata      = *LOGICAL(get_value(fp, "ancsitedata"));
      ancrt            = *INTEGER(get_value(fp, "ancrtInt")); // char to int
    }
    if (has_value(fp, "logitiota")) 
      logitiota        = *LOGICAL(get_value(fp, "logitiota"));
    if (has_value(fp, "rw_start")) 
      rw_start         = *REAL(get_value(fp, "rw_start"));
    incidmod           = *INTEGER(get_value(fp, "incidmodInt")); // char to int
    if (has_value(fp, "rvec"))
      rvec             = REAL(get_value(fp, "rvec"));
    if (has_value(fp, "iota"))
      iota             = *REAL(get_value(fp, "iota"));
    if (has_value(fp, "ancbias")) {// double
      ancbias          = *REAL(get_value(fp, "ancbias"));      
      v_infl           = *REAL(get_value(fp, "v_infl"));      
      log_frr_adjust   = *REAL(get_value(fp, "log_frr_adjust"));      
      ancrtcens_vinfl  = *REAL(get_value(fp, "ancrtcens_vinfl"));      
      ancrtsite_beta   = *REAL(get_value(fp, "ancrtsite_beta"));
    }
    // rt parameters
    if (has_value(fp, "rt")) // direct incidence model does not have rt
      rt.init_me(fp);
}

CeppFP::CeppFP(SEXP fp) : 
// Boost array class init
  fp_ss(get_value(fp, "ss")),
  p_fert_idx(INTEGER(get_value(fp_ss, "p_fert_idx")),
      extents[GET_LENGTH(get_value(fp_ss, "p_fert_idx"))]),
  p_age15to49_idx(INTEGER(get_value(fp_ss, "p_age15to49_idx")),
      extents[GET_LENGTH(get_value(fp_ss, "p_age15to49_idx"))]),
  p_age15plus_idx(INTEGER(get_value(fp_ss, "p_age15plus_idx")),
      extents[GET_LENGTH(get_value(fp_ss, "p_age15plus_idx"))]),
  h_ag_span(REAL(get_value(fp_ss, "h_ag_span")),
      extents[GET_LENGTH(get_value(fp_ss, "h_ag_span"))]),
  ag_idx(INTEGER(get_value(fp_ss, "ag_idx")),
      extents[GET_LENGTH(get_value(fp_ss, "ag_idx"))]),
  agfirst_idx(INTEGER(get_value(fp_ss, "agfirst_idx")),
      extents[GET_LENGTH(get_value(fp_ss, "agfirst_idx"))]),
  aglast_idx(INTEGER(get_value(fp_ss, "aglast_idx")),
      extents[GET_LENGTH(get_value(fp_ss, "aglast_idx"))]),
  h_fert_idx(INTEGER(get_value(fp_ss, "h_fert_idx")),
      extents[GET_LENGTH(get_value(fp_ss, "h_fert_idx"))]),
  h_age15to49_idx(INTEGER(get_value(fp_ss, "h_age15to49_idx")),
      extents[GET_LENGTH(get_value(fp_ss, "h_age15to49_idx"))]),
  h_age15plus_idx(INTEGER(get_value(fp_ss, "h_age15plus_idx")),
      extents[GET_LENGTH(get_value(fp_ss, "h_age15plus_idx"))]),
  p(fp)
  { // all .names --> _names
// Non class model state space
  proj_start         = *INTEGER(get_value(fp_ss, "proj_start"));
  PROJ_YEARS         = *INTEGER(get_value(fp_ss, "PROJ_YEARS"));
  AGE_START          = *INTEGER(get_value(fp_ss, "AGE_START"));
  hiv_steps_per_year = *INTEGER(get_value(fp_ss, "hiv_steps_per_year"));
  if (has_value(fp_ss, "time_epi_start"))
    time_epi_start     = *REAL(get_value(fp_ss, "time_epi_start"));
  NG                 = (int) *REAL(get_value(fp_ss, "NG"));
  pDS                = (int) *REAL(get_value(fp_ss, "pDS"));
  m_idx              = (int) *REAL(get_value(fp_ss, "m_idx")) - 1;
  f_idx              = (int) *REAL(get_value(fp_ss, "f_idx")) - 1;
  hivn_idx           = (int) *REAL(get_value(fp_ss, "hivn_idx")) - 1;
  hivp_idx           = (int) *REAL(get_value(fp_ss, "hivp_idx")) - 1;
  pAG                = (int) *REAL(get_value(fp_ss, "pAG"));
  ag_rate            = *REAL(get_value(fp_ss, "ag_rate"));
  hAG                = *INTEGER(get_value(fp_ss, "hAG"));
  hDS                = (int) *REAL(get_value(fp_ss, "hDS"));
  hTS                = (int) *REAL(get_value(fp_ss, "hTS"));
  DT                 = *REAL(get_value(fp_ss, "DT"));
  if (has_value(fp_ss, "pDB")) {
    pDB            = *INTEGER(get_value(fp_ss, "pDB"));
    hDB            = pDB; // single-year sexual debut pop
  }
  n_steps          = (PROJ_YEARS-1) * hiv_steps_per_year;
  pAG_FERT         = (p_fert_idx[0] - 1) + p_fert_idx.num_elements();
  hAG_FERT         = (h_fert_idx[0] - 1) + h_fert_idx.num_elements();
  pAG_1549         = (p_age15to49_idx[0] - 1) + p_age15to49_idx.num_elements();
  hAG_1549         = (h_age15to49_idx[0] - 1) + h_age15to49_idx.num_elements();
  pAG_15plus       = (p_age15plus_idx[0] - 1) + p_age15plus_idx.num_elements();
  hAG_15plus       = (h_age15plus_idx[0] - 1) + h_age15plus_idx.num_elements();
}