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
#include "fpClass.h"

void fp_rt::init_me(SEXP fp) {
  SEXP fp_rt = get_value(fp, "rt");
  proj_steps = arma_vec(get_value(fp_rt, "proj_steps")); // vec
  rw_start = *REAL(get_value(fp_rt, "rw_start")); // double
  rw_trans = *REAL(get_value(fp_rt, "rw_trans")); // double
  rlogistic_steps = arma_vec(get_value(fp_rt, "rlogistic_steps")); // vec
  rw_steps = arma_vec(get_value(fp_rt, "rw_steps")); // vec
  n_rw = *REAL(get_value(fp_rt, "n_rw")); // double
  rw_dk = *REAL(get_value(fp_rt, "rw_dk")); // double
  rw_knots = arma_vec(get_value(fp_rt, "rw_knots")); // vec
  rw_idx = arma_uvec(get_value(fp_rt, "rw_idx")); // int vec
  n_param = *REAL(get_value(fp_rt, "n_param")); // double
  rw_transition = arma_vec(get_value(fp_rt, "rw_transition")); // vec
};

void fp_main::init_me(SEXP fp) {
  SIM_YEARS = *INTEGER(get_value(fp, "SIM_YEARS")); // int
  proj_steps = arma_vec(get_value(fp, "proj_steps")); // vec
  basepop = arma_mat(get_value(fp, "basepop")); // mat
  Sx = arma_cube(get_value(fp, "Sx")); // cube
  asfr = arma_mat(get_value(fp, "asfr")); // mat
  srb = arma_mat(get_value(fp, "srb")); // mat         
  netmigr = arma_cube(get_value(fp, "netmigr")); // cube
  birthslag = arma_mat(get_value(fp, "birthslag")); // mat         
  cumsurv = arma_mat(get_value(fp, "cumsurv")); // mat         
  cumnetmigr = arma_mat(get_value(fp, "cumnetmigr")); // mat     
  popadjust = *LOGICAL(get_value(fp, "popadjust")); // bool
  entrantpop = arma_mat(get_value(fp, "entrantpop")); // mat         
  targetpop = arma_cube(get_value(fp, "targetpop")); // cube        
  births = arma_vec(get_value(fp, "births")); // vec
  relinfectART = *REAL(get_value(fp, "relinfectART")); // double
  incrr_sex = arma_vec(get_value(fp, "incrr_sex")); // vec         
  incrr_age = arma_cube(get_value(fp, "incrr_age")); // cube        
  cd4_initdist = arma_cube(get_value(fp, "cd4_initdist")); // cube        
  cd4_prog = arma_cube(get_value(fp, "cd4_prog")); // cube        
  cd4_mort = arma_cube(get_value(fp, "cd4_mort")); // cube        
  art_mort = arma_4D(get_value(fp, "art_mort")); // field<cube>
  artmx_timerr = arma_mat(get_value(fp, "artmx_timerr")); // mat         
  frr_cd4 = arma_cube(get_value(fp, "frr_cd4")); // cube        
  frr_art = arma_4D(get_value(fp, "frr_art")); // field<cube>
  art15plus_num = arma_mat(get_value(fp, "art15plus_num")); // mat         
  art15plus_isperc = arma_mat(get_value(fp, "art15plus_isperc")); // umat
  specpop_percelig = arma_vec(get_value(fp, "specpop_percelig")); // vec         
  artcd4elig_idx = arma_uvec(get_value(fp, "artcd4elig_idx")); // vec         
  pw_artelig = arma_vec(get_value(fp, "pw_artelig")); // vec         
  who34percelig = *REAL(get_value(fp, "who34percelig")); // double      
  art_dropout = arma_vec(get_value(fp, "art_dropout")); // vec         
  median_cd4init = arma_vec(get_value(fp, "median_cd4init")); // vec         
  med_cd4init_input = arma_uvec(get_value(fp, "med_cd4init_input")); // vec
  med_cd4init_cat = arma_uvec(get_value(fp, "med_cd4init_cat")); // vec
  tARTstart = *INTEGER(get_value(fp, "tARTstart")); // int
  art_alloc_method = *INTEGER(get_value(fp, "art_alloc_method")); // int
  art_alloc_mxweight = *REAL(get_value(fp, "art_alloc_mxweight")); // double
  scale_cd4_mort = *INTEGER(get_value(fp, "scale_cd4_mort")); // int
  verttrans_lag = arma_vec(get_value(fp, "verttrans_lag")); // vec
  paedsurv_lag = arma_vec(get_value(fp, "paedsurv_lag")); // vec
  entrantprev = arma_mat(get_value(fp, "entrantprev")); // mat     
  entrantartcov = arma_mat(get_value(fp, "entrantartcov")); // mat 
  paedsurv_cd4dist = arma_cube(get_value(fp, "paedsurv_cd4dist")); // cube   
  paedsurv_artcd4dist = arma_4D(get_value(fp, "paedsurv_artcd4dist")); // field<cube>
  netmig_hivprob = *REAL(get_value(fp, "netmig_hivprob")); // double
  netmighivsurv = *REAL(get_value(fp, "netmighivsurv")); // double
  circ_incid_rr = *REAL(get_value(fp, "circ_incid_rr")); // double
  circ_prop = arma_mat(get_value(fp, "circ_prop")); // mat
  tsEpidemicStart = *REAL(get_value(fp, "tsEpidemicStart")); // double
  eppmod = *INTEGER(get_value(fp, "eppmodInt")); // char to int
  ancsitedata = *LOGICAL(get_value(fp, "ancsitedata")); // bool
  ancrt = *INTEGER(get_value(fp, "ancrtInt")); // char to int
  logitiota = *LOGICAL(get_value(fp, "logitiota")); // bool
  rw_start = *REAL(get_value(fp, "rw_start")); // double
  incidmod = *INTEGER(get_value(fp, "incidmodInt")); // char to int
  rvec = arma_vec(get_value(fp, "rvec")); // vec    
  iota = *REAL(get_value(fp, "iota")); // double      
  ancbias = *REAL(get_value(fp, "ancbias")); // double      
  v_infl = *REAL(get_value(fp, "v_infl")); // double      
  log_frr_adjust = *REAL(get_value(fp, "log_frr_adjust")); // double      
  ancrtcens_vinfl = *REAL(get_value(fp, "ancrtcens_vinfl")); // double      
  ancrtsite_beta = *REAL(get_value(fp, "ancrtsite_beta")); // double      
  mat_m = arma_mat(get_value(fp, "mat_m")); // mat
  mat_f = arma_mat(get_value(fp, "mat_f")); // mat
  db_pr = arma_mat(get_value(fp, "db_pr")); // mat

  // rt parameters
  rt.init_me(fp);
};

CeppFP::CeppFP(SEXP fp) { // all .names --> _names
  SEXP fp_ss = get_value(fp, "ss");
  // Model state space
  proj_start         = *INTEGER(get_value(fp_ss, "proj_start"));
  PROJ_YEARS         = *INTEGER(get_value(fp_ss, "PROJ_YEARS"));
  AGE_START          = *INTEGER(get_value(fp_ss, "AGE_START"));
  hiv_steps_per_year = *INTEGER(get_value(fp_ss, "hiv_steps_per_year"));
  time_epi_start     = *REAL(get_value(fp_ss, "time_epi_start"));
  NG                 = *REAL(get_value(fp_ss, "NG"));
  pDS                = *REAL(get_value(fp_ss, "pDS"));
  m_idx              = *REAL(get_value(fp_ss, "m_idx")) - 1;
  f_idx              = *REAL(get_value(fp_ss, "f_idx")) - 1;
  hivn_idx           = *REAL(get_value(fp_ss, "hivn_idx")) - 1;
  hivp_idx           = *REAL(get_value(fp_ss, "hivp_idx")) - 1;
  pAG                = *REAL(get_value(fp_ss, "pAG"));
  ag_rate            = *REAL(get_value(fp_ss, "ag_rate"));
  p_fert_idx         = arma_uvec(get_value(fp_ss, "p_fert_idx"));
  p_age15to49_idx    = arma_uvec(get_value(fp_ss, "p_age15to49_idx"));
  p_age15plus_idx    = arma_uvec(get_value(fp_ss, "p_age15plus_idx"));
  h_ag_span          = arma_vec(get_value(fp_ss, "h_ag_span"));
  hAG                = *INTEGER(get_value(fp_ss, "hAG"));
  hDS                = *REAL(get_value(fp_ss, "hDS"));
  hTS                = *REAL(get_value(fp_ss, "hTS"));
  ag_idx             = arma_uvec(get_value(fp_ss, "ag_idx"));
  agfirst_idx        = arma_uvec(get_value(fp_ss, "agfirst_idx"));
  aglast_idx         = arma_uvec(get_value(fp_ss, "aglast_idx"));
  h_fert_idx         = arma_uvec(get_value(fp_ss, "h_fert_idx"));
  h_age15to49_idx    = arma_uvec(get_value(fp_ss, "h_age15to49_idx"));
  h_age15plus_idx    = arma_uvec(get_value(fp_ss, "h_age15plus_idx"));
  DT                 = *REAL(get_value(fp_ss, "DT"));
  pDB                = *REAL(get_value(fp_ss, "pDB"));
  // Fix parameters
  p.init_me(fp);
};