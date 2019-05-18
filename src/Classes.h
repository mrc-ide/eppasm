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
#include "utils.h"
#include "Rdefines.h"
class artC;
class hivC;

// Pop class
class popC : public CeppFP {
public: // Pop inits
  popC(SEXP fp, int inMODEL, bool inMIX) : CeppFP(fp) {
    int np = 0;
    MODEL = inMODEL;
    MIX = inMIX;
    year = 1;
    
    // SEXP inits and create pointers
    data_sexp = PROTECT(NEW_NUMERIC(pAG * NG * pDS * PROJ_YEARS)); ++np;
    SEXP data_sexp_dim = PROTECT(NEW_INTEGER(4)); ++np;
    INTEGER(data_sexp_dim)[0] = pAG;
    INTEGER(data_sexp_dim)[1] = NG;
    INTEGER(data_sexp_dim)[2] = pDS;
    INTEGER(data_sexp_dim)[3] = PROJ_YEARS;
    SET_DIM(data_sexp, data_sexp_dim);
    data = fourD2ptr(data_sexp);
    data(0).slice(hivn_idx) = p.basepop;
    // 
    SEXP AgeSexYear_dim = PROTECT(NEW_INTEGER(3)); ++np;
    INTEGER(AgeSexYear_dim)[0] = pAG;
    INTEGER(AgeSexYear_dim)[1] = NG;
    INTEGER(AgeSexYear_dim)[2] = PROJ_YEARS;

    infections_sexp = PROTECT(NEW_NUMERIC(pAG * NG * PROJ_YEARS)); ++np;
    hivdeaths_sexp  = PROTECT(NEW_NUMERIC(pAG * NG * PROJ_YEARS)); ++np;
    natdeaths_sexp  = PROTECT(NEW_NUMERIC(pAG * NG * PROJ_YEARS)); ++np;
    popadjust_sexp  = PROTECT(NEW_NUMERIC(pAG * NG * PROJ_YEARS)); ++np;
    SET_DIM(infections_sexp, AgeSexYear_dim);
    SET_DIM(hivdeaths_sexp, AgeSexYear_dim);
    SET_DIM(natdeaths_sexp, AgeSexYear_dim);
    SET_DIM(popadjust_sexp, AgeSexYear_dim);
    infections = cube2ptr(infections_sexp);
    hivdeaths  = cube2ptr(hivdeaths_sexp);
    natdeaths  = cube2ptr(natdeaths_sexp);
    popadjust  = cube2ptr(popadjust_sexp);

    pregprevlag_sexp = PROTECT(NEW_NUMERIC(PROJ_YEARS)); ++np;
    prev15to49_sexp = PROTECT(NEW_NUMERIC(PROJ_YEARS)); ++np;
    incid15to49_sexp = PROTECT(NEW_NUMERIC(PROJ_YEARS)); ++np;
    entrantprev_sexp = PROTECT(NEW_NUMERIC(PROJ_YEARS)); ++np;
    pregprevlag = vec2ptr(pregprevlag_sexp);
    prev15to49  = vec2ptr(prev15to49_sexp);
    incid15to49 = vec2ptr(incid15to49_sexp);
    entrantprev = vec2ptr(entrantprev_sexp);
    // 
    int n_steps = (PROJ_YEARS-1) * hiv_steps_per_year;
    incrate15to49_ts_sexp = PROTECT(NEW_NUMERIC(n_steps)); ++np;
    prev15to49_ts_sexp    = PROTECT(NEW_NUMERIC(n_steps)); ++np;
    rvec_sexp             = PROTECT(NEW_NUMERIC(n_steps)); ++np;
    incrate15to49_ts = vec2ptr(incrate15to49_ts_sexp);
    prev15to49_ts    = vec2ptr(prev15to49_ts_sexp);
    rvec             = vec2ptr(rvec_sexp);
    if (p.eppmod == 1)
      rvec = p.rvec;

    // for external use
    hivp_entrants_out = zeros(NG, PROJ_YEARS);
    birth_agrp        = zeros(h_fert_idx.n_elem);
    birth_age         = zeros(p_fert_idx.n_elem);
    hiv_sx_prob       = zeros(agfirst_idx.n_elem, NG);
    hiv_mr_prob       = hiv_sx_prob;
    adj_prob          = zeros(pAG, NG);
    birthslag         = p.birthslag;

    data_all = zeros(pAG, NG, pDS); // 1 year only
    if (MODEL==2) { // @debut empty pop, all inactive starts from new entrants
      if (pDB==1) 
        Rf_warning("Debut model state-space seem not exist >> ?update_fp_debut");
      data_db = field<cube>(PROJ_YEARS);
      data_db.for_each( [&] (cube& X) { X.zeros(pDB, NG, pDS); }); // debut ages
    }
    incrate15to49_ts_m = cube(pAG, NG, p.rvec.n_elem);
    prev15to49_ts_m = incrate15to49_ts_m;
    if (MIX) {
      prev15to49_ts_m.zeros(); 
      incrate15to49_ts_m.zeros();
    }
    UNPROTECT(np);
  }
// Pop methods 
  void my_all (int when) ;
  void aging () ;
  void add_entrants () ;
  void sexual_debut () ;
  mat hiv_aging_prob () ;
  vec entrant_art () ;
  void deaths () ;
  void migration () ;
  void update_fertile () ;
  void adjust_pop () ;
  void cal_prev_pregant (hivC& hivpop, artC& artpop); // only on active pop
  void save_prev_n_inc () ;
  mat infect_mix (uword ii);
  mat infect_spec (hivC& hivpop, artC& artpop, uword time_step);
  // epp_disease_model_direct (hivpop, artpop)
  double calc_rtrend_rt (uword ts, double time_step) ;
  void update_rvec (double time_step) ;
  void update_infection (mat infect) ;
  void remove_hiv_death (cube cd4_mx, hivC& hivpop, artC& artpop);
  cube update_preg (cube art_elig, hivC& hivpop, artC& artpop);
  vec artInit (vec art_curr, cube art_elig, int time_step);
  cube artDist (cube art_elig, vec art_need);
  cube scale_cd4_mort (hivC& hivpop, artC& artpop);
  void epp_art_init (hivC& hivpop, artC& artpop, int time_step);
  void finalize (hivC& hivpop, artC& artpop);
public: // Pop fields
  int         MODEL;
  bool        MIX;
  SEXP        data_sexp;
  field<cube> data; // pointer to data_sexp, the same for others
  uword       year;
  vec         birth_age;
  vec         birth_agrp;
  SEXP        prev15to49_sexp;
  vec         prev15to49;  
  SEXP        incid15to49_sexp;
  vec         incid15to49;
  double      prev = 0.0;
  double      incid = 0.0;
  SEXP        entrantprev_sexp;
  vec         entrantprev;
  SEXP        pregprevlag_sexp;
  vec         pregprevlag;
  mat         birthslag;
  SEXP        infections_sexp;
  cube        infections;
  SEXP        hivdeaths_sexp;
  cube        hivdeaths;
  SEXP        natdeaths_sexp;
  cube        natdeaths;
  SEXP        popadjust_sexp;
  cube        popadjust;
  mat         hivp_entrants_out;
  SEXP        incrate15to49_ts_sexp;  
  vec         incrate15to49_ts;  
  SEXP        prev15to49_ts_sexp;
  vec         prev15to49_ts;
  cube        incrate15to49_ts_m; // for storing in mixing model
  cube        prev15to49_ts_m; // for storing in mixing model
  
  field<cube> data_db; // debut only population
  cube        data_all; // all populations in the year requested
  vec         artcov = zeros(2); //numeric(2), // initially no one on treatment
  double      prev_last = 0.0; // = 0 last time step prevalence
  mat         hiv_sx_prob;
  mat         hiv_mr_prob;
  mat         adj_prob;
  SEXP        rvec_sexp;
  vec         rvec;
};

// HIV class
class hivC : public CeppFP {
public: // inits
  hivC(SEXP fp, int inMODEL) : CeppFP(fp) {
    int np = 0;
    MODEL = inMODEL;
    year = 1;

    data_sexp = PROTECT(NEW_NUMERIC(hDS * hAG * NG * PROJ_YEARS)); ++np;
    SEXP data_dim = PROTECT(NEW_INTEGER(4)); ++np;
    INTEGER(data_dim)[0] = hDS;
    INTEGER(data_dim)[1] = hAG;
    INTEGER(data_dim)[2] = NG;
    INTEGER(data_dim)[3] = PROJ_YEARS;
    SET_DIM(data_sexp, data_dim);
    data = fourD2ptr(data_sexp);
    // 
    grad = zeros(hDS, hAG, NG);
    data_all = zeros(pAG, NG, pDS); // 1 year only
    if (MODEL==2) {
      data_db = data;
      data_all = grad_db = grad;
    }
    UNPROTECT(np);
  };
// methods
  void aging(mat ag_prob);
  void add_entrants(vec artYesNo) ;
  void sexual_debut() ;
  void deaths (mat survival_pr) ;
  void migration (mat migration_pr) ;
  void update_infection (mat new_infect) ;
  void grad_progress (cube mortality_rate) ;
  vec eligible_for_art () ;
  cube distribute_artinit (cube artinit, artC& artpop);
  void add_grad_to_pop () ;
  void adjust_pop (mat adj_prob) ;
  // set_data = function(FUN="+", x, DS=T, AG=T, NG=T, YEAR=NULL) {
  // get = function(YEAR=NULL, DS=T, AG=T, NG=T) {
  // sweep_sex = function(FUN="*", x, year) {
public: // fields
  uword       year;
  int         MODEL;
  SEXP        data_sexp;
  field<cube> data;
  field<cube> data_db; // debut only population
  cube        grad;
  cube        grad_db;
  cube        data_all; // all populations in the year requested
};

// ART class
class artC : public CeppFP {
public: // Inits
  artC(SEXP fp, int inMODEL) : CeppFP(fp) {
    int np = 0;
    MODEL = inMODEL;
    year = 1;

    data_sexp = PROTECT(NEW_NUMERIC(hTS * hDS * hAG * NG * PROJ_YEARS)); ++np;
    SEXP data_dim = PROTECT(NEW_INTEGER(5)); ++np;
    INTEGER(data_dim)[0] = hTS;
    INTEGER(data_dim)[1] = hDS;
    INTEGER(data_dim)[2] = hAG;
    INTEGER(data_dim)[3] = NG;
    INTEGER(data_dim)[4] = PROJ_YEARS;
    SET_DIM(data_sexp, data_dim);

    data = fiveD2ptr(data_sexp);
    // data = field<cube>(PROJ_YEARS, NG); // 52 X 2
    // data.for_each( [&](cube& X) { X.zeros(hTS, hDS, hAG); } );
    gradART = field<cube>(NG);
    gradART.for_each( [&](cube& X) { X.zeros(hTS, hDS, hAG); } );
    if (MODEL==2) {
      data_db = data;
      gradART_db = gradART;
    }
    UNPROTECT(np);
  }
// Methods
  void aging (mat ag_prob) ;
  void add_entrants (vec artYesNo) ;
  void sexual_debut () ;
  void deaths (mat survival_pr) ;
  void migration (mat migration_pr) ;
  void grad_progress () ;
  void art_dropout (hivC& hivpop) ;
  vec current_on_art () ;
  void grad_init (cube artinit) ;
  void grad_db_init (cube artinit_db) ;
  void adjust_pop (mat adj_prob) ;
  // set_data (FUN="+", x, TS=T, DS=T, AG=T, NG=T, YEAR=NULL) {
  // get (YEAR=NULL, TS=T, DS=T, AG=T, NG=T) {
  // sweep_sex (FUN="*", x, year) {
public: // fields
  uword       year;
  int         MODEL;
  SEXP        data_sexp;
  field<cube> data;
  field<cube> data_db; // debut only population
  field<cube> gradART;
  field<cube> gradART_db;
};