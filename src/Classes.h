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
  popC(SEXP fp, int inMODEL, bool inMIX) : CeppFP(fp),
// boost array class inits
    data_sexp(PROTECT(NEW_NUMERIC(pAG * NG * pDS * PROJ_YEARS))),
    data(REAL(data_sexp), extents[PROJ_YEARS][pDS][NG][pAG]),
    
    birth_age(extents[pAG_FERT]),
    birth_agrp(extents[hAG_FERT]),
    
    prev15to49_sexp(PROTECT(NEW_NUMERIC(PROJ_YEARS))),
    incid15to49_sexp(PROTECT(NEW_NUMERIC(PROJ_YEARS))),
    prev15to49(REAL(prev15to49_sexp), extents[PROJ_YEARS]),
    incid15to49(REAL(incid15to49_sexp), extents[PROJ_YEARS]),
    
    pregprevlag_sexp(PROTECT(NEW_NUMERIC(PROJ_YEARS))),
    entrantprev_sexp(PROTECT(NEW_NUMERIC(PROJ_YEARS))),
    entrantprev(REAL(entrantprev_sexp), extents[PROJ_YEARS]),
    pregprevlag(REAL(pregprevlag_sexp), extents[PROJ_YEARS]),
    
    inci15to49_ts_sexp(PROTECT(NEW_NUMERIC(n_steps))),
    prev15to49_ts_sexp(PROTECT(NEW_NUMERIC(n_steps))),
    incrate15to49_ts(REAL(inci15to49_ts_sexp), extents[n_steps]),
    prev15to49_ts(REAL(prev15to49_ts_sexp), extents[n_steps]),
    
    rvec_sexp(PROTECT(NEW_NUMERIC(n_steps))),
    rvec(REAL(rvec_sexp), extents[n_steps]),
    
    birthslag(p.birthslag),
    hivp_entrants_out(extents[PROJ_YEARS][NG]),
    hiv_sx_prob(extents[NG][hAG]),
    hiv_mr_prob(extents[NG][hAG]),
    adj_prob(extents[NG][hAG]),
    
    infections_sexp ( PROTECT(NEW_NUMERIC(pAG * NG * PROJ_YEARS)) ),
    hivdeaths_sexp  ( PROTECT(NEW_NUMERIC(pAG * NG * PROJ_YEARS)) ),
    natdeaths_sexp  ( PROTECT(NEW_NUMERIC(pAG * NG * PROJ_YEARS)) ),
    popadjust_sexp  ( PROTECT(NEW_NUMERIC(pAG * NG * PROJ_YEARS)) ),
    infections (REAL(infections_sexp), extents[PROJ_YEARS][NG][pAG]),
    hivdeaths  (REAL(hivdeaths_sexp ), extents[PROJ_YEARS][NG][pAG]),
    natdeaths  (REAL(natdeaths_sexp ), extents[PROJ_YEARS][NG][pAG]),
    popadjust  (REAL(popadjust_sexp ), extents[PROJ_YEARS][NG][pAG]),
    
    data_all(extents[pDS][NG][pAG]), // 1 year only
    data_db(extents[PROJ_YEARS][pDS][NG][pDB])
  {
// Non class init
    UNPROTECT(12);
    MODEL = inMODEL;
    MIX = inMIX;
    year = 1;

    memset(REAL(data_sexp         ), 0, pAG * NG * pDS * PROJ_YEARS * sizeof(double));
    memset(REAL(prev15to49_sexp   ), 0, PROJ_YEARS                  * sizeof(double));
    memset(REAL(incid15to49_sexp  ), 0, PROJ_YEARS                  * sizeof(double));
    memset(REAL(pregprevlag_sexp  ), 0, PROJ_YEARS                  * sizeof(double));
    memset(REAL(entrantprev_sexp  ), 0, PROJ_YEARS                  * sizeof(double));
    memset(REAL(inci15to49_ts_sexp), 0, n_steps                     * sizeof(double));
    memset(REAL(prev15to49_ts_sexp), 0, n_steps                     * sizeof(double));
    memset(REAL(rvec_sexp         ), 0, n_steps                     * sizeof(double));
    memset(REAL(infections_sexp   ), 0, pAG * NG * PROJ_YEARS       * sizeof(double));
    memset(REAL(hivdeaths_sexp    ), 0, pAG * NG * PROJ_YEARS       * sizeof(double));
    memset(REAL(natdeaths_sexp    ), 0, pAG * NG * PROJ_YEARS       * sizeof(double));
    memset(REAL(popadjust_sexp    ), 0, pAG * NG * PROJ_YEARS       * sizeof(double));

    for (int sex = 0; sex < NG; sex++)
      for (int age = 0; age < pAG; age++)
        data[0][hivn_idx][sex][age] = p.basepop[sex][age];

    if (p.eppmod == 0)
      rvec = p.rvec;

    if ( MODEL==2 && pDB==1 )
      Rf_warning("Debut model state-space not exist, see update_fp_debut()");
  }
// Pop methods 
  void my_all (int when) ;
  void aging () ;
  void add_entrants () ;
  void sexual_debut () ;
  boost2D hiv_aging_prob () ;
  boost1D entrant_art () ;
  void deaths () ;
  void migration () ;
  void update_fertile () ;
  void adjust_pop () ;
  void cal_prev_pregant (hivC& hivpop, artC& artpop); // only on active pop
  void save_prev_n_inc () ;
  boost2D infect_mix (int ii);
  boost2D infect_spec (hivC& hivpop, artC& artpop, int time_step);
  void epp_disease_model_direct (hivC& hivpop, artC& artpop) ;
  double calc_rtrend_rt (int ts, double time_step) ;
  void update_rvec (double time_step) ;
  void update_infection (boost2D infect) ;
  void remove_hiv_death (boost3D cd4_mx, hivC& hivpop, artC& artpop);
  boost3D update_preg (boost3D art_elig, hivC& hivpop, artC& artpop);
  boost1D artInit (boost1D art_curr, boost3D art_elig, int time_step);
  boost3D artDist (boost3D art_elig, boost1D art_need);
  boost3D scale_cd4_mort (hivC& hivpop, artC& artpop);
  void epp_art_init (hivC& hivpop, artC& artpop, int time_step);
  void finalize (hivC& hivpop, artC& artpop);
public: // Pop fields
  int         MODEL;
  bool        MIX;
  int         year;
  SEXP        data_sexp;
  boost4D_ptr data; // pointer to pop_data_sexp, the same for others
  boost1D     birth_age;
  boost1D     birth_agrp;
  SEXP        prev15to49_sexp;
  SEXP        incid15to49_sexp;
  boost1D_ptr prev15to49;  
  boost1D_ptr incid15to49;
  double      prev = 0.0;
  double      incid = 0.0;
  SEXP        pregprevlag_sexp;
  SEXP        entrantprev_sexp;
  boost1D_ptr entrantprev;
  boost1D_ptr pregprevlag;
  SEXP        inci15to49_ts_sexp;
  SEXP        prev15to49_ts_sexp;
  boost1D_ptr incrate15to49_ts;  
  boost1D_ptr prev15to49_ts;
  SEXP        rvec_sexp;
  boost1D_ptr rvec;
  boost2D     birthslag;
  boost2D     hivp_entrants_out;
  boost2D     hiv_sx_prob;
  boost2D     hiv_mr_prob;
  boost2D     adj_prob;
  SEXP        infections_sexp;
  SEXP        hivdeaths_sexp;
  SEXP        natdeaths_sexp;
  SEXP        popadjust_sexp;
  boost3D_ptr infections;
  boost3D_ptr hivdeaths;
  boost3D_ptr natdeaths;
  boost3D_ptr popadjust;
  // boost3D     incrate15to49_ts_m; // for storing in mixing model
  // boost3D     prev15to49_ts_m; // for storing in mixing model
  boost3D     data_all; // all populations in the year requested
  boost4D     data_db; // debut only population
  double      artcov[2] = {0, 0}; //numeric(2), // initially no one on treatment
  double      prev_last = 0.0; // = 0 last time step prevalence
};

// HIV class
class hivC : public CeppFP {
public: // inits
  hivC(SEXP fp, int inMODEL) : CeppFP(fp),
// Boost array class init
  data_sexp(PROTECT(NEW_NUMERIC(hDS * hAG * NG * PROJ_YEARS))),
  data(REAL(data_sexp), extents[PROJ_YEARS][NG][hAG][hDS]),
  data_db(extents[PROJ_YEARS][NG][hAG][hDS]), // later return this as well
  grad(extents[NG][hAG][hDS]),
  grad_db(extents[NG][hAG][hDS]),
  data_all(extents[NG][hAG][hDS])
  {
    UNPROTECT(1);
    MODEL = inMODEL;
    memset(REAL(data_sexp), 0.0, hDS * hAG * NG * PROJ_YEARS * sizeof(double));
  };
// methods
  void aging(boost2D ag_prob);
  void add_entrants(boost1D artYesNo) ;
  void sexual_debut() ;
  void deaths (boost2D survival_pr) ;
  void migration (boost2D migration_pr) ;
  void update_infection (boost2D new_infect) ;
  void grad_progress (boost3D mortality_rate) ;
  boost1D eligible_for_art () ;
  boost3D distribute_artinit (boost3D artinit, artC& artpop);
  void add_grad_to_pop () ;
  void adjust_pop (boost2D adj_prob) ;
public: // fields
  int         year = 1;
  int         MODEL;
  SEXP        data_sexp;
  boost4D_ptr data;
  boost4D     data_db; // debut only population
  boost3D     grad;
  boost3D     grad_db;
  boost3D     data_all; // all populations in the year requested
};

// ART class
class artC : public CeppFP {
public: // Inits
  artC(SEXP fp, int inMODEL) : CeppFP(fp), 
    data_sexp(PROTECT(NEW_NUMERIC(hTS * hDS * hAG * NG * PROJ_YEARS))),
    data(REAL(data_sexp), extents[PROJ_YEARS][NG][hAG][hDS][hTS]),
    data_db(extents[PROJ_YEARS][NG][hAG][hDS][hTS]),
    gradART(extents[NG][hAG][hDS][hTS]),
    gradART_db(extents[NG][hAG][hDS][hTS])
  {
    UNPROTECT(1);
    MODEL = inMODEL;
    memset(REAL(data_sexp), 0, hTS * hDS * hAG * NG * PROJ_YEARS * sizeof(double));
  }
// Methods
  void aging (boost2D ag_prob) ;
  void add_entrants (boost1D artYesNo) ;
  void sexual_debut () ;
  void deaths (boost2D survival_pr) ;
  void migration (boost2D migration_pr) ;
  void grad_progress () ;
  void art_dropout (hivC& hivpop) ;
  boost1D current_on_art () ;
  void grad_init (boost3D artinit) ;
  void grad_db_init (boost3D artinit_db) ;
  void adjust_pop (boost2D adj_prob) ;
public: // fields
  int         year = 1;
  int         MODEL;
  SEXP        data_sexp;
  boost5D_ptr data;
  boost5D     data_db; // debut only population
  boost4D     gradART;
  boost4D     gradART_db;
};