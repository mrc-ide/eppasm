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
#include "utils.hpp"
#include "Rdefines.h"
class artC;
class hivC;

// Pop class
class popC : public CeppFP {
public: // Pop inits
  popC(SEXP fp, int inMODEL, bool inMIX) : CeppFP(fp),
// boost array class inits
    pop_sexp(PROTECT(NEW_NUMERIC(pAG * NG * pDS * PROJ_YEARS))),
    data(REAL(pop_sexp), extents[PROJ_YEARS][pDS][NG][pAG]),
    
    birth_age(pAG_FERT),
    birth_agrp(hAG_FERT),
    
    prev15to49_sexp(PROTECT(NEW_NUMERIC(PROJ_YEARS))),
    prev15to49(REAL(prev15to49_sexp), extents[PROJ_YEARS]),
    incid15to49_sexp(PROTECT(NEW_NUMERIC(PROJ_YEARS))),
    incid15to49(REAL(incid15to49_sexp), extents[PROJ_YEARS]),
    
    entrantprev_sexp(PROTECT(NEW_NUMERIC(PROJ_YEARS))),
    entrantprev(REAL(entrantprev_sexp), extents[PROJ_YEARS]),
    pregprevlag_sexp(PROTECT(NEW_NUMERIC(PROJ_YEARS))),
    pregprevlag(REAL(pregprevlag_sexp), extents[PROJ_YEARS]),
    
    inci15to49_ts_sexp(PROTECT(NEW_NUMERIC(n_steps))),
    incrate15to49_ts(REAL(inci15to49_ts_sexp), extents[n_steps]),
    prev15to49_ts_sexp(PROTECT(NEW_NUMERIC(n_steps))),
    prev15to49_ts(REAL(prev15to49_ts_sexp), extents[n_steps]),
    
    rvec_sexp(PROTECT(NEW_NUMERIC(n_steps))),
    rvec(REAL(rvec_sexp), extents[n_steps]),
    
    birthslag(p.birthslag),
    hivp_entrants_out(extents[PROJ_YEARS][NG]),
    hiv_sx_prob(extents[NG][hAG]),
    hiv_mr_prob(extents[NG][hAG]),
    adj_prob(extents[NG][hAG]),
    
    infections_sexp(PROTECT(NEW_NUMERIC(pAG * NG * PROJ_YEARS))),
    infections(REAL(infections_sexp), extents[PROJ_YEARS][NG][pAG]),
    hivdeaths_sexp(PROTECT(NEW_NUMERIC(pAG * NG * PROJ_YEARS))),
    hivdeaths(REAL(hivdeaths_sexp ), extents[PROJ_YEARS][NG][pAG]),
    natdeaths_sexp(PROTECT(NEW_NUMERIC(pAG * NG * PROJ_YEARS))),
    natdeaths(REAL(natdeaths_sexp ), extents[PROJ_YEARS][NG][pAG]),
    popadjust_sexp(PROTECT(NEW_NUMERIC(pAG * NG * PROJ_YEARS))),
    popadjust(REAL(popadjust_sexp ), extents[PROJ_YEARS][NG][pAG]),
    
    data_all(extents[pDS][NG][pAG]), // 1 year only
    data_db_sexp(PROTECT(NEW_NUMERIC(pDB * NG * pDS * PROJ_YEARS))),
    data_db(REAL(data_db_sexp), extents[PROJ_YEARS][pDS][NG][pDB]),
    data_active(extents[pDS][NG][pAG]) // 1 year only
  {
// Non class init
    MODEL = inMODEL;
    MIX = inMIX;

    memset(REAL(pop_sexp          ), 0, pAG * NG * pDS * PROJ_YEARS * sizeof(double));
    memset(REAL(data_db_sexp      ), 0, pDB * NG * pDS * PROJ_YEARS * sizeof(double));
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
      for (int i = 0; i < n_steps; ++i)
        rvec[i] = p.rvec[i];

    if ( MODEL==2 && pDB==1 )
      Rf_warning("Debut model state-space not exist, see update_fp_debut()");
  }
// Pop methods 
  void my_all (int when) ;
  void update_active_pop_to (int year) ;
  boost3D get_active_pop_in (int year) ;
  void aging () ;
  void add_entrants () ;
  void sexual_debut () ;
  boost2D hiv_aging_prob () ;
  dvec entrant_art () ;
  void deaths () ;
  void migration () ;
  void update_fertile () ;
  void adjust_pop () ;
  void cal_prev_pregant (const hivC& hivpop, const artC& artpop); // only on active pop
  void save_prev_n_inc () ;
  boost2D infect_mix (int ii);
  boost2D infect_spec (const hivC& hivpop, const artC& artpop, int time_step);
  void epp_disease_model_direct (hivC& hivpop, artC& artpop) ;
  double calc_rtrend_rt (int ts, double time_step) ;
  void update_rvec (double time_step) ;
  void update_infection (const boost2D& infect) ;
  void remove_hiv_death (const boost3D& cd4_mx,
                         const hivC& hivpop, const artC& artpop);
  void update_preg (boost3D& art_elig,
                    const hivC& hivpop, const artC& artpop);
  dvec art_initiate (const dvec& art_curr, const boost3D& art_elig,
                   int time_step);
  boost3D art_distribute (const boost3D& art_elig, const dvec& art_need);
  boost3D scale_cd4_mort (hivC& hivpop, artC& artpop);
  void epp_art_init (hivC& hivpop, artC& artpop, int time_step);
  void finalize (hivC& hivpop, artC& artpop);
public: // Pop fields
  int         MODEL;
  bool        MIX;
  int         year = 1;
  SEXP        pop_sexp;
  boost4D_ptr data; // pointer to pop_data_sexp, the same for others
  dvec        birth_age;
  dvec        birth_agrp;
  SEXP        prev15to49_sexp;
  boost1D_ptr prev15to49;  
  SEXP        incid15to49_sexp;
  boost1D_ptr incid15to49;
  double      prev = 0.0;
  double      incid = 0.0;
  SEXP        entrantprev_sexp;
  boost1D_ptr entrantprev;
  SEXP        pregprevlag_sexp;
  boost1D_ptr pregprevlag;
  SEXP        inci15to49_ts_sexp;
  boost1D_ptr incrate15to49_ts;  
  SEXP        prev15to49_ts_sexp;
  boost1D_ptr prev15to49_ts;
  SEXP        rvec_sexp;
  boost1D_ptr rvec;
  boost2D     birthslag;
  boost2D     hivp_entrants_out;
  boost2D     hiv_sx_prob;
  boost2D     hiv_mr_prob;
  boost2D     adj_prob;
  SEXP        infections_sexp;
  boost3D_ptr infections;
  SEXP        hivdeaths_sexp;
  boost3D_ptr hivdeaths;
  SEXP        natdeaths_sexp;
  boost3D_ptr natdeaths;
  SEXP        popadjust_sexp;
  boost3D_ptr popadjust;
  // boost3D     incrate15to49_ts_m; // for storing in mixing model
  // boost3D     prev15to49_ts_m; // for storing in mixing model
  boost3D     data_all; // all populations in the year requested
  SEXP        data_db_sexp; // debut only population
  boost4D_ptr data_db; // debut only population
  boost3D     data_active;
  double      artcov[2] = {0.0, 0.0}; // initially no one on treatment
  double      prev_last = 0.0; // = 0 last time step prevalence
};

// HIV class
class hivC : public CeppFP {
public: // inits
  hivC(SEXP fp, int inMODEL) : CeppFP(fp),
// Boost array class init
  hiv_sexp(PROTECT(NEW_NUMERIC(hDS * hAG * NG * PROJ_YEARS))),
  data(REAL(hiv_sexp), extents[PROJ_YEARS][NG][hAG][hDS]),
  data_db(extents[PROJ_YEARS][NG][hAG][hDS]), // later return this as well
  grad(extents[NG][hAG][hDS]),
  grad_db(extents[NG][hAG][hDS]),
  data_all(extents[NG][hAG][hDS]),
  _death(extents[NG][hAG][hDS]),
  _death_db(extents[NG][hAG][hDS])
  {
    MODEL = inMODEL;
    memset(REAL(hiv_sexp), 0, hDS * hAG * NG * PROJ_YEARS * sizeof(double));
    zeroing(data_db); zeroing(grad); zeroing(grad_db);
  };
// methods
  void aging(const boost2D& ag_prob);
  void add_entrants(const dvec& artYesNo) ;
  void sexual_debut() ;
  void deaths (const boost2D& survival_pr) ;
  void migration (const boost2D& migration_pr) ;
  void update_infection (const boost2D& new_infect) ;
  void grad_progress (const boost3D& mortality_rate) ;
  dvec eligible_for_art () ;
  void distribute_artinit (boost3D& artinit, artC& artpop);
  void add_grad_to_pop () ;
  void adjust_pop (const boost2D& adj_prob) ;
public: // fields
  int         year = 1;
  int         MODEL;
  SEXP        hiv_sexp;
  boost4D_ptr data;
  boost4D     data_db; // debut only population
  boost3D     grad;
  boost3D     grad_db;
  boost3D     data_all; // all populations in the year requested
  boost3D     _death; // death in this year
  boost3D     _death_db; // death in this year
};

// ART class
class artC : public CeppFP {
public: // Inits
  artC(SEXP fp, int inMODEL) : CeppFP(fp), 
    art_sexp(PROTECT(NEW_NUMERIC(hTS * hDS * hAG * NG * PROJ_YEARS))),
    data(REAL(art_sexp), extents[PROJ_YEARS][NG][hAG][hDS][hTS]),
    data_db(extents[PROJ_YEARS][NG][hAG][hDS][hTS]),
    gradART(extents[NG][hAG][hDS][hTS]),
    gradART_db(extents[NG][hAG][hDS][hTS]),
    _death(extents[NG][hAG][hDS][hTS]),
    _death_db(extents[NG][hAG][hDS][hTS])
  {
    MODEL = inMODEL;
    memset(REAL(art_sexp), 0, hTS * hDS * hAG * NG * PROJ_YEARS * sizeof(double));
    zeroing(data_db); zeroing(gradART); zeroing(gradART_db);
  }
// Methods
  void aging (const boost2D& ag_prob) ;
  void add_entrants (const dvec& artYesNo) ;
  void sexual_debut () ;
  void deaths (const boost2D& survival_pr) ;
  void migration (const boost2D& migration_pr) ;
  void grad_progress () ;
  void art_dropout (hivC& hivpop) ;
  dvec current_on_art () ;
  void grad_init (const boost3D& artinit) ;
  void grad_db_init (const boost3D& artinit_db) ;
  void adjust_pop (const boost2D& adj_prob) ;
  void count_death () ;
public: // fields
  int         year = 1;
  int         MODEL;
  SEXP        art_sexp;
  boost5D_ptr data;
  boost5D     data_db; // debut only population
  boost4D     gradART;
  boost4D     gradART_db;
  boost4D     _death;
  boost4D     _death_db;
};
