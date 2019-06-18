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
class oSEXP : public CeppFP { // outputs for R
public:
  oSEXP(SEXP fp) : CeppFP(fp) {
    pop           = PROTECT(NEW_NUMERIC(pAG * NG * pDS * PROJ_YEARS)); ++np;
    prev15to49    = PROTECT(NEW_NUMERIC(PROJ_YEARS)); ++np;
    incid15to49   = PROTECT(NEW_NUMERIC(PROJ_YEARS)); ++np;
    entrantprev   = PROTECT(NEW_NUMERIC(PROJ_YEARS)); ++np;
    pregprevlag   = PROTECT(NEW_NUMERIC(PROJ_YEARS)); ++np;
    inci15to49_ts = PROTECT(NEW_NUMERIC(n_steps)); ++np;
    prev15to49_ts = PROTECT(NEW_NUMERIC(n_steps)); ++np;
    rvec          = PROTECT(NEW_NUMERIC(n_steps)); ++np;
    infections    = PROTECT(NEW_NUMERIC(pAG * NG * PROJ_YEARS)); ++np;
    hivdeaths     = PROTECT(NEW_NUMERIC(pAG * NG * PROJ_YEARS)); ++np;
    natdeaths     = PROTECT(NEW_NUMERIC(pAG * NG * PROJ_YEARS)); ++np;
    popadjust     = PROTECT(NEW_NUMERIC(pAG * NG * PROJ_YEARS)); ++np;
    data_db       = PROTECT(NEW_NUMERIC(pDB * NG * pDS * PROJ_YEARS)); ++np;

    memset(REAL(pop          ), 0, pAG * NG * pDS * PROJ_YEARS * sizeof(double));
    memset(REAL(prev15to49   ), 0, PROJ_YEARS                  * sizeof(double));
    memset(REAL(incid15to49  ), 0, PROJ_YEARS                  * sizeof(double));
    memset(REAL(pregprevlag  ), 0, PROJ_YEARS                  * sizeof(double));
    memset(REAL(entrantprev  ), 0, PROJ_YEARS                  * sizeof(double));
    memset(REAL(inci15to49_ts), 0, n_steps                     * sizeof(double));
    memset(REAL(prev15to49_ts), 0, n_steps                     * sizeof(double));
    memset(REAL(rvec         ), 0, n_steps                     * sizeof(double));
    memset(REAL(infections   ), 0, pAG * NG * PROJ_YEARS       * sizeof(double));
    memset(REAL(hivdeaths    ), 0, pAG * NG * PROJ_YEARS       * sizeof(double));
    memset(REAL(natdeaths    ), 0, pAG * NG * PROJ_YEARS       * sizeof(double));
    memset(REAL(popadjust    ), 0, pAG * NG * PROJ_YEARS       * sizeof(double));
    memset(REAL(data_db      ), 0, pDB * NG * pDS * PROJ_YEARS * sizeof(double));
    
    hivpop = PROTECT(NEW_NUMERIC(hDS * hAG * NG * PROJ_YEARS)); ++np;
    memset(REAL(hivpop), 0, hDS * hAG * NG * PROJ_YEARS * sizeof(double));
    
    artpop = PROTECT(NEW_NUMERIC(hTS * hDS * hAG * NG * PROJ_YEARS)); ++np;
    memset(REAL(artpop), 0, hTS * hDS * hAG * NG * PROJ_YEARS * sizeof(double));
  }
  // output fields
  void finalize();
public:
  int np = 0;
  SEXP pop;
  SEXP prev15to49;
  SEXP incid15to49;
  SEXP entrantprev;
  SEXP pregprevlag;
  SEXP inci15to49_ts;
  SEXP prev15to49_ts;
  SEXP rvec;
  SEXP infections;
  SEXP hivdeaths;
  SEXP natdeaths;
  SEXP popadjust;
  SEXP data_db;
  SEXP hivpop;
  SEXP artpop;
};

// Pop class
class popC : public CeppFP {
public: // Pop inits
  popC(oSEXP& O, SEXP fp, bool inMIX) : CeppFP(fp),
// boost array class inits
    data(REAL(O.pop), extents[PROJ_YEARS][pDS][NG][pAG]),
    birth_age(pAG_FERT),
    birth_agrp(hAG_FERT),
    prev15to49(REAL(O.prev15to49), extents[PROJ_YEARS]),
    incid15to49(REAL(O.incid15to49), extents[PROJ_YEARS]),
    entrantprev(REAL(O.entrantprev), extents[PROJ_YEARS]),
    pregprevlag(REAL(O.pregprevlag), extents[PROJ_YEARS]),
    incrate15to49_ts(REAL(O.inci15to49_ts), extents[n_steps]),
    prev15to49_ts(REAL(O.prev15to49_ts), extents[n_steps]),
    rvec(REAL(O.rvec), extents[n_steps]),
    birthslag(p.birthslag),
    hivp_entrants_out(extents[PROJ_YEARS][NG]),
    hiv_sx_prob(extents[NG][hAG]),
    hiv_mr_prob(extents[NG][hAG]),
    adj_prob(extents[NG][hAG]),
    infections(REAL(O.infections), extents[PROJ_YEARS][NG][pAG]),
    hivdeaths(REAL(O.hivdeaths ), extents[PROJ_YEARS][NG][pAG]),
    natdeaths(REAL(O.natdeaths ), extents[PROJ_YEARS][NG][pAG]),
    popadjust(REAL(O.popadjust ), extents[PROJ_YEARS][NG][pAG]),
    data_all(extents[pDS][NG][pAG]), // 1 year only
    data_db(REAL(O.data_db), extents[PROJ_YEARS][pDS][NG][pDB]),
    data_active(extents[pDS][NG][pAG]) // 1 year only
  {
// Non class init
    MIX = inMIX;

    for (int sex = 0; sex < NG; sex++)
      for (int age = 0; age < pAG; age++)
        data[0][hivn_idx][sex][age] = p.basepop[sex][age];

    if (p.eppmod == 0)
      for (int i = 0; i < n_steps; ++i)
        rvec[i] = p.rvec[i];

    if (MODEL == 2 && pDB == 1)
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
  bool        MIX;
  int         year = 1;
  boost4D_ptr data; // pointer to pop_data_sexp, the same for others
  dvec        birth_age;
  dvec        birth_agrp;
  boost1D_ptr prev15to49;  
  boost1D_ptr incid15to49;
  double      prev = 0.0;
  double      incid = 0.0;
  boost1D_ptr entrantprev;
  boost1D_ptr pregprevlag;
  boost1D_ptr incrate15to49_ts;  
  boost1D_ptr prev15to49_ts;
  boost1D_ptr rvec;
  boost2D     birthslag;
  boost2D     hivp_entrants_out;
  boost2D     hiv_sx_prob;
  boost2D     hiv_mr_prob;
  boost2D     adj_prob;
  boost3D_ptr infections;
  boost3D_ptr hivdeaths;
  boost3D_ptr natdeaths;
  boost3D_ptr popadjust;
  // boost3D     incrate15to49_ts_m; // for storing in mixing model
  // boost3D     prev15to49_ts_m; // for storing in mixing model
  boost3D     data_all; // all populations in the year requested
  boost4D_ptr data_db; // debut only population
  boost3D     data_active;
  double      artcov[2] = {0.0, 0.0}; // initially no one on treatment
  double      prev_last = 0.0; // = 0 last time step prevalence
};

// HIV class
class hivC : public CeppFP {
public: // inits
  hivC(oSEXP& O, SEXP fp) : CeppFP(fp),
// Boost array class init
  data(REAL(O.hivpop), extents[PROJ_YEARS][NG][hAG][hDS]),
  data_db(extents[PROJ_YEARS][NG][hAG][hDS]), // later return this as well
  grad(extents[NG][hAG][hDS]),
  grad_db(extents[NG][hAG][hDS]),
  data_all(extents[NG][hAG][hDS]),
  _death(extents[NG][hAG][hDS]),
  _death_db(extents[NG][hAG][hDS]) {};
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
  artC(oSEXP& O, SEXP fp) : CeppFP(fp),
    data(REAL(O.artpop), extents[PROJ_YEARS][NG][hAG][hDS][hTS]),
    data_db(extents[PROJ_YEARS][NG][hAG][hDS][hTS]),
    gradART(extents[NG][hAG][hDS][hTS]),
    gradART_db(extents[NG][hAG][hDS][hTS]),
    _death(extents[NG][hAG][hDS][hTS]),
    _death_db(extents[NG][hAG][hDS][hTS]) {}
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
  boost5D_ptr data;
  boost5D     data_db; // debut only population
  boost4D     gradART;
  boost4D     gradART_db;
  boost4D     _death;
  boost4D     _death_db;
};
