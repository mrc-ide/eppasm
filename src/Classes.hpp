#include "fpClass.hpp"
#include "utils.hpp"
#include "Rdefines.h"
class artC;
class hivC;

class outputSEXP { // outputs for R
public:
  outputSEXP(StateSpace& s) {
    pop = PROTECT(NEW_NUMERIC(s.pAG * s.NG * s.pDS * s.PROJ_YEARS)); ++np;
    memset(REAL(pop), 0, s.pAG * s.NG * s.pDS * s.PROJ_YEARS * sizeof(double));

    prev15to49    = PROTECT(NEW_NUMERIC(s.PROJ_YEARS)); ++np;
    incid15to49   = PROTECT(NEW_NUMERIC(s.PROJ_YEARS)); ++np;
    entrantprev   = PROTECT(NEW_NUMERIC(s.PROJ_YEARS)); ++np;
    pregprevlag   = PROTECT(NEW_NUMERIC(s.PROJ_YEARS)); ++np;
    inci15to49_ts = PROTECT(NEW_NUMERIC(s.n_steps)); ++np;
    prev15to49_ts = PROTECT(NEW_NUMERIC(s.n_steps)); ++np;
    rvec          = PROTECT(NEW_NUMERIC(s.n_steps)); ++np;
    infections    = PROTECT(NEW_NUMERIC(s.pAG * s.NG * s.PROJ_YEARS)); ++np;
    hivdeaths     = PROTECT(NEW_NUMERIC(s.pAG * s.NG * s.PROJ_YEARS)); ++np;
    natdeaths     = PROTECT(NEW_NUMERIC(s.pAG * s.NG * s.PROJ_YEARS)); ++np;
    popadjust     = PROTECT(NEW_NUMERIC(s.pAG * s.NG * s.PROJ_YEARS)); ++np;

    memset(REAL(prev15to49   ), 0, s.PROJ_YEARS * sizeof(double));
    memset(REAL(incid15to49  ), 0, s.PROJ_YEARS * sizeof(double));
    memset(REAL(pregprevlag  ), 0, s.PROJ_YEARS * sizeof(double));
    memset(REAL(entrantprev  ), 0, s.PROJ_YEARS * sizeof(double));
    
    memset(REAL(inci15to49_ts), 0, s.n_steps    * sizeof(double));
    memset(REAL(prev15to49_ts), 0, s.n_steps    * sizeof(double));
    memset(REAL(rvec         ), 0, s.n_steps    * sizeof(double));

    memset(REAL(infections), 0, s.pAG * s.NG * s.PROJ_YEARS * sizeof(double));
    memset(REAL(hivdeaths ), 0, s.pAG * s.NG * s.PROJ_YEARS * sizeof(double));
    memset(REAL(natdeaths ), 0, s.pAG * s.NG * s.PROJ_YEARS * sizeof(double));
    memset(REAL(popadjust ), 0, s.pAG * s.NG * s.PROJ_YEARS * sizeof(double));
    
    data_db = PROTECT(NEW_NUMERIC(s.pDB * s.NG * s.pDS * s.PROJ_YEARS)); ++np;
    memset(REAL(data_db), 0, s.pDB * s.NG * s.pDS * s.PROJ_YEARS * sizeof(double));
    
    hivpop = PROTECT(NEW_NUMERIC(s.hDS * s.hAG * s.NG * s.PROJ_YEARS)); ++np;
    memset(REAL(hivpop), 0, s.hDS * s.hAG * s.NG * s.PROJ_YEARS * sizeof(double));
    
    artpop = PROTECT(NEW_NUMERIC(s.hTS * s.hDS * s.hAG * s.NG * s.PROJ_YEARS)); ++np;
    memset(REAL(artpop), 0, s.hTS * s.hDS * s.hAG * s.NG * s.PROJ_YEARS * sizeof(double));
  }
  // output fields
  void finalize(const StateSpace& s);
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
class popC {
public: // Pop inits
  popC(outputSEXP& O, const StateSpace& s) :
// boost array class inits
    data(REAL(O.pop), extents[s.PROJ_YEARS][s.pDS][s.NG][s.pAG]),
    birth_age          (s.pAG_FERT),
    birth_agrp         (s.hAG_FERT),
    prev15to49         (REAL(O.prev15to49), extents[s.PROJ_YEARS]),
    incid15to49        (REAL(O.incid15to49), extents[s.PROJ_YEARS]),
    entrantprev        (REAL(O.entrantprev), extents[s.PROJ_YEARS]),
    pregprevlag        (REAL(O.pregprevlag), extents[s.PROJ_YEARS]),
    incrate15to49_ts   (REAL(O.inci15to49_ts), extents[s.n_steps]),
    prev15to49_ts      (REAL(O.prev15to49_ts), extents[s.n_steps]),
    rvec               (REAL(O.rvec), extents[s.n_steps]),
    birthslag          (extents[s.PROJ_YEARS][s.NG]),
    hivp_entrants_out  (extents[s.PROJ_YEARS][s.NG]),
    hiv_sx_prob        (extents[s.NG][s.hAG]),
    hiv_mr_prob        (extents[s.NG][s.hAG]),
    adj_prob           (extents[s.NG][s.hAG]),
    infections         (REAL(O.infections), extents[s.PROJ_YEARS][s.NG][s.pAG]),
    hivdeaths          (REAL(O.hivdeaths ), extents[s.PROJ_YEARS][s.NG][s.pAG]),
    natdeaths          (REAL(O.natdeaths ), extents[s.PROJ_YEARS][s.NG][s.pAG]),
    popadjust          (REAL(O.popadjust ), extents[s.PROJ_YEARS][s.NG][s.pAG]),
    data_all           (extents[s.pDS][s.NG][s.pAG]), // 1 year only
    data_db            (REAL(O.data_db), extents[s.PROJ_YEARS][s.pDS][s.NG][s.pDB]),
    data_active        (extents[s.pDS][s.NG][s.pAG]), // 1 year only
    active_last_year_  (extents[s.pDS][s.NG][s.pAG]), // 1 year only
    infections_        (extents[s.NG][s.pAG]), // avoid repeated allocations
    art_elig_          (extents[s.NG][s.hAG][s.hDS]),
    art_init_          (extents[s.NG][s.hAG][s.hDS]),
    hiv_by_agrp_       (extents[s.NG][s.hAG]),
    num_death_         (extents[s.pDS][s.NG][s.pAG]),
    death_by_agrp_     (extents[s.NG][s.hAG]),
    migrant_by_agrp_   (extents[s.NG][s.hAG]),
    hiv_aging_prob_    (extents[s.NG][s.hAG]),
    num_migrate_       (extents[s.NG][s.pAG]),
    migrate_prob_      (extents[s.NG][s.pAG]),
    num_adjust_        (extents[s.NG][s.pAG]),
    num_adjust_by_agrp_(extents[s.NG][s.hAG])
  {
// Non class init
    if (s.MODEL == 2 && s.pDB == 1)
      Rf_error("Debut model state-space not exist, see update_fp_debut()");
  }
// Pop methods 
  void initiate(const Parameters& p, const StateSpace& s);
  void update_active_pop_to(int year, const StateSpace& s);
  void update_active_last_year(const StateSpace& s);
  void aging(const StateSpace& s);
  void add_entrants(const Parameters& p, const StateSpace& s);
  void sexual_debut(const Parameters& p, const StateSpace& s);
  void update_hiv_aging_prob(const StateSpace& s);
  dvec entrant_art(const Parameters& p, const StateSpace& s);
  void deaths(const Parameters& p, const StateSpace& s);
  void migration(const Parameters& p, const StateSpace& s);
  void update_fertile(const Parameters& p, const StateSpace& s);
  void adjust_pop(const Parameters& p, const StateSpace& s);
  void cal_prev_pregant (const hivC& hivpop, const artC& artpop,
                         const Parameters& p, const StateSpace& s);
  void save_prev_n_inc (const StateSpace& s);
  void infect_mix (int ii, const Parameters& p, const StateSpace& s);
  void infect_spec (const hivC& hivpop, const artC& artpop, int time_step,
                    const Parameters& p, const StateSpace& s);
  void epp_disease_model_direct (hivC& hivpop, artC& artpop,
                                 const Parameters& p, const StateSpace& s);
  double calc_rtrend_rt (int ts, double time_step, const StateSpace& s);
  void update_rvec (double time_step, const Parameters& p, const StateSpace& s);
  void update_infection (const StateSpace& s);
  void remove_hiv_death (const hivC& hivpop, const artC& artpop,
                         const Parameters& p, const StateSpace& s);
  void update_preg (const hivC& hivpop, const artC& artpop,
                    const Parameters& p, const StateSpace& s);
  dvec art_initiate (const dvec& art_curr, int time_step,
                     const Parameters& p, const StateSpace& s);
  void art_distribute (const dvec& art_need, const Parameters& p,
                       const StateSpace& s);
  void epp_art_init (hivC& hivpop, artC& artpop, int time_step,
                     const Parameters& p, const StateSpace& s);
public: // Pop fields
  boost4D_ptr data; // pointer to pop_data_sexp, the same for others
  dvec        birth_age;
  dvec        birth_agrp;
  boost1D_ptr prev15to49;  
  boost1D_ptr incid15to49;
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
  boost3D     data_all; // all populations in the year requested
  boost4D_ptr data_db; // debut only population
  boost3D     data_active;
  boost3D     active_last_year_;
  double      artcov[2] = {0.0, 0.0}; // initially no one on treatment
  double      prev_last = 0.0; // = 0 last time step prevalence
  boost2D     infections_;
  boost3D     art_elig_;
  boost3D     art_init_;
  boost2D     hiv_by_agrp_;
  boost3D     num_death_;
  boost2D     death_by_agrp_;
  boost2D     migrant_by_agrp_;
  boost2D     hiv_aging_prob_;
  boost2D     num_migrate_;
  boost2D     migrate_prob_;
  boost2D     num_adjust_;
  boost2D     num_adjust_by_agrp_;
  dvec        entrant_art_;
};

// HIV class
class hivC {
public: // inits
  hivC(outputSEXP& O, const StateSpace& s) :
// Boost array class init
  data(REAL(O.hivpop), extents[s.PROJ_YEARS][s.NG][s.hAG][s.hDS]),
  data_db        (extents[s.PROJ_YEARS][s.NG][s.hAG][s.hDS]), // later return this as well
  grad           (extents[s.NG][s.hAG][s.hDS]),
  grad_db        (extents[s.NG][s.hAG][s.hDS]),
  data_all       (extents[s.NG][s.hAG][s.hDS]),
  death_         (extents[s.NG][s.hAG][s.hDS]),
  death_db_      (extents[s.NG][s.hAG][s.hDS]),
  cd4_mort_      (extents[s.NG][s.hAG][s.hDS]),
  infect_by_agrp_(extents[s.NG][s.hAG]) {};
// methods
  void aging(const boost2D& ag_prob, const StateSpace& s);
  void add_entrants(const dvec& artYesNo, const Parameters& p, const StateSpace& s);
  void sexual_debut(const Parameters& p, const StateSpace& s);
  void deaths(const boost2D& survival_pr, const StateSpace& s);
  void migration(const boost2D& migration_pr, const StateSpace& s);
  void update_infection(const boost2D& new_infect, const Parameters& p, const StateSpace& s);
  void scale_cd4_mort(artC& artpop, const Parameters& p, const StateSpace& s);
  void grad_progress(const Parameters& p, const StateSpace& s);
  dvec eligible_for_art(const Parameters& p, const StateSpace& s);
  void distribute_artinit(boost3D& artinit, artC& artpop, const StateSpace& s);
  void add_grad_to_pop(const StateSpace& s);
  void adjust_pop(const boost2D& adj_prob, const StateSpace& s);
public: // fields
  boost4D_ptr data;
  boost4D     data_db; // debut only population
  boost3D     grad;
  boost3D     grad_db;
  boost3D     data_all; // all populations in the year requested
  boost3D     death_; // death in this year
  boost3D     death_db_; // death in this year
  boost3D     cd4_mort_;
  boost2D     infect_by_agrp_;
};

// ART class
class artC {
public: // Inits
  artC(outputSEXP& O, const StateSpace& s) :
    data(REAL(O.artpop), extents[s.PROJ_YEARS][s.NG][s.hAG][s.hDS][s.hTS]),
    data_db   (extents[s.PROJ_YEARS][s.NG][s.hAG][s.hDS][s.hTS]),
    gradART   (extents[s.NG][s.hAG][s.hDS][s.hTS]),
    gradART_db(extents[s.NG][s.hAG][s.hDS][s.hTS]),
    death_    (extents[s.NG][s.hAG][s.hDS][s.hTS]),
    death_db_ (extents[s.NG][s.hAG][s.hDS][s.hTS]) {}
// Methods
  void aging(const boost2D& ag_prob, const StateSpace& s);
  void add_entrants(const dvec& artYesNo, const Parameters& p, const StateSpace& s);
  void sexual_debut(const Parameters& p, const StateSpace& s);
  void deaths(const boost2D& survival_pr, const StateSpace& s);
  void migration(const boost2D& migration_pr, const StateSpace& s);
  void grad_progress(const StateSpace& s);
  void art_dropout(hivC& hivpop, const Parameters& p, const StateSpace& s);
  void update_current_on_art(const StateSpace& s);
  void grad_init(const boost3D& artinit, const StateSpace& s);
  void grad_db_init(const boost3D& artinit_db, const StateSpace& s);
  void adjust_pop(const boost2D& adj_prob, const StateSpace& s);
  void count_death(const Parameters& p, const StateSpace& s);
public: // fields
  boost5D_ptr data;
  boost5D     data_db; // debut only population
  boost4D     gradART;
  boost4D     gradART_db;
  boost4D     death_;
  boost4D     death_db_;
  dvec        art_by_sex_ = {0.0, 0.0};
};

class Model
{
public:
  Model(outputSEXP& O, const StateSpace& s, const Parameters& p) :
    s(s),
    p(p),
    pop(O, s),
    hivpop(O, s),
    artpop(O, s)
    {}
  void initiate();
  void run(int t);
  void aging();
  void death();
  void migration();
  void adjust_pop();
  void infection_process();
  void save_outputs();
public:
  StateSpace s;
  Parameters p;
  popC pop;
  hivC hivpop;
  artC artpop;
};