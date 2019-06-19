#ifndef _EPPASM_MODEL_H_
#define _EPPASM_MODEL_H_

#include <boost/multi_array.hpp>

#include <iostream>

struct SsDim {

  SsDim(int ng,
	int p_ag,
	int p_ds,
	int h_ag,
	int h_ds,
	int h_ts,
	int proj_start,
	int proj_years,
	int sim_years,
	int hiv_steps_per_year,
	int t_art_start,
	int *h_ag_span) :
    ng(ng),
    p_ag(p_ag),
    p_ds(p_ds),
    h_ag(h_ag),
    h_ds(h_ds),
    h_ts(h_ts),
    proj_start(proj_start),
    proj_years(proj_years),
    sim_years(sim_years),
    hiv_steps_per_year(hiv_steps_per_year),
    proj_steps((sim_years - 1) * hiv_steps_per_year + 1),
    dt(1.0 / hiv_steps_per_year),
    t_art_start(t_art_start),
    h_ag_span(h_ag_span),
    h_ag_start(h_ag)
  {
    h_ag_start[0] = 0;
    for(int ha = 1; ha < h_ag; ha++)
      h_ag_start[ha] = h_ag_start[ha-1] + h_ag_span[ha-1];
  }

  // State space dimensions
  
  const int ng;
  const int p_ag;
  const int p_ds;
  const int h_ag;
  const int h_ds;
  const int h_ts;

  // Projection metadata 
  const int proj_start;
  const int proj_years;         // maximum projection years defined by inputs
  const int sim_years;          // number of years to simulate (<= proj_years)
  const int hiv_steps_per_year; 
  const int proj_steps;         // number of HIV time steps
  const double dt;

  // const int t_hiv_start;        // first year of HIV infections
  // const int hts_hiv_start;      // time step of HIV start
  const int t_art_start;        // first year ART available
  

  // Indices for accessing data

  const size_t i_hivn {0};
  const size_t i_hivp {1};
  const size_t i_male {0};
  const size_t i_female {1};

  const size_t pop_age_start {15};

  // Indices for iterating over age subsets
  const size_t i_p_ag_fert {0};  // first age index for fertility
  const size_t r_p_ag_fert {35}; // number of indices for fertility

  const size_t i_p_ag_15to49 {0};  // first age index for age 15-49
  const size_t r_p_ag_15to49 {35}; // number of indices for age 15-49

  const size_t i_p_ag_15plus {0};  // first age index for age 15+
  const size_t r_p_ag_15plus {66}; // number of indices for age 15+

  const size_t i_h_ag_fert {0};  // first age index for fertility
  const size_t r_h_ag_fert {8}; // number of indices for fertility

  const size_t i_h_ag_15to49 {0};  // first age index for age 15-49
  const size_t r_h_ag_15to49 {8}; // number of indices for age 15-49

  const size_t i_h_ag_15plus {0};  // first age index for age 15+
  const size_t r_h_ag_15plus {9};  // numer of indices for age 15+


  const int* h_ag_span;
  std::vector<int> h_ag_start;
    
};


struct DemogParam {

  boost::const_multi_array_ref<double, 2> basepop;
  boost::const_multi_array_ref<double, 3> Sx;
  boost::const_multi_array_ref<double, 3> netmigr;
  boost::const_multi_array_ref<double, 2> asfr;
  boost::const_multi_array_ref<double, 2> srb;
  boost::const_multi_array_ref<double, 2> birthslag;
  boost::const_multi_array_ref<double, 2> cumsurv;
  boost::const_multi_array_ref<double, 2> cumnetmigr;
  const int flag_popadjust;
  boost::const_multi_array_ref<double, 3> targetpop;
  boost::const_multi_array_ref<double, 2> entrantpop;

  DemogParam(const SsDim &ss,
	     double *ptr_basepop,
	     double *ptr_Sx,
	     double *ptr_netmigr,
	     double *ptr_asfr,
	     double *ptr_srb,
	     double *ptr_birthslag,
	     double *ptr_cumsurv,
	     double *ptr_cumnetmigr,
	     int bin_popadjust,
	     double *ptr_targetpop,
	     double *ptr_entrantpop) :
    basepop(   ptr_basepop,    boost::extents[ss.ng][ss.p_ag]),
    Sx(        ptr_Sx,         boost::extents[ss.proj_years][ss.ng][ss.p_ag]),
    netmigr(   ptr_netmigr,    boost::extents[ss.proj_years][ss.ng][ss.p_ag]),
    asfr(      ptr_asfr,       boost::extents[ss.proj_years][ss.r_p_ag_fert]),
    srb(       ptr_srb,        boost::extents[ss.proj_years][ss.ng]),
    birthslag( ptr_birthslag,  boost::extents[ss.proj_years][ss.ng]),
    cumsurv(   ptr_cumsurv,    boost::extents[ss.proj_years][ss.ng]),
    cumnetmigr(ptr_cumnetmigr, boost::extents[ss.proj_years][ss.ng]) ,
    flag_popadjust(bin_popadjust),
    targetpop( ptr_targetpop,  boost::extents[ss.proj_years][ss.ng][ss.p_ag]),
    entrantpop(ptr_entrantpop, boost::extents[ss.proj_years][ss.ng])
  {}

};

struct NaturalHistoryParam {

  const int art0mos {0};
  const int art6mos {1};
  const int art1yr {2};

  const int i_h_cd4_350 {2};

  const double art_stage_prog_rate {2.0}; // HARD CODED: ART stage progression rate

  boost::const_multi_array_ref<double, 3> cd4_initdist;
  boost::const_multi_array_ref<double, 3> cd4_prog;
  boost::const_multi_array_ref<double, 3> cd4_mort;
  boost::const_multi_array_ref<double, 4> art_mort;
  boost::const_multi_array_ref<double, 2> artmx_timerr;
  boost::const_multi_array_ref<double, 3> frr_cd4;
  boost::const_multi_array_ref<double, 4> frr_art;
  
  NaturalHistoryParam(const SsDim &ss,
		      double *ptr_cd4_initdist,
		      double *ptr_cd4_prog,     
		      double *ptr_cd4_mort,     
		      double *ptr_art_mort,     
		      double *ptr_artmx_timerr,
		      double *ptr_frr_cd4,
		      double *ptr_frr_art) :
    cd4_initdist( ptr_cd4_initdist,  boost::extents[ss.ng][ss.h_ag][ss.h_ds]),
    cd4_prog(     ptr_cd4_prog,      boost::extents[ss.ng][ss.h_ag][ss.h_ds-1]),
    cd4_mort(     ptr_cd4_mort,      boost::extents[ss.ng][ss.h_ag][ss.h_ds]),
    art_mort(     ptr_art_mort,      boost::extents[ss.ng][ss.h_ag][ss.h_ds][ss.h_ts]),
    artmx_timerr( ptr_artmx_timerr,  boost::extents[ss.proj_years][ss.h_ts]),
    frr_cd4(      ptr_frr_cd4,       boost::extents[ss.proj_years][ss.r_h_ag_fert][ss.h_ds]),
    frr_art(      ptr_frr_art,       boost::extents[ss.proj_years][ss.r_h_ag_fert][ss.h_ds][ss.h_ts])
  {}
};

struct ArtData {

  boost::const_multi_array_ref<double, 2> artnum15plus;
  boost::const_multi_array_ref<int, 2> art15plus_isperc;
  
  const int *artcd4elig_idx;  // NOTE: 1-based indexing
  const double *special_pop_percelig;
  const double *preg_women_artelig;
  const double who34percelig;
  
  const double *art_dropout;
  const double *median_cd4init;
  
  const int *med_cd4init_cat;
  const int *med_cd4init_input;
  
  const int art_alloc_method;
  const double art_alloc_mxweight;
  const int scale_cd4_mort;

  ArtData(const SsDim &ss,
	  const double *ptr_artnum15plus,
	  const int *ptr_art15plus_isperc,
	  const int *t_artcd4elig_idx,
	  const double *t_special_pop_percelig,
	  const double *t_preg_women_artelig,
	  const double t_who34percelig,
	  const double *t_art_dropout,
	  const double *t_median_cd4init,
	  const int *t_med_cd4init_cat,
	  const int *t_med_cd4init_input,
	  const int t_art_alloc_method,
	  const double t_art_alloc_mxweight,
	  const int t_scale_cd4_mort) :
    artnum15plus(        ptr_artnum15plus,     boost::extents[ss.proj_years][ss.ng]),
    art15plus_isperc(    ptr_art15plus_isperc, boost::extents[ss.proj_years][ss.ng]),
    artcd4elig_idx(t_artcd4elig_idx),
    special_pop_percelig(t_special_pop_percelig),
    preg_women_artelig(t_preg_women_artelig),
    who34percelig(t_who34percelig),
    art_dropout(t_art_dropout),
    median_cd4init(t_median_cd4init),
    med_cd4init_cat(t_med_cd4init_cat),
    med_cd4init_input(t_med_cd4init_input),
    art_alloc_method(t_art_alloc_method),   
    art_alloc_mxweight(t_art_alloc_mxweight),
    scale_cd4_mort(t_scale_cd4_mort)
  {}
};


struct IncidenceParam {

private:

  // eppmod consts
  const static int EPP_RSPLINE {0};
  const static int EPP_RTREND {1};
  const static int EPP_DIRECTINCID {2};  // annual direct incidence inputs (as Spectrum)

  // reference age group for incidence input
  const static int INCIDPOP_15TO49 {0};
  const static int INCIDPOP_15PLUS {1};

public:

  const int eppmod; 

  const int i_p_ag_incidpop;  // start of incidence reference population
  const int r_p_ag_incidpop;  // span of incidence reference population

  // direct incidence model inputs
  const double *incidinput;
  
  // EPP model inputs 
  const int ts_epidemic_start;
  const double iota;
  const double rel_infect_art;
  const double *rspline_rvec;
  const double *rtrend_beta;
  const double rtrend_tstab;
  const double rtrend_r0;

  // incidence rate ratios for allocating incidence input
  const double *incrr_sex;
  boost::const_multi_array_ref<double, 3> incrr_age;

  IncidenceParam(const SsDim &ss,
		 const int t_eppmod,
		 const int t_incidpop_age,
		 const int t_ts_epidemic_start,
		 const double t_iota,
		 const double t_rel_infect_art,
		 const double *ptr_rspline_rvec,
		 const double *ptr_rtrend_beta,
		 const double t_rtrend_tstab,
		 const double t_rtrend_r0,
		 const double *ptr_incidinput,
		 const double *ptr_incrr_sex,
		 const double *ptr_incrr_age) :
    eppmod(t_eppmod),
    i_p_ag_incidpop(t_incidpop_age == INCIDPOP_15TO49 ? ss.i_p_ag_15to49 : ss.i_p_ag_15plus),
    r_p_ag_incidpop(t_incidpop_age == INCIDPOP_15TO49 ? ss.r_p_ag_15to49 : ss.r_p_ag_15plus),
    incidinput(ptr_incidinput),
    ts_epidemic_start(t_ts_epidemic_start),
    iota(t_iota),
    rel_infect_art(t_rel_infect_art),
    rspline_rvec(ptr_rspline_rvec),
    rtrend_beta(ptr_rtrend_beta),
    rtrend_tstab(t_rtrend_tstab),
    rtrend_r0(t_rtrend_r0),
    incrr_sex(ptr_incrr_sex),
    incrr_age(ptr_incrr_age, boost::extents[ss.proj_years][ss.ng][ss.p_ag])
  {}
};

struct PaediatricHivParam {

  const int use_entrantprev; 
  const double *verttrans_lag;
  const double *paedsurv_lag;
  const double netmig_hivprob;
  boost::const_multi_array_ref<double, 2> entrantprev;
  boost::const_multi_array_ref<double, 2> entrantartcov;
  boost::const_multi_array_ref<double, 3> paedsurv_cd4dist;
  boost::const_multi_array_ref<double, 4> paedsurv_artcd4dist;
  
  PaediatricHivParam(const SsDim &ss,
		     const int t_use_entrantprev,
		     const double *ptr_verttrans_lag,
		     const double *ptr_paedsurv_lag,
		     const double netmig_hivprob,
		     const double *ptr_entrantprev,
		     const double *ptr_entrantartcov,
		     const double *ptr_paedsurv_cd4dist,
		     const double *ptr_paedsurv_artcd4dist) :
    use_entrantprev(t_use_entrantprev),
    verttrans_lag(ptr_verttrans_lag),
    paedsurv_lag(ptr_paedsurv_lag),
    netmig_hivprob(netmig_hivprob),
    entrantprev(ptr_entrantprev, boost::extents[ss.proj_years][ss.ng]),
    entrantartcov(ptr_entrantartcov, boost::extents[ss.proj_years][ss.ng]),
    paedsurv_cd4dist(ptr_paedsurv_cd4dist, boost::extents[ss.proj_years][ss.ng][ss.h_ds]),
    paedsurv_artcd4dist(ptr_paedsurv_artcd4dist, boost::extents[ss.proj_years][ss.ng][ss.h_ds][ss.h_ts])
  {}

};
    

class Model {

  // data
private:
  const SsDim ss;
  const DemogParam demp;
  const PaediatricHivParam paedhp;
  const NaturalHistoryParam nhp;
  const ArtData artp;
  
public:
  boost::multi_array_ref<double, 4> pop;        // population by single-year age and HIV status
  boost::multi_array_ref<double, 4> hivpop;     // HIV population with stage stratification
  boost::multi_array_ref<double, 5> artpop;     // ART population with stage stratification
  boost::multi_array_ref<double, 3> infections;  // count of infections
  boost::multi_array_ref<double, 3> hivdeaths;   // count of HIV deaths
  boost::multi_array_ref<double, 3> natdeaths;   // count of non-HIV deaths
  boost::multi_array_ref<double, 3> popadjust;   // proportion population adjustment

  double *prev15to49;
  double *incid15to49;
  double *pregprev;
  double *entrantprev;

  double *prev15to49_ts;
  double *incid15to49_ts;
  double *rvec_ts;

  // internal calculations
  std::vector<double> pregprevlag;
  boost::multi_array<double, 2> birthslag;
  double hivn15to49;
  double prevcurr {0.0};
  double prevlast {0.0};

  int everARTelig_idx;
  int cd4elig_idx;
  int anyelig_idx;

  std::vector<double> births_by_ha;
  double births;

  boost::multi_array<double, 3> grad;
  boost::multi_array<double, 4> gradART;
  boost::multi_array<double, 2> hivdeaths_ha;
  

public:

  // constructors
  Model(const SsDim &t_ss,
	const DemogParam &t_demp,
	const PaediatricHivParam &t_paedhp,
	const NaturalHistoryParam &t_nhp,
	const ArtData &t_artp,
	double* ptr_pop,
	double* ptr_hivpop,
	double* ptr_artpop,
	double* ptr_infections,
	double* ptr_hivdeaths,
	double* ptr_natdeaths,
	double* ptr_popadjust,
	double* ptr_prev15to49,
	double* ptr_incid15to49,
	double* ptr_pregprev,
	double* ptr_entrantprev,
	double* ptr_prev15to49_ts,
	double* ptr_incid15to49_ts,
	double* ptr_rvec_ts) :
    ss(t_ss),
    demp(t_demp),
    paedhp(t_paedhp),
    nhp(t_nhp),
    artp(t_artp),
    pop(ptr_pop,       boost::extents[ss.proj_years][ss.p_ds][ss.ng][ss.p_ag]),  
    hivpop(ptr_hivpop, boost::extents[ss.proj_years][ss.ng][ss.h_ag][ss.h_ds]),
    artpop(ptr_artpop, boost::extents[ss.proj_years][ss.ng][ss.h_ag][ss.h_ds][ss.h_ts]),
    infections(ptr_infections, boost::extents[ss.proj_years][ss.ng][ss.p_ag]),
    hivdeaths(ptr_hivdeaths, boost::extents[ss.proj_years][ss.ng][ss.p_ag]),
    natdeaths(ptr_natdeaths, boost::extents[ss.proj_years][ss.ng][ss.p_ag]),
    popadjust(ptr_popadjust, boost::extents[ss.proj_years][ss.ng][ss.p_ag]),
    prev15to49(ptr_prev15to49),
    incid15to49(ptr_incid15to49),
    pregprev(ptr_pregprev),
    entrantprev(ptr_entrantprev),
    prev15to49_ts(ptr_prev15to49_ts),
    incid15to49_ts(ptr_incid15to49_ts),
    rvec_ts(ptr_rvec_ts),
    pregprevlag(ss.proj_years, 0),
    birthslag(demp.birthslag),
    everARTelig_idx(ss.h_ds),
    births_by_ha(ss.r_h_ag_fert, 0),
    births(0),
    grad(boost::extents[ss.ng][ss.h_ag][ss.h_ds]),
    gradART(boost::extents[ss.ng][ss.h_ag][ss.h_ds][ss.h_ts]),
    hivdeaths_ha(boost::extents[ss.ng][ss.h_ag])
  {

    // initialise memory to zero
    // some of these are required, some not required but conservative
    memset(pop.data(),    0, pop.num_elements() * sizeof(double));
    memset(hivpop.data(), 0, hivpop.num_elements() * sizeof(double));
    memset(artpop.data(), 0, artpop.num_elements() * sizeof(double));
    memset(infections.data(), 0, infections.num_elements() * sizeof(double));
    memset(hivdeaths.data(),  0, hivdeaths.num_elements() * sizeof(double));
    memset(natdeaths.data(),  0, natdeaths.num_elements() * sizeof(double));
    memset(popadjust.data(),  0, popadjust.num_elements() * sizeof(double));

    memset(prev15to49, 0, ss.proj_years * sizeof(double));
    memset(incid15to49, 0, ss.proj_years * sizeof(double));
    memset(pregprev, 0, ss.proj_years * sizeof(double));
    memset(entrantprev, 0, ss.proj_years * sizeof(double));
  }

  // methods
  void initialise();
  void ageing(int t);
  void pop_entrants(int t);
  void death_and_migration(int t);
  void fertility(int t);
  void adjust_population(int t);

  void initialise_hiv_ts(int t);
  void hiv_progression_mortality(int t);
  void art_initiation(int t, int hts);
  void do_hts_step(int t);
  
  void remove_hiv_deaths(int t);

};


#endif
