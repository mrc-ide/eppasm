#include <unsupported/Eigen/CXX11/Tensor>

#include <R.h>
#include <Rinternals.h>


#define AGE_START 15

#define NG 2
#define pAG 66
#define pDS 2

#define pIDX_FERT 0
#define pAG_FERT 35
#define pIDX_15TO49 0
#define pAG_15TO49  35
#define pIDX_15PLUS 0
#define pAG_15PLUS  66

#define hAG 9
#define hDS 7
#define hTS 3

#define hIDX_FERT 0
#define hAG_FERT 8
#define hIDX_15TO49 0
#define hAG_15TO49  8
#define hIDX_15PLUS 0
#define hAG_15PLUS  9

#define hIDX_CD4_350 2

#define MALE 0
#define FEMALE 1

#define HIVN 0
#define HIVP 1

#define ART0MOS 0

#define EPP_RSPLINE 0
#define EPP_RTREND 1
#define EPP_DIRECTINCID 2  // annual direct incidence inputs (as Spectrum)

#define INCIDMOD_EPPSPEC 0
#define INCIDMOD_TRANSM 1

#define INCIDPOP_15TO49 0 // age range corresponding to incidence input
#define INCIDPOP_15PLUS 1

using namespace Eigen;


// Function declarations
SEXP getListElement(SEXP list, const char *str);
int checkListElement(SEXP list, const char *str);

void calc_infections_eppspectrum_tensor(const TensorMap<Tensor<double, 4>> pop_t,
					const TensorMap<Tensor<double, 4>> hivpop_t,
					const TensorMap<Tensor<double, 5>> artpop_t,					
					double r_ts,
					double relinfectART,
					double iota,
					double *incrr_sex,
					const TensorMap<Tensor<double, 3>> incrr_age,
					int t_ART_start,
					double DT,
					int t,
					int hts,
					int *hAG_START,
					int *hAG_SPAN,
					double *prevcurr,
					double *incrate15to49_ts,
					TensorFixedSize<double, Sizes<pAG, NG>>& infections_ts);

extern "C" {

  SEXP eppasm_tensorC(SEXP s_fp){

    ////////////////////////////////
    ////  set parameter values  ////
    ////////////////////////////////

    // state space dimensions
    SEXP s_ss = getListElement(s_fp, "ss");
    int PROJ_YEARS = *INTEGER(getListElement(s_ss, "PROJ_YEARS"));
    int HIVSTEPS_PER_YEAR = *INTEGER(getListElement(s_ss, "hiv_steps_per_year"));
    double DT = 1.0/HIVSTEPS_PER_YEAR;
    int *hAG_SPAN = INTEGER(getListElement(s_ss, "h.ag.span"));
    double *h_art_stage_dur = REAL(getListElement(s_ss, "h_art_stage_dur"));

    int hAG_START[hAG];
    hAG_START[0] = 0;
    for(int ha = 1; ha < hAG; ha++)
      hAG_START[ha] = hAG_START[ha-1] + hAG_SPAN[ha-1];

    int SIM_YEARS = *INTEGER(getListElement(s_fp, "SIM_YEARS"));
    double *projsteps = REAL(getListElement(s_fp, "proj.steps"));

    // demographic projection
    TensorMap<Tensor<double, 2>> basepop(REAL(getListElement(s_fp, "basepop")), pAG, NG);
    TensorMap<Tensor<double, 3>> Sx(REAL(getListElement(s_fp, "Sx")), pAG, NG, PROJ_YEARS);
    TensorMap<Tensor<double, 3>> netmigr(REAL(getListElement(s_fp, "netmigr")), pAG, NG, PROJ_YEARS);
    TensorMap<Tensor<double, 2>> asfr(REAL(getListElement(s_fp, "asfr")), pAG_FERT, PROJ_YEARS);
    TensorMap<Tensor<double, 2>> srb(REAL(getListElement(s_fp, "srb")), NG, PROJ_YEARS);
    TensorMap<Tensor<double, 2>> birthslag(REAL(getListElement(s_fp, "birthslag")), NG, PROJ_YEARS);
    TensorMap<Tensor<double, 2>> cumsurv(REAL(getListElement(s_fp, "cumsurv")), NG, PROJ_YEARS);
    TensorMap<Tensor<double, 2>> cumnetmigr(REAL(getListElement(s_fp, "cumnetmigr")), NG, PROJ_YEARS);

    int bin_popadjust = *INTEGER(getListElement(s_fp, "popadjust"));
    double *ptr_targetpop;
    double *ptr_entrantpop;
    if(bin_popadjust){
      ptr_targetpop = REAL(getListElement(s_fp, "targetpop"));
      ptr_entrantpop = REAL(getListElement(s_fp, "entrantpop"));
    }
    TensorMap<Tensor<double, 3>> targetpop(ptr_targetpop, pAG, NG, PROJ_YEARS);
    TensorMap<Tensor<double, 2>> entrantpop(ptr_entrantpop, NG, PROJ_YEARS);

    // disease progression
    TensorMap<Tensor<double, 3>> cd4_initdist(REAL(getListElement(s_fp, "cd4_initdist")), hDS, hAG, NG);
    TensorMap<Tensor<double, 3>> cd4_prog(REAL(getListElement(s_fp, "cd4_prog")), hDS-1, hAG, NG);
    TensorMap<Tensor<double, 3>> cd4_mort(REAL(getListElement(s_fp, "cd4_mort")), hDS, hAG, NG);
    TensorMap<Tensor<double, 4>> art_mort(REAL(getListElement(s_fp, "art_mort")), hTS, hDS, hAG, NG);
    TensorMap<Tensor<double, 2>> artmx_timerr(REAL(getListElement(s_fp, "artmx_timerr")), hTS, PROJ_YEARS);

    // sub-fertility
    TensorMap<Tensor<double, 3>> frr_cd4(REAL(getListElement(s_fp, "frr_cd4")), hDS, hAG_FERT, PROJ_YEARS);
    TensorMap<Tensor<double, 4>> frr_art(REAL(getListElement(s_fp, "frr_art")), hTS, hDS, hAG_FERT, PROJ_YEARS);

    // ART inputs
    int t_ART_start = *INTEGER(getListElement(s_fp, "tARTstart")) - 1; // -1 for 0-based indexing in C vs. 1-based in R
    TensorMap<Tensor<double, 2>> artnum15plus(REAL(getListElement(s_fp, "art15plus_num")), NG, PROJ_YEARS);
    TensorMap<Tensor<int, 2>> art15plus_isperc(LOGICAL(getListElement(s_fp, "art15plus_isperc")), NG, PROJ_YEARS);;

    int *artcd4elig_idx = INTEGER(getListElement(s_fp, "artcd4elig_idx"));  // NOTE: 1-based indexing
    double *specpop_percelig = REAL(getListElement(s_fp, "specpop_percelig"));
    double *pw_artelig = REAL(getListElement(s_fp, "pw_artelig"));
    double who34percelig = *REAL(getListElement(s_fp, "who34percelig"));

    double *art_dropout = REAL(getListElement(s_fp, "art_dropout"));
    double *median_cd4init = REAL(getListElement(s_fp, "median_cd4init"));

    int *med_cd4init_cat = INTEGER(getListElement(s_fp, "med_cd4init_cat"));
    int *med_cd4init_input = INTEGER(getListElement(s_fp, "med_cd4init_input"));

    int art_alloc_method = *INTEGER(getListElement(s_fp, "art_alloc_method"));
    double art_alloc_mxweight = *REAL(getListElement(s_fp, "art_alloc_mxweight"));
    int scale_cd4_mort = *INTEGER(getListElement(s_fp, "scale_cd4_mort"));


    // incidence model
    // double *prev15to49 = REAL(getListElement(s_fp, "prev15to49"));
    int incidmod = *INTEGER(getListElement(s_fp, "incidmodInt"));
    double *incrr_sex;
    double *mf_transm_rr;
    double *relsexact_cd4cat;
    if(incidmod == INCIDMOD_EPPSPEC)
      incrr_sex = REAL(getListElement(s_fp, "incrr_sex"));
    else {
      mf_transm_rr = REAL(getListElement(s_fp, "mf_transm_rr"));
      relsexact_cd4cat = REAL(getListElement(s_fp, "relsexact_cd4cat"));
    }
    
    TensorMap<Tensor<double, 3>> incrr_age(REAL(getListElement(s_fp, "incrr_age")), pAG, NG, PROJ_YEARS);

    int eppmod = *INTEGER(getListElement(s_fp, "eppmodInt"));

    double *incidinput;
    int pIDX_INCIDPOP, pAG_INCIDPOP;
    double tsEpidemicStart, iota, relinfectART;
    double *rspline_rvec;
    double *rtrend_beta, rtrend_tstab, rtrend_r0;
    if(eppmod == EPP_DIRECTINCID){
      incidinput = REAL(getListElement(s_fp, "incidinput"));
      pIDX_INCIDPOP = 0;
      if(*INTEGER(getListElement(s_fp, "incidpopage")) == INCIDPOP_15TO49)
	pAG_INCIDPOP = pAG_15TO49;
      else
	pAG_INCIDPOP = pAG_15PLUS;
    } else {
      relinfectART = *REAL(getListElement(s_fp, "relinfectART"));
      tsEpidemicStart = *REAL(getListElement(s_fp, "tsEpidemicStart"));
      iota = *REAL(getListElement(s_fp, "iota"));
      
      if(eppmod == EPP_RSPLINE)
	rspline_rvec = REAL(getListElement(s_fp, "rvec"));
      else if(eppmod == EPP_RTREND){
	SEXP s_rtrend = getListElement(s_fp, "rtrend");
	rtrend_beta = REAL(getListElement(s_rtrend, "beta"));
	rtrend_tstab = *REAL(getListElement(s_rtrend, "tStabilize"));
	rtrend_r0 = *REAL(getListElement(s_rtrend, "r0"));
      }
    }

    // vertical transmission and survival
    double *verttrans_lag = REAL(getListElement(s_fp, "verttrans_lag"));
    double *paedsurv_lag = REAL(getListElement(s_fp, "paedsurv_lag"));
    double netmig_hivprob = *REAL(getListElement(s_fp, "netmig_hivprob"));
    double netmighivsurv = *REAL(getListElement(s_fp, "netmighivsurv"));

    double *a_entrantprev;
    int use_entrantprev = checkListElement(s_fp, "entrantprev");
    if(use_entrantprev)
      a_entrantprev = REAL(getListElement(s_fp, "entrantprev"));
    TensorMap<Tensor<double, 2>> entrantprev(a_entrantprev, NG, PROJ_YEARS);

    double *a_entrantartcov;
    if(checkListElement(s_fp, "entrantartcov"))
      a_entrantartcov = REAL(getListElement(s_fp, "entrantartcov"));
    else {
      a_entrantartcov = (double*) R_alloc(PROJ_YEARS, sizeof(double));
      memset(a_entrantartcov, 0, PROJ_YEARS*sizeof(double));
    }
    TensorMap<Tensor<double, 2>> entrantartcov(a_entrantartcov, NG, PROJ_YEARS);

    TensorMap<Tensor<double, 3>> paedsurv_cd4dist(REAL(getListElement(s_fp, "paedsurv_cd4dist")), hDS, NG, PROJ_YEARS);
    TensorMap<Tensor<double, 4>> paedsurv_artcd4dist(REAL(getListElement(s_fp, "paedsurv_artcd4dist")), hTS, hDS, NG, PROJ_YEARS);

    // initialize output
    SEXP s_pop = PROTECT(allocVector(REALSXP, pAG * NG * pDS * PROJ_YEARS));
    SEXP s_pop_dim = PROTECT(allocVector(INTSXP, 4));
    INTEGER(s_pop_dim)[0] = pAG;
    INTEGER(s_pop_dim)[1] = NG;
    INTEGER(s_pop_dim)[2] = pDS;
    INTEGER(s_pop_dim)[3] = PROJ_YEARS;
    setAttrib(s_pop, R_DimSymbol, s_pop_dim);

    SEXP s_hivpop = PROTECT(allocVector(REALSXP, hDS * hAG * NG * PROJ_YEARS));
    SEXP s_hivpop_dim = PROTECT(allocVector(INTSXP, 4));
    INTEGER(s_hivpop_dim)[0] = hDS;
    INTEGER(s_hivpop_dim)[1] = hAG;
    INTEGER(s_hivpop_dim)[2] = NG;
    INTEGER(s_hivpop_dim)[3] = PROJ_YEARS;
    setAttrib(s_hivpop, R_DimSymbol, s_hivpop_dim);
    setAttrib(s_pop, install("hivpop"), s_hivpop);

    SEXP s_artpop = PROTECT(allocVector(REALSXP, hTS * hDS * hAG * NG * PROJ_YEARS));
    SEXP s_artpop_dim = PROTECT(allocVector(INTSXP, 5));
    INTEGER(s_artpop_dim)[0] = hTS;
    INTEGER(s_artpop_dim)[1] = hDS;
    INTEGER(s_artpop_dim)[2] = hAG;
    INTEGER(s_artpop_dim)[3] = NG;
    INTEGER(s_artpop_dim)[4] = PROJ_YEARS;
    setAttrib(s_artpop, R_DimSymbol, s_artpop_dim);
    setAttrib(s_pop, install("artpop"), s_artpop);

    SEXP s_infections = PROTECT(allocVector(REALSXP, pAG * NG * PROJ_YEARS));
    SEXP s_infections_dim = PROTECT(allocVector(INTSXP, 3));
    INTEGER(s_infections_dim)[0] = pAG;
    INTEGER(s_infections_dim)[1] = NG;
    INTEGER(s_infections_dim)[2] = PROJ_YEARS;
    setAttrib(s_infections, R_DimSymbol, s_infections_dim);
    setAttrib(s_pop, install("infections"), s_infections);

    SEXP s_hivdeaths = PROTECT(allocVector(REALSXP, pAG * NG * PROJ_YEARS));
    SEXP s_hivdeaths_dim = PROTECT(allocVector(INTSXP, 3));
    INTEGER(s_hivdeaths_dim)[0] = pAG;
    INTEGER(s_hivdeaths_dim)[1] = NG;
    INTEGER(s_hivdeaths_dim)[2] = PROJ_YEARS;
    setAttrib(s_hivdeaths, R_DimSymbol, s_hivdeaths_dim);
    setAttrib(s_pop, install("hivdeaths"), s_hivdeaths);

    SEXP s_natdeaths = PROTECT(allocVector(REALSXP, pAG * NG * PROJ_YEARS));
    SEXP s_natdeaths_dim = PROTECT(allocVector(INTSXP, 3));
    INTEGER(s_natdeaths_dim)[0] = pAG;
    INTEGER(s_natdeaths_dim)[1] = NG;
    INTEGER(s_natdeaths_dim)[2] = PROJ_YEARS;
    setAttrib(s_natdeaths, R_DimSymbol, s_natdeaths_dim);
    setAttrib(s_pop, install("natdeaths"), s_natdeaths);
    
    SEXP s_aidsdeaths_noart = PROTECT(allocVector(REALSXP, hDS * hAG * NG * PROJ_YEARS));
    SEXP s_aidsdeaths_noart_dim = PROTECT(allocVector(INTSXP, 4));
    INTEGER(s_aidsdeaths_noart_dim)[0] = hDS;
    INTEGER(s_aidsdeaths_noart_dim)[1] = hAG;
    INTEGER(s_aidsdeaths_noart_dim)[2] = NG;
    INTEGER(s_aidsdeaths_noart_dim)[3] = PROJ_YEARS;
    setAttrib(s_aidsdeaths_noart, R_DimSymbol, s_aidsdeaths_noart_dim);
    setAttrib(s_pop, install("aidsdeaths_noart"), s_aidsdeaths_noart);

    SEXP s_aidsdeaths_art = PROTECT(allocVector(REALSXP, hTS * hDS * hAG * NG * PROJ_YEARS));
    SEXP s_aidsdeaths_art_dim = PROTECT(allocVector(INTSXP, 5));
    INTEGER(s_aidsdeaths_art_dim)[0] = hTS;
    INTEGER(s_aidsdeaths_art_dim)[1] = hDS;
    INTEGER(s_aidsdeaths_art_dim)[2] = hAG;
    INTEGER(s_aidsdeaths_art_dim)[3] = NG;
    INTEGER(s_aidsdeaths_art_dim)[4] = PROJ_YEARS;
    setAttrib(s_aidsdeaths_art, R_DimSymbol, s_aidsdeaths_art_dim);
    setAttrib(s_pop, install("aidsdeaths_art"), s_aidsdeaths_art);

    SEXP s_popadjust = PROTECT(allocVector(REALSXP, pAG * NG * PROJ_YEARS));
    SEXP s_popadjust_dim = PROTECT(allocVector(INTSXP, 3));
    INTEGER(s_popadjust_dim)[0] = pAG;
    INTEGER(s_popadjust_dim)[1] = NG;
    INTEGER(s_popadjust_dim)[2] = PROJ_YEARS;
    setAttrib(s_popadjust, R_DimSymbol, s_popadjust_dim);
    setAttrib(s_pop, install("popadjust"), s_popadjust);

    SEXP s_artinit = PROTECT(allocVector(REALSXP, hDS * hAG * NG * PROJ_YEARS));
    SEXP s_artinit_dim = PROTECT(allocVector(INTSXP, 4));
    INTEGER(s_artinit_dim)[0] = hDS;
    INTEGER(s_artinit_dim)[1] = hAG;
    INTEGER(s_artinit_dim)[2] = NG;
    INTEGER(s_artinit_dim)[3] = PROJ_YEARS;
    setAttrib(s_artinit, R_DimSymbol, s_artinit_dim);
    setAttrib(s_pop, install("artinit"), s_artinit);

    SEXP s_pregprevlag = PROTECT(allocVector(REALSXP, PROJ_YEARS));
    setAttrib(s_pop, install("pregprevlag"), s_pregprevlag);

    SEXP s_incrate15to49_ts = PROTECT(allocVector(REALSXP, (PROJ_YEARS-1) * HIVSTEPS_PER_YEAR));
    setAttrib(s_pop, install("incrate15to49_ts"), s_incrate15to49_ts);

    SEXP s_prev15to49_ts = PROTECT(allocVector(REALSXP, (PROJ_YEARS-1) * HIVSTEPS_PER_YEAR));
    setAttrib(s_pop, install("prev15to49_ts"), s_prev15to49_ts);

    SEXP s_rvec_ts = PROTECT(allocVector(REALSXP, (PROJ_YEARS-1) * HIVSTEPS_PER_YEAR));
    setAttrib(s_pop, install("rvec_ts"), s_rvec_ts);

    SEXP s_prev15to49 = PROTECT(allocVector(REALSXP, PROJ_YEARS));
    setAttrib(s_pop, install("prev15to49"), s_prev15to49);

    SEXP s_pregprev = PROTECT(allocVector(REALSXP, PROJ_YEARS));
    setAttrib(s_pop, install("pregprev"), s_pregprev);

    SEXP s_incid15to49 = PROTECT(allocVector(REALSXP, PROJ_YEARS));
    setAttrib(s_pop, install("incid15to49"), s_incid15to49);

    SEXP s_entrantprev_out = PROTECT(allocVector(REALSXP, PROJ_YEARS));
    setAttrib(s_pop, install("entrantprev"), s_entrantprev_out);

    double *hivn15to49 = (double*) R_alloc(PROJ_YEARS, sizeof(double));
    double *hivp15to49 = (double*) R_alloc(PROJ_YEARS, sizeof(double));
    memset(hivn15to49, 0, PROJ_YEARS*sizeof(double));
    memset(hivp15to49, 0, PROJ_YEARS*sizeof(double));

    
    // initialize population

    // population by single-year age
    TensorMap<Tensor<double, 4>> pop_t(REAL(s_pop), pAG, NG, pDS, PROJ_YEARS);
    pop_t.setZero();
    
    for(int g = 0; g < NG; g++)
      for(int a = 0; a < pAG; a++){
        pop_t(a, g, HIVN, 0) = basepop(a, g);
        if(a >= pIDX_15TO49 & a < pIDX_15TO49+pAG_15TO49)
          hivn15to49[0] += basepop(a, g);
      }

    // HIV population with stage stratification
    TensorMap<Tensor<double, 4>> hivpop_t(REAL(s_hivpop), hDS, hAG, NG, PROJ_YEARS);
    hivpop_t.setZero();

    // ART population with stage stratification
    TensorMap<Tensor<double, 5>> artpop_t(REAL(s_artpop), hTS, hDS, hAG, NG, PROJ_YEARS);
    artpop_t.setZero();


    // assign references to output arrays

    TensorMap<Tensor<double, 3>> infections(REAL(s_infections), pAG, NG, PROJ_YEARS);
    infections.setZero();

    TensorMap<Tensor<double, 3>> hivdeaths(REAL(s_hivdeaths), pAG, NG, PROJ_YEARS);
    hivdeaths.setZero();

    TensorMap<Tensor<double, 3>> natdeaths(REAL(s_natdeaths), pAG, NG, PROJ_YEARS);
    natdeaths.setZero();

    TensorMap<Tensor<double, 4>> aidsdeaths_noart(REAL(s_aidsdeaths_noart), hDS, hAG, NG, PROJ_YEARS);
    aidsdeaths_noart.setZero();

    TensorMap<Tensor<double, 5>> aidsdeaths_art(REAL(s_aidsdeaths_art), hTS, hDS, hAG, NG, PROJ_YEARS);
    aidsdeaths_art.setZero();

    TensorMap<Tensor<double, 3>> popadjust(REAL(s_popadjust), pAG, NG, PROJ_YEARS);
    popadjust.setZero();

    TensorMap<Tensor<double, 4>> artinit(REAL(s_artinit), hDS, hAG, NG, PROJ_YEARS);
    artinit.setZero();

    
    double *incrate15to49_ts_out = REAL(s_incrate15to49_ts);
    memset(incrate15to49_ts_out, 0, length(s_incrate15to49_ts)*sizeof(double));
    
    double *prev15to49_ts_out = REAL(s_prev15to49_ts);
    memset(prev15to49_ts_out, 0, length(s_prev15to49_ts)*sizeof(double));
    
    double *rvec = REAL(s_rvec_ts);
    
    double *prev15to49 = REAL(s_prev15to49);
    prev15to49[0] = 0.0;
    
    double *pregprev = REAL(s_pregprev);
    pregprev[0] = 0.0;
    
    double *incid15to49 = REAL(s_incid15to49);
    memset(incid15to49, 0, length(s_incid15to49)*sizeof(double));
    
    double *entrantprev_out = REAL(s_entrantprev_out);
    memset(entrantprev_out, 0, length(s_entrantprev_out)*sizeof(double));

    

    // array to store lagged prevalence among pregnant women
    double *pregprevlag = REAL(s_pregprevlag); // (double*) R_alloc(PROJ_YEARS, sizeof(double));
    memset(pregprevlag, 0, AGE_START*sizeof(double));


    double prevcurr = 0.0, prevlast; // store prevalence at last time step for r-trend model

    int everARTelig_idx = hDS;

    ////////////////////////////////////
    ////  do population projection  ////
    ////////////////////////////////////

    for(int t = 1; t < SIM_YEARS; t++){

      // age the population one year
      for(int m = 0; m < pDS; m++)
        for(int g = 0; g < NG; g++){
          for(int a = 1; a < pAG; a++)
            pop_t(a, g, m, t) = pop_t(a-1, g, m, t-1);
	  pop_t(pAG-1, g, m, t) += pop_t(pAG-1, g, m, t-1); // open age group
        }

      TensorFixedSize<double, Sizes<hAG, NG>> hiv_ag_prob;
      // Tensor<double, 2> hiv_ag_prob(hAG, NG);      
      hiv_ag_prob.setZero();
      
      for(int g = 0; g < NG; g++){
        int a = 0;
        for(int ha = 0; ha < (hAG-1); ha++){
          for(int i = 0; i < hAG_SPAN[ha]; i++){
            hiv_ag_prob(ha, g) += pop_t(a, g, HIVP, t-1);
            a++;
          }
          hiv_ag_prob(ha, g) = (hiv_ag_prob(ha, g) > 0) ? pop_t(a-1, g, HIVP, t-1) / hiv_ag_prob(ha, g) : 0;
        }
	// Note: loop stops at hAG-1; no one ages out of the open-ended age group
      }

      for(int g = 0; g < NG; g++)
        for(int ha = 1; ha < hAG; ha++)
          for(int hm = 0; hm < hDS; hm++){
            hivpop_t(hm, ha, g, t) = (1-hiv_ag_prob(ha, g)) * hivpop_t(hm, ha, g, t-1) + hiv_ag_prob(ha-1, g)*hivpop_t(hm, ha-1, g, t-1);
            if(t > t_ART_start)
              for(int hu = 0; hu < hTS; hu++)
                artpop_t(hu, hm, ha, g, t) = (1-hiv_ag_prob(ha, g)) * artpop_t(hu, hm, ha, g, t-1) + hiv_ag_prob(ha-1, g)*artpop_t(hu, hm, ha-1, g, t-1);
          }

      // add lagged births to youngest age group
      for(int g = 0; g < NG; g++){

        double paedsurv_g;
        double entrant_prev;
	
	if(use_entrantprev)
	  entrant_prev = entrantprev(g, t);
	else
	  entrant_prev = pregprevlag[t-1] * verttrans_lag[t-1] * paedsurv_lag[t-1];

        if(bin_popadjust){
          pop_t(0, g, HIVN, t) =  entrantpop(g, t-1) * (1.0-entrant_prev);
          paedsurv_g = entrantpop(g, t-1) * entrant_prev;
        } else {
          pop_t(0, g, HIVN, t) = birthslag(g, t-1) * cumsurv(g, t-1) * (1.0-entrant_prev / paedsurv_lag[t-1]) + cumnetmigr(g, t-1) * (1.0-pregprevlag[t-1] * netmig_hivprob);
          paedsurv_g = birthslag(g, t-1) * cumsurv(g, t-1) * entrant_prev + cumnetmigr(g, t-1) * entrant_prev;
        }

        pop_t(0, g, HIVP, t) = paedsurv_g;

        entrantprev_out[t] = (pop_t(0, MALE, HIVP, t) + pop_t(0, FEMALE, HIVP, t)) / (pop_t(0, MALE, HIVN, t) + pop_t(0, FEMALE, HIVN, t) + pop_t(0, MALE, HIVP, t) + pop_t(0, FEMALE, HIVP, t));

        for(int hm = 0; hm < hDS; hm++){
          hivpop_t(hm, 0, g, t) = (1-hiv_ag_prob(0, g)) * hivpop_t(hm, 0, g, t-1) + paedsurv_g * paedsurv_cd4dist(hm, g, t) * (1.0 - entrantartcov(g, t));
          if(t > t_ART_start){
            for(int hu = 0; hu < hTS; hu++){
              artpop_t(hu, hm, 0, g, t) = (1-hiv_ag_prob(0, g)) * artpop_t(hu, hm, 0, g, t-1);
	      artpop_t(hu, hm, 0, g, t) += paedsurv_g * paedsurv_artcd4dist(hu, hm, g, t) * entrantartcov(g, t);
	    }
	  }
        }
      }


      // non-HIV mortality and netmigration
      for(int g = 0; g < NG; g++){
        int a = 0;
        for(int ha = 0; ha < hAG; ha++){
          double deathsmig_ha = 0, hivpop_ha = 0;
          for(int i = 0; i < hAG_SPAN[ha]; i++){

            hivpop_ha += pop_t(a, g, HIVP, t);

            // non-HIV mortality
            double qx = 1.0 - Sx(a, g, t);
            double ndeaths_a = pop_t(a, g, HIVN, t) * qx;
            pop_t(a, g, HIVN, t) -= ndeaths_a; // survival HIV- population
            double hdeaths_a = pop_t(a, g, HIVP, t) * qx;
            deathsmig_ha -= hdeaths_a;
            pop_t(a, g, HIVP, t) -= hdeaths_a;   // survival HIV+ population
            natdeaths(a, g, t) = ndeaths_a + hdeaths_a;

            // net migration
            double migrate_a = netmigr(a, g, t) * (1+Sx(a, g, t))/2.0 / (pop_t(a, g, HIVN, t) + pop_t(a, g, HIVP, t));
            pop_t(a, g, HIVN, t) *= 1+migrate_a;
            double hmig_a = migrate_a * pop_t(a, g, HIVP, t);
            deathsmig_ha += hmig_a;
            pop_t(a, g, HIVP, t) += hmig_a;

            a++;
          }

          // migration and deaths for hivpop
          double deathmigrate_ha = hivpop_ha > 0 ? deathsmig_ha / hivpop_ha : 0.0;
          for(int hm = 0; hm < hDS; hm++){
            hivpop_t(hm, ha, g, t) *= 1+deathmigrate_ha;
            if(t > t_ART_start)
              for(int hu = 0; hu < hTS; hu++)
                artpop_t(hu, hm, ha, g , t) *= 1+deathmigrate_ha;
          } // loop over hm
        } // loop over ha
      } // loop over g


      // fertility
      double births = 0.0, births_by_ha[hAG_FERT];
      memset(births_by_ha, 0, hAG_FERT*sizeof(double));
      for(int m = 0; m < pDS; m++){
        int a = pIDX_FERT;
        for(int ha = hIDX_FERT; ha < hIDX_FERT+hAG_FERT; ha++){
          for(int i = 0; i < hAG_SPAN[ha]; i++){
            births_by_ha[ha-hIDX_FERT] += (pop_t(a, FEMALE, m, t-1) + pop_t(a, FEMALE, m, t))/2 * asfr(a, t);
            a++;
          }
        }
      }
      for(int ha = hIDX_FERT; ha < hAG_FERT; ha++)
        births += births_by_ha[ha-hIDX_FERT];

      if(t + AGE_START < PROJ_YEARS)
        for(int g = 0; g < NG; g++)
          birthslag(g, t + AGE_START-1) = srb(g, t) * births;

      ////////////////////////////////
      ////  HIV model simulation  ////
      ////////////////////////////////

      int cd4elig_idx = artcd4elig_idx[t] - 1; // -1 for 0-based indexing vs. 1-based in R
      int anyelig_idx = (specpop_percelig[t] > 0 | pw_artelig[t] > 0) ? 0 : (who34percelig > 0) ? hIDX_CD4_350 : cd4elig_idx;
      everARTelig_idx = anyelig_idx < everARTelig_idx ? anyelig_idx : everARTelig_idx;
      
      for(int hts = 0; hts < HIVSTEPS_PER_YEAR; hts++){

        int ts = (t-1)*HIVSTEPS_PER_YEAR + hts;

	TensorFixedSize<double, Sizes<hAG, NG>> hivdeaths_ha;
	hivdeaths_ha.setZero();

        // untreated population

        // disease progression and mortality
	TensorFixedSize<double, Sizes<hDS, hAG, NG>> grad_t;
	// Tensor<double, 3> grad_t(hDS, hAG, NG);
	
        for(int g = 0; g < NG; g++)
          for(int ha = 0; ha < hAG; ha++){
            for(int hm = 0; hm < hDS; hm++){

	      double cd4mx_scale = 1.0;
	      if(scale_cd4_mort & t >= t_ART_start & hm >= everARTelig_idx){
		double artpop_hahm = 0.0;
		for(int hu = 0; hu < hTS; hu++)
		  artpop_hahm += artpop_t(hu, hm, ha, g, t);
		cd4mx_scale = hivpop_t(hm, ha, g, t) / (hivpop_t(hm, ha, g, t) + artpop_hahm);
	      }
	      
              double deaths = cd4mx_scale * cd4_mort(hm, ha, g) * hivpop_t(hm, ha, g, t);
              hivdeaths_ha(ha, g) += DT*deaths;
	      aidsdeaths_noart(hm, ha, g, t) += DT*deaths;
	      grad_t(hm, ha, g) = -deaths;
            }
            for(int hm = 1; hm < hDS; hm++){
              grad_t(hm-1, ha, g) -= cd4_prog(hm-1, ha, g) * hivpop_t(hm-1, ha, g, t);
              grad_t(hm, ha, g) += cd4_prog(hm-1, ha, g) * hivpop_t(hm-1, ha, g, t);
            }
          }

        if(eppmod != EPP_DIRECTINCID){

          // incidence

          // calculate r(t)
          if(eppmod == EPP_RSPLINE) {
            rvec[ts] = rspline_rvec[ts];
          } else {
	    Rprintf("R-trend model not implemented.\n");
            // rvec[ts] = calc_rtrend_rt(pop, rtrend_tstab, rtrend_beta, rtrend_r0,
            //                           projsteps[ts], tsEpidemicStart, DT, t, hts,
            //                           rvec[ts-1], &prevlast, &prevcurr);
	  }

          // calculate new infections by sex and age
	  TensorFixedSize<double, Sizes<pAG, NG>> infections_ts;
          if(incidmod == INCIDMOD_EPPSPEC)
            calc_infections_eppspectrum_tensor(pop_t,
					       hivpop_t,
					       artpop_t,
                                        rvec[ts], relinfectART, (projsteps[ts] == tsEpidemicStart) ? iota : 0.0,
                                        incrr_sex,
					       incrr_age, 
					t_ART_start, DT, t, hts, hAG_START, hAG_SPAN,
					       &prevcurr, &incrate15to49_ts_out[ts], infections_ts);

          prev15to49_ts_out[ts] = prevcurr;

          // add new infections to HIV population
          for(int g = 0; g < NG; g++){
            int a = 0;
            for(int ha = 0; ha < hAG; ha++){
              double infections_a, infections_ha = 0.0;
              for(int i = 0; i < hAG_SPAN[ha]; i++){
                infections_ha += infections_a = infections_ts(a, g);
                infections(a, g, t) += DT*infections_a;
                pop_t(a, g, HIVN, t) -= DT*infections_a;
                pop_t(a, g, HIVP, t) += DT*infections_a;
                a++;
              }
              if(ha < hIDX_15TO49+hAG_15TO49 )
                incid15to49[t] += DT*infections_ha;

              // add infections to grad hivpop
              for(int hm = 0; hm < hDS; hm++) {
                grad_t(hm, ha, g) += infections_ha * cd4_initdist(hm, ha, g);
	      }
            }
          }
        }

        // ART progression, mortality, and initiation
        if(t >= t_ART_start){

	  TensorFixedSize<double, Sizes<hTS, hDS, hAG, NG>> gradART_t;
	  // Tensor<double, 4> gradART_t(hTS, hDS, hAG, NG);
	  
          // progression and mortality
          for(int g = 0; g < NG; g++)
            for(int ha = 0; ha < hAG; ha++)
              for(int hm = everARTelig_idx; hm < hDS; hm++){

                for(int hu = 0; hu < hTS; hu++){
                  double deaths = art_mort(hu, hm, ha, g) * artmx_timerr(hu, t) * artpop_t(hu, hm, ha, g, t);
                  hivdeaths_ha(ha, g) += DT*deaths;
		  aidsdeaths_art(hu, hm, ha, g, t) += DT*deaths;
                  gradART_t(hu, hm, ha, g) = -deaths;
                }

		for(int hu = 0; hu < (hTS - 1); hu++) {
		  gradART_t(hu, hm, ha, g) += -artpop_t(hu, hm, ha, g, t) / h_art_stage_dur[hu];
		  gradART_t(hu+1, hm, ha, g) += artpop_t(hu, hm, ha, g, t) / h_art_stage_dur[hu];
		}

		// ART dropout
		if(art_dropout[t] > 0)
		  for(int hu = 0; hu < hTS; hu++){
		    grad_t(hm, ha, g) += art_dropout[t] * artpop_t(hu, hm, ha, g, t);
                    gradART_t(hu, hm, ha, g) -= art_dropout[t] * artpop_t(hu, hm, ha, g, t);
		  }

              }


          // ART initiation
          for(int g = 0; g < NG; g++){

	    TensorFixedSize<double, Sizes<hDS, hAG_15PLUS>> artelig_hahm;
	      
	    double Xart_15plus = 0.0, Xartelig_15plus = 0.0, expect_mort_artelig15plus = 0.0;
            for(int ha = hIDX_15PLUS; ha < hAG; ha++){
              for(int hm = everARTelig_idx; hm < hDS; hm++){
		if(hm >= anyelig_idx){
		  double prop_elig = (hm >= cd4elig_idx) ? 1.0 : (hm >= hIDX_CD4_350) ? 1.0 - (1.0-specpop_percelig[t])*(1.0-who34percelig) : specpop_percelig[t];
		  Xartelig_15plus += artelig_hahm(hm, ha-hIDX_15PLUS) = prop_elig * hivpop_t(hm, ha, g, t);
		  expect_mort_artelig15plus += cd4_mort(hm, ha, g) * artelig_hahm(hm, ha-hIDX_15PLUS);
		}
                for(int hu = 0; hu < hTS; hu++)
                  Xart_15plus += artpop_t(hu, hm, ha, g, t) + DT * gradART_t(hu, hm, ha, g);
              }

              // if pw_artelig, add pregnant women to artelig_hahm population
              if(g == FEMALE & pw_artelig[t] > 0 & ha < hAG_FERT){
                double frr_pop_ha = 0;
                for(int a =  hAG_START[ha]; a < hAG_START[ha]+hAG_SPAN[ha]; a++)
                  frr_pop_ha += pop_t(a, g, HIVN, t); // add HIV- population
                for(int hm = 0; hm < hDS; hm++){
                  frr_pop_ha += frr_cd4(hm, ha-hIDX_FERT, t) * hivpop_t(hm, ha, g, t);
                  for(int hu = 0; hu < hTS; hu++)
                    frr_pop_ha += frr_art(hu, hm, ha-hIDX_FERT, t) * artpop_t(hu, hm, ha, g, t);
                }
                for(int hm = anyelig_idx; hm < cd4elig_idx; hm++){
                  double pw_elig_hahm = births_by_ha[ha-hIDX_FERT] * frr_cd4(hm, ha-hIDX_FERT, t) * hivpop_t(hm, ha, g, t) / frr_pop_ha;
                  artelig_hahm(hm, ha-hIDX_15PLUS) += pw_elig_hahm;
                  Xartelig_15plus += pw_elig_hahm;
                  expect_mort_artelig15plus += cd4_mort(hm, ha, g) * pw_elig_hahm;
                }
              }
            } // loop over ha

            // calculate number on ART at end of ts, based on number or percent
            double artnum_hts = 0.0;
            if(DT*(hts+1) < 0.5){
              if(!art15plus_isperc(g, t-2) & !art15plus_isperc(g, t-1)){ // both numbers
                artnum_hts = (0.5-DT*(hts+1))*artnum15plus(g, t-2) + (DT*(hts+1)+0.5)*artnum15plus(g, t-1);
              } else if(art15plus_isperc(g, t-2) & art15plus_isperc(g, t-1)){ // both percentages
                double artcov_hts = (0.5-DT*(hts+1))*artnum15plus(g, t-2) + (DT*(hts+1)+0.5)*artnum15plus(g, t-1);
                artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
              } else if(!art15plus_isperc(g, t-2) & art15plus_isperc(g, t-1)){ // transition from number to percentage
                double curr_coverage = Xart_15plus / (Xart_15plus + Xartelig_15plus);
                double artcov_hts = curr_coverage + (artnum15plus(g, t-1) - curr_coverage) * DT / (0.5-DT*hts);
                artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
              }
            } else {
              if(!art15plus_isperc(g, t-1) & !art15plus_isperc(g, t)){ // both numbers
                artnum_hts = (1.5-DT*(hts+1))*artnum15plus(g, t-1) + (DT*(hts+1)-0.5)*artnum15plus(g, t);
              } else if(art15plus_isperc(g, t-1) & art15plus_isperc(g, t)){ // both percentages
                double artcov_hts = (1.5-DT*(hts+1))*artnum15plus(g, t-1) + (DT*(hts+1)-0.5)*artnum15plus(g, t);
                artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
              } else if(!art15plus_isperc(g, t-1) & art15plus_isperc(g, t)){ // transition from number to percentage
                double curr_coverage = Xart_15plus / (Xart_15plus + Xartelig_15plus);
                double artcov_hts = curr_coverage + (artnum15plus(g, t) - curr_coverage) * DT / (1.5-DT*hts);
                artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
              }
            }

            double artinit_hts = artnum_hts > Xart_15plus ? artnum_hts - Xart_15plus : 0;

            // median CD4 at initiation inputs
            if(med_cd4init_input[t]){

              const int CD4_LOW_LIM[hDS] = {500, 350, 250, 200, 100, 50, 0};
              const int CD4_UPP_LIM[hDS] = {1000, 500, 350, 250, 200, 100, 50};

              int medcd4_idx = med_cd4init_cat[t] - 1; // -1 for 0-based indexing vs. 1-based in R
              double medcat_propbelow = (median_cd4init[t] - CD4_LOW_LIM[medcd4_idx]) / (CD4_UPP_LIM[medcd4_idx] - CD4_LOW_LIM[medcd4_idx]);

              double elig_below = 0.0, elig_above = 0.0;
              for(int ha = hIDX_15PLUS; ha < hAG; ha++){
                for(int hm = anyelig_idx; hm < medcd4_idx; hm++)
                  elig_above += artelig_hahm(hm, ha-hIDX_15PLUS);
                elig_above += (1.0 - medcat_propbelow) * artelig_hahm(medcd4_idx, ha-hIDX_15PLUS);
                elig_below += medcat_propbelow * artelig_hahm(medcd4_idx, ha-hIDX_15PLUS);
                for(int hm = medcd4_idx+1; hm < hDS; hm++)
                  elig_below += artelig_hahm(hm, ha-hIDX_15PLUS);
              }

              double initprob_below = (elig_below > artinit_hts * 0.5) ? artinit_hts * 0.5 / elig_below : 1.0;
              double initprob_above = (elig_above > artinit_hts * 0.5) ? artinit_hts * 0.5 / elig_above : 1.0;
              double initprob_medcat = initprob_below * medcat_propbelow + initprob_above * (1.0-medcat_propbelow);

              for(int ha = hIDX_15PLUS; ha < hAG; ha++)
                for(int hm = anyelig_idx; hm < hDS; hm++){
                  double artinit_hahm;
                  if(hm < medcd4_idx)
                    artinit_hahm = artelig_hahm(hm, ha-hIDX_15PLUS) * initprob_above;
                  else if(hm == medcd4_idx)
                    artinit_hahm = artelig_hahm(hm, ha-hIDX_15PLUS) * initprob_medcat;
                  if(hm > medcd4_idx)
                    artinit_hahm = artelig_hahm(hm, ha-hIDX_15PLUS) * initprob_below;
                  if (artinit_hahm > hivpop_t(hm, ha, g, t) + DT * grad_t(hm, ha, g)) {
		    artinit_hahm = hivpop_t(hm, ha, g, t) + DT * grad_t(hm, ha, g);
		  }
		  grad_t(hm, ha, g) -= artinit_hahm / DT;
                  gradART_t(ART0MOS, hm, ha, g) += artinit_hahm / DT;
		  artinit(hm, ha, g, t) += artinit_hahm; 
                }

            } else if(art_alloc_method == 4) {  // lowest CD4 first

	      for(int hm = hDS-1; hm >= anyelig_idx; hm--){
		double artelig_hm = 0;
		for(int ha = hIDX_15PLUS; ha < hAG; ha++)
		  artelig_hm += artelig_hahm(hm, ha-hIDX_15PLUS);
		double init_prop = (artelig_hm == 0 | artinit_hts > artelig_hm) ? 1.0 : artinit_hts / artelig_hm;

		for(int ha = hIDX_15PLUS; ha < hAG; ha++){
		  double artinit_hahm = init_prop * artelig_hahm(hm, ha-hIDX_15PLUS);

		  if (artinit_hahm > hivpop_t(hm, ha, g, t) + DT * grad_t(hm, ha, g)) {
		    artinit_hahm = hivpop_t(hm, ha, g, t) + DT * grad_t(hm, ha, g);
		  }

                  grad_t(hm, ha, g) -= artinit_hahm / DT;
                  gradART_t(ART0MOS, hm, ha, g) += artinit_hahm / DT;
		  artinit(hm, ha, g, t) += artinit_hahm; 
		}
		if(init_prop < 1.0)
		  break;
		artinit_hts -= init_prop * artelig_hm;
	      }
	      
	    } else { // Use mixture of eligibility and expected mortality for initiation distribution

              for(int ha = hIDX_15PLUS; ha < hAG; ha++)
                for(int hm = anyelig_idx; hm < hDS; hm++){
                  double artinit_hahm = artinit_hts * artelig_hahm(hm, ha-hIDX_15PLUS) * ((1.0 - art_alloc_mxweight)/Xartelig_15plus + art_alloc_mxweight * cd4_mort(hm, ha, g) / expect_mort_artelig15plus);
                  if (artinit_hahm > artelig_hahm(hm, ha-hIDX_15PLUS)) {
		    artinit_hahm = artelig_hahm(hm, ha-hIDX_15PLUS);
		  } if (artinit_hahm > hivpop_t(hm, ha, g, t) + DT * grad_t(hm, ha, g)) {
		    artinit_hahm = hivpop_t(hm, ha, g, t) + DT * grad_t(hm, ha, g);
		  }
                  grad_t(hm, ha, g) -= artinit_hahm / DT;
                  gradART_t(ART0MOS, hm, ha, g) += artinit_hahm / DT;
		  artinit(hm, ha, g, t) += artinit_hahm; 
                }
            }
          }

	  for(int g = 0; g < NG; g++)
	    for(int ha = 0; ha < hAG; ha++)
	      for(int hm = everARTelig_idx; hm < hDS; hm++)
		for(int hu = 0; hu < hTS; hu++)
		  artpop_t(hu, hm, ha, g, t) += DT*gradART_t(hu, hm, ha, g);
	  
	} // if(t >= t_ART_start)

	for(int g = 0; g < NG; g++)
          for(int ha = 0; ha < hAG; ha++)
            for(int hm = 0; hm < hDS; hm++)
              hivpop_t(hm, ha, g, t) += DT*grad_t(hm, ha, g);


	// remove hivdeaths from pop
	for(int g = 0; g < NG; g++){
	  
	  // sum HIV+ population size in each hivpop age group
	  double hivpop_ha[hAG];
	  int a = 0;
	  for(int ha = 0; ha < hAG; ha++){
	    hivpop_ha[ha] = 0.0;
	    for(int i = 0; i < hAG_SPAN[ha]; i++){
	      hivpop_ha[ha] += pop_t(a, g, HIVP, t);
	      a++;
	    }
	  }
	  
	  // remove hivdeaths proportionally to age-distribution within each age group
	  a = 0;
	  for(int ha = 0; ha < hAG; ha++){
	    if(hivpop_ha[ha] > 0){
	      double hivqx_ha = hivdeaths_ha(ha, g) / hivpop_ha[ha];
	      for(int i = 0; i < hAG_SPAN[ha]; i++){
		hivdeaths(a, g, t) += pop_t(a, g, HIVP, t) * hivqx_ha;
		pop_t(a, g, HIVP, t) *= (1.0-hivqx_ha);
		a++;
	      }
	    } else {
	      a += hAG_SPAN[ha];
	    }  // end if(pop_ha[ha] > 0)
	  }
	}



      } // loop HIVSTEPS_PER_YEAR


      
      if(eppmod == EPP_DIRECTINCID){
	// Calculating new infections once per year (like Spectrum)
	
	double Xhivp = 0.0, Xhivn[NG], Xhivn_incagerr[NG];
	
	for(int g = 0; g < NG; g++){
	  Xhivn[g] = 0.0;
	  Xhivn_incagerr[g] = 0.0;
	  for(int a = pIDX_INCIDPOP; a < pIDX_INCIDPOP+pAG_INCIDPOP; a++){
	    Xhivp += pop_t(a, g, HIVP, t-1);
	    Xhivn[g] += pop_t(a, g, HIVN, t-1);
	    Xhivn_incagerr[g] += incrr_age(a, g, t) * pop_t(a, g, HIVN, t-1);
	  }
	}
	// double prev_i = Xhivp / (Xhivn[MALE] + Xhivn[FEMALE] + Xhivp);
	// double incrate15to49_i = (prev15to49[t] - prev_i)/(1.0 - prev_i);
	double incrate_i = incidinput[t];
	double incrate_g[NG];
	incrate_g[MALE] = incrate_i * (Xhivn[MALE]+Xhivn[FEMALE]) / (Xhivn[MALE] + incrr_sex[t]*Xhivn[FEMALE]);
	incrate_g[FEMALE] = incrate_i * incrr_sex[t]*(Xhivn[MALE]+Xhivn[FEMALE]) / (Xhivn[MALE] + incrr_sex[t]*Xhivn[FEMALE]);
	
	for(int g = 0; g < NG; g++){
	  int a = 0;
	  for(int ha = 0; ha < hAG; ha++){
	    double infections_a, infections_ha = 0.0;
	    for(int i = 0; i < hAG_SPAN[ha]; i++){
	      infections_ha += infections_a = pop_t(a, g, HIVN, t-1) * incrate_g[g] * incrr_age(a, g, t) * Xhivn[g] / Xhivn_incagerr[g];
	      infections(a, g, t) += infections_a;
	      pop_t(a, g, HIVN, t) -= infections_a;
	      pop_t(a, g, HIVP, t) += infections_a;
	      a++;
	    }
	    if(ha < hIDX_15TO49+hAG_15TO49)
	      incid15to49[t] += infections_ha;
	    
	    // add infections to hivpop
	    for(int hm = 0; hm < hDS; hm++)
	      hivpop_t(hm, ha, g, t) += infections_ha * cd4_initdist(hm, ha, g);
	  }
	}
      }

      // adjust population to match target population
      if(bin_popadjust){
        for(int g = 0; g < NG; g++){
          int a = 0;
          for(int ha = 0; ha < hAG; ha++){
            double popadj_ha = 0, hivpop_ha = 0;
            for(int i = 0; i < hAG_SPAN[ha]; i++){

              hivpop_ha += pop_t(a, g, HIVP, t);

              double popadjrate_a = popadjust(a, g, t) = targetpop(a, g, t) / (pop_t(a, g, HIVN, t) + pop_t(a, g, HIVP, t));
              pop_t(a, g, HIVN, t) *= popadjrate_a;
              double hpopadj_a = (popadjrate_a-1.0) * pop_t(a, g, HIVP, t);
              popadj_ha += hpopadj_a;
              pop_t(a, g, HIVP, t) += hpopadj_a;
              a++;
            }

            // population adjustment for hivpop
            double popadjrate_ha = hivpop_ha > 0 ? popadj_ha / hivpop_ha : 0.0;
            for(int hm = 0; hm < hDS; hm++){
              hivpop_t(hm, ha, g, t) *= 1+popadjrate_ha;
              if(t >= t_ART_start)
                for(int hu = 0; hu < hTS; hu++)
                  artpop_t(hu, hm, ha, g, t) *= 1+popadjrate_ha;
            } // loop over hm
          } // loop over ha
        } // loop over g
      } // if(bin_popadjust)



      // prevalence among pregnant women

      double hivbirths = 0;
      for(int ha = hIDX_FERT; ha < hIDX_FERT+hAG_FERT; ha++){
        double hivn_ha = 0, frr_hivpop_ha = 0;
        for(int a =  hAG_START[ha]; a < hAG_START[ha]+hAG_SPAN[ha]; a++)
          hivn_ha += (pop_t(a, FEMALE, HIVN, t-1) + pop_t(a, FEMALE, HIVN, t))/2;
        for(int hm = 0; hm < hDS; hm++){
          frr_hivpop_ha += frr_cd4(hm, ha-hIDX_FERT, t) * (hivpop_t(hm, ha, FEMALE, t-1) + hivpop_t(hm, ha, FEMALE, t))/2;
          if(t == t_ART_start)
            for(int hu = 0; hu < hTS; hu++)
              frr_hivpop_ha += frr_art(hu, hm, ha-hIDX_FERT, t) * artpop_t(hu, hm, ha, FEMALE, t)/2;
          else if(t > t_ART_start)
            for(int hu = 0; hu < hTS; hu++)
              frr_hivpop_ha += frr_art(hu, hm, ha-hIDX_FERT, t) * (artpop_t(hu, hm, ha, FEMALE, t-1) + artpop_t(hu, hm, ha, FEMALE, t))/2;
        }
        hivbirths += births_by_ha[ha-hIDX_FERT] * frr_hivpop_ha / (hivn_ha + frr_hivpop_ha);
      }

      pregprev[t] = hivbirths/births;
      if(t + AGE_START < PROJ_YEARS)
        pregprevlag[t + AGE_START-1] = pregprev[t];

      // prevalence 15 to 49
      for(int g = 0; g < NG; g++)
        for(int a = pIDX_15TO49; a < pIDX_15TO49 + pAG_15TO49; a++){
          hivn15to49[t] += pop_t(a, g, HIVN, t);
          hivp15to49[t] += pop_t(a, g, HIVP, t);
        }
      prev15to49[t] = hivp15to49[t]/(hivn15to49[t] + hivp15to49[t]);
      incid15to49[t] /= hivn15to49[t-1];
    }

    UNPROTECT(28);
    return s_pop;
  }
}


void calc_infections_eppspectrum_tensor(const TensorMap<Tensor<double, 4>> pop_t,
					const TensorMap<Tensor<double, 4>> hivpop_t,
					const TensorMap<Tensor<double, 5>> artpop_t,
					double r_ts,
					double relinfectART,
					double iota,
					double *incrr_sex,
					const TensorMap<Tensor<double, 3>> incrr_age,
					int t_ART_start,
					double DT,
					int t,
					int hts,
					int *hAG_START,
					int *hAG_SPAN,
					double *prevcurr,
					double *incrate15to49_ts,
					TensorFixedSize<double, Sizes<pAG, NG>>& infections_ts)
{

  // sum population sizes
  double Xhivn_g[NG], Xhivn_incagerr[NG], Xhivp_noart = 0.0, Xart = 0.0;
  for(int g = 0; g < NG; g++){
    Xhivn_g[g] = 0.0;
    Xhivn_incagerr[g] = 0.0;
    for(int a = pIDX_15TO49; a < pIDX_15TO49+pAG_15TO49; a++){
      Xhivn_g[g] += pop_t(a, g, HIVN, t);
      Xhivn_incagerr[g] += incrr_age(a, g, t) * pop_t(a, g, HIVN, t);
    }

    for(int ha = hIDX_15TO49; ha < hIDX_15TO49+hAG_15TO49+1; ha++){

      // adjustment to first and last age group for partial year time step
      // calculation proportion of HIV population to include / exclude based on hivpop in single-year ages.
      double prop_include;
      if(ha == hIDX_15TO49){
        double hivp_ha = 0.0;
        for(int a = hAG_START[ha]; a < hAG_START[ha]+hAG_SPAN[ha]; a++)
          hivp_ha += pop_t(a, g, HIVP, t);
        prop_include = (hivp_ha > 0) ? 1.0 - pop_t(hAG_START[ha], g, HIVP, t) / hivp_ha * (1.0 - DT*hts) : 1.0;
      } else if(ha == hIDX_15TO49+hAG_15TO49) {
        double hivp_ha = 0.0;
        for(int a = hAG_START[ha]; a < hAG_START[ha]+hAG_SPAN[ha]; a++)
          hivp_ha += pop_t(a, g, HIVP, t);
        prop_include = (hivp_ha > 0) ? pop_t(hAG_START[ha], g, HIVP, t) / hivp_ha * (1.0 - DT*hts) : 1.0;
      } else
        prop_include = 1.0;

      for(int hm = 0; hm < hDS; hm++){
        Xhivp_noart += hivpop_t(hm, ha, g, t) * prop_include;
        if(t >= t_ART_start)
          for(int hu = 0; hu < hTS; hu++)
            Xart += artpop_t(hu, hm, ha, g, t) * prop_include;
      }
    }
  } // end loop over g
  double Xhivn = Xhivn_g[MALE] + Xhivn_g[FEMALE];

  // adjust HIV negative population for partial year time step
  for(int g = 0; g < NG; g++){
    Xhivn -= pop_t(pIDX_15TO49, g, HIVN, t) * (1.0 - DT*hts);
    Xhivn += pop_t(pIDX_15TO49+pAG_15TO49, g, HIVN, t) * (1.0 - DT*hts);
  }

  double Xtot = Xhivn + Xhivp_noart + Xart;
  *prevcurr = (Xhivp_noart + Xart) / Xtot;

  *incrate15to49_ts = r_ts * (Xhivp_noart + relinfectART * Xart) / Xtot + iota;

  // incidence by sex
  double incrate15to49_g[NG];
  incrate15to49_g[MALE] = *incrate15to49_ts * (Xhivn_g[MALE]+Xhivn_g[FEMALE]) / (Xhivn_g[MALE] + incrr_sex[t]*Xhivn_g[FEMALE]);
  incrate15to49_g[FEMALE] = *incrate15to49_ts * incrr_sex[t]*(Xhivn_g[MALE]+Xhivn_g[FEMALE]) / (Xhivn_g[MALE] + incrr_sex[t]*Xhivn_g[FEMALE]);

  // annualized infections by age and sex
  for(int g = 0; g < NG; g++)
    for(int a = 0; a < pAG; a++){
      infections_ts(a, g) = pop_t(a, g, HIVN, t) * incrate15to49_g[g] * incrr_age(a, g, t) * Xhivn_g[g] / Xhivn_incagerr[g];
    }

  return;
}
