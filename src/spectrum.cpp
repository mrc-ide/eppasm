#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>

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
#define ART6MOS 1
#define ART1YR 2

#define ART_STAGE_PROG_RATE 2.0 // HARD CODED: ART stage progression rate


#define EPP_RSPLINE 0
#define EPP_RTREND 1


SEXP getListElement(SEXP list, const char *str);

extern "C" {

  SEXP checkBoostAsserts(){
#ifndef BOOST_DISABLE_ASSERTS
    Rprintf("BOOST ASSERTS ENABLED\n");
#endif
    return R_NilValue;
  }

  SEXP spectrumC(SEXP s_fp){

    ////////////////////////////////
    ////  set parameter values  ////
    ////////////////////////////////

    using namespace boost;

    // state space dimensions
    SEXP s_ss = getListElement(s_fp, "ss");
    int PROJ_YEARS = *INTEGER(getListElement(s_ss, "PROJ_YEARS"));
    int HIVSTEPS_PER_YEAR = *INTEGER(getListElement(s_ss, "hiv_steps_per_year"));
    double DT = 1.0/HIVSTEPS_PER_YEAR;
    int *hAG_SPAN = INTEGER(getListElement(s_ss, "h.ag.span"));

    int hAG_START[hAG];
    hAG_START[0] = 0;
    for(int ha = 1; ha < hAG; ha++)
      hAG_START[ha] = hAG_START[ha-1] + hAG_SPAN[ha-1];

    double *projsteps = REAL(getListElement(s_fp, "proj.steps"));

    // demographic projection
    multi_array_ref<double, 2> basepop(REAL(getListElement(s_fp, "basepop")), extents[NG][pAG]);
    multi_array_ref<double, 3> Sx(REAL(getListElement(s_fp, "Sx")), extents[PROJ_YEARS][NG][pAG]);
    multi_array_ref<double, 3> netmigr(REAL(getListElement(s_fp, "netmigr")), extents[PROJ_YEARS][NG][pAG]);
    multi_array_ref<double, 2> asfr(REAL(getListElement(s_fp, "asfr")), extents[PROJ_YEARS][pAG_FERT]);
    multi_array_ref<double, 2> srb(REAL(getListElement(s_fp, "srb")), extents[PROJ_YEARS][NG]);
    multi_array_ref<double, 2> birthslag(REAL(getListElement(s_fp, "birthslag")), extents[PROJ_YEARS][NG]);
    multi_array_ref<double, 2> cumsurv(REAL(getListElement(s_fp, "cumsurv")), extents[PROJ_YEARS][NG]);
    multi_array_ref<double, 2> cumnetmigr(REAL(getListElement(s_fp, "cumnetmigr")), extents[PROJ_YEARS][NG]);

    int bin_popadjust = *INTEGER(getListElement(s_fp, "popadjust"));
    double *ptr_targetpop;
    double *ptr_entrantpop;
    if(bin_popadjust){
      ptr_targetpop = REAL(getListElement(s_fp, "targetpop"));
      ptr_entrantpop = REAL(getListElement(s_fp, "entrantpop"));
    }
    multi_array_ref<double, 3> targetpop(ptr_targetpop, extents[PROJ_YEARS][NG][pAG]);
    multi_array_ref<double, 2> entrantpop(ptr_entrantpop, extents[PROJ_YEARS][NG]);

    // disease progression
    multi_array_ref<double, 3> cd4_initdist(REAL(getListElement(s_fp, "cd4_initdist")), extents[NG][hAG][hDS]);
    multi_array_ref<double, 3> cd4_prog(REAL(getListElement(s_fp, "cd4_prog")), extents[NG][hAG][hDS-1]);
    multi_array_ref<double, 3> cd4_mort(REAL(getListElement(s_fp, "cd4_mort")), extents[NG][hAG][hDS]);
    multi_array_ref<double, 4> art_mort(REAL(getListElement(s_fp, "art_mort")), extents[NG][hAG][hDS][hTS]);

    // sub-fertility
    multi_array_ref<double, 3> frr_cd4(REAL(getListElement(s_fp, "frr_cd4")), extents[PROJ_YEARS][hAG_FERT][hDS]);
    multi_array_ref<double, 4> frr_art(REAL(getListElement(s_fp, "frr_art")), extents[PROJ_YEARS][hAG_FERT][hDS][hTS]);

    // ART inputs
    int t_ART_start = *INTEGER(getListElement(s_fp, "tARTstart")) - 1; // -1 for 0-based indexing in C vs. 1-based in R
    multi_array_ref<double, 2> artnum15plus(REAL(getListElement(s_fp, "art15plus_num")), extents[PROJ_YEARS][NG]);
    multi_array_ref<int, 2> art15plus_isperc(LOGICAL(getListElement(s_fp, "art15plus_isperc")), extents[PROJ_YEARS][NG]);

    int *artcd4elig_idx = INTEGER(getListElement(s_fp, "artcd4elig_idx"));  // NOTE: 1-based indexing
    double *specpop_percelig = REAL(getListElement(s_fp, "specpop_percelig"));
    double *pw_artelig = REAL(getListElement(s_fp, "pw_artelig"));
    double who34percelig = *REAL(getListElement(s_fp, "who34percelig"));

    double *art_dropout = REAL(getListElement(s_fp, "art_dropout"));
    double *median_cd4init = REAL(getListElement(s_fp, "median_cd4init"));

    int *med_cd4init_cat = INTEGER(getListElement(s_fp, "med_cd4init_cat"));
    int *med_cd4init_input = INTEGER(getListElement(s_fp, "med_cd4init_input"));
      

    // incidence model
    // double *prev15to49 = REAL(getListElement(s_fp, "prev15to49"));
    double *incrr_sex = REAL(getListElement(s_fp, "incrr_sex"));
    multi_array_ref<double, 3> incrr_age(REAL(getListElement(s_fp, "incrr_age")), extents[PROJ_YEARS][NG][pAG]);

    double relinfectART = *REAL(getListElement(s_fp, "relinfectART"));
    // double ts_epidemic_start = *INTEGER(getListElement(s_fp, "ts_epi_start")) - 1; // -1 for 0-based indexing in C vs. 1-based in R
    double tsEpidemicStart = *REAL(getListElement(s_fp, "tsEpidemicStart")); // -1 for 0-based indexing in C vs. 1-based in R
    double iota = *REAL(getListElement(s_fp, "iota"));

    int eppmod = *INTEGER(getListElement(s_fp, "eppmodInt"));
    double *rspline_rvec;
    double *rtrend_beta, rtrend_tstab, rtrend_r0;
    if(eppmod == EPP_RSPLINE)
      rspline_rvec = REAL(getListElement(s_fp, "rvec"));
    else {
      SEXP s_rtrend = getListElement(s_fp, "rtrend");
      rtrend_beta = REAL(getListElement(s_rtrend, "beta"));
      rtrend_tstab = *REAL(getListElement(s_rtrend, "tStabilize"));
      rtrend_r0 = *REAL(getListElement(s_rtrend, "r0"));
    }
    

    // vertical transmission and survival
    double verttrans = *REAL(getListElement(s_fp, "verttrans"));
    double paedsurv = *REAL(getListElement(s_fp, "paedsurv"));
    double netmig_hivprob = *REAL(getListElement(s_fp, "netmig_hivprob"));
    double netmighivsurv = *REAL(getListElement(s_fp, "netmighivsurv"));
    double *paedsurv_cd4dist = REAL(getListElement(s_fp, "paedsurv_cd4dist"));


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
    multi_array_ref<double, 3> infections(REAL(s_infections), extents[PROJ_YEARS][NG][pAG]);
    memset(REAL(s_infections), 0, length(s_infections)*sizeof(double));
    
    SEXP s_hivdeaths = PROTECT(allocVector(REALSXP, pAG * NG * PROJ_YEARS));
    SEXP s_hivdeaths_dim = PROTECT(allocVector(INTSXP, 3));
    INTEGER(s_hivdeaths_dim)[0] = pAG;
    INTEGER(s_hivdeaths_dim)[1] = NG;
    INTEGER(s_hivdeaths_dim)[2] = PROJ_YEARS;
    setAttrib(s_hivdeaths, R_DimSymbol, s_hivdeaths_dim);
    setAttrib(s_pop, install("hivdeaths"), s_hivdeaths);
    multi_array_ref<double, 3> hivdeaths(REAL(s_hivdeaths), extents[PROJ_YEARS][NG][pAG]);
    memset(REAL(s_hivdeaths), 0, length(s_hivdeaths)*sizeof(double));

    SEXP s_natdeaths = PROTECT(allocVector(REALSXP, pAG * NG * PROJ_YEARS));
    SEXP s_natdeaths_dim = PROTECT(allocVector(INTSXP, 3));
    INTEGER(s_natdeaths_dim)[0] = pAG;
    INTEGER(s_natdeaths_dim)[1] = NG;
    INTEGER(s_natdeaths_dim)[2] = PROJ_YEARS;
    setAttrib(s_natdeaths, R_DimSymbol, s_natdeaths_dim);
    setAttrib(s_pop, install("natdeaths"), s_natdeaths);
    multi_array_ref<double, 3> natdeaths(REAL(s_natdeaths), extents[PROJ_YEARS][NG][pAG]);
    memset(REAL(s_natdeaths), 0, length(s_natdeaths)*sizeof(double));

    SEXP s_popadjust = PROTECT(allocVector(REALSXP, pAG * NG * PROJ_YEARS));
    SEXP s_popadjust_dim = PROTECT(allocVector(INTSXP, 3));
    INTEGER(s_popadjust_dim)[0] = pAG;
    INTEGER(s_popadjust_dim)[1] = NG;
    INTEGER(s_popadjust_dim)[2] = PROJ_YEARS;
    setAttrib(s_popadjust, R_DimSymbol, s_popadjust_dim);
    setAttrib(s_pop, install("popadjust"), s_popadjust);
    multi_array_ref<double, 3> popadjust(REAL(s_popadjust), extents[PROJ_YEARS][NG][pAG]);
    memset(REAL(s_popadjust), 0, length(s_popadjust)*sizeof(double));

    SEXP s_pregprevlag = PROTECT(allocVector(REALSXP, PROJ_YEARS));
    setAttrib(s_pop, install("pregprevlag"), s_pregprevlag);

    SEXP s_incrate15to49_ts = PROTECT(allocVector(REALSXP, (PROJ_YEARS-1) * HIVSTEPS_PER_YEAR));
    setAttrib(s_pop, install("incrate15to49_ts"), s_incrate15to49_ts);
    double *incrate15to49_ts_out = REAL(s_incrate15to49_ts);
    memset(incrate15to49_ts_out, 0, length(s_incrate15to49_ts)*sizeof(double));

    SEXP s_prev15to49_ts = PROTECT(allocVector(REALSXP, (PROJ_YEARS-1) * HIVSTEPS_PER_YEAR));
    setAttrib(s_pop, install("prev15to49_ts"), s_prev15to49_ts);
    double *prev15to49_ts_out = REAL(s_prev15to49_ts);
    memset(prev15to49_ts_out, 0, length(s_prev15to49_ts)*sizeof(double));

    SEXP s_rvec_ts = PROTECT(allocVector(REALSXP, (PROJ_YEARS-1) * HIVSTEPS_PER_YEAR));
    setAttrib(s_pop, install("rvec_ts"), s_rvec_ts);
    double *rvec = REAL(s_rvec_ts);

    SEXP s_prev15to49 = PROTECT(allocVector(REALSXP, PROJ_YEARS));
    setAttrib(s_pop, install("prev15to49"), s_prev15to49);
    double *prev15to49 = REAL(s_prev15to49);
    prev15to49[0] = 0.0;

    SEXP s_pregprev = PROTECT(allocVector(REALSXP, PROJ_YEARS));
    setAttrib(s_pop, install("pregprev"), s_pregprev);
    double *pregprev = REAL(s_pregprev);
    pregprev[0] = 0.0;

    SEXP s_incid15to49 = PROTECT(allocVector(REALSXP, PROJ_YEARS));
    setAttrib(s_pop, install("incid15to49"), s_incid15to49);
    double *incid15to49 = REAL(s_incid15to49);
    memset(incid15to49, 0, length(s_incid15to49)*sizeof(double));

    double *hivn15to49 = (double*) R_alloc(PROJ_YEARS, sizeof(double));
    double *hivp15to49 = (double*) R_alloc(PROJ_YEARS, sizeof(double));
    memset(hivn15to49, 0, PROJ_YEARS*sizeof(double));
    memset(hivp15to49, 0, PROJ_YEARS*sizeof(double));

    
    // initialize population

    // population by single-year age
    // double pop[PROJ_YEARS][pDS][NG][pAG];
    multi_array_ref<double, 4> pop(REAL(s_pop), extents[PROJ_YEARS][pDS][NG][pAG]);
    for(int g = 0; g < NG; g++)
      for(int a = 0; a < pAG; a++){
        pop[0][HIVN][g][a] = basepop[g][a];
        pop[0][HIVP][g][a] = 0.0;
	if(a >= pIDX_15TO49 & a < pIDX_15TO49+pAG_15TO49)
	  hivn15to49[0] += basepop[g][a];
      }

    // HIV population with stage stratification
    // double hivpop[PROJ_YEARS][NG][hAG][hDS];
    multi_array_ref<double, 4> hivpop(REAL(s_hivpop), extents[PROJ_YEARS][NG][hAG][hDS]);
    for(int g = 0; g < NG; g++)
      for(int ha = 0; ha < hAG; ha++)
        for(int hm = 0; hm < hDS; hm++)
          hivpop[0][g][ha][hm] = 0.0;

    // ART population with stage stratification
    // double artpop[PROJ_YEARS][NG][hAG][hDS][hTS];
    multi_array_ref<double, 5> artpop(REAL(s_artpop), extents[PROJ_YEARS][NG][hAG][hDS][hTS]);
    // memset(REAL(s_artpop), 0, length(s_artpop) * sizeof(double)); // initialize artpop to 0
    if(t_ART_start < PROJ_YEARS)
      for(int g = 0; g < NG; g++)
	for(int ha = 0; ha < hAG; ha++)
	  for(int hm = 0; hm < hDS; hm++)
	    for(int hu = 0; hu < hTS; hu++)
	      artpop[t_ART_start][g][ha][hm][hu] = 0.0;  // initialize to zero in year of ART start

          
    // array to store lagged prevalence among pregnant women
    double *pregprevlag = REAL(s_pregprevlag); // (double*) R_alloc(PROJ_YEARS, sizeof(double));
    memset(pregprevlag, 0, AGE_START*sizeof(double));


    double prevcurr = 0.0, prevlast; // store prevalence at last time step for r-trend model

    ////////////////////////////////////
    ////  do population projection  ////
    ////////////////////////////////////

    for(int t = 1; t < PROJ_YEARS; t++){

      // age the population one year
      for(int m = 0; m < pDS; m++)
        for(int g = 0; g < NG; g++){
          for(int a = 1; a < pAG; a++)
            pop[t][m][g][a] = pop[t-1][m][g][a-1];
          pop[t][m][g][pAG-1] += pop[t-1][m][g][pAG-1]; // open age group
        }

      double hiv_ag_prob[NG][hAG];
      for(int g = 0; g < NG; g++){
        int a = 0;
        for(int ha = 0; ha < (hAG-1); ha++){
          hiv_ag_prob[g][ha] = 0;
          for(int i = 0; i < hAG_SPAN[ha]; i++){
            hiv_ag_prob[g][ha] += pop[t-1][HIVP][g][a];
            a++;
          }
          hiv_ag_prob[g][ha] = (hiv_ag_prob[g][ha] > 0) ? pop[t-1][HIVP][g][a-1] / hiv_ag_prob[g][ha] : 0;
        }
        hiv_ag_prob[g][hAG-1] = 0.0; // no one ages out of the open-ended age group
      }

      for(int g = 0; g < NG; g++)
        for(int ha = 1; ha < hAG; ha++)
          for(int hm = 0; hm < hDS; hm++){
            hivpop[t][g][ha][hm] = (1-hiv_ag_prob[g][ha]) * hivpop[t-1][g][ha][hm] + hiv_ag_prob[g][ha-1]*hivpop[t-1][g][ha-1][hm];
            if(t > t_ART_start)
              for(int hu = 0; hu < hTS; hu++)
                artpop[t][g][ha][hm][hu] = (1-hiv_ag_prob[g][ha]) * artpop[t-1][g][ha][hm][hu] + hiv_ag_prob[g][ha-1]*artpop[t-1][g][ha-1][hm][hu];
          }

      // add lagged births to youngest age group
      for(int g = 0; g < NG; g++){

	double paedsurv_g;
	if(bin_popadjust){
	  double entrant_prev = pregprevlag[t-1] * verttrans * paedsurv;
	  pop[t][HIVN][g][0] =  entrantpop[t-1][g] * (1.0-entrant_prev);
	  paedsurv_g = entrantpop[t-1][g] * entrant_prev;
	} else {
	  pop[t][HIVN][g][0] = birthslag[t-1][g] * cumsurv[t-1][g] * (1.0-pregprevlag[t-1] * verttrans) + cumnetmigr[t-1][g] * (1.0-pregprevlag[t-1] * netmig_hivprob);
	  paedsurv_g = birthslag[t-1][g] * cumsurv[t-1][g] * pregprevlag[t-1] * verttrans * paedsurv + cumnetmigr[t-1][g] * pregprevlag[t-1] * netmig_hivprob * netmighivsurv;
	}

	pop[t][HIVP][g][0] = paedsurv_g;

        for(int hm = 0; hm < hDS; hm++){
          hivpop[t][g][0][hm] = (1-hiv_ag_prob[g][0]) * hivpop[t-1][g][0][hm] + paedsurv_g * paedsurv_cd4dist[hm]; // !!! UPDATE FOR VERTICAL TRANSMISSION
          if(t > t_ART_start)
            for(int hu = 0; hu < hTS; hu++)
              artpop[t][g][0][hm][hu] = (1-hiv_ag_prob[g][0]) * artpop[t-1][g][0][hm][hu];
        }
      }


      // non-HIV mortality and netmigration
      for(int g = 0; g < NG; g++){
        int a = 0;
        for(int ha = 0; ha < hAG; ha++){
          double deathsmig_ha = 0, hivpop_ha = 0;
          for(int i = 0; i < hAG_SPAN[ha]; i++){

            hivpop_ha += pop[t][HIVP][g][a];

            // non-HIV mortality
	    double qx = 1.0 - Sx[t][g][a];
	    double ndeaths_a = pop[t][HIVN][g][a] * qx;
            pop[t][HIVN][g][a] -= ndeaths_a; // survival HIV- population
	    double hdeaths_a = pop[t][HIVP][g][a] * qx;
            deathsmig_ha -= hdeaths_a;
            pop[t][HIVP][g][a] -= hdeaths_a;   // survival HIV+ population
	    natdeaths[t][g][a] = ndeaths_a + hdeaths_a;

            // net migration
            double migrate_a = netmigr[t][g][a] * (1+Sx[t][g][a])/2.0 / (pop[t][HIVN][g][a] + pop[t][HIVP][g][a]);
            pop[t][HIVN][g][a] *= 1+migrate_a;
	    double hmig_a = migrate_a * pop[t][HIVP][g][a];
            deathsmig_ha += hmig_a;
            pop[t][HIVP][g][a] += hmig_a;

            a++;
          }

          // migration and deaths for hivpop
          double deathmigrate_ha = hivpop_ha > 0 ? deathsmig_ha / hivpop_ha : 0.0;
          for(int hm = 0; hm < hDS; hm++){
            hivpop[t][g][ha][hm] *= 1+deathmigrate_ha;
            if(t > t_ART_start)
              for(int hu = 0; hu < hTS; hu++)
                artpop[t][g][ha][hm][hu] *= 1+deathmigrate_ha;
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
            births_by_ha[ha-hIDX_FERT] += (pop[t-1][m][FEMALE][a] + pop[t][m][FEMALE][a])/2 * asfr[t][a];
            a++;
          }
        }
      }
      for(int ha = hIDX_FERT; ha < hAG_FERT; ha++)
        births += births_by_ha[ha-hIDX_FERT];
      
      if(t + AGE_START < PROJ_YEARS)
        for(int g = 0; g < NG; g++)
          birthslag[t + AGE_START-1][g] = srb[t][g] * births;

      
      ////////////////////////////////
      ////  HIV model simulation  ////
      ////////////////////////////////
      
      for(int hts = 0; hts < HIVSTEPS_PER_YEAR; hts++){

	int ts = (t-1)*HIVSTEPS_PER_YEAR + hts;
	
	double hivdeaths_ha[NG][hAG];
	memset(hivdeaths_ha, 0, sizeof(double)*NG*hAG);

	// untreated population

	// disease progression and mortality
	double grad[NG][hAG][hDS];
        for(int g = 0; g < NG; g++)
          for(int ha = 0; ha < hAG; ha++){
            for(int hm = 0; hm < hDS; hm++){
              double deaths = cd4_mort[g][ha][hm] * hivpop[t][g][ha][hm];
              hivdeaths_ha[g][ha] += DT*deaths;
              grad[g][ha][hm] = -deaths;
            }
            for(int hm = 1; hm < hDS; hm++){
              grad[g][ha][hm-1] -= cd4_prog[g][ha][hm-1] * hivpop[t][g][ha][hm-1];
              grad[g][ha][hm] += cd4_prog[g][ha][hm-1] * hivpop[t][g][ha][hm-1];
	    }
	  }
	
	// incidence

	// sum population sizes
	double Xhivn_g[NG], Xhivn_incagerr[NG], Xhivp_noart = 0.0, Xart = 0.0;
	for(int g = 0; g < NG; g++){
	  Xhivn_g[g] = 0.0;
	  Xhivn_incagerr[g] = 0.0;
	  for(int a = pIDX_15TO49; a < pIDX_15TO49+pAG_15TO49; a++){
	    Xhivn_g[g] += pop[t][HIVN][g][a];
	    Xhivn_incagerr[g] += incrr_age[t][g][a] * pop[t][HIVN][g][a];
	  }
	  for(int ha = hIDX_15TO49; ha < hIDX_15TO49+hAG_15TO49; ha++)
	    for(int hm = 0; hm < hDS; hm++){
	      Xhivp_noart += hivpop[t][g][ha][hm];
	      if(t >= t_ART_start)
		for(int hu = 0; hu < hTS; hu++)
		  Xart += artpop[t][g][ha][hm][hu];
	    }
	}
	double Xhivn = Xhivn_g[MALE] + Xhivn_g[FEMALE];
	double Xhivp = Xhivp_noart + Xart;
	double prop_art_ts = Xhivp > 0 ? Xart / Xhivp : 0.0;
	
	// adjust HIV population for partial year time step
	for(int g = 0; g < NG; g++){
	  Xhivn -= pop[t][HIVN][g][pIDX_15TO49] * (1.0 - DT*hts);
	  Xhivp -= pop[t][HIVP][g][pIDX_15TO49] * (1.0 - DT*hts);
	  Xhivn += pop[t][HIVN][g][pIDX_15TO49+pAG_15TO49] * (1.0 - DT*hts);
	  Xhivp += pop[t][HIVP][g][pIDX_15TO49+pAG_15TO49] * (1.0 - DT*hts);
	}
	    
	double Xtot = Xhivn + Xhivp;

	prevlast = prevcurr;
	prevcurr = Xhivp / Xtot;
	prev15to49_ts_out[ts] = prevcurr;

	// calculate r(t)
	if(eppmod == EPP_RSPLINE)
	  rvec[ts] = rspline_rvec[ts];
	else {
	  if(projsteps[ts] > tsEpidemicStart){
	    double gamma_ts = (projsteps[ts] < rtrend_tstab)?0.0:(prevcurr-prevlast) * (projsteps[ts] - rtrend_tstab) / (DT * prevlast);
	    double logr_diff = rtrend_beta[1]*(rtrend_beta[0] - rvec[ts-1]) + rtrend_beta[2]*prevlast + rtrend_beta[3]*gamma_ts;
	    rvec[ts] = exp(log(rvec[ts-1]) + logr_diff);
	  } else {
	    rvec[ts] = rtrend_r0;
	  }
	}

	double incrate15to49_ts = rvec[ts] * Xhivp * (1.0 - (1.0 - relinfectART) * prop_art_ts) / Xtot + ((projsteps[ts] == tsEpidemicStart) ? iota : 0.0);
	incrate15to49_ts_out[ts] = incrate15to49_ts;


	// incidence by sex
	double incrate15to49_g[NG];
	incrate15to49_g[MALE] = incrate15to49_ts * (Xhivn_g[MALE]+Xhivn_g[FEMALE]) / (Xhivn_g[MALE] + incrr_sex[t]*Xhivn_g[FEMALE]);
	incrate15to49_g[FEMALE] = incrate15to49_ts * incrr_sex[t]*(Xhivn_g[MALE]+Xhivn_g[FEMALE]) / (Xhivn_g[MALE] + incrr_sex[t]*Xhivn_g[FEMALE]);

	for(int g = 0; g < NG; g++){
	  int a = 0;
	  for(int ha = 0; ha < hAG; ha++){
	    double infections_a, infections_ha = 0.0;
	    for(int i = 0; i < hAG_SPAN[ha]; i++){
	      infections_ha += infections_a = pop[t][HIVN][g][a] * incrate15to49_g[g] * incrr_age[t][g][a] * Xhivn_g[g] / Xhivn_incagerr[g];
	      infections[t][g][a] += DT*infections_a;
	      pop[t][HIVN][g][a] -= DT*infections_a;
	      pop[t][HIVP][g][a] += DT*infections_a;
	      a++;
	    }
	    if(ha < hIDX_15TO49+hAG_15TO49 )
	      incid15to49[t] += DT*infections_ha;
	    
	    // add infections to grad hivpop
	    for(int hm = 0; hm < hDS; hm++)
	      grad[g][ha][hm] += infections_ha * cd4_initdist[g][ha][hm];
	  }
	}
	
	
	for(int g = 0; g < NG; g++)
	  for(int ha = 0; ha < hAG; ha++)
	    for(int hm = 0; hm < hDS; hm++)
	      hivpop[t][g][ha][hm] += DT*grad[g][ha][hm];

	// ART progression, mortality, and initiation
	if(t >= t_ART_start){
	  int cd4elig_idx = artcd4elig_idx[t] - 1; // -1 for 0-based indexing vs. 1-based in R
	  int anyelig_idx = (specpop_percelig[t] > 0 | pw_artelig[t] > 0) ? 0 : (who34percelig > 0) ? hIDX_CD4_350 : cd4elig_idx;

	  // progression and mortality
	  for(int g = 0; g < NG; g++)
	    for(int ha = 0; ha < hAG; ha++)
		for(int hm = anyelig_idx; hm < hDS; hm++){
		  double gradART[hTS];
		  
		  for(int hu = 0; hu < hTS; hu++){
		    double deaths = art_mort[g][ha][hm][hu] * artpop[t][g][ha][hm][hu];
		    hivdeaths_ha[g][ha] += DT*deaths;
		    gradART[hu] = -deaths;
		  }
		  
		  gradART[ART0MOS] += -ART_STAGE_PROG_RATE * artpop[t][g][ha][hm][ART0MOS];
		  gradART[ART6MOS] += ART_STAGE_PROG_RATE * artpop[t][g][ha][hm][ART0MOS] - ART_STAGE_PROG_RATE * artpop[t][g][ha][hm][ART6MOS];
		  gradART[ART1YR] += ART_STAGE_PROG_RATE * artpop[t][g][ha][hm][ART6MOS];
		  
		  for(int hu = 0; hu < hTS; hu++)
		    artpop[t][g][ha][hm][hu] += DT*gradART[hu];
		}

	  // ART dropout
	  if(art_dropout[t] > 0){
	  for(int g = 0; g < NG; g++)
	    for(int ha = 0; ha < hAG; ha++)
	      for(int hm = anyelig_idx; hm < hDS; hm++)
		for(int hu = 0; hu < hTS; hu++){
		  hivpop[t][g][ha][hm] += DT * art_dropout[t] * artpop[t][g][ha][hm][hu];
		  artpop[t][g][ha][hm][hu] -= DT * art_dropout[t] * artpop[t][g][ha][hm][hu];
		}
	  }

	  // ART initiation
          for(int g = 0; g < NG; g++){
	    
            double artelig_hahm[hAG_15PLUS][hDS], Xart_15plus = 0.0, Xartelig_15plus = 0.0, expect_mort_artelig15plus = 0.0;
            for(int ha = hIDX_15PLUS; ha < hAG; ha++){
              for(int hm = anyelig_idx; hm < hDS; hm++){
		double prop_elig = (hm >= cd4elig_idx) ? 1.0 : (hm >= hIDX_CD4_350) ? 1.0 - (1.0-specpop_percelig[t])*(1.0-who34percelig) : specpop_percelig[t];
                Xartelig_15plus += artelig_hahm[ha-hIDX_15PLUS][hm] = prop_elig * hivpop[t][g][ha][hm] ;
                expect_mort_artelig15plus += cd4_mort[g][ha][hm] * artelig_hahm[ha-hIDX_15PLUS][hm];
		  for(int hu = 0; hu < hTS; hu++)
                  Xart_15plus += artpop[t][g][ha][hm][hu];
              }

              // if pw_artelig, add pregnant women to artelig_hahm population
              if(g == FEMALE & pw_artelig[t] > 0 & ha < hAG_FERT){
                double frr_pop_ha = 0;
                for(int a =  hAG_START[ha]; a < hAG_START[ha]+hAG_SPAN[ha]; a++)
                  frr_pop_ha += pop[t][HIVN][g][a]; // add HIV- population
                for(int hm = 0; hm < hDS; hm++){
                  frr_pop_ha += frr_cd4[t][ha-hIDX_FERT][hm] * hivpop[t][g][ha][hm];
                  for(int hu = 0; hu < hTS; hu++)
                    frr_pop_ha += frr_art[t][ha-hIDX_FERT][hm][hu] * artpop[t][g][ha][hm][hu];
                }

                for(int hm = anyelig_idx; hm < cd4elig_idx; hm++){
                  double pw_elig_hahm = DT * births_by_ha[ha-hIDX_FERT] * frr_cd4[t][ha-hIDX_FERT][hm] * hivpop[t][g][ha][hm] / frr_pop_ha;
                  artelig_hahm[ha-hIDX_15PLUS][hm] += pw_elig_hahm;
                  Xartelig_15plus += pw_elig_hahm;
                  expect_mort_artelig15plus += cd4_mort[g][ha][hm] * pw_elig_hahm;
                }
	      }
            } // loop over ha

	    // calculate number on ART at end of ts, based on number or percent
            double artnum_hts = 0.0;
            if(DT*(hts+1) < 0.5){
	      if(!art15plus_isperc[t-2][g] & !art15plus_isperc[t-1][g]){ // both numbers
		artnum_hts = (0.5-DT*(hts+1))*artnum15plus[t-2][g] + (DT*(hts+1)+0.5)*artnum15plus[t-1][g];
	      } else if(art15plus_isperc[t-2][g] & art15plus_isperc[t-1][g]){ // both percentages
		double artcov_hts = (0.5-DT*(hts+1))*artnum15plus[t-2][g] + (DT*(hts+1)+0.5)*artnum15plus[t-1][g];
		artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
	      } else if(!art15plus_isperc[t-2][g] & art15plus_isperc[t-1][g]){ // transition from number to percentage
		double curr_coverage = Xart_15plus / (Xart_15plus + Xartelig_15plus);
		double artcov_hts = curr_coverage + (artnum15plus[t-1][g] - curr_coverage) * DT / (0.5-DT*hts);
		artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
	      } 
            } else {
	      if(!art15plus_isperc[t-1][g] & !art15plus_isperc[t][g]){ // both numbers
		artnum_hts = (1.5-DT*(hts+1))*artnum15plus[t-1][g] + (DT*(hts+1)-0.5)*artnum15plus[t][g];
	      } else if(art15plus_isperc[t-1][g] & art15plus_isperc[t][g]){ // both percentages
		double artcov_hts = (1.5-DT*(hts+1))*artnum15plus[t-1][g] + (DT*(hts+1)-0.5)*artnum15plus[t][g];
		artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
	      } else if(!art15plus_isperc[t-1][g] & art15plus_isperc[t][g]){ // transition from number to percentage
		double curr_coverage = Xart_15plus / (Xart_15plus + Xartelig_15plus);
		double artcov_hts = curr_coverage + (artnum15plus[t][g] - curr_coverage) * DT / (1.5-DT*hts);
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
		   elig_above += artelig_hahm[ha-hIDX_15PLUS][hm];
		 elig_above += (1.0 - medcat_propbelow) * artelig_hahm[ha-hIDX_15PLUS][medcd4_idx];
		 elig_below += medcat_propbelow * artelig_hahm[ha-hIDX_15PLUS][medcd4_idx];
		 for(int hm = medcd4_idx+1; hm < hDS; hm++)
		   elig_below += artelig_hahm[ha-hIDX_15PLUS][hm];
	       }

	       double initprob_below = artinit_hts * 0.5 / elig_below;
	       double initprob_above = artinit_hts * 0.5 / elig_above;
	       double initprob_medcat = initprob_below * medcat_propbelow + initprob_above * (1.0-medcat_propbelow);

	       for(int ha = hIDX_15PLUS; ha < hAG; ha++)
		 for(int hm = anyelig_idx; hm < hDS; hm++){
		   double artinit_hahm;
		   if(hm < medcd4_idx)
		     artinit_hahm = artelig_hahm[ha-hIDX_15PLUS][hm] * initprob_above;
		   else if(hm == medcd4_idx)
		     artinit_hahm = artelig_hahm[ha-hIDX_15PLUS][hm] * initprob_medcat;
		   if(hm > medcd4_idx)
		     artinit_hahm = artelig_hahm[ha-hIDX_15PLUS][hm] * initprob_below;
		   if(artinit_hahm > hivpop[t][g][ha][hm]) artinit_hahm = hivpop[t][g][ha][hm];
		   hivpop[t][g][ha][hm] -= artinit_hahm;
		   artpop[t][g][ha][hm][ART0MOS] += artinit_hahm;
		 }

	    } else { // Use mixture of eligibility and expected mortality for initiation distribution
	    
            for(int ha = hIDX_15PLUS; ha < hAG; ha++)
              for(int hm = anyelig_idx; hm < hDS; hm++){
                double artinit_hahm = artinit_hts * artelig_hahm[ha-hIDX_15PLUS][hm] * 0.5 * (1.0/Xartelig_15plus + cd4_mort[g][ha][hm] / expect_mort_artelig15plus);
                if(artinit_hahm > artelig_hahm[ha-hIDX_15PLUS][hm]) artinit_hahm = artelig_hahm[ha-hIDX_15PLUS][hm];
                hivpop[t][g][ha][hm] -= artinit_hahm;
                artpop[t][g][ha][hm][ART0MOS] += artinit_hahm;
              }
	    }

          }
        }


	// remove hivdeaths from pop
	for(int g = 0; g < NG; g++){
	  
	  // sum HIV+ population size in each hivpop age group
	  double hivpop_ha[hAG];
	  int a = 0;
	  for(int ha = 0; ha < hAG; ha++){
	    hivpop_ha[ha] = 0.0;
	    for(int i = 0; i < hAG_SPAN[ha]; i++){
	      hivpop_ha[ha] += pop[t][HIVP][g][a];
	      a++;
	    }
	  }
	  
	  // remove hivdeaths proportionally to age-distribution within each age group
	  a = 0;
	  for(int ha = 0; ha < hAG; ha++){
	    if(hivpop_ha[ha] > 0){
	      double hivqx_ha = hivdeaths_ha[g][ha] / hivpop_ha[ha];
	      for(int i = 0; i < hAG_SPAN[ha]; i++){
		hivdeaths[t][g][a] += pop[t][HIVP][g][a] * hivqx_ha;
		pop[t][HIVP][g][a] *= (1.0-hivqx_ha);
		a++;
	      }
	    }  // end if(pop_ha[ha] > 0)
	  }
	}

      } // loop HIVSTEPS_PER_YEAR


      /*
      // Code for calculating new infections once per year to match prevalence (like Spectrum)

      // (1) incidence from prevalence input
      double Xhivp = 0.0, Xhivn[NG], Xhivn_incagerr[NG];

      for(int g = 0; g < NG; g++){
        Xhivn[g] = 0.0;
        Xhivn_incagerr[g] = 0.0;
        for(int a = pIDX_15TO49; a < pIDX_15TO49+pAG_15TO49; a++){
          Xhivp += pop[t][HIVP][g][a];
          Xhivn[g] += pop[t][HIVN][g][a];
          Xhivn_incagerr[g] += incrr_age[t][g][a] * pop[t][HIVN][g][a];
        }
      }
      double prev_i = Xhivp / (Xhivn[MALE] + Xhivn[FEMALE] + Xhivp);
      double incrate15to49_i = (prev15to49[t] - prev_i)/(1.0 - prev_i);
      double incrate15to49_g[NG];
      incrate15to49_g[MALE] = incrate15to49_i * (Xhivn[MALE]+Xhivn[FEMALE]) / (Xhivn[MALE] + incrr_sex[t]*Xhivn[FEMALE]);
      incrate15to49_g[FEMALE] = incrate15to49_i * incrr_sex[t]*(Xhivn[MALE]+Xhivn[FEMALE]) / (Xhivn[MALE] + incrr_sex[t]*Xhivn[FEMALE]);

      for(int g = 0; g < NG; g++){
        int a = 0;
        for(int ha = 0; ha < hAG; ha++){
          double infections_a, infections_ha = 0.0;
          for(int i = 0; i < hAG_SPAN[ha]; i++){
            infections_ha += infections_a = pop[t][HIVN][g][a] * incrate15to49_g[g] * incrr_age[t][g][a] * Xhivn[g] / Xhivn_incagerr[g];
            pop[t][HIVN][g][a] -= infections_a;
            pop[t][HIVP][g][a] += infections_a;
            a++;
          }

          // add infections to hivpop
          for(int hm = 0; hm < hDS; hm++)
            hivpop[t][g][ha][hm] += infections_ha * cd4_initdist[g][ha][hm];
        }
      }
      */


      // adjust population to match target population
      if(bin_popadjust){
	for(int g = 0; g < NG; g++){
	  int a = 0;
	  for(int ha = 0; ha < hAG; ha++){
	    double popadj_ha = 0, hivpop_ha = 0;
	    for(int i = 0; i < hAG_SPAN[ha]; i++){
	      
	      hivpop_ha += pop[t][HIVP][g][a];
	      
	      double popadjrate_a = popadjust[t][g][a] = targetpop[t][g][a] / (pop[t][HIVN][g][a] + pop[t][HIVP][g][a]);
	      pop[t][HIVN][g][a] *= popadjrate_a;
	      double hpopadj_a = (popadjrate_a-1.0) * pop[t][HIVP][g][a];
	      popadj_ha += hpopadj_a;
	      pop[t][HIVP][g][a] += hpopadj_a;
	      a++;
	    }
	    
	    // population adjustment for hivpop
	    double popadjrate_ha = hivpop_ha > 0 ? popadj_ha / hivpop_ha : 0.0;
	    for(int hm = 0; hm < hDS; hm++){
	      hivpop[t][g][ha][hm] *= 1+popadjrate_ha;
	      if(t >= t_ART_start)
		for(int hu = 0; hu < hTS; hu++)
		  artpop[t][g][ha][hm][hu] *= 1+popadjrate_ha;
	    } // loop over hm
	  } // loop over ha
	} // loop over g
      } // if(bin_popadjust)

      

      // prevalence among pregnant women
      
        double hivbirths = 0;
	for(int ha = hIDX_FERT; ha < hIDX_FERT+hAG_FERT; ha++){
	  double hivn_ha = 0, frr_hivpop_ha = 0;
	  for(int a =  hAG_START[ha]; a < hAG_START[ha]+hAG_SPAN[ha]; a++)
            hivn_ha += (pop[t-1][HIVN][FEMALE][a] + pop[t][HIVN][FEMALE][a])/2;
	  for(int hm = 0; hm < hDS; hm++){
	    frr_hivpop_ha += frr_cd4[t][ha-hIDX_FERT][hm] * (hivpop[t-1][FEMALE][ha][hm]+hivpop[t][FEMALE][ha][hm])/2;
	    if(t == t_ART_start)
	      for(int hu = 0; hu < hTS; hu++)
		frr_hivpop_ha += frr_art[t][ha-hIDX_FERT][hm][hu] * artpop[t][FEMALE][ha][hm][hu]/2;
	    else if(t > t_ART_start)
	      for(int hu = 0; hu < hTS; hu++)
		frr_hivpop_ha += frr_art[t][ha-hIDX_FERT][hm][hu] * (artpop[t-1][FEMALE][ha][hm][hu]+artpop[t][FEMALE][ha][hm][hu])/2;
	  }
	  hivbirths += births_by_ha[ha-hIDX_FERT] * frr_hivpop_ha / (hivn_ha + frr_hivpop_ha);
	}

	pregprev[t] = hivbirths/births;
	if(t + AGE_START < PROJ_YEARS)
	  pregprevlag[t + AGE_START-1] = pregprev[t];

	// prevalence 15 to 49
	for(int g = 0; g < NG; g++)
	  for(int a = pIDX_15TO49; a < pIDX_15TO49 + pAG_15TO49; a++){
	    hivn15to49[t] += pop[t][HIVN][g][a];
	    hivp15to49[t] += pop[t][HIVP][g][a];
	  }
	prev15to49[t] = hivp15to49[t]/(hivn15to49[t] + hivp15to49[t]);
	incid15to49[t] /= hivn15to49[t-1];
    }

    UNPROTECT(21);
    return s_pop;
  }
}



SEXP getListElement(SEXP list, const char *str)
{
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;
  for ( i = 0; i < length(list); i++ )
    if ( strcmp(CHAR(STRING_ELT(names, i)), str) == 0 ) {
      elmt = VECTOR_ELT(list, i);
      break;
    }

    if ( elmt == R_NilValue )
    error("%s missing from list", str);

  return elmt;
}
