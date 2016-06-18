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

#define MALE 0
#define FEMALE 1

#define HIVN 0
#define HIVP 1

#define ART0MOS 0
#define ART6MOS 1
#define ART1YR 2

#define ART_STAGE_PROG_RATE 2.0 // HARD CODED: ART stage progression rate



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

    // demographic projection
    multi_array_ref<double, 2> basepop(REAL(getListElement(s_fp, "basepop")), extents[NG][pAG]);
    multi_array_ref<double, 3> Sx(REAL(getListElement(s_fp, "Sx")), extents[PROJ_YEARS][NG][pAG]);
    multi_array_ref<double, 3> netmigr(REAL(getListElement(s_fp, "netmigr")), extents[PROJ_YEARS][NG][pAG]);
    multi_array_ref<double, 2> asfr(REAL(getListElement(s_fp, "asfr")), extents[PROJ_YEARS][pAG_FERT]);
    multi_array_ref<double, 2> srb(REAL(getListElement(s_fp, "srb")), extents[PROJ_YEARS][NG]);
    multi_array_ref<double, 2> birthslag(REAL(getListElement(s_fp, "birthslag")), extents[PROJ_YEARS][NG]);
    multi_array_ref<double, 2> cumsurv(REAL(getListElement(s_fp, "cumsurv")), extents[PROJ_YEARS][NG]);
    multi_array_ref<double, 2> cumnetmigr(REAL(getListElement(s_fp, "cumnetmigr")), extents[PROJ_YEARS][NG]);

    // disease progression
    multi_array_ref<double, 3> cd4_initdist(REAL(getListElement(s_fp, "cd4_initdist")), extents[NG][hAG][hDS]);
    multi_array_ref<double, 3> cd4_prog(REAL(getListElement(s_fp, "cd4_prog")), extents[NG][hAG][hDS-1]);
    multi_array_ref<double, 3> cd4_mort(REAL(getListElement(s_fp, "cd4_mort")), extents[NG][hAG][hDS]);
    multi_array_ref<double, 4> art_mort(REAL(getListElement(s_fp, "art_mort")), extents[NG][hAG][hDS][hTS]);

    // sub-fertility
    multi_array_ref<double, 2> frr_cd4(REAL(getListElement(s_fp, "frr_cd4")), extents[hAG_FERT][hDS]);
    multi_array_ref<double, 3> frr_art(REAL(getListElement(s_fp, "frr_art")), extents[hAG_FERT][hDS][hTS]);

    // ART inputs
    int t_ART_start = *INTEGER(getListElement(s_fp, "tARTstart")) - 1; // -1 for 0-based indexing in C vs. 1-based in R
    multi_array_ref<double, 2> artnum15plus(REAL(getListElement(s_fp, "art15plus_num")), extents[PROJ_YEARS][NG]);
    multi_array_ref<double, 2> artcd4elig(REAL(getListElement(s_fp, "artcd4elig")), extents[PROJ_YEARS][hDS]);
    multi_array_ref<double, 2> pw_artelig(REAL(getListElement(s_fp, "pw_artelig")), extents[PROJ_YEARS][hDS]);

    // incidence model
    // double *prev15to49 = REAL(getListElement(s_fp, "prev15to49"));
    double *incrr_sex = REAL(getListElement(s_fp, "incrr_sex"));
    multi_array_ref<double, 3> incrr_age(REAL(getListElement(s_fp, "incrr_age")), extents[PROJ_YEARS][NG][pAG]);

    double *rvec = REAL(getListElement(s_fp, "rvec"));
    double iota = *REAL(getListElement(s_fp, "iota"));
    double relinfectART = *REAL(getListElement(s_fp, "relinfectART"));
    double ts_epidemic_start = *INTEGER(getListElement(s_fp, "ts_epi_start")) - 1; // -1 for 0-based indexing in C vs. 1-based in R
    

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
    
    SEXP s_artpop = PROTECT(allocVector(REALSXP, hTS * hDS * hAG * NG * PROJ_YEARS));
    SEXP s_artpop_dim = PROTECT(allocVector(INTSXP, 5));
    INTEGER(s_artpop_dim)[0] = hTS;
    INTEGER(s_artpop_dim)[1] = hDS;
    INTEGER(s_artpop_dim)[2] = hAG;
    INTEGER(s_artpop_dim)[3] = NG;
    INTEGER(s_artpop_dim)[4] = PROJ_YEARS;
    setAttrib(s_artpop, R_DimSymbol, s_artpop_dim);

    SEXP s_pregprevlag = PROTECT(allocVector(REALSXP, PROJ_YEARS));
    
    setAttrib(s_pop, install("hivpop"), s_hivpop);
    setAttrib(s_pop, install("artpop"), s_artpop);
    setAttrib(s_pop, install("pregprevlag"), s_pregprevlag);

    SEXP s_incrate15to49_ts = PROTECT(allocVector(REALSXP, (PROJ_YEARS-1) * HIVSTEPS_PER_YEAR));
    setAttrib(s_pop, install("incrate15to49_ts"), s_incrate15to49_ts);
    double *incrate15to49_ts_out = REAL(s_incrate15to49_ts);
    memset(incrate15to49_ts_out, 0, length(s_incrate15to49_ts)*sizeof(double));
    
    // initialize population

    // population by single-year age
    // double pop[PROJ_YEARS][pDS][NG][pAG];
    multi_array_ref<double, 4> pop(REAL(s_pop), extents[PROJ_YEARS][pDS][NG][pAG]);
    for(int g = 0; g < NG; g++)
      for(int a = 0; a < pAG; a++){
        pop[0][HIVN][g][a] = basepop[g][a];
        pop[0][HIVP][g][a] = 0.0;
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
    for(int g = 0; g < NG; g++)
      for(int ha = 0; ha < hAG; ha++)
        for(int hm = 0; hm < hDS; hm++)
          for(int hu = 0; hu < hTS; hu++)
            artpop[t_ART_start][g][ha][hm][hu] = 0.0;  // initialize to zero in year of ART start


    // array to store lagged prevalence among pregnant women
    double *pregprevlag = REAL(s_pregprevlag); // (double*) R_alloc(PROJ_YEARS, sizeof(double));
    memset(pregprevlag, 0, AGE_START*sizeof(double));

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
	pop[t][HIVN][g][0] = birthslag[t-1][g] * cumsurv[t-1][g] * (1.0-pregprevlag[t-1] * verttrans) + cumnetmigr[t-1][g] * (1.0-pregprevlag[t-1] * netmig_hivprob);

	double paedsurv_g = birthslag[t-1][g] * cumsurv[t-1][g] * pregprevlag[t-1] * verttrans * paedsurv + cumnetmigr[t-1][g] * pregprevlag[t-1] * netmig_hivprob * netmighivsurv;
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
          double deathsmig_ha = 0, hivpop_ha = 0, hdeaths_a, hmig_a;
          for(int i = 0; i < hAG_SPAN[ha]; i++){

            hivpop_ha += pop[t][HIVP][g][a];

            // non-HIV mortality
            pop[t][HIVN][g][a] *= Sx[t][g][a]; // survival HIV- population
            deathsmig_ha -= hdeaths_a = pop[t][HIVP][g][a] * (1.0 - Sx[t][g][a]);
            pop[t][HIVP][g][a] -= hdeaths_a;   // survival HIV+ population

            // net migration
            double migrate_a = netmigr[t][g][a] * (1+Sx[t][g][a])/2.0 / (pop[t][HIVN][g][a] + pop[t][HIVP][g][a]);
            pop[t][HIVN][g][a] *= 1+migrate_a;
            deathsmig_ha += hmig_a = migrate_a * pop[t][HIVP][g][a];
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
	double hivdeaths[NG][hAG];
	memset(hivdeaths, 0, sizeof(double)*NG*hAG);

	/*
	// incidence

	*/

        // CD4 progression and mortality

	/*
        for(int g = 0; g < NG; g++)
          for(int ha = 0; ha < hAG; ha++){
            double grad[hDS];

            // untreated population
            for(int hm = 0; hm < hDS; hm++){
              double deaths = cd4_mort[g][ha][hm] * hivpop[t][g][ha][hm];
              hivdeaths[g][ha] += DT*deaths;
              grad[hm] = -deaths;
            }
            for(int hm = 1; hm < hDS; hm++){
              grad[hm-1] -= cd4_prog[g][ha][hm-1] * hivpop[t][g][ha][hm-1];
              grad[hm] += cd4_prog[g][ha][hm-1] * hivpop[t][g][ha][hm-1];
            }
            for(int hm = 0; hm < hDS; hm++)
              hivpop[t][g][ha][hm] += DT*grad[hm];

            // ART population
            if(t >= t_ART_start)
              for(int hm = 0; hm < hDS; hm++){
                double gradART[hTS];

                for(int hu = 0; hu < hTS; hu++){
                  double deaths = art_mort[g][ha][hm][hu] * artpop[t][g][ha][hm][hu];
                  hivdeaths[g][ha] += DT*deaths;
                  gradART[hu] = -deaths;
                }

                gradART[ART0MOS] += -ART_STAGE_PROG_RATE * artpop[t][g][ha][hm][ART0MOS];
                gradART[ART6MOS] += ART_STAGE_PROG_RATE * artpop[t][g][ha][hm][ART0MOS] - ART_STAGE_PROG_RATE * artpop[t][g][ha][hm][ART6MOS];
                gradART[ART1YR] += ART_STAGE_PROG_RATE * artpop[t][g][ha][hm][ART6MOS];

                for(int hu = 0; hu < hTS; hu++)
                  artpop[t][g][ha][hm][hu] += DT*gradART[hu];
              } // loop over hm
          } // loop over ha

	*/

	// untreated population
	double grad[NG][hAG][hDS];
        for(int g = 0; g < NG; g++)
          for(int ha = 0; ha < hAG; ha++){
            for(int hm = 0; hm < hDS; hm++){
              double deaths = cd4_mort[g][ha][hm] * hivpop[t][g][ha][hm];
              hivdeaths[g][ha] += DT*deaths;
              grad[g][ha][hm] = -deaths;
            }
            for(int hm = 1; hm < hDS; hm++){
              grad[g][ha][hm-1] -= cd4_prog[g][ha][hm-1] * hivpop[t][g][ha][hm-1];
              grad[g][ha][hm] += cd4_prog[g][ha][hm-1] * hivpop[t][g][ha][hm-1];
	    }
	  }
	
	// incidence

	double Xhivn[NG], Xhivn_incagerr[NG], Xhivp_noart = 0.0, Xart = 0.0;
	for(int g = 0; g < NG; g++){
	  Xhivn[g] = 0.0;
	  Xhivn_incagerr[g] = 0.0;
	  for(int a = pIDX_15TO49; a < pIDX_15TO49+pAG_15TO49; a++){
	    Xhivn[g] += pop[t][HIVN][g][a];
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
	double Xtot = Xhivn[MALE] + Xhivn[FEMALE] + Xhivp_noart + Xart;

	int ts = (t-1)*HIVSTEPS_PER_YEAR + hts;
	double incrate15to49_ts = rvec[ts] * (Xhivp_noart + relinfectART * Xart)/Xtot + ((ts == ts_epidemic_start) ? iota : 0.0);
	double incrate15to49_g[NG];
	incrate15to49_g[MALE] = incrate15to49_ts * (Xhivn[MALE]+Xhivn[FEMALE]) / (Xhivn[MALE] + incrr_sex[t]*Xhivn[FEMALE]);
	incrate15to49_g[FEMALE] = incrate15to49_ts * incrr_sex[t]*(Xhivn[MALE]+Xhivn[FEMALE]) / (Xhivn[MALE] + incrr_sex[t]*Xhivn[FEMALE]);


	incrate15to49_ts_out[ts] = incrate15to49_ts; // !!TEMP CODE
	
	for(int g = 0; g < NG; g++){
	  int a = 0;
	  for(int ha = 0; ha < hAG; ha++){
	    double infections_a, infections_ha = 0.0;
	    for(int i = 0; i < hAG_SPAN[ha]; i++){
	      infections_ha += infections_a = pop[t][HIVN][g][a] * incrate15to49_g[g] * incrr_age[t][g][a] * Xhivn[g] / Xhivn_incagerr[g];
	      pop[t][HIVN][g][a] -= DT*infections_a;
	      pop[t][HIVP][g][a] += DT*infections_a;
	      a++;
	    }
	    // add infections to grad hivpop
	    for(int hm = 0; hm < hDS; hm++)
	      grad[g][ha][hm] += infections_ha * cd4_initdist[g][ha][hm];
	  }
	}
	
	
	for(int g = 0; g < NG; g++)
	  for(int ha = 0; ha < hAG; ha++)
	    for(int hm = 0; hm < hDS; hm++)
	      hivpop[t][g][ha][hm] += DT*grad[g][ha][hm];

	// Progression and mortality for ART population
	if(t >= t_ART_start)
	  for(int g = 0; g < NG; g++)
	    for(int ha = 0; ha < hAG; ha++)
		for(int hm = 0; hm < hDS; hm++){
		  double gradART[hTS];
		  
		  for(int hu = 0; hu < hTS; hu++){
		    double deaths = art_mort[g][ha][hm][hu] * artpop[t][g][ha][hm][hu];
		    hivdeaths[g][ha] += DT*deaths;
		    gradART[hu] = -deaths;
		  }
		  
		  gradART[ART0MOS] += -ART_STAGE_PROG_RATE * artpop[t][g][ha][hm][ART0MOS];
		  gradART[ART6MOS] += ART_STAGE_PROG_RATE * artpop[t][g][ha][hm][ART0MOS] - ART_STAGE_PROG_RATE * artpop[t][g][ha][hm][ART6MOS];
		  gradART[ART1YR] += ART_STAGE_PROG_RATE * artpop[t][g][ha][hm][ART6MOS];
		  
		  for(int hu = 0; hu < hTS; hu++)
		    artpop[t][g][ha][hm][hu] += DT*gradART[hu];
		}

	// remove hivdeaths from pop
	for(int g = 0; g < NG; g++){
	  
	  // sum HIV+ population size in each hivpop age group
	  double pop_ha[hAG];
	  int a = 0;
	  for(int ha = 0; ha < hAG; ha++){
	    pop_ha[ha] = 0.0;
	    for(int i = 0; i < hAG_SPAN[ha]; i++){
	      pop_ha[ha] += pop[t][HIVP][g][a];
	      a++;
	    }
	  }
	  
	  // remove hivdeaths proportionally to age-distribution within each age group
	  a = 0;
	  for(int ha = 0; ha < hAG; ha++){
	    if(pop_ha[ha] > 0){
	      double hivsurv_ha = 1 - hivdeaths[g][ha] / pop_ha[ha];
	      for(int i = 0; i < hAG_SPAN[ha]; i++){
              pop[t][HIVP][g][a] *= hivsurv_ha;
              a++;
	      }
	    }  // end if(pop_ha[ha] > 0)
	  }
	}
	

        // ART initiation
        // calculated separately for men and women
        if(t >= t_ART_start){
          for(int g = 0; g < NG; g++){
	    
            double artelig_hahm[hAG_15PLUS][hDS], Xart_15plus = 0.0, Xartelig_15plus = 0.0, expect_mort_artelig15plus = 0.0;
            for(int ha = hIDX_15PLUS; ha < hAG; ha++){
              for(int hm = 0; hm < hDS; hm++){
                Xartelig_15plus += artelig_hahm[ha-hIDX_15PLUS][hm] = artcd4elig[t][hm] * hivpop[t][g][ha][hm];
                expect_mort_artelig15plus += cd4_mort[g][ha][hm] * artelig_hahm[ha-hIDX_15PLUS][hm];
                for(int hu = 0; hu < hTS; hu++)
                  Xart_15plus += artpop[t][g][ha][hm][hu];
              }

              // if pw_artelig, add pregnant women to artelig_hahm population
              if(g == FEMALE & pw_artelig[t][0] > 0 & ha < hAG_FERT){
                double frr_pop_ha = 0;
                for(int a =  hAG_START[ha]; a < hAG_START[ha]+hAG_SPAN[ha]; a++)
                  frr_pop_ha += pop[t][HIVN][g][a]; // add HIV- population
                for(int hm = 0; hm < hDS; hm++){
                  frr_pop_ha += frr_cd4[ha-hIDX_FERT][hm] * hivpop[t][g][ha][hm];
                  for(int hu = 0; hu < hTS; hu++)
                    frr_pop_ha += frr_art[ha-hIDX_FERT][hm][hu] * artpop[t][g][ha][hm][hu];
                }

                for(int hm = 0; hm < hDS; hm++){
                  double pw_elig_hahm = DT * births_by_ha[ha-hIDX_FERT] * pw_artelig[t][hm] * frr_cd4[ha-hIDX_FERT][hm] * hivpop[t][g][ha][hm] / frr_pop_ha;
                  artelig_hahm[ha-hIDX_15PLUS][hm] += pw_elig_hahm;
                  Xartelig_15plus += pw_elig_hahm;
                  expect_mort_artelig15plus += cd4_mort[g][ha][hm] * pw_elig_hahm;
                }
	      }
            } // loop over ha


            double artnum_hts;
            if(DT*(hts+1) < 0.5)
              artnum_hts = (0.5-DT*(hts+1))*artnum15plus[t-2][g] + (DT*(hts+1)+0.5)*artnum15plus[t-1][g];
            else
              artnum_hts = (1.5-DT*(hts+1))*artnum15plus[t-1][g] + (DT*(hts+1)-0.5)*artnum15plus[t][g];

            double artinit_hts = artnum_hts > Xart_15plus ? artnum_hts - Xart_15plus : 0;

            for(int ha = hIDX_15PLUS; ha < hAG; ha++)
              for(int hm = 0; hm < hDS; hm++){
                double artinit_hahm = artinit_hts * artelig_hahm[ha-hIDX_15PLUS][hm] * 0.5 * (1.0/Xartelig_15plus + cd4_mort[g][ha][hm] / expect_mort_artelig15plus);
                if(artinit_hahm > hivpop[t][g][ha][hm]) artinit_hahm = hivpop[t][g][ha][hm];
                hivpop[t][g][ha][hm] -= artinit_hahm;
                artpop[t][g][ha][hm][ART0MOS] += artinit_hahm;
              }

          }
        }


      } // loop HIVSTEPS_PER_YEAR

      /*
      // remove hivdeaths from pop
      for(int g = 0; g < NG; g++){

        // sum HIV+ population size in each hivpop age group
        double pop_ha[hAG];
        int a = 0;
        for(int ha = 0; ha < hAG; ha++){
          pop_ha[ha] = 0.0;
          for(int i = 0; i < hAG_SPAN[ha]; i++){
            pop_ha[ha] += pop[t][HIVP][g][a];
            a++;
          }
        }

        // remove hivdeaths proportionally to age-distribution within each age group
        a = 0;
        for(int ha = 0; ha < hAG; ha++){
          if(pop_ha[ha] > 0){
            double hivsurv_ha = 1 - hivdeaths[g][ha] / pop_ha[ha];
            for(int i = 0; i < hAG_SPAN[ha]; i++){
              pop[t][HIVP][g][a] *= hivsurv_ha;
              a++;
            }
          }  // end if(pop_ha[ha] > 0)
        }
      }
      */

      /*
      // incidence

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
      // prevalence among pregnant women

      if(t + AGE_START < PROJ_YEARS){
	double hivbirths = 0;
        for(int ha = hIDX_FERT; ha < hIDX_FERT+hAG_FERT; ha++){
	  double hivn_ha = 0, frr_hivpop_ha = 0;
	  for(int a =  hAG_START[ha]; a < hAG_START[ha]+hAG_SPAN[ha]; a++)
            hivn_ha += (pop[t-1][HIVN][FEMALE][a] + pop[t][HIVN][FEMALE][a])/2;
	  for(int hm = 0; hm < hDS; hm++){
	    frr_hivpop_ha += frr_cd4[ha-hIDX_FERT][hm] * (hivpop[t-1][FEMALE][ha][hm]+hivpop[t][FEMALE][ha][hm])/2;
	    if(t == t_ART_start)
	      for(int hu = 0; hu < hTS; hu++)
		frr_hivpop_ha += frr_art[ha-hIDX_FERT][hm][hu] * artpop[t][FEMALE][ha][hm][hu]/2;
	    else if(t > t_ART_start)
	      for(int hu = 0; hu < hTS; hu++)
		frr_hivpop_ha += frr_art[ha-hIDX_FERT][hm][hu] * (artpop[t-1][FEMALE][ha][hm][hu]+artpop[t][FEMALE][ha][hm][hu])/2;
	  }
	  hivbirths += births_by_ha[ha-hIDX_FERT] * frr_hivpop_ha / (hivn_ha + frr_hivpop_ha);
	}
	pregprevlag[t + AGE_START-1] = hivbirths/births;
      }
	
      }

    //////////////////////////
    ////  Record outputs  ////
    //////////////////////////
    /*
      multi_array_ref<double, 4> ma_pop(REAL(s_pop), extents[PROJ_YEARS][pDS][NG][pAG]);
      for(int t = 0; t < PROJ_YEARS; t++)
      for(int m = 0; m < pDS; m++)
      for(int g = 0; g < NG; g++)
      for(int a = 0; a < pAG; a++)
      ma_pop[t][m][g][a] = pop[t][m][g][a];
    */
    /*
      multi_array_ref<double, 4> ma_hivpop(REAL(s_hivpop), extents[PROJ_YEARS][NG][hAG][hDS]);
      for(int t = 0; t < PROJ_YEARS; t++)
      for(int g = 0; g < NG; g++)
      for(int ha = 0; ha < hAG; ha++)
      for(int hm = 0; hm < hDS; hm++)
      ma_hivpop[t][g][ha][hm] = hivpop[t][g][ha][hm];


      multi_array_ref<double, 5> ma_artpop(REAL(s_artpop), extents[PROJ_YEARS][NG][hAG][hDS][hTS]);
      for(int t = 0; t < PROJ_YEARS; t++)
      for(int g = 0; g < NG; g++)
      for(int ha = 0; ha < hAG; ha++)
      for(int hm = 0; hm < hDS; hm++)
      for(int hu = 0; hu < hTS; hu++)
      ma_artpop[t][g][ha][hm][hu] = t < t_ART_start ? 0 : artpop[t][g][ha][hm][hu];
    */
    UNPROTECT(8);
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
