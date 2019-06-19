#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>

#include <R.h>
#include <Rinternals.h>

#include "rutils.h"
#include "model.hpp"

#define EPP_RSPLINE 0
#define EPP_RTREND 1
#define EPP_DIRECTINCID 2  // annual direct incidence inputs (as Spectrum)

void simulate_eppasm2(Model &md,
		      const SsDim &ss,
		      const DemogParam &demp,
		      const PaediatricHivParam &paedhp,
		      const NaturalHistoryParam &nhp,
		      const ArtData &artp,
		      const IncidenceParam &incp);

void simulate_eppasm3(Model &md,
		      const SsDim &ss,
		      const DemogParam &demp,
		      const PaediatricHivParam &paedhp,
		      const NaturalHistoryParam &nhp,
		      const ArtData &artp,
		      const IncidenceParam &incp);

void simulate_eppasm4(Model &md,
		      const SsDim &ss,
		      const DemogParam &demp,
		      const PaediatricHivParam &paedhp,
		      const NaturalHistoryParam &nhp,
		      const ArtData &artp,
		      const IncidenceParam &incp);

void simulate_eppasm5(Model &md,
		      const SsDim &ss,
		      const DemogParam &demp,
		      const PaediatricHivParam &paedhp,
		      const NaturalHistoryParam &nhp,
		      const ArtData &artp,
		      const IncidenceParam &incp);

extern "C" {

  SEXP eppasmR(SEXP s_fp, SEXP s_version){

    ////////////////////////////////
    ////  set parameter values  ////
    ////////////////////////////////

    // state space dimensions
    SEXP s_ss = getListElement(s_fp, "ss");
    SsDim ss(*INTEGER(getListElement(s_ss, "NG")),
             *INTEGER(getListElement(s_ss, "pAG")),
             *INTEGER(getListElement(s_ss, "pDS")),
             *INTEGER(getListElement(s_ss, "hAG")),
             *INTEGER(getListElement(s_ss, "hDS")),
             *INTEGER(getListElement(s_ss, "hTS")),
             *INTEGER(getListElement(s_ss, "proj_start")),
             *INTEGER(getListElement(s_ss, "PROJ_YEARS")),
             *INTEGER(getListElement(s_fp, "SIM_YEARS")),
             *INTEGER(getListElement(s_ss, "hiv_steps_per_year")),
             *INTEGER(getListElement(s_fp, "tARTstart")) - 1, // -1 for 0-based indexing in C vs. 1-based in R
             INTEGER(getListElement(s_ss, "h.ag.span")));

    double *projsteps = REAL(getListElement(s_fp, "proj.steps"));

    // demographic parameters

    int bin_popadjust = *INTEGER(getListElement(s_fp, "popadjust"));
    double *ptr_targetpop = nullptr;
    double *ptr_entrantpop = nullptr;
    if(bin_popadjust){
      ptr_targetpop = REAL(getListElement(s_fp, "targetpop"));
      ptr_entrantpop = REAL(getListElement(s_fp, "entrantpop"));
    }

    DemogParam demp(ss,
                    REAL(getListElement(s_fp, "basepop")),
                    REAL(getListElement(s_fp, "Sx")),
                    REAL(getListElement(s_fp, "netmigr")),
                    REAL(getListElement(s_fp, "asfr")),
                    REAL(getListElement(s_fp, "srb")),
                    REAL(getListElement(s_fp, "birthslag")),
                    REAL(getListElement(s_fp, "cumsurv")),
                    REAL(getListElement(s_fp, "cumnetmigr")),
                    bin_popadjust,
                    ptr_targetpop,
                    ptr_entrantpop);

    // disease progression

    NaturalHistoryParam nhp(ss,
                            REAL(getListElement(s_fp, "cd4_initdist")),
                            REAL(getListElement(s_fp, "cd4_prog")),
                            REAL(getListElement(s_fp, "cd4_mort")),
                            REAL(getListElement(s_fp, "art_mort")),
                            REAL(getListElement(s_fp, "artmx_timerr")),
                            REAL(getListElement(s_fp, "frr_cd4")),
                            REAL(getListElement(s_fp, "frr_art")));

    // ART inputs

    ArtData artp(ss,
                 REAL(getListElement(s_fp, "art15plus_num")),
                 LOGICAL(getListElement(s_fp, "art15plus_isperc")),
                 INTEGER(getListElement(s_fp, "artcd4elig_idx")),
                 REAL(getListElement(s_fp, "specpop_percelig")),
                 REAL(getListElement(s_fp, "pw_artelig")),
                 *REAL(getListElement(s_fp, "who34percelig")),
                 REAL(getListElement(s_fp, "art_dropout")),
                 REAL(getListElement(s_fp, "median_cd4init")),
                 INTEGER(getListElement(s_fp, "med_cd4init_cat")),
                 INTEGER(getListElement(s_fp, "med_cd4init_input")),
                 *INTEGER(getListElement(s_fp, "art_alloc_method")),
                 *REAL(getListElement(s_fp, "art_alloc_mxweight")),
                 *INTEGER(getListElement(s_fp, "scale_cd4_mort")));

    // incidence model

    int eppmod = *INTEGER(getListElement(s_fp, "eppmodInt"));

    double *incidinput = nullptr;
    double tsEpidemicStart = 0, iota = 0, rel_infect_art = 0;

    double *rspline_rvec = nullptr;
    double *rtrend_beta = nullptr, rtrend_tstab = 0, rtrend_r0 = 0;
    double relinfectART;

    int incidpopage = 0; // default INCIDPOP_15TO49

    if(eppmod == EPP_DIRECTINCID){
      incidinput = REAL(getListElement(s_fp, "incidinput"));
      incidpopage = *INTEGER(getListElement(s_fp, "incidpopage"));
    } else {
      rel_infect_art = *REAL(getListElement(s_fp, "relinfectART"));
      tsEpidemicStart = *REAL(getListElement(s_fp, "tsEpidemicStart"));
      iota = *REAL(getListElement(s_fp, "iota"));
      relinfectART = rel_infect_art;

      if(eppmod == EPP_RSPLINE)
        rspline_rvec = REAL(getListElement(s_fp, "rvec"));
      else if(eppmod == EPP_RTREND){
        SEXP s_rtrend = getListElement(s_fp, "rtrend");
        rtrend_beta = REAL(getListElement(s_rtrend, "beta"));
        rtrend_tstab = *REAL(getListElement(s_rtrend, "tStabilize"));
        rtrend_r0 = *REAL(getListElement(s_rtrend, "r0"));
      }
    }

    int ts_epidemic_start;
    for(ts_epidemic_start = 0; ts_epidemic_start < ss.proj_steps; ts_epidemic_start++)
      if(projsteps[ts_epidemic_start] >= tsEpidemicStart)
        break;

    IncidenceParam incp(ss,
                        eppmod,
                        incidpopage,
                        ts_epidemic_start,
                        iota,
                        rel_infect_art,
                        rspline_rvec,
                        rtrend_beta,
                        rtrend_tstab,
                        rtrend_r0,
                        incidinput,
                        REAL(getListElement(s_fp, "incrr_sex")),
                        REAL(getListElement(s_fp, "incrr_age")));

    // vertical transmission and survival
    double *a_entrantprev = nullptr;
    int use_entrantprev = checkListElement(s_fp, "entrantprev");
    if(use_entrantprev)
      a_entrantprev = REAL(getListElement(s_fp, "entrantprev"));

    double *a_entrantartcov;
    if(checkListElement(s_fp, "entrantartcov"))
      a_entrantartcov = REAL(getListElement(s_fp, "entrantartcov"));
    else {
      a_entrantartcov = (double*) R_alloc(ss.proj_years, sizeof(double));
      memset(a_entrantartcov, 0, ss.proj_years*sizeof(double));
    }

    PaediatricHivParam paedhp(ss,
                              use_entrantprev,
                              REAL(getListElement(s_fp, "verttrans_lag")),
                              REAL(getListElement(s_fp, "paedsurv_lag")),
                              *REAL(getListElement(s_fp, "netmig_hivprob")),
                              a_entrantprev,
                              a_entrantartcov,
                              REAL(getListElement(s_fp, "paedsurv_cd4dist")),
                              REAL(getListElement(s_fp, "paedsurv_artcd4dist")));

    // initialize output

    int pop_dim[4] = {ss.p_ag, ss.ng, ss.p_ds, ss.proj_years};
    int hivpop_dim[4] = {ss.h_ds, ss.h_ag, ss.ng, ss.proj_years};
    int artpop_dim[5] = {ss.h_ts, ss.h_ds, ss.h_ag, ss.ng, ss.proj_years};
    int pop1_dim[3] = {ss.p_ag, ss.ng, ss.proj_years};

    SEXP s_pop = sexp_numeric_array(4, pop_dim);
    SEXP s_hivpop = sexp_numeric_array(4, hivpop_dim);
    SEXP s_artpop = sexp_numeric_array(5, artpop_dim);
    SEXP s_infections = sexp_numeric_array(3, pop1_dim);
    SEXP s_hivdeaths = sexp_numeric_array(3, pop1_dim);
    SEXP s_natdeaths = sexp_numeric_array(3, pop1_dim);
    SEXP s_popadjust = sexp_numeric_array(3, pop1_dim);

    SEXP s_incid15to49_ts = PROTECT(allocVector(REALSXP, (ss.proj_years-1) * ss.hiv_steps_per_year));
    SEXP s_prev15to49_ts = PROTECT(allocVector(REALSXP, (ss.proj_years-1) * ss.hiv_steps_per_year));
    SEXP s_rvec_ts = PROTECT(allocVector(REALSXP, (ss.proj_years-1) * ss.hiv_steps_per_year));
    SEXP s_prev15to49 = PROTECT(allocVector(REALSXP, ss.proj_years));
    SEXP s_incid15to49 = PROTECT(allocVector(REALSXP, ss.proj_years));
    SEXP s_pregprev = PROTECT(allocVector(REALSXP, ss.proj_years));
    SEXP s_entrantprev_out = PROTECT(allocVector(REALSXP, ss.proj_years));

    setAttrib(s_pop, install("hivpop"), s_hivpop);
    setAttrib(s_pop, install("artpop"), s_artpop);
    setAttrib(s_pop, install("infections"), s_infections);
    setAttrib(s_pop, install("hivdeaths"), s_hivdeaths);
    setAttrib(s_pop, install("natdeaths"), s_natdeaths);
    setAttrib(s_pop, install("popadjust"), s_popadjust);

    setAttrib(s_pop, install("prev15to49"), s_prev15to49);
    setAttrib(s_pop, install("pregprev"), s_pregprev);
    setAttrib(s_pop, install("incid15to49"), s_incid15to49);
    setAttrib(s_pop, install("entrantprev"), s_entrantprev_out);

    setAttrib(s_pop, install("prev15to49_ts"), s_prev15to49_ts);
    setAttrib(s_pop, install("incid15to49_ts"), s_incid15to49_ts);
    setAttrib(s_pop, install("rvec_ts"), s_rvec_ts);

    // create model object
    Model md(ss,
             demp,
             paedhp,
	     nhp,
	     artp,
             REAL(s_pop),
             REAL(s_hivpop),
             REAL(s_artpop),
             REAL(s_infections),
             REAL(s_hivdeaths),
             REAL(s_natdeaths),
             REAL(s_popadjust),
             REAL(s_prev15to49),
             REAL(s_incid15to49),
             REAL(s_pregprev),
             REAL(s_entrantprev_out),
             REAL(s_prev15to49_ts),
             REAL(s_incid15to49_ts),
             REAL(s_rvec_ts));

    if(*INTEGER(s_version) == 2)
      simulate_eppasm2(md, ss, demp, paedhp, nhp, artp, incp);
    else if(*INTEGER(s_version) == 3)
      simulate_eppasm3(md, ss, demp, paedhp, nhp, artp, incp);
    else if(*INTEGER(s_version) == 4)
      simulate_eppasm4(md, ss, demp, paedhp, nhp, artp, incp);
    else if(*INTEGER(s_version) == 5)
      simulate_eppasm5(md, ss, demp, paedhp, nhp, artp, incp);

    UNPROTECT(21);
    return s_pop;
  }
}
