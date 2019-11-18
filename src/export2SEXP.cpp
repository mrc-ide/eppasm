#include "Classes.hpp"

void outputSEXP::finalize(const StateSpace& s) {

  SEXP pop_sexp_dim = PROTECT(NEW_INTEGER(4)); ++np;
  INTEGER(pop_sexp_dim)[0] = s.pAG;
  INTEGER(pop_sexp_dim)[1] = s.NG;
  INTEGER(pop_sexp_dim)[2] = s.pDS;
  INTEGER(pop_sexp_dim)[3] = s.PROJ_YEARS;
  SET_DIM(pop, pop_sexp_dim);

  if (s.MODEL!=0) {
    SEXP age_sex_year_dim = PROTECT(NEW_INTEGER(3)); ++np;
    INTEGER(age_sex_year_dim)[0] = s.pAG;
    INTEGER(age_sex_year_dim)[1] = s.NG;
    INTEGER(age_sex_year_dim)[2] = s.PROJ_YEARS;
    SET_DIM(infections, age_sex_year_dim);
    SET_DIM(hivdeaths, age_sex_year_dim);
    SET_DIM(natdeaths, age_sex_year_dim);
    SET_DIM(popadjust, age_sex_year_dim);

    SEXP hiv_sexp_dim = PROTECT(NEW_INTEGER(4)); ++np;
    INTEGER(hiv_sexp_dim)[0] = s.hDS;
    INTEGER(hiv_sexp_dim)[1] = s.hAG;
    INTEGER(hiv_sexp_dim)[2] = s.NG;
    INTEGER(hiv_sexp_dim)[3] = s.PROJ_YEARS;
    SET_DIM(hivpop, hiv_sexp_dim);

    SEXP art_sexp_dim = PROTECT(NEW_INTEGER(5)); ++np;
    INTEGER(art_sexp_dim)[0] = s.hTS;
    INTEGER(art_sexp_dim)[1] = s.hDS;
    INTEGER(art_sexp_dim)[2] = s.hAG;
    INTEGER(art_sexp_dim)[3] = s.NG;
    INTEGER(art_sexp_dim)[4] = s.PROJ_YEARS;
    SET_DIM(artpop, art_sexp_dim);

    SET_ATTR(pop, Rf_install("infections"), infections);
    SET_ATTR(pop, Rf_install("hivdeaths"), hivdeaths);
    SET_ATTR(pop, Rf_install("natdeaths"), natdeaths);
    SET_ATTR(pop, Rf_install("popadjust"), popadjust);
    SET_ATTR(pop, Rf_install("pregprevlag"), pregprevlag);
    SET_ATTR(pop, Rf_install("rvec_ts"), rvec);
    SET_ATTR(pop, Rf_install("prev15to49"), prev15to49);
    SET_ATTR(pop, Rf_install("incid15to49"), incid15to49);
    SET_ATTR(pop, Rf_install("entrantprev"), entrantprev);
    SET_ATTR(pop, Rf_install("incrate15to49_ts"), inci15to49_ts);
    SET_ATTR(pop, Rf_install("prev15to49_ts_sexp"), prev15to49_ts);
    SET_ATTR(pop, Rf_install("artpop"), artpop);
    SET_ATTR(pop, Rf_install("hivpop"), hivpop);
    SET_CLASS(pop, Rf_mkString("spec"));
    if (s.MODEL == 2) {
      SEXP pop_db_sexp_dim = PROTECT(NEW_INTEGER(4)); ++np;
      INTEGER(pop_db_sexp_dim)[0] = s.pDB;
      INTEGER(pop_db_sexp_dim)[1] = s.NG;
      INTEGER(pop_db_sexp_dim)[2] = s.pDS;
      INTEGER(pop_db_sexp_dim)[3] = s.PROJ_YEARS;
      SET_DIM(data_db, pop_db_sexp_dim);

      SET_ATTR(pop, Rf_install("debut_pop"), data_db);
    }
  }
  else 
    SET_CLASS(pop, Rf_mkString("dempp"));
  UNPROTECT(np); // pop 
}