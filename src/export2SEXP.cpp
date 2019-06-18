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
#include "Classes.hpp"

void popC::finalize (hivC& hivpop, artC& artpop) {
  int np = 0;

  SEXP pop_sexp_dim = PROTECT(NEW_INTEGER(4)); ++np;
  INTEGER(pop_sexp_dim)[0] = pAG;
  INTEGER(pop_sexp_dim)[1] = NG;
  INTEGER(pop_sexp_dim)[2] = pDS;
  INTEGER(pop_sexp_dim)[3] = PROJ_YEARS;
  SET_DIM(pop_sexp, pop_sexp_dim);

  if (MODEL!=0) {
    SEXP age_sex_year_dim = PROTECT(NEW_INTEGER(3)); ++np;
    INTEGER(age_sex_year_dim)[0] = pAG;
    INTEGER(age_sex_year_dim)[1] = NG;
    INTEGER(age_sex_year_dim)[2] = PROJ_YEARS;
    SET_DIM(infections_sexp, age_sex_year_dim);
    SET_DIM(hivdeaths_sexp, age_sex_year_dim);
    SET_DIM(natdeaths_sexp, age_sex_year_dim);
    SET_DIM(popadjust_sexp, age_sex_year_dim);

    SEXP hiv_sexp_dim = PROTECT(NEW_INTEGER(4)); ++np;
    INTEGER(hiv_sexp_dim)[0] = hDS;
    INTEGER(hiv_sexp_dim)[1] = hAG;
    INTEGER(hiv_sexp_dim)[2] = NG;
    INTEGER(hiv_sexp_dim)[3] = PROJ_YEARS;
    SET_DIM(hivpop.hiv_sexp, hiv_sexp_dim);

    SEXP art_sexp_dim = PROTECT(NEW_INTEGER(5)); ++np;
    INTEGER(art_sexp_dim)[0] = hTS;
    INTEGER(art_sexp_dim)[1] = hDS;
    INTEGER(art_sexp_dim)[2] = hAG;
    INTEGER(art_sexp_dim)[3] = NG;
    INTEGER(art_sexp_dim)[4] = PROJ_YEARS;
    SET_DIM(artpop.art_sexp, art_sexp_dim);

    SET_ATTR(pop_sexp, Rf_install("infections"), infections_sexp);
    SET_ATTR(pop_sexp, Rf_install("hivdeaths"), hivdeaths_sexp);
    SET_ATTR(pop_sexp, Rf_install("natdeaths"), natdeaths_sexp);
    SET_ATTR(pop_sexp, Rf_install("popadjust"), popadjust_sexp);
    SET_ATTR(pop_sexp, Rf_install("pregprevlag"), pregprevlag_sexp);
    SET_ATTR(pop_sexp, Rf_install("rvec_ts"), rvec_sexp);
    SET_ATTR(pop_sexp, Rf_install("prev15to49"), prev15to49_sexp);
    SET_ATTR(pop_sexp, Rf_install("incid15to49"), incid15to49_sexp);
    SET_ATTR(pop_sexp, Rf_install("entrantprev"), entrantprev_sexp);
    SET_ATTR(pop_sexp, Rf_install("incrate15to49_ts"), inci15to49_ts_sexp);
    SET_ATTR(pop_sexp, Rf_install("prev15to49_ts_sexp"), prev15to49_ts_sexp);
    SET_ATTR(pop_sexp, Rf_install("artpop"), artpop.art_sexp);
    SET_ATTR(pop_sexp, Rf_install("hivpop"), hivpop.hiv_sexp);
    SET_CLASS(pop_sexp, Rf_mkString("spec"));
    if (MODEL==2) {
      SEXP pop_db_sexp_dim = PROTECT(NEW_INTEGER(4)); ++np;
      INTEGER(pop_db_sexp_dim)[0] = pDB;
      INTEGER(pop_db_sexp_dim)[1] = NG;
      INTEGER(pop_db_sexp_dim)[2] = pDS;
      INTEGER(pop_db_sexp_dim)[3] = PROJ_YEARS;
      SET_DIM(data_db_sexp, pop_db_sexp_dim);
      SET_ATTR(pop_sexp, Rf_install("debut_pop"), data_db_sexp);
    }
  }
  else 
    SET_CLASS(pop_sexp, Rf_mkString("dempp"));
  UNPROTECT(np + 15); // pop 
}