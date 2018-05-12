#include <R.h>
#include <Rinternals.h>

SEXP ageprevC(SEXP s_mod, SEXP s_aidx, SEXP s_sidx, SEXP s_yidx, SEXP s_agspan) {

  double *mod = REAL(s_mod);
  int *dim = INTEGER(getAttrib(s_mod, R_DimSymbol));
  int *aidx = INTEGER(s_aidx);
  int *sidx = INTEGER(s_sidx);
  int *yidx = INTEGER(s_yidx);
  int *agspan = INTEGER(s_agspan);
  
  int n = length(s_aidx);
  SEXP s_out = PROTECT(allocVector(REALSXP, n));
  double *out = REAL(s_out);

  for(int i = 0; i < n; i++){

    double hivn = 0;
    double hivp = 0;

    int s1 = sidx[i] == 0 ? 1 : sidx[i];
    int s2 = sidx[i] == 0 ? 2 : sidx[i];

    for(int s = s1; s <= s2; s++){
      int hivnidx = (aidx[i] - 1) + (s - 1) * dim[0] + 0 * dim[0]*dim[1] + (yidx[i] - 1) * dim[0]*dim[1]*dim[2];
      int hivpidx = (aidx[i] - 1) + (s - 1) * dim[0] + 1 * dim[0]*dim[1] + (yidx[i] - 1) * dim[0]*dim[1]*dim[2];

      for(int j = 0; j < agspan[i]; j++){
        hivn += mod[hivnidx + j];
        hivp += mod[hivpidx + j];
      }
    }

    out[i] = hivp / (hivn + hivp);
  }

  UNPROTECT(1);
  return s_out;
}
