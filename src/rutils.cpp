#include "rutils.h"
#include <R.h>
#include <Rinternals.h>


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

int checkListElement(SEXP list, const char *str)
{
  SEXP names = getAttrib(list, R_NamesSymbol);
  for (int i = 0; i < length(list); i++ )
    if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0 )
      return 1;

  return 0;
}

SEXP sexp_numeric_array(const int n, const int* dim)
{
  int len = 1;
  SEXP s_dim = PROTECT(allocVector(INTSXP, n));
  for(int i = 0; i < n; i++)
    len *= INTEGER(s_dim)[i] = dim[i];
  SEXP s_val = PROTECT(allocVector(REALSXP, len));
  setAttrib(s_val, R_DimSymbol, s_dim);
  return s_val;
}
