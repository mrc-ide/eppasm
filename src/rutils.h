#ifndef _EPPASM_RUTILS_H_
#define _EPPASM_RUTILS_H_

#include <Rinternals.h>

SEXP getListElement(SEXP list, const char *str);
int checkListElement(SEXP list, const char *str);
SEXP sexp_numeric_array(const int n, const int* dim);

#endif
