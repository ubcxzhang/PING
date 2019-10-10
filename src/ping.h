#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>

SEXP fitPING(SEXP segReadsList, SEXP paraEM, SEXP paraPrior, SEXP minReads, SEXP detailS, SEXP rescaleS, SEXP calphaS, SEXP PES);

/*SEXP getDensity(SEXP ping, SEXP strand, SEXP step, SEXP filter, SEXP sum, SEXP scale);*/
/*SEXP getDensityList(SEXP pingList, SEXP strand, SEXP step, SEXP filter, SEXP sum, SEXP scale);*/
SEXP getListElement(SEXP list, const char *str);
