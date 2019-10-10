#include <R.h>
#include <R_ext/Rdynload.h>
#include "ping.h"


static const R_CallMethodDef CallEntries[] = {
    {"fitPING", (DL_FUNC)&fitPING, 8},
//    {"getDensity", (DL_FUNC)&getDensity, 6},
//    {"getDensityList", (DL_FUNC)&getDensityList, 6},

    {NULL, NULL, 0}
};

static const R_CMethodDef CEntries[] = {
    {NULL, NULL, 0}
};

void R_init_PING(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);

}





