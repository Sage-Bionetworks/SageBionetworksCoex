#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include "dsyrk.h"

static const R_CallMethodDef callMethods[] = {
    {"aat", (DL_FUNC) &aat, 1},
    {NULL, NULL, 0}
};

void
R_init_SageBionetworksCoex(DllInfo *info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

void
R_unload_SageBionetworksCoex(DllInfo *info)
{
    /* any clean-up when package unloaded */
}
