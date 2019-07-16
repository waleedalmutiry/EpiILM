

#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Fortran calls */
extern void F77_NAME(conmcmc)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(datacon)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(dataconsir)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(dataxy)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(dataxysir)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(like)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(likecon)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(likeconsir)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(likesir)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mcmc)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(rconsir)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(rxysir)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"conmcmc",    (DL_FUNC) &F77_NAME(conmcmc),    41},
    {"datacon",    (DL_FUNC) &F77_NAME(datacon),    12},
    {"dataconsir", (DL_FUNC) &F77_NAME(dataconsir), 14},
    {"dataxy",     (DL_FUNC) &F77_NAME(dataxy),     13},
    {"dataxysir",  (DL_FUNC) &F77_NAME(dataxysir),  15},
    {"like",       (DL_FUNC) &F77_NAME(like),       13},
    {"likecon",    (DL_FUNC) &F77_NAME(likecon),    12},
    {"likeconsir", (DL_FUNC) &F77_NAME(likeconsir), 13},
    {"likesir",    (DL_FUNC) &F77_NAME(likesir),    14},
    {"mcmc",       (DL_FUNC) &F77_NAME(mcmc),       42},
    {"rconsir",    (DL_FUNC) &F77_NAME(rconsir),    14},
    {"rxysir",     (DL_FUNC) &F77_NAME(rxysir),     15},
    {NULL, NULL, 0}
};

void R_init_EpiILM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
