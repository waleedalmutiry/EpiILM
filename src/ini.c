// RegisteringDynamic Symbols

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R_ext/RS.h>


/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Fortran calls */
extern void F77_NAME(datacon)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(dataconsir)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(dataxy)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(dataxysir)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(like)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(likecon)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(likeconsir)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(likesir)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(rconsir)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void * );
extern void F77_NAME(rxysir)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

void F77_SUB(seedin)(void)
{
   GetRNGstate();
}


void F77_SUB(seedout)(void)
{
   PutRNGstate();
}


void F77_SUB(randomnumber)(double *x)
{
        *x = unif_rand();
}


static const R_FortranMethodDef FortranEntries[] = {
    {"datacon",    (DL_FUNC) &F77_NAME(datacon),    14},
    {"dataconsir", (DL_FUNC) &F77_NAME(dataconsir), 16},
    {"dataxy",     (DL_FUNC) &F77_NAME(dataxy),     15},
    {"dataxysir",  (DL_FUNC) &F77_NAME(dataxysir),  17},
    {"like",       (DL_FUNC) &F77_NAME(like),       16},
    {"likecon",    (DL_FUNC) &F77_NAME(likecon),    15},
    {"likeconsir", (DL_FUNC) &F77_NAME(likeconsir), 16},
    {"likesir",    (DL_FUNC) &F77_NAME(likesir),    17},
    {"rconsir",    (DL_FUNC) &F77_NAME(rconsir),    16},
    {"rxysir",     (DL_FUNC) &F77_NAME(rxysir),     17},
    {NULL, NULL, 0}
};

void R_init_EpiILM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
