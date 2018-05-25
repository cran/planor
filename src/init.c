#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP PLANORdesignkeynonsym(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP PLANORdesignkeysym(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP PLANORinversesbasep(SEXP);
extern SEXP PLANORkernelcheck(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP PLANORprintLib(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP PLANORsubgroup(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP PLANORweightorder(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"PLANORdesignkeynonsym", (DL_FUNC) &PLANORdesignkeynonsym, 12},
    {"PLANORdesignkeysym",    (DL_FUNC) &PLANORdesignkeysym,    12},
    {"PLANORinversesbasep",   (DL_FUNC) &PLANORinversesbasep,    1},
    {"PLANORkernelcheck",     (DL_FUNC) &PLANORkernelcheck,     10},
    {"PLANORprintLib",        (DL_FUNC) &PLANORprintLib,         5},
    {"PLANORsubgroup",        (DL_FUNC) &PLANORsubgroup,         5},
    {"PLANORweightorder",     (DL_FUNC) &PLANORweightorder,      5},
    {NULL, NULL, 0}
};

void R_init_planor(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
