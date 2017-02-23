#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;

/* +++  includes from R +++ */
#include <R.h>
#include <Rmath.h>// useful for pow
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Utils.h> // to allow user interrupt
/* +++ system includes +++ */
#include <float.h>
#include <math.h> 

/* +++ local includes  +++ */
#include "util.h"

/* ++++++++++++++++++++++++++++++++++++++
inverses_basep
  # Raw calculation of the inverses modulo p
  # ARGUMENTS:
  #  p: a prime (not checked)
  # OUTPUT:
  #  vector of length p-1
++++++++++++++++++++++++++++++++++++++++ */
IntegerVector  inverses_basep( int p) {
  int un = (int)1;
  IntegerVector outvec(p-un);
  int i, dim, fone;

  if (p <= (int)3) {
    // Return (1), if p=2, or (1,2), if p=3
    // outvec <- seq(1, p-1)
    outvec = sequence(un, p-un);
  } else {
    dim=p-(int)3;
    IntegerVector vec(dim);
    vec  = sequence((int)2, p-(int)2);
    for (i=0; i<dim; i++) {
      IntegerVector onerow = vec[i] * vec;
      onerow =  moduloVec(onerow,p);
      // index of the first value equal to 1
      fone = findFOneIndex(onerow);
      outvec[i+1] = fone + un;
    } // end i
    outvec[0]=un;
    outvec[p-2] = p-un;
  } // end p>3

return(outvec);
} // end inverses_basep

