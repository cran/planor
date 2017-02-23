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
#include "inverses_basep.h"

/* ++++++++++++++++++++++++++++++++++++++
findFactnum
 Return the elements of factnum in the positions where x is not zero.
factnum and x are supposed of same size
++++++++++++++++++++++++++++++++++++++++ */
IntegerVector findFactnum (IntegerVector factnum, IntegerVector x) {
  int i, nfind = 0;
  int moinsun =  (int)(-1);
  int zero = (int)0;
  IntegerVector outvec(x.size(), moinsun);

  for (i=0; i< x.size(); i++) {
    if (x[i] != zero) {
      outvec[nfind++] =factnum[i];
    }
  }

  if (nfind == zero) {
    // All values in x are zeros
    Rcpp::stop("INTERNAL ERROR findFactnum: all values in x are equal to zero");
    return(IntegerVector(0));
  } else {
    return(outvec[Range(0, nfind-1)]);
  }

} // end findFactnum


/* ++++++++++++++++++++++++++++++++++++++
attrWeight
Calculate pseudoweight, binrank, modrank
Return a vector with these three values in this order.
++++++++++++++++++++++++++++++++++++++++ */
 IntegerVector attrWeight( IntegerVector vec, int p) {
   int i;
   int pseudoweight = (int)0;
   int binrank = (int)0;
   int modrank = (int)0;
  int zero = (int)0;

   for (i=0; i< vec.size(); i++) {
     if (vec[i] != zero) {
       pseudoweight +=1;
       binrank  += (int) R_pow_di( (double)(p), i);
       modrank  += ( (int)(vec[i]) * 	
		       R_pow_di( (double)(p), (vec.size()-i-1)));
     }
   }
   IntegerVector outvec(3);
   outvec[0] = pseudoweight;
   outvec[1] =  binrank;
   outvec[2] = modrank;

   return(outvec);
 } // end attrWeight


