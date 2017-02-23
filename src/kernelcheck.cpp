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
kernelcheck_basep
RETURN
a logical vector of length admissible.ncol
CALLED BY
PLANORkernelcheck, designkeysym, designkeynonsym
++++++++++++++++++++++++++++++++++++++++ */
LogicalVector kernelcheck_basep(IntegerMatrix PhiStar, 
				IntegerMatrix admissible, 
				IntegerMatrix IneligibleSet, 
				int p) {



  // b.ImagesIS <- (- b.PhiStar %*% b.IneligibleSet)%%p 
  int k,j,l;
  bool trouve;
  int nb_admissible = admissible.ncol();
  int nb_ineligible = IneligibleSet.ncol();

  /* Output */
 LogicalVector test(nb_admissible, true);



 // b.ImagesIS <- (- b.PhiStar %*% b.IneligibleSet)%%p

// Note: IntegerMatrix may induce conversion arma error 
 NumericMatrix ImagesIS= wrap(as<arma::mat>(PhiStar) * as<arma::mat>(IneligibleSet));
 
ImagesIS= moduloMatNum(ImagesIS, p);

R_CheckUserInterrupt(); // allow user interrupt


  for (k=0; k <  nb_admissible; k++) {
    for (j=0; j < nb_ineligible; j++) {
      trouve = true;
      for (l=0; l < admissible.nrow(); l++) {
	// VOIR	if ( (-ImagesIS(l,j)) !=admissible(l,k)) {
         if ( int(ImagesIS(l,j)) !=admissible(l,k)) {
	  trouve = false;
	  break ; // go out the loop on l
	}
      } // end for l
      if (trouve == true) {
	test[k] = false;
	break; // go out the loop on j
      }
    } // end for j
  } // end for k


  return(test);
} // end kernelcheck_basep


