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


/* ++++++++++++++++++++++++++++++++++++++++++++++++ 
crossing
 Generates all n1 x n2 x ... x ns combinations of size nbg 
 with n1,...,ns integers
 Here, the values n1,...,ns are all identical and equal to p
INPUT
  p: the value of the series of integers
  nbg
  N=p^nbg
  start: integer from where to start the series of integers
OUTPUT
  crosses: an integer matrix with p^nbg rows and nbg columns 
  giving all combinations along the rows, in lexicographic order.
CALLED BY
subgroup
 ++++++++++++++++++++++++++++++++++++++++++++++++ */

void crossing( int p,  int nbg, int N,
			short int start, 
		     IntegerMatrix crosses) 
{

  int  i,l,k,ik, il, nbgi, imotif, iunmotif;

  Rcpp::checkUserInterrupt(); 

  IntegerVector motif= int(start) - 1 +seq_len(p);

 imotif = N/p; // number of times the pattern is repeated
 iunmotif = (int)1 ; // number of times each element of the pattern is repeated
 
 for(i=1; i<=nbg; i++) {
   l=(int)0;
   nbgi=nbg-i; // indice of the current column
     for (ik=0; ik<imotif; ik++) {
       for (k=0; k<p; k++) {
	 for (il=0; il< iunmotif; il++) {
	   if (l>=N) {
	     Rcpp::stop("INTERNAL ERROR crossing: l >=N \n");
	     return;
	   }
	   crosses(l, nbgi) = motif[k];
	   l++;
	 } // end il (end repetition of one element of the pattern)
       } // end k (end one element of the pattern)
     }   // end ik (end repetition of the pattern)
     // next column
     imotif = imotif/p;
     iunmotif = iunmotif*p;
 } // end i (ned of the comumn)

} // end crossing



/* ++++++++++++++++++++++++++++++++++++++
subgroup
++++++++++++++++++++++++++++++++++++++++ */

IntegerMatrix subgroup(IntegerMatrix matrice,
	      int p,
	      bool all = true) {

int nbg=matrice.ncol();
int N=R_pow_di(p, nbg);
 int zero = (int)0;

 IntegerMatrix coeffs(N, nbg);
 crossing(p, nbg, N, zero, coeffs);

 //  invert the columns
 IntegerVector ind=rev(seq_len(nbg)) - (int)1 ;
 coeffs=takeCol(coeffs, ind);

// remove the first line
 coeffs = coeffs(Range(1,  coeffs.nrow()-1), Range(0, coeffs.ncol()-1));

 if (all== false) {
   // we kept only the rows where the first non-zero element is  1
   LogicalVector indices = findFnzOne(coeffs);
   coeffs=takeRowLog(coeffs, indices);
 } // end all

 //  outmat <- mat %*% t(coeffs)%%p
 IntegerMatrix outmat = multMatMod(matrice, coeffs,p );

 return(outmat);

} // end subgroup
