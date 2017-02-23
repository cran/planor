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


/* ++++++++++++++++++++++++++++++++++++++
sequence
Generate first, first+1, first+2, ... last
(simplified version of iota(vec.begin(), vec.end(), p)
++++++++++++++++++++++++++++++++++++++++ */
  IntegerVector sequence( int first, int last) {
  IntegerVector outvec(last-first+1);
  int i, a;

  i= (int)0;
  for (a=first; a<=last; a++) {
    outvec[i++] = a;
  }
  return(outvec);
 } // end sequence
/* ++++++++++++++++++++++++++++++++++++++
Rmodulo
Return x % y, calculated as R does
x%y=x - y*(the greatest integer inf or equal to x/y)
When x or y is negative, the result is different from the C operator %
++++++++++++++++++++++++++++++++++++++++ */
int Rmodulo(int x, int y) {
double ra= double(x)/double(y);
 int res = x - y * floor(ra);
 return(res);
}


/* ++++++++++++++++++++++++++++++++++++++
moduloVec
Return x % p, x is a vector
++++++++++++++++++++++++++++++++++++++++ */
IntegerVector moduloVec (IntegerVector x, const int p) {
  int i;
  IntegerVector outvec(x.size());

  for (i=0; i< x.size(); i++) {
    //    outvec[i] = x[i] % p;
    outvec[i] = Rmodulo( x[i], p);
  }

  return (outvec);
}

/* ++++++++++++++++++++++++++++++++++++++
modulo
NOTE : le modulo des nbres <0 est différent en R et C
    ### (-15)%%300 = -15 en C, et est 285 en R
    ### On adopte la solution R
++++++++++++++++++++++++++++++++++++++++ */
/* DOES NOT WORK template ne marche pas utilisé avec les
types IntegerMatrix ou Numericmatrix
template <typename T>
T moduloMatT( T mat, const int p) {
int i,j;


  for (i=0; i< mat.nrow(); i++) {
    for (j=0; j<mat.ncol(); j++) {
      //     mat(i,j) = int(mat(i,j)) % p;
          mat(i,j) = Rmodulo(int(mat(i,j)), p);
    }
  }

  return(mat);
} // end moduloMatT
**************** */

IntegerMatrix moduloMat( IntegerMatrix mat, const int p) {
int i,j;


  for (i=0; i< mat.nrow(); i++) {
    for (j=0; j<mat.ncol(); j++) {
      //     mat(i,j) = int(mat(i,j)) % p;
          mat(i,j) = Rmodulo(int(mat(i,j)), p);
    }
  }

  return(mat);
} // end moduloMat

NumericMatrix moduloMatNum (NumericMatrix mat, const int p) {
int i,j;


  for (i=0; i< mat.nrow(); i++) {
    for (j=0; j<mat.ncol(); j++) {
      //     mat(i,j) = int(mat(i,j)) % p;
          mat(i,j) = Rmodulo(int(mat(i,j)), p);
    }
  }

  return(mat);
} // end moduloMatNum

/* ++++++++++++++++++++++++++++++++++++++
multTrans
 return matrix= mat %*% trans(coeff)
using Armadillo (matricial multiplication)
++++++++++++++++++++++++++++++++++++++++ */
arma::mat multTrans(arma::mat matrice, arma::mat coeff) {
  arma::mat result = matrice * arma::trans(coeff);
  return(result);

}
/* ++++++++++++++++++++++++++++++++++++++
cbind
++++++++++++++++++++++++++++++++++++++++ */
IntegerMatrix cbind(IntegerMatrix mat1, IntegerMatrix mat2) {
  IntegerMatrix output(mat1.nrow(), mat1.ncol()+mat2.ncol());
  int j, l;

  for ( j=0; j<mat1.ncol(); j++) {
    output(_,j) = mat1(_,j);
  }
  l=mat1.ncol();
  for ( j=0; j<mat2.ncol(); j++) {
    output(_,l++) = mat2(_,j);
  }

  return(output);
} //end cbind


/* ++++++++++++++++++++++++++++++++++++++
takeElt
 return vec with the elements whose numbers are stored in indices
Careful: indices are supposed numbers from 0
++++++++++++++++++++++++++++++++++++++++ */

IntegerVector takeElt(IntegerVector vec,
	      NumericVector indices) {
IntegerVector output(vec.size());


for (int j=0; j<indices.size(); j++) {
  output[j]= vec[(int)indices[j]];
    }
return(output);
} // end takeElt

/* ++++++++++++++++++++++++++++++++++++++
takeCol
 return mat with only the columns whose numbers are stored in indices
Careful: indices are supposed numbers from 0
++++++++++++++++++++++++++++++++++++++++ */

IntegerMatrix takeCol(IntegerMatrix mat,
	      IntegerVector indices) {
IntegerMatrix output(mat.nrow(),indices.size());

for (int j=0; j<indices.size(); j++) {
    output(_, j)= mat(_, indices[j]);
    }


return(output);
} // end takeCol

/* ++++++++++++++++++++++++++++++++++++++
takeColLog
 return mat with only the columns for which indices is true
++++++++++++++++++++++++++++++++++++++++ */

IntegerMatrix takeColLog(IntegerMatrix mat,
	      LogicalVector indices) {
  IntegerMatrix output(mat.nrow(),indices.size());
  int zero = (int)0;
  int l = zero;

for (int j=0; j<indices.size(); j++) {
  if (indices[j]== true) {
    output(_, l++)= mat(_, j);
  }
 } // end j


  if (l==zero) {
    return(IntegerMatrix(zero, zero));
  } else {
    return(output(_,  Range(0, l-1)));
  }
} // end takeColLog


/* ++++++++++++++++++++++++++++++++++++++
takeRowLog
 return mat with only the rows for which indices is true
++++++++++++++++++++++++++++++++++++++++ */

IntegerMatrix takeRowLog(IntegerMatrix mat,
	      LogicalVector indices) {
  IntegerMatrix output(indices.size(), mat.ncol());
  int l = (int)0;
  int zero = (int)0;

 
for (int j=0; j<indices.size(); j++) {
  if (indices[j]) {
    output(l++, _)= mat( j, _);
  }
 }
  if (l==zero) {
    return(IntegerMatrix(0,0));
  } else {
    return(output(Range(0,l-1), _));
  }
} // end takeRowLog

/* ++++++++++++++++++++++++++++++++++++++
findIndexLog
Return the indices of the true or false elements in vec
Careful: indices are  numbers from 0
++++++++++++++++++++++++++++++++++++++++ */

IntegerVector findIndexLog(LogicalVector vec, bool val) {
  IntegerVector output(vec.size());
  int l = (int)0;
  int zero = (int)0;

  for (int j=0; j<vec.size(); j++) {
    if (vec[j] == val) {
      output[l++] = j;
    }
  }
  if (l==zero) {
     return(IntegerVector(0));
  } else {
    return(output[Range(0,l-1)]);
  }
} // end findIndexLog


/* ++++++++++++++++++++++++++++++++++++++
findIndexVal
Return the indices of the elements in vec equal to s
Careful: indices are  numbers from 0
Called by designkeynonsym, designkeysym
++++++++++++++++++++++++++++++++++++++++ */

IntegerVector findIndexVal(IntegerVector vec, int val) {
  IntegerVector output(vec.size());
  int l = (int)0;
  int zero = (int)0;

  for (int j=0; j<vec.size(); j++) {
    if (vec[j] == val) {
      output[l++] = j;
    }
  }
  if (l==zero) {
    return(IntegerVector(0));
  } else {
    return(output[Range(0,l-1)]);
  }
} // end findIndexVal


/* ++++++++++++++++++++++++++++++++++++++
findFnzOne
Search for the first non-zero element on each row
of a matrix and return  true if this first non-zero element is 1,
false if not. 
Called by subgroup
++++++++++++++++++++++++++++++++++++++++ */
LogicalVector findFnzOne(IntegerMatrix mat) {
  LogicalVector indices(mat.nrow(),  NA_INTEGER);
  int i,j;

  for (i=0; i< mat.nrow(); i++) {
    for (j=0;j<mat.ncol(); j++) {
      if (mat(i,j) != (int)0) {
	if (mat(i,j) == (int)1) {
	  indices[i]= true;
	} else {
	  // the first nz is not 1
	  indices[i]= false;
	}
	break; // exit from j-loop
      } 
    } // end j
  } // end i

  if (any(is_na(indices))) {
      Rcpp::stop("ERROR findfnz: at least a row is nul\n");
    }
  return(indices);
} // end findFnzOne


/* ++++++++++++++++++++++++++++++++++++++
findFnzIndex 
 Return the index of the  first non-zero element in vec
Index from 0;
Return -1 if none element is equal to 0
++++++++++++++++++++++++++++++++++++++++ */
int findFnzIndex (IntegerVector vec) {
  int i, fnz = (int)-1;
  int zero = (int)0;

  for (i=0; i< vec.size(); i++) {
    if (vec[i] != zero) {
      fnz =i;
      break;
    }
  }
  return(fnz);
} // end findFnzIndex


/* ++++++++++++++++++++++++++++++++++++++
findFOneIndex 
 Return the index of the first element equal to 1 in vec
Index from 0; -1 if none element is equal to 1
++++++++++++++++++++++++++++++++++++++++ */
int findFOneIndex (IntegerVector vec) {
  int i, fnz = (int)-1;
  for (i=0; i< vec.size(); i++) {
    if (vec[i] == (int)1) {
      fnz =i;
      break;
    }
  }
  return(fnz);
} // end findFOneIndex

/* ++++++++++++++++++++++++++++++++++++++
multMatMod
Compute   outmat <- mat %*% t(coeffs)%%p
++++++++++++++++++++++++++++++++++++++++ */
IntegerMatrix multMatMod(IntegerMatrix matrice,
			 IntegerMatrix coeffs,
			 int p) {

 // 1/  outmat <- mat %*% t(coeffs)
 IntegerMatrix outmat= wrap(multTrans(as<arma::mat>(matrice), as<arma::mat>(coeffs)));

 // 2/ outmat <- outmat modulo p
 outmat= moduloMat(outmat,  p);
 return(outmat);
} // end multMatMod



/* ++++++++++++++++++++++++++++++++++++++
fillMat
Fill a matrix with the values of a vector 
The vector is for a R-matrix  argument 
++++++++++++++++++++++++++++++++++++++++ */
IntegerMatrix fillMat(IntegerVector mat, int nr, int nc) {
  int i,j, l = (int)0;
 IntegerMatrix outmat(nr, nc);
 for ( j=0; j< nc;j++) {
   for (i=0; i< nr; i++) {
     outmat(i,j)= mat[l++];
   }
 }
 return(outmat);
} // end fillMat




