/* +++  includes from package bigmemory +++ */
/* For access to bigmatrices  */
// NOTE: these includes should be before the others
#include "bigmemory/BigMatrix.h"
#include "bigmemory/MatrixAccessor.hpp"
/* +++  includes from R +++ */
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>
/* +++ system includes +++ */
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <math.h>
/* TYPE OF THE big matrices */
#define  TYPEOFBIG short int 

/* ++++++++++++++++++++++++++++++++++++++++++++++++
Functions to access and manipulate a  big matrix object
of type integer
+++++++++++++++++++++++++++++++++++++++++++++++++ */
/* -----------------------------------------------
 LIST OF FUNCTIONS
PLANORmodBig: modulo
PLANORcol0Big: 1 if a column entirely 0

PLANORmultBigmod:   Calculate (bigmatrix1 %*% bigmatrix2) %%mp
               or (-bigmatrix1 %*% bigmatrix2) %%mp
PLANORmult: multiplication
PLANORcbindBig, PLANORrbindBig: cbind and rbind of two big matrices
PLANORrbind1Big: rbind of a big matrix and a no-big one
PLANOR0cumsum: cumsum of a big.matrix, preceeded by one zero
--------------------------------------------------- */
/* ++++++++++++++++++++++++++++++++++++++++++++++++
MACRO
  Return 1 if a non negative  integer is zero
  Use the fonction imax2 of R.h to ensure
  exact arithmetic
 ++++++++++++++++++++++++++++++++++++++++++++++++ */
#define ISZERO(a) ( (imax2(a, 0) == 0)? 1 : 0 )


/* ++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION
  Return the modulo p of a big matrix mat
INPUT
   addressofmat: address of the big matrix mat
   p
INPUT/OUTPUT
   addressofout: address of the resulting big matrix
+++++++++++++++++++++++++++++++++++++++++++++++++ */
extern "C" {
SEXP PLANORmodBig(SEXP addressofmat, SEXP gp,
			 SEXP addressofout)
 {
  int *p = INTEGER_POINTER(gp);
  int i,j,a, pp=*p;
  // Access to the big matrices
  BigMatrix *ptrin = (BigMatrix *) (R_ExternalPtrAddr(addressofmat));
  MatrixAccessor<TYPEOFBIG> bigin(*ptrin);
  BigMatrix *ptrout = (BigMatrix *) (R_ExternalPtrAddr(addressofout));
  MatrixAccessor<TYPEOFBIG> bigout(*ptrout);
  int nrow= ptrin->nrow();
  int ncol=ptrin->ncol();

  for (i=0; i< nrow; i++) {
    for (j=0; j< ncol; j++) {
//NOTE: big matrix indexes: first the column index
      // a= ( (int) bigin[i, j]% pp);
      a= ( (int) bigin[j][i]% pp);
      /* VU difference du calcul du mod en C et en R quand <0 */
      if (a<0) {
	// On fait comme en R:
	//	bigout[i, j]= (TYPEOFBIG) (pp-abs(a));
	bigout[j][i]= (TYPEOFBIG) (pp-abs(a));
      } // fin a<0
      else
	// bigout[i, j]= (TYPEOFBIG) a; 
	bigout[j][i]= (TYPEOFBIG) a; 
    } // fin j
  } // fin i
 return(addressofout);
 } // fin PLANORmodBig

/* ++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION
 return 0 if no column of a matrix mat is entirely equal to zero
 return 1 if one column at least is entirely equal to zero
INPUT
 nrow, ncol: dimension of mat
 addressofmat: address of mat (big matrix)
OUTPUT
 res: the returned value (0 or 1)
+++++++++++++++++++++++++++++++++++++++++++++++++ */
SEXP PLANORcol0Big( SEXP gnrow, SEXP gncol,
			 SEXP addressofmat, 
			 SEXP gres)
 {
   int i,j;
  int *nrow = INTEGER_POINTER(gnrow);
  int *ncol = INTEGER_POINTER(gncol);
  int *res = INTEGER_POINTER(gres);
  BigMatrix *ptrin = (BigMatrix *) (R_ExternalPtrAddr(addressofmat));
  MatrixAccessor<TYPEOFBIG> bigin(*ptrin);

  for (j=0; j<*ncol; j++) {
    *res=1;
    for (i=0; i<*nrow;i++) {
      //      if (ISZERO(ptr[i,j]) 
  //NOTE: big matrix indexes: first the column index
     if (!ISZERO(bigin[j][i]) ) {
	*res=0;
	break; // passer a la col suivante
      }
    } // fin i
    if (*res==1) {
    //une colonne est entierement nulle
      return(gres);
    }
  } // fin j
 return(gres);
 } // fin PLANORcol0Big


/* ++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION
Calculate (bigmatrix1 %*% bigmatrix2) %%mp
or (-bigmatrix1 %*% bigmatrix2) %%mp
INPUT
 addressofbig1: address of bigmatrix1 (n,p)
 addressofbig2:  address of bigmatrix2 (p,q)
 sign:  -1 or +1 according to the sign of the first member
INPUT/OUTPUT
 addressofout: address of the resulting bigmatrix (n,q)
NOTE
 The dimensions are not checked
 ++++++++++++++++++++++++++++++++++++++++++++++++ */
  SEXP PLANORmultBigmod(SEXP addressofbig1, SEXP addressofbig2,
			SEXP gn, SEXP gp, SEXP gq,
			SEXP gmp, SEXP gsigne,
			SEXP addressofout)
  {
    int a, i, j, l;
  int *n = INTEGER_POINTER(gn);
  int *p = INTEGER_POINTER(gp);
  int *q = INTEGER_POINTER(gq);
  int *mp = INTEGER_POINTER(gmp);
  int *signe=INTEGER_POINTER(gsigne);
  // Access to the big matrices
  BigMatrix *ptrin1 = (BigMatrix *) (R_ExternalPtrAddr(addressofbig1));
  MatrixAccessor<TYPEOFBIG> big1(*ptrin1);
  BigMatrix *ptrin2 = (BigMatrix *) (R_ExternalPtrAddr(addressofbig2));
  MatrixAccessor<TYPEOFBIG> big2(*ptrin2);
  BigMatrix *ptrout = (BigMatrix *) (R_ExternalPtrAddr(addressofout));
  MatrixAccessor<TYPEOFBIG> bigout(*ptrout);

R_CheckUserInterrupt(); // permettre a l'utilisateur d'interrompre

for (i=0; i< *n; i++) {
  for (j=0; j< *q; j++) {
    //    ires= (*n) * j +i; // indice courant (=[i,j]) dans res
  //NOTE: big matrix indexes: first the column index
    //   bigout[i, j] =0;
   bigout[j][i] =0;
    for (l=0; l< *p; l++) {
	 //      res[i,j] += (big1[i,l] * big2[l, j]);
      bigout[j][i] += ((*signe)*big1[l][i] * big2[ j][l]);
    } // fin for l
    //    a= ( (int)bigout[i,j]%(*mp));
    a= ( (int)bigout[j][i]%(*mp));
    /* VU difference du calcul du mod en C et en R quand <0 */
    if (a<0) {
      // Version R
      //bigout[i,j] =  (TYPEOFBIG) (*mp-abs(a));
      bigout[j][i] =  (TYPEOFBIG) (*mp-abs(a));
    } // fin a<0
    else
      // bigout[i, j] = (TYPEOFBIG)a;
      bigout[j][i] = (TYPEOFBIG)a;
  } // fin j
 } //fin i
  return(addressofout);

} // fin PLANORmultBigmod


/* ++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION
Calculate (bigmatrix1 %*% matrix2) %%mp
or (-bigmatrix1 %*% matrix2) %%mp
INPUT
 addressofbig1: address of bigmatrix1 (n,p)
 matrix2:  matrix2 (p,q)
 sign:  -1 or +1 according to the sign of the first member
INPUT/OUTPUT
 addressofout: address of the resulting bigmatrix (n,q)
NOTE
 The dimensions are not checked
 ++++++++++++++++++++++++++++++++++++++++++++++++ */
  SEXP PLANORmultBigSmod(SEXP addressofbig1, SEXP m2,
			SEXP gn, SEXP gp, SEXP gq,
			SEXP gmp, SEXP gsigne,
			SEXP addressofout)
  {
    int a, i, j, l;
    double *matrix2;
  int *n = INTEGER_POINTER(gn);
  int *p = INTEGER_POINTER(gp);
  int *q = INTEGER_POINTER(gq);
  int *mp = INTEGER_POINTER(gmp);
  int *signe=INTEGER_POINTER(gsigne);
  PROTECT(m2 = coerceVector(m2, REALSXP));
  matrix2 = REAL(m2);

  // Access to the big matrices
  BigMatrix *ptrin1 = (BigMatrix *) (R_ExternalPtrAddr(addressofbig1));
  MatrixAccessor<TYPEOFBIG> big1(*ptrin1);
  BigMatrix *ptrout = (BigMatrix *) (R_ExternalPtrAddr(addressofout));
  MatrixAccessor<TYPEOFBIG> bigout(*ptrout);

R_CheckUserInterrupt(); // permettre a l'utilisateur d'interrompre

for (i=0; i< *n; i++) {
  for (j=0; j< *q; j++) {
    //    ires= (*n) * j +i; // indice courant (=[i,j]) dans res
  //NOTE: big matrix indexes: first the column index
    //   bigout[i, j] =0;
   bigout[j][i] =0;
    for (l=0; l< *p; l++) {
	 //      res[i,j] += (big1[i,l] * matrix2[l, j]);
      bigout[j][i] += ((*signe)*big1[l][i] * (short int)matrix2[ j * (*p) +l]);

    } // fin for l
    //    a= ( (int)bigout[i,j]%(*mp));
    a= ( (int)bigout[j][i]%(*mp));

    /* VU difference du calcul du mod en C et en R quand <0 */
    if (a<0) {
      // Version R
      //bigout[i,j] =  (TYPEOFBIG) (*mp-abs(a));
      bigout[j][i] =  (TYPEOFBIG) (*mp-abs(a));
    } // fin a<0
    else
      // bigout[i, j] = (TYPEOFBIG)a;
      bigout[j][i] = (TYPEOFBIG)a;
  } // fin j
 } //fin i
 UNPROTECT(1);

  return(addressofout);

} // fin PLANORmultBigSmod


/* ++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION
Calculate (bigmatrix1 %*% bigmatrix2)
INPUT
 addressofbig1: address of bigmatrix1 (n,p)
 addressofbig2:  address of bigmatrix2 (p,q)
OUTPUT
 addressofout: address of the resulting bigmatrix (n,q)
NOTE
 The dimensions are not checked
 ++++++++++++++++++++++++++++++++++++++++++++++++ */
  SEXP PLANORmultBig(SEXP addressofbig1, SEXP addressofbig2,
			SEXP gn, SEXP gp, SEXP gq,
			SEXP addressofout)
  {
    int  i, j, l;
  int *n = INTEGER_POINTER(gn);
  int *p = INTEGER_POINTER(gp);
  int *q = INTEGER_POINTER(gq);
 // Access to the big matrices
  BigMatrix *ptrin1 = (BigMatrix *) (R_ExternalPtrAddr(addressofbig1));
  MatrixAccessor<TYPEOFBIG> big1(*ptrin1);
  BigMatrix *ptrin2 = (BigMatrix *) (R_ExternalPtrAddr(addressofbig2));
  MatrixAccessor<TYPEOFBIG> big2(*ptrin2);
  BigMatrix *ptrout = (BigMatrix *) (R_ExternalPtrAddr(addressofout));
  MatrixAccessor<TYPEOFBIG> bigout(*ptrout);

R_CheckUserInterrupt(); // permettre a l'utilisateur d'interrompre

for (i=0; i< *n; i++) {
  for (j=0; j< *q; j++) {
    //    ires= (*n) * j +i; // indice courant (=[i,j]) dans res
//NOTE: big matrix indexes: first the column index
//     bigout[i, j] =0;
    bigout[j][i] =0;
    for (l=0; l< *p; l++) {
	 //      res[i,j] += (big1[i,l] * big2[l, j]);
      bigout[j][i] +=  (big1[l][i] * big2[j][l]);
    } // fin for l
  } // fin j
 } //fin i
  return(addressofout);

} // fin PLANORmultBig
/* ++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION
  cbind of two big.matrices
INPUT
  addressofbig1, addressofbig2: address of the two big.matrices
  nrow: number of rows of the big.matrices
  ncol1, ncol2: number of columns  of the two big.matrices
INPUT/OUTPUT
  addressofout: address of the returned matrix (previously allocated)
    cbind(big1, big2)
+++++++++++++++++++++++++++++++++++++++++++++++++++ */

  SEXP PLANORcbindBig(SEXP addressofbig1, SEXP addressofbig2,
		      SEXP addressofout, SEXP gnrow,
		      SEXP gncol1, SEXP gncol2)
  {
    int i,j;
    int *nr = INTEGER_POINTER(gnrow);
    int *nc1 = INTEGER_POINTER(gncol1);
    int *nc2 = INTEGER_POINTER(gncol2);
// Access to the big matrices
  BigMatrix *ptrin1 = (BigMatrix *) (R_ExternalPtrAddr(addressofbig1));
  MatrixAccessor<TYPEOFBIG> big1(*ptrin1);
  BigMatrix *ptrin2 = (BigMatrix *) (R_ExternalPtrAddr(addressofbig2));
  MatrixAccessor<TYPEOFBIG> big2(*ptrin2);
  BigMatrix *ptrout = (BigMatrix *) (R_ExternalPtrAddr(addressofout));
  MatrixAccessor<TYPEOFBIG> bigout(*ptrout);

  for (i=0; i<*nr; i++) {
      for (j=0;j< *nc1; j++) {
	//	res[i,j]=big1[i,j];
	//NOTE: big matrix indexes: first the column index
	bigout[j][i]=big1[j][i];
      } // fin j

      for (j=0;j< *nc2; j++) {
	//	res[i,j+nc1]=big2[i,j];
	bigout[j+ *nc1][i]= big2[j][i];
      } // fin j
    } // fin i
    return(addressofout);
  } // fin PLANORcbindBig

/* ++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION
  rbind of two big.matrices
INPUT
  addressofbig1, addressofbig2: address of the two big.matrices
  nrow1, nrow2: number of rows  of the two big.matrices
  ncol: number of columns of the big.matrices
INPUT/OUTPUT
  addressofout: address of the returned matrix (previously allocated)
    =rbind(big1, big2)
+++++++++++++++++++++++++++++++++++++++++++++++++++ */
  SEXP PLANORrbindBig(SEXP addressofbig1, SEXP addressofbig2,
		      SEXP addressofout, 
		      SEXP gnrow1, SEXP gnrow2, SEXP gncol)
  {
    int i,j;
    int *nc = INTEGER_POINTER(gncol);
    int *nr1 = INTEGER_POINTER(gnrow1);
    int *nr2 = INTEGER_POINTER(gnrow2);

  // Access to the big matrices
  BigMatrix *ptrin1 = (BigMatrix *) (R_ExternalPtrAddr(addressofbig1));
  MatrixAccessor<TYPEOFBIG> big1(*ptrin1);
  BigMatrix *ptrin2 = (BigMatrix *) (R_ExternalPtrAddr(addressofbig2));
  MatrixAccessor<TYPEOFBIG> big2(*ptrin2);
  BigMatrix *ptrout = (BigMatrix *) (R_ExternalPtrAddr(addressofout));
  MatrixAccessor<TYPEOFBIG> bigout(*ptrout);

    for (j=0; j< *nc; j++) {
      for (i=0; i< *nr1; i++) {
	//	res[i,j]=big1[i,j];
//NOTE: big matrix indexes: first the column index
	bigout[j][i] = big1[j][i];
      } // fin j
      for (i=0;i< *nr2; i++) {
	//	res[i+nr1,j]=big2[i,j];
	bigout[j][i+ *nr1] = big2[j][i];
      } // fin j
    } // fin i
 
    return(addressofout);
  } // fin PLANORrbindBig

/* ++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION
  rbind of a big.matrix and a non big matrix
INPUT
  addressofbig1: address of a big.matrix
  mat2: second, non-big matrix
  nrow1, nrow2: number of rows  of the two matrices
  ncol: number of columns of the matrices
INPUT/OUTPUT
  addressofout: address of the returned matrix (previously allocated)
  a big matrix =rbind(big1, mat2)
+++++++++++++++++++++++++++++++++++++++++++++++++++ */
  SEXP PLANORrbind1Big(SEXP addressofbig1, SEXP gmat2,
		      SEXP addressofout, 
		      SEXP gnrow1, SEXP gnrow2, SEXP gncol)
  {
    int i,j;
    int *nc = INTEGER_POINTER(gncol);
    int *nr1 = INTEGER_POINTER(gnrow1);
    int *nr2 = INTEGER_POINTER(gnrow2);
  // Access to the big matrices
  BigMatrix *ptrin1 = (BigMatrix *) (R_ExternalPtrAddr(addressofbig1));
  MatrixAccessor<TYPEOFBIG> big1(*ptrin1);
  BigMatrix *ptrout = (BigMatrix *) (R_ExternalPtrAddr(addressofout));
  MatrixAccessor<TYPEOFBIG> bigout(*ptrout);
    double *mat2 = NUMERIC_POINTER(gmat2);

    for (j=0; j< *nc; j++) {
      for (i=0; i< *nr1; i++) {
	//	res[i,j]=big1[i,j];
      //NOTE: big matrix indexes: first the column index
	bigout[j][i] = big1[j][i];
      } // fin j
      for (i=0;i< *nr2; i++) {
	//	res[i+nr1,j]=mat2[i,j];
	bigout[j][i+ *nr1] =(TYPEOFBIG)mat2[j * (*nr2)+i];
      } // fin j
    } // fin i
 
    return(addressofout);
  } // fin PLANORrbind1Big




} // fin extern C
  
