/* +++  includes from package bigmemory +++ */
/* For access to bigmatrices  */
// NOTE: these includes should be before the others
#include "bigmemory/BigMatrix.h"
#include "bigmemory/MatrixAccessor.hpp"
/* +++  includes from R +++ */
#include <R.h>
#include <Rmath.h>// pour pow
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Utils.h> // pour permettre d'interrompre
/* +++ system includes +++ */
#include <float.h>
#include <math.h> // pour pow
/* TYPE OF THE big matrices */
#define  TYPEOFBIG short int 
extern "C" {
/* ++++++++++++++++++++++++++++++++++++++++++++++++
MACRO
  Return 1 if two integers are equal
  Use the fonctions imax2 and imin2 of R.h to ensure
  exact arithmetic
 ++++++++++++++++++++++++++++++++++++++++++++++++ */
#define EQUAL(a, b) (((imax2(a,b) == a) && (imin2(a,b) == a)) ? 1 : 0 )
/* ++++++++++++++++++++++++++++++++++++++++++++++++
MACRO
  Return 1 if a non negative  integer is zero
  Use the fonction imax2 of R.h to ensure
  exact arithmetic
 ++++++++++++++++++++++++++++++++++++++++++++++++ */
#define ISZERO(a) ( (imax2(a, 0) == 0)? 1 : 0 )
/* ++++++++++++++++++++++++++++++++++++++++++++++++
MACRO
  Return 1 if a non negative  integer is 1
  Use the fonction imax2 of R.h to ensure
  exact arithmetic
 ++++++++++++++++++++++++++++++++++++++++++++++++ */
#define ISUN(a) ( (imax2(a, 1) == 1)? 1 : 0 )
/* ++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION
 Access to a component by name from a R-list
INPUT
 list: a R-list
 str: name of a component of list
RETURN VALUE
 the component "str" of list
 ++++++++++++++++++++++++++++++++++++++++++++++++ */
SEXP getListElement(SEXP list, const char *str)
{
  SEXP elmt = R_NilValue, names=getAttrib(list, R_NamesSymbol);
  int i;
  for (i=0; i< length(list); i++)
    if (strcmp(CHAR(STRING_ELT(names, i)), str) ==0) {
      elmt=VECTOR_ELT(list, i);
      break;
    }
  return(elmt);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION
 Multiplication of a given column, nc, of the big matrix A
 by the column nc of a matrix b. This column is stored in B.
 If iz[i]==pp1, B[i] is ignored.
 iz[i]=pp1, if the first nonzero value of the row i of t(b)
 is not 1
INPUT
 addressofA: address of big matrix A (nrow, nbg)
 B: vector (N1)
 iz: vector (N1)
 all: see subgroup.base
 nc: current row of A
 nrow, N1: dimensions
 pp1: marker of values to ignored 
OUTPUT
 addressofres:  address of resulting big matrix (nrow,N1)
CALLED BY
 PLANORsubgroup
NOTE
 The dimensions are not checked
NOTE AB:
 En fait, on doit faire mat %*% t(coeffs)
 Mais, comme on calcule coeffs, colonne par colonne,
 on ne connait qu'une ligne de t(coeffs)
 Les elements de cette ligne, la ligne nc, doivent etre multipli�s
 par la colonne correspondante de mat, soit la colonne nc
 Tous les elements de res sont modifi�s. Au prochain appel,
 B correspondra � la ligne suivante de  t(coeffs): tous les
  elements de res seront � nouveau modifi�s.
 ++++++++++++++++++++++++++++++++++++++++++++++++ */

  void multpp1col (SEXP addressofa, short int *B,
		   short int *iz,
		   int all, int nc, int nrow,  int N1,
		   short int pp1,
		   SEXP addressofres)
{
  int l, r;
  /* Access to the values of the big matrices */
  BigMatrix *ptrin = (BigMatrix *) (R_ExternalPtrAddr(addressofa));
  MatrixAccessor<TYPEOFBIG> bigA(*ptrin);
  BigMatrix *ptrout = (BigMatrix *) (R_ExternalPtrAddr(addressofres));
  MatrixAccessor<TYPEOFBIG> bigres(*ptrout);

R_CheckUserInterrupt(); // permettre a l'utilisateur d'interrompre

 for (r=0; r< N1;r++) {
     // calculer iz
     /* on ne d�tecte les lignes � �liminer que lorsqu'on
rencontre un elt >1: comme on balaie colonne par colonne,
en commen�ant par la premi�re,
on a d�j� vus les elts pr�c�dents sur la ligne.
Ils sont soit nuls, et c'est comme si on les ignorait, 
soit �gaux � 1, et alors, on a bien fait de les consid�rer */
   if (all==0) {
     // calculer iz
     /* on ne d�tecte les lignes � �liminer que lorsqu'on
rencontre un elt >1 : les elts pr�c�dents sur la ligne ont donc �t�
pris en compte; mais ce n'est pas grave: soit ils sont nuls,
et c'est comme si on les ignorait, soit ils sont �gaux � 1,
et alors, il faut consid�rerla ligne */
     if ((iz[r]==0) && (B[r]==1))
	 iz[r]=1;
       if ((iz[r]!=1) && (B[r]!=0) && (B[r]!=1))
	   iz[r]=pp1;
   } // fin all


   if (iz[r] != pp1) {
     // ligne pas a ignorer
     for (l=0; l <nrow; l++) {
      //      res[l,r]+= A[l,nc]*B[r]
       //NOTE: big matrix indexes: first the column index
       bigres[r][l] += (bigA[nc][l] * B[r]);
     } // fin l
   } // fin if iz
     else {
       // ligne de coeffs a ignorer donc colonne de B
       // res[0,r] =-1;
       bigres[r ][0] =-1;
     }
 } //fin r

} // fin multcolpp1col


/* ++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION
 Repeat the value motif, prod1 times in crosses, 
 from the position l included.
 Indexes begin from zero.
INPUT
  motif: the value to repeat;
  prod1: the number of times, motif should be repeated;
  l: index of the beginning 
INPUT-OUTPUT
  crosses: in output, crosses[l:(l+prod1-1)] contain motif
CALLED BY
 gcrossingcol
 ++++++++++++++++++++++++++++++++++++++++++++++++ */
void     repetcol( int motif, int prod1, int l, 
	        short int *crosses)
{
  int i;
R_CheckUserInterrupt(); // permettre a l'utilisateur d'interrompre

  for (i=l; i<(l+prod1); i++) {
    //    crosses[i, nc] = motif;
        crosses[i ] = motif;
  }
} //fin repetcol




/* ++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION
 Generates all n1 x n2 x ... x ns combinations of size s 
 with n1,...,ns integers
 Here, the values n1,...,ns are all identical and equal to p
 (NOTE AB: Remplace la fonction R crossing dans subgroup,
o� les valeurs de n sont rep(p,nbg) donc toutes egales a p)
INPUT
  nrow, ncol: dimension of the output matrix crosses
  p: the value of the series of integers
  start: integer from where to start the series of integers
OUTPUT
  crosses: an integer matrix with nrows rows and ncol columns 
  giving all combinations along the rows, in lexicographic order.
  (NOTE AB: l'ordre des colonnes est invers� par rapport a l'original)
  This program is called for each column, so the ouput is not
  the complete matrix, but one column
CALLED BY
 PLANORsubgroup

 ++++++++++++++++++++++++++++++++++++++++++++++++ */

void gcrossingcol( int nrow, int p,  int prod1,
		      short int start,
		     short int *crosses) 
{
  int  k, l, motif;
R_CheckUserInterrupt(); // permettre a l'utilisateur d'interrompre

    motif= start;

    for (k=0; k < p; k++) {
      l=k*prod1;
      while ( (l+prod1) <= nrow) {
      /* R�p�ter motif, prod1 fois dans les elements
	 de  crosses, a partir du li�me  */
	  repetcol(motif, prod1, l, crosses);
	l += ( p * prod1);
      }
      motif++;
    } // fin k
} // fin gcrossingcol




/* ++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION
  Calculate the non null elements of the subgroup H generated 
  by the columns of mat, considered as vectors in (Zp)^s
INPUT
  nrow, nbg: dimension of mat
  N: number of columns of res (= p^nbg)
  addressofmat: address of mat:
    a big matrix of integers modulo p whose columns are assumed to
       be independent vectors in (Zp)^s
  p: a prime
  all: if TRUE all elements in H are given, if FALSE only elements up
       to multiplication by an integer modulo p are given
OUTPUT:
 addressofres:  resulting big matrix:
   a matrix of integers modulo p whose columns are the subgroup 
   elements
   dimension (nrow, N-1)
DETAILS:
  it is not checked whether the column vectors of mat are independent, so there
  may be several times 
CALLED BY
 subgroup.basep
 ++++++++++++++++++++++++++++++++++++++++++++++++ */


SEXP PLANORsubgroup( SEXP gnrow, SEXP gnbg, SEXP gN,
		     SEXP addressofmat, SEXP gp, SEXP gall,
		     SEXP addressofres)
 {
   short int  *coeffs, start, pp1, *iz;
   int prod1, nc,i, j;
  int *nrow = INTEGER_POINTER(gnrow);
  int *nbg = INTEGER_POINTER(gnbg);
  int  *N  = INTEGER_POINTER(gN);
  int *p = INTEGER_POINTER(gp);
  int *all = INTEGER_POINTER(gall);


  start=0;  
  pp1= (short int) (*p); // marqueur des lignes a ignorer

   // vn <- rep(p,nbg)
   // N <- prod(vn)
   //N= (int) R_pow_di( (double)(*p), (*nbg));

   /* coeffs <- matrix(0, nrow=N, ncol=nbg),
      mais, ici on la calcule colonne par colonne */

  coeffs = (short int *) R_alloc((*N), sizeof(short int));
  iz = (short int *) S_alloc((long) (*N - 1), sizeof(short int));
//iz: initialis� � z�ro

  prod1=1;

  // Calcul de la colonne nc de coeffs
  for (nc=0; nc< *nbg; nc++) {
    gcrossingcol((*N),  (*p), prod1, start, coeffs);
    prod1 *= (*p);
 
   /* Calcul matriciel en ne consid�rant que les elts
      de coeffs/iz!= pp1   */
   /* La 1iere ligne de coeffs est tjrs ignor�e
c'est pourquoi on d�marre � coeffs+1, et on ne consid�re
que N-1 elts */
   multpp1col (addressofmat, (coeffs+1), iz, *all,
	       nc, (*nrow), (*N-1), pp1, addressofres);
  } // fin nc

    /* Faire modulo p sur le r�sultat */
  // NOTE AB: on ne le fait pas dans la foul�e lorsqu'on calcule
  // res car celui-ci est une accumulation progressive de valeurs
  // dont le r�sultat final n'est d�termin� qu'en fin de boucle 
  //Access to the resulting big matrix
   BigMatrix *ptrout = (BigMatrix *) (R_ExternalPtrAddr(addressofres));
  MatrixAccessor<TYPEOFBIG> bigres(*ptrout);

  for (i=0; i< ptrout->nrow(); i++) {
    for (j=0; j<  ptrout->ncol() ;j++) {
       //NOTE: big matrix indexes: first the column index
       //      a= ( (int) bigres[i,j]% (*p));
      if (bigres[j][i] >0) {
	// Les valeurs nulles doivent etre ignor�es
	//	bigres[i,j]= bigres[i,j]%p
	bigres[j][i]= (TYPEOFBIG) ( (int) bigres[j][i]% (*p)); 
      }
    } // fin j
  } // fin i

  return(addressofres);
} // fin PLANORsubgroup




/* ++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION
 Ouput function.
 Replace the time consuming R commands:
    for(j in 1:ncol(H)){
      select.kj <- H[,j] != 0
      coeff.char <- paste("^",as.character(H[select.kj,j]), sep="")
      coeff.char[H[select.kj,j] == 1] <- ""
      LIBS.k <- paste(LIBtpf.k[select.kj], coeff.char, sep="")
      cat("1 = ") ; cat(LIBS.k,sep=" ") ; cat("\n")
    }

INPUT
  nrow, ncol: dimension of H
  addressofH: address of H: a big matrix
  LIBtpf.k: character vector
  MAXPRINT: maximum number of rows and columns to print
OUTPUT
  Nothing, but print on standard output
CALLED BY
  The R function:  summary.designkey
 ++++++++++++++++++++++++++++++++++++++++++++++++ */

SEXP PLANORlibsk(SEXP gnrow, SEXP gncol,
		 SEXP addressofH, SEXP gLIBtpf,
		 SEXP gMAXPRINT)
{
  int i,j;
  int *nrow = INTEGER_POINTER(gnrow);
  int *ncol = INTEGER_POINTER(gncol);
  int *maxprint = INTEGER_POINTER(gMAXPRINT);
  // Nbre max de colonnes a ip, pour limiter les ip
  // imin2: fonction de R.h
  int br = imin2(*nrow, *maxprint);
  int bc = imin2(*ncol, *maxprint);

R_CheckUserInterrupt(); // permettre a l'utilisateur d'interrompre
// Access to the big matrix H
  BigMatrix *ptrin = (BigMatrix *) (R_ExternalPtrAddr(addressofH));
  MatrixAccessor<TYPEOFBIG> bigH(*ptrin);

  for (j=0; j< bc; j++) {
    Rprintf("1 =") ; 
    for (i=0; i < br; i++) {

      //      if (H[i,j] ==0) continue;
      //      if (H[j * (*nrow) + i] ==0) continue;
      // que l'on remplace, pour eviter les erreurs d'arrondis par:
      // (NB: H contient des entiers positifs)
       //NOTE: big matrix indexes: first the column index
      //      if (ISZERO(bigH[i,j])) continue;
      if (ISZERO(bigH[j][i])) continue;
   Rprintf(" %s", CHAR(STRING_ELT(gLIBtpf,i)));

      //      if (H[i,j] !=1) 
      if (!ISUN(bigH[j][i])) {
	Rprintf("^%d ", bigH[j][i]);
      } 

    } //fin i
    Rprintf("\n");
  } // fin j
  return(addressofH);

} // fin PLANORlibsk


/* ++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION
 Return the index of the multiple of k
 such as the remainder by p is 1
INPUT
  p: integer
  k: integer [2, p-2]
OUTPUT
  inv: integer
CALLED BY
  PLANORinv 
EXAMPLE
p=11
The values successively returned when k=1...10 are:
 1  6  4  3  9  2  8  7  5 10
 ++++++++++++++++++++++++++++++++++++++++++++++++ */
int prodk(int p, int k)
{
int l, a, prodkv;
div_t d;
  a=k;

  for (l=2; l <= (p-2); l++)
    {
      prodkv = a *l;
      d= div(prodkv, p);
      if (d.rem ==1) 
	return(l);
    } // fin l
error("Internal error: prodk\n");
 return(-1);
}


/* ++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION
 Raw calculation of the inverses of the value J modulo p
NOTE AB: remplace inverses.basep
   PLANORinv(p,J) = inverses.basep(p)[J]
INPUT
  p: integer = a prime (not checked)
  J: integer = index of the value to calculate (>0 and <p)
OUTPUT
  inv: integer
CALLED BY
  kernelmatrix.basep, PLANORweightorder
 ++++++++++++++++++++++++++++++++++++++++++++++++ */
void PLANORinv (int *p, int *J,
		int *inv)
{
  //  if ( (*J==1) || (*J==(*p-1))) {
  if ( EQUAL(*J, 1) || EQUAL(*J, (*p-1))) {
    *inv = *J;
  } else {
    *inv =  prodk(*p, *J);
  }
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION
  Return 1 if a value X is present in the vector factvu
  and 0 otherwise. In the case it is not present,
  insert X into factvu at the last position.
  This position, ifv, is then incremented.
INPUT
  X: value to search; integer scalar
INPUT-OUTPUT
  ifv: current index in factvu; integer scalar (init by 0)
  factvu: the values already seen
CALLED BY
  Called by PLANORcalcweight  
 ++++++++++++++++++++++++++++++++++++++++++++++++ */
int pasvu(int X, int *ifv, int *factvu)
{
  int i;
  for (i=0; i <= *ifv; i++) {
// egalit� stricte assur�e car factvu est rempli par des valeurs X
    if (factvu[i] == X) 
      return 0;
  } // fin i
  factvu[(*ifv++)]=X;
  return 1;
} // fin pasvu

/* ++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION
  Replace the R following commands:
# fnz <- apply(mat, 2, function(x){min(seq(along=x)[x!=0])})
# for(j in seq(nc)){
#   mat[,j] <- (inverses.basep(p)[mat[fnz[j],j]] * mat[,j]) %% p
#   } # fin j
#  pseudoweight <- apply(mat, 2, function(x){sum(x!=0)})
# weight <- apply(mat, 2, function(x){length(unique(factnum[x!=0]))})
INPUT
  nrow, ncol: mat dimensions
  p: a prime
  factnum: a numeric factor to identify rows of mat associated
           with the same factor; length=nrow
INPUT-OUTPUT
  addressofmat: a design key matrix in base p (big matrix)
OUTPUT
  weight, pseudoweight, binrank, modrank : vectors of length ncol
CALLED BY
  Called by the R function  weightorder.basep 
 ++++++++++++++++++++++++++++++++++++++++++++++++ */
SEXP PLANORweightorder(SEXP gnrow, SEXP gncol, SEXP gp,
		  SEXP gfactnum,
		       SEXP addressofmat, SEXP retour)

{

  int i, j, inv,a, fnz=0, ifv, *factvu;
  div_t d;
  /* Acceder aux arguments */
  SEXP gweight = getListElement(retour, "weight");
  SEXP gpseudoweight = getListElement(retour, "pseudoweight");
  SEXP gbinrank = getListElement(retour, "binrank");
  SEXP gmodrank = getListElement(retour, "modrank");

  int *nrow = INTEGER_POINTER(gnrow);
  int *ncol = INTEGER_POINTER(gncol);
  int *p = INTEGER_POINTER(gp);
  int *factnum = INTEGER_POINTER(gfactnum);
  double *weight = NUMERIC_POINTER(gweight);
  double *pseudoweight = NUMERIC_POINTER(gpseudoweight);
  double *binrank  = NUMERIC_POINTER(gbinrank);
  double *modrank = NUMERIC_POINTER(gmodrank);
  //Access to the big matrix mat
  BigMatrix *ptrin = (BigMatrix *) (R_ExternalPtrAddr(addressofmat));
  MatrixAccessor<TYPEOFBIG> bigmat(*ptrin);


R_CheckUserInterrupt(); // permettre a l'utilisateur d'interrompre

  factvu = (int *) R_alloc((*nrow), sizeof(int));

  for (j=0; j<*ncol; j++) {
    binrank[j]=0;
    modrank[(*ncol)-j-1]=0; // on inverse l'ordre des colonnes
    weight[j] =0;
    pseudoweight[j] =0;
  //    fnz est l' indice du 1ier non-zero dans la colonne j ou -1
    fnz =-1;
    for (i= 0; i< (*nrow); i++) {
      //      if (mat[i, j] != 0) {
//NOTE: big matrix indexes: first the column index
      if ( !ISZERO(bigmat[j ][i])) {
	  fnz=(i+1);
	break; // sortir de la boucle i
      }
    } //fin i
       //NOTE: big matrix indexes: first is the column index
    //    a=(int)bigmat[ fnz-1, j];
    // AB, 8/6/11 Ajout du test sur fnz:
    // la matrice restreinte aux facteurs traitements ou
    // blocs  peut  contenir  z�ro valeurs non nulles
    if (fnz==-1) {
      // pas de non nul sur la colonne j
      for (i=0; i< *nrow; i++) {
	bigmat[j][i] = 0;
      } 
    } else {
	
    a=(int)bigmat[j][ fnz-1];
    PLANORinv(p, &a, &inv);

    // factvu: les differentes valeurs deja rencontrees
    // dans la colonne j, sans double
    // ifv: indice du 1ier emplacement libre dans factvu
    ifv =0;
    for (i=0; i< *nrow; i++)
      factvu[i]=-1;
    for (i=0; i< *nrow; i++) {
      //      d= div(inv * mat[i,j], p);
      d= div((inv * (int)bigmat[j][i]), *p);
      //      mat[i,j] = d.rem;
      bigmat[j][i] = (short int) d.rem;

      //      if (mat[i,j] != 0) {
      if (ISZERO(bigmat[j][i]) == 0) {
	//	pseudoweight[j] += 1;
	pseudoweight[j] += 1;
	weight[j] += pasvu(factnum[i], &ifv, factvu);
	binrank[j] += (double) R_pow_di( (double)(*p), i);
	//	modrank[j] += (mat[i,j] * pow(p, (nr-i-1)));
	modrank[(*ncol)-j-1] += ( (double)(bigmat[j][i]) * 
			R_pow_di( (double)(*p), ((*nrow)-i-1)));

      } // fin if

    } // fin i
    } 


  } // fin j
  return(retour);

} // fin PLANORweightorder


/* ++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION
 Replace the time consuming R commands:
   for(k in seq(nb.admissible)){
    test.mat <- ImagesIS == matrix(admissible[,k,drop=F], nrow=r, ncol=nb.ineligible)
    test[k] <- sum( apply(test.mat, 2, prod) )==0
  }
  in the R function planor.kernelcheck.basep
INPUT
  r: integer
  nbadmissible: integer
  nbineligible: integer
  addressofImagesIS: (r x  nbineligible) matrix of integers
     (address of a short big.matrix)
  addressofadmissible: (r x nbadmissible) matrix of integers 
     (address of short big.matrix)
OUTPUT
  test:  a logical vector of length nbadmissible
CALLED BY
  planor.kernelcheck.basep
 ++++++++++++++++++++++++++++++++++++++++++++++++ */
SEXP PLANORloopkernelcheck( SEXP gr, SEXP gnbadmissible,
			    SEXP gnbineligible, SEXP addressofImagesIS,
			    SEXP addressofadmissible, SEXP gtest)
{
  int k,j,l, trouve;
  int *r = INTEGER_POINTER(gr);
  int *nbadmissible = INTEGER_POINTER(gnbadmissible);
  int *nbineligible = INTEGER_POINTER(gnbineligible);
  int *test = INTEGER_POINTER(gtest);
  // Access to the big matrices
  BigMatrix *ptrImagesIS = (BigMatrix *) (R_ExternalPtrAddr(addressofImagesIS));
  MatrixAccessor<TYPEOFBIG> bigImagesIS(*ptrImagesIS);

  BigMatrix *ptradmissible = (BigMatrix *) (R_ExternalPtrAddr(addressofadmissible));
  MatrixAccessor<TYPEOFBIG> bigadmissible(*ptradmissible);


R_CheckUserInterrupt(); // permettre a l'utilisateur d'interrompre

  for (k=0; k< *nbadmissible; k++) {
    test[k] =1;
    for (j=0; j < *nbineligible; j++) {
      trouve =1; // 1 si toute la colonne est egale
      for (l=0; l< *r; l++) {
	//	if (ImagesIS[l,j] !=admissible [l,k]) {
       //NOTE: big matrix indexes: first the column index
	if (!EQUAL(bigImagesIS[j][l], bigadmissible[k][l])) {
	  trouve=0;
	  break; // sortir de la boucle sur l
	} // fin if
      } // fin l
      if (trouve == 1) {
	test[k] = 0;
	break; // sortir de la boucle sur j
      } 
    } // fin j
  } // fin k
  return(gtest);

} // fin function

} // fin extern C
