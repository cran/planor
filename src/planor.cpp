/* TYPE OF THE matrices */
#define  TYPEOFMAT int

/* +++  includes from R +++ */
#include <R.h>
#include <Rmath.h>// useful for pow
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Utils.h> // to allow user interrupt
/* +++ system includes +++ */
#include <float.h>
#include <math.h> 
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
/* ++++++++++++++++++++++++++++++++++++++++++++++++ */
SEXP getListElement(SEXP list, const char *str)
{
  /* FUNCTION
 Access to a component by name from a R-list
INPUT
 list: a R-list
 str: name of a component of list
RETURN VALUE
 the component "str" of list
 ++++++++++++++++++++++++++++++++++++++++++++++++ */
  SEXP elmt = R_NilValue, names=getAttrib(list, R_NamesSymbol);
  int i;
  for (i=0; i< length(list); i++)
    if (strcmp(CHAR(STRING_ELT(names, i)), str) ==0) {
      elmt=VECTOR_ELT(list, i);
      break;
    }
  return(elmt);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++ */

  void multpp1col (SEXP addressofa, short int *B,
		   short int *iz,
		   int all, int nc, int nrow,  int N1,
		   short int pp1,
		   SEXP addressofres)
{
/* FUNCTION
 Multiplication of a given column, nc, of the  matrix A
 by the column nc of a matrix b. This column is stored in B.
 If iz[i]==pp1, B[i] is ignored.
 iz[i]=pp1, if the first nonzero value of the row i of t(b)
 is not 1
INPUT
 addressofA: address of  matrix A (nrow, N1)
 B: vector (N1)
 iz: vector (N1)
 all: see subgroup.base
 nc: current row of A
 nrow, N1: dimensions
 pp1: marker of values to ignored 
OUTPUT
 addressofres:  address of resulting  matrix (nrow,N1)
CALLED BY
 PLANORsubgroup
NOTE
 The dimensions are not checked
 ++++++++++++++++++++++++++++++++++++++++++++++++ */
  int l, r;
  /* Access to the values of the matrices */
  TYPEOFMAT *A = (TYPEOFMAT *) INTEGER_POINTER(addressofa);
  TYPEOFMAT *res = (TYPEOFMAT *) INTEGER_POINTER(addressofres);


R_CheckUserInterrupt(); // to allow user interrupt

 for (r=0; r< N1;r++) {
     // Compute iz
     /* we detect the lines to eliminate, only when a element 
	greater than 1  is encountered. As we progress column
	by column, beginning by the first one, the preceeding
	elements on the line have already been seen. They are,
	either nuls, -- and then, it is as they are ignored --,
	or equal to 1 -- and then, the line should be considered. */
   if (all==0) {
     if ((iz[r]==0) && (B[r]==1))
	 iz[r]=1;
       if ((iz[r]!=1) && (B[r]!=0) && (B[r]!=1))
	   iz[r]=pp1;
   } // end all


   if (iz[r] != pp1) {
     // line to be looked at
     for (l=0; l <nrow; l++) {
      //      res[l,r]+= A[l,nc]*B[r]
//NOTE: in C, the values are stored by column
       res[nrow *r + l] += (A[nrow* nc +l] * B[r]);

     } // end l
   } // end if iz
     else {
       // line of coeffs to be ignored, so the: res[0,r] =-1;
res[nrow * r ] =-1;
     }
 } // end r

} // end multcolpp1col


/* ++++++++++++++++++++++++++++++++++++++++++++++++ */
void     repetcol( int motif, int prod1, int l, 
	        short int *crosses)
{
  /* FUNCTION
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
  int i;
R_CheckUserInterrupt(); // for the user to interrupt

  for (i=l; i<(l+prod1); i++) {
    //    crosses[i, nc] = motif;
        crosses[i ] = motif;
  }
} // end repetcol




/* ++++++++++++++++++++++++++++++++++++++++++++++++ */

void gcrossingcol( int nrow, int p,  int prod1,
		      short int start,
		     short int *crosses) 
{
   /* FUNCTION
 Generates all n1 x n2 x ... x ns combinations of size s 
 with n1,...,ns integers
 Here, the values n1,...,ns are all identical and equal to p
 (NOTE AB: Replace the R function crossing of subgroup,
where all the n values are rep(p,nbg), so all equal to p)
INPUT
  nrow, ncol: dimension of the output matrix crosses
  p: the value of the series of integers
  start: integer from where to start the series of integers
OUTPUT
  crosses: an integer matrix with nrows rows and ncol columns 
  giving all combinations along the rows, in lexicographic order.
  (NOTE AB: the column order is reversed from the original)
  This program is called for each column, so the ouput is not
  the complete matrix, but one column
CALLED BY
 PLANORsubgroup

 ++++++++++++++++++++++++++++++++++++++++++++++++ */
  int  k, l, motif;
R_CheckUserInterrupt(); 

    motif= start;

    for (k=0; k < p; k++) {
      l=k*prod1;
      while ( (l+prod1) <= nrow) {
      /* Repeat motif, prod1 times in the elements
	 of crosses, from the l-st  */
	  repetcol(motif, prod1, l, crosses);
	l += ( p * prod1);
      }
      motif++;
    } // end k
} // end gcrossingcol


/* ++++++++++++++++++++++++++++++++++++++++++++++++ */

SEXP PLANORsubgroup( SEXP gnrow, SEXP gnbg, SEXP gN,
		     SEXP addressofmat, SEXP gp, SEXP gall,
		     SEXP addressofres)
 {
  /* FUNCTION
  Calculate the non null elements of the subgroup H generated 
  by the columns of mat, considered as vectors in (Zp)^s
INPUT
  nrow, nbg: dimension of mat
  N: number of columns of res (= p^nbg)
  addressofmat: address of mat:
    a  matrix of integers modulo p whose columns are assumed to
       be independent vectors in (Zp)^s
  p: a prime
  all: if TRUE all elements in H are given, if FALSE only elements up
       to multiplication by an integer modulo p are given
OUTPUT:
 addressofres:  resulting  matrix:
   a matrix of integers modulo p whose columns are the subgroup 
   elements
   dimension (nrow, N-1)
DETAILS:
  it is not checked whether the column vectors of mat are independent, so there
  may be several times 
CALLED BY
 subgroup.basep
 ++++++++++++++++++++++++++++++++++++++++++++++++ */
   short int  *coeffs, start, pp1, *iz;
   int prod1, nc,i, j;
   int aux;

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
      but, here, it is calculated column by column */

  coeffs = (short int *) R_alloc((*N), sizeof(short int));
  iz = (short int *) S_alloc((long) (*N - 1), sizeof(short int));
//iz: initialized at zero

  prod1=1;

  // Compute column nc of coeffs
  for (nc=0; nc< *nbg; nc++) {
    gcrossingcol((*N),  (*p), prod1, start, coeffs);
    prod1 *= (*p);
 
   /* Matrix calculation; only, the elts de coefs
      such as coeffs/iz!= pp1 are considered.
      The first line of coeffs is always ignored.
      It is why we begin at coeffs+1 and consider
      N-1 elts, only */
   multpp1col (addressofmat, (coeffs+1), iz, *all,
	       nc, (*nrow), (*N-1), pp1, addressofres);
  } // end nc

  /* Take modulo p of the result */
  // NOTE AB: we don't do it in the run while res is calculated
  // because this one is a progressive cumulation of values
  // whose final result is determined at the end of loops only
  //Access to the resulting  matrix
TYPEOFMAT *res = (TYPEOFMAT *) INTEGER_POINTER(addressofres);


  for (i=0; i< (*nrow); i++) {
    for (j=0; j< (*N)-1 ;j++) {
     //NOTE:  matrix strored by column
       //      a= ( (int) res[i,j]% (*p));
      if (res[(*nrow) *j +i] >0) {
	// The null values should be ignored
	//	res[i,j]= res[i,j]%p
	aux= (int) res[(*nrow) *j +i]% (*p); 
	res[(*nrow) *j +i]= (TYPEOFMAT)aux; 
      }
    } // end j
  } // end i

  return(addressofres);
} // end PLANORsubgroup




  /* ++++++++++++++++++++++++++++++++++++++++++++++++ */

SEXP PLANORlibsk(SEXP gnrow, SEXP gncol,
		 SEXP addressofH, SEXP gLIBtpf,
		 SEXP gMAXPRINT)
{
/* FUNCTION
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
  addressofH: address of H: a  matrix
  LIBtpf.k: character vector
  MAXPRINT: maximum number of rows and columns to print
OUTPUT
  Nothing, but print on standard output
CALLED BY
  The R function:  summary.designkey
 ++++++++++++++++++++++++++++++++++++++++++++++++ */
  int i,j;
  int *nrow = INTEGER_POINTER(gnrow);
  int *ncol = INTEGER_POINTER(gncol);
  int *maxprint = INTEGER_POINTER(gMAXPRINT);
  // Nbre max de colonnes a ip, pour limiter les ip
  // imin2: fonction de R.h
  int br = imin2(*nrow, *maxprint);
  int bc = imin2(*ncol, *maxprint);

R_CheckUserInterrupt(); // check User interrupt
// Access to the  matrix H
TYPEOFMAT *H = (TYPEOFMAT *) INTEGER_POINTER(addressofH);


  for (j=0; j< bc; j++) {
    Rprintf("1 =") ; 
    for (i=0; i < br; i++) {

      //      if (H[i,j] ==0) continue;
      //      if (H[j * (*nrow) + i] ==0) continue;
      // are replaced, to avoid rounding errors.
      // (NB: H contains positive integers)
if (ISZERO((double)H[  (*nrow)*j+i])) continue;

   Rprintf(" %s", CHAR(STRING_ELT(gLIBtpf,i)));

      //      if (H[i,j] !=1) 
if (!ISUN((double)H[(*nrow)*j+i])) {
Rprintf("^%d ", H[(*nrow)*j+i]);
   }

    } //fin i
    Rprintf("\n");
  } // end j

  /* AB, July 2015. Warn if all is not printed */
  if ((br <*nrow) || (bc <*ncol)) {
    Rprintf("\n");
  }
  if (br <*nrow) {
    Rprintf("The first %d rows on a total of %d\n", br, *nrow);
  }

  if (bc <*ncol) {
    Rprintf("The first %d columns on a total of %d\n", bc, *ncol);
  }

  return(addressofH);

} // end PLANORlibsk


/* ++++++++++++++++++++++++++++++++++++++++++++++++ */
int prodk(int p, int k)
{
/* FUNCTION
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
int l, a, prodkv;
div_t d;
  a=k;

  for (l=2; l <= (p-2); l++)
    {
      prodkv = a *l;
      d= div(prodkv, p);
      if (d.rem ==1) 
	return(l);
    } // end l
error("Internal error: prodk\n");
 return(-1);
}


/* ++++++++++++++++++++++++++++++++++++++++++++++++ */
void PLANORinv (int *p, int *J,
		int *inv)
{
  /* FUNCTION
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
  //  if ( (*J==1) || (*J==(*p-1))) {
  if ( EQUAL((double)(*J), 1) || EQUAL((double)(*J), (double)(*p-1))) {
    *inv = *J;
  } else {
    *inv =  prodk(*p, *J);
  }
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++ */
int pasvu(int X, int *ifv, int *factvu)
{
/* FUNCTION
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
  int i;
  for (i=0; i <= *ifv; i++) {
// strict equality ensured because factvu is filled in by values from X
    if (factvu[i] == X) 
      return 0;
  } // end i
  factvu[(*ifv++)]=X;
  return 1;
} // end pasvu

/* ++++++++++++++++++++++++++++++++++++++++++++++++ */
SEXP PLANORweightorder(SEXP gnrow, SEXP gncol, SEXP gp,
		  SEXP gfactnum,
		       SEXP addressofmat, SEXP retour)

{

/* FUNCTION
  Replace the R following commands:
# fnz <- apply(mat, 2, function(x){min(seq(along=x)[x!=0])})
# for(j in seq(nc)){
#   mat[,j] <- (inverses.basep(p)[mat[fnz[j],j]] * mat[,j]) %% p
#   } # end j
#  pseudoweight <- apply(mat, 2, function(x){sum(x!=0)})
# weight <- apply(mat, 2, function(x){length(unique(factnum[x!=0]))})
INPUT
  nrow, ncol: mat dimensions
  p: a prime
  factnum: a numeric factor to identify rows of mat associated
           with the same factor; length=nrow
INPUT-OUTPUT
  addressofmat: a design key matrix in base p
OUTPUT
  weight, pseudoweight, binrank, modrank : vectors of length ncol
CALLED BY
  Called by the R function  weightorder.basep 
 ++++++++++++++++++++++++++++++++++++++++++++++++ */
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
  int *weight = INTEGER_POINTER(gweight);
  int *pseudoweight = INTEGER_POINTER(gpseudoweight);
  int *binrank  = INTEGER_POINTER(gbinrank);
  int *modrank = INTEGER_POINTER(gmodrank);
  //Access to the  matrix mat

TYPEOFMAT *mat =  INTEGER_POINTER(addressofmat);


R_CheckUserInterrupt(); 

  factvu = (int *) R_alloc((*nrow), sizeof(int));

  for (j=0; j<*ncol; j++) {
    binrank[j]=0;
    modrank[(*ncol)-j-1]=0; // we invert columns order
    weight[j] =0;
    pseudoweight[j] =0;
  //    fnz is the index of the first no-zero in the column j ou -1
    fnz =-1;
    for (i= 0; i< (*nrow); i++) {
      //      if (mat[i, j] != 0) {
 if ( !ISZERO((double)mat[(*nrow)* j +i])) {

	  fnz=(i+1);
	break; // sortir de la boucle i
      }
    } //fin i
    // AB, 8/6/11 Add the test on fnz:
    // the matrix restricted to the treatment or block factors
    // may contain no non-nul values
    if (fnz==-1) {
      // no non-nul value in the column j
      for (i=0; i< *nrow; i++) {
mat[(*nrow)* j +i] = 0;

      } 
    } else {
a=(int)mat[(*nrow)*j +( fnz-1)];

    PLANORinv(p, &a, &inv);

    // factvu: the different values already encountered in column j
    // without double
    // ifv: index of the first empty place in factvu
    ifv =0;
    for (i=0; i< *nrow; i++)
      factvu[i]=-1;
    for (i=0; i< *nrow; i++) {
      //      d= div(inv * mat[i,j], p);
d= div((inv * (int)mat[(*nrow)* j +i]), *p);
      //      mat[i,j] = d.rem;
      mat[(*nrow)* j +i] = (short int) d.rem;

      //      if (mat[i,j] != 0) {
      if (ISZERO((double)mat[(*nrow)* j +i]) == 0) {

	//	pseudoweight[j] += 1;
	pseudoweight[j] += 1;
	weight[j] += pasvu(factnum[i], &ifv, factvu);
	binrank[j] += (int) R_pow_di( (double)(*p), i);
	//	modrank[j] += (mat[i,j] * pow(p, (nr-i-1)));
modrank[(*ncol)-j-1] += ( (int)(mat[(*nrow)* j +i]) * 
		R_pow_di( (double)(*p), ((*nrow)-i-1)));


      } // end if

    } // end i
    } 


  } // end j
  return(retour);

} // end PLANORweightorder


/* ++++++++++++++++++++++++++++++++++++++++++++++++ */
SEXP PLANORloopkernelcheck( SEXP gr, SEXP gnbadmissible,
			    SEXP gnbineligible, SEXP addressofImagesIS,
			    SEXP addressofadmissible, SEXP gtest)
{
/* FUNCTION
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
  addressofImagesIS: (r x  nbineligible)  matrix of integers
  addressofadmissible: (r x nbadmissible)  matrix of integers 
OUTPUT
  test:  a logical vector of length nbadmissible
CALLED BY
  planor.kernelcheck.basep
 ++++++++++++++++++++++++++++++++++++++++++++++++ */
  int k,j,l, trouve;
  int *r = INTEGER_POINTER(gr);
  int *nbadmissible = INTEGER_POINTER(gnbadmissible);
  int *nbineligible = INTEGER_POINTER(gnbineligible);
  int *test = INTEGER_POINTER(gtest);

  // Access to the matrices
  TYPEOFMAT *ImagesIS = (TYPEOFMAT *) INTEGER_POINTER(addressofImagesIS);

  TYPEOFMAT *admissible = (TYPEOFMAT *)  INTEGER_POINTER(addressofadmissible);



R_CheckUserInterrupt(); // allow user interrupt

  for (k=0; k< *nbadmissible; k++) {
    test[k] =1;
    for (j=0; j < *nbineligible; j++) {
      trouve =1; // 1 if the whole column is equal
      for (l=0; l< *r; l++) {
	//	if (ImagesIS[l,j] !=admissible [l,k]) {
//NOTE: values stored by column
	if (!EQUAL((double)ImagesIS[(*r)*j +l], (double)admissible[(*r)*k+l])) {

	  trouve=0;
	  break; // go out the loop on l
	} // end if
      } // end l
      if (trouve == 1) {
	test[k] = 0;
	break; // go out the loop on j
      } 
    } // end j
  } // end k
  return(gtest);

} // end function

} // end extern C
