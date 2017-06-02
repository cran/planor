#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;

/* +++  includes from R +++ */
#include <R.h>
#include <Rmath.h>// useful for pow
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Utils.h> // to allow user interrupt
#include <R_ext/Rdynload.h>

/* +++ system includes +++ */
#include <float.h>
#include <math.h> 
#include <stdlib.h>
#include <string.h>

/* +++ local includes  +++ */
#include "util.h"
#include "subgroup.h"
#include "inverses_basep.h"
#include "kernelcheck.h"
#include "weightorder.h"




/* ++++++++++++++++++++++++++++++++++++++
PLANORweightorder
Interface from R
INPUT
- matrice: a matrix (nrow X ncol)
- p: a prime
- factnum: a numeric factor to identify rows of mat associated
           with the same factor; length=nrow 
OUTPUT
 weight, pseudoweight, binrank, modrank : vectors of length ncol
CALLED BY
 Called by the R function weightorder.basep
 which is called by summary.keymatrix 
++++++++++++++++++++++++++++++++++++++++ */

  
RcppExport SEXP PLANORweightorder(SEXP Smatrice,
				  SEXP nrow,
				  SEXP ncol,
				  SEXP Sp,
				     SEXP Sfactnum) {

  /* Arguments */
  IntegerVector matrice(Smatrice);
  int * nr = INTEGER_POINTER(nrow);
  int * nc = INTEGER_POINTER(ncol);
  IntegerMatrix mat = fillMat(matrice, *nr, *nc);

  int *pp = INTEGER_POINTER(Sp);
  int p=*pp;

  IntegerVector factnum(Sfactnum);
  /* Working */
  int j, fnz, zero= (int)0;
  IntegerVector onecol;

  /* Output */
  IntegerVector  pseudoweight(*nc), weight(*nc), 
    binrank(*nc), modrank(*nc);


  /* Raw calculation of the inverses modulo p */
  IntegerVector invp= inverses_basep(p);

  R_CheckUserInterrupt(); 

for (j=0; j <*nc; j++) {
    onecol = mat(_,j); 

    // fnz=position of the first non-zero
    fnz = findFnzIndex(onecol);

    if (fnz == -1) {
      // the column is entirely zero
      std::fill(onecol.begin(), onecol.end(), zero);
    } else {
      onecol =invp[mat(fnz,j)-1] * mat(_,j);
      onecol = moduloVec(onecol, p);
    } // end else (fnz ==  onecol.end())

    /* weight = the  number of different values in factnum[onecol!=0] */ 
    IntegerVector w = findFactnum(factnum,onecol);
    sort(w.begin(), w.end());
    IntegerVector::iterator new_endw = unique(w.begin(), w.end());
    weight[j] = ( new_endw - w.begin()) ;

    IntegerVector attr=attrWeight( onecol, p);
    pseudoweight[j]  = attr[0];
    binrank[j]  = attr[1];
    modrank[j] = attr[2];

    mat(_,j) =onecol;
 } // end j

 std::reverse(modrank.begin(), modrank.end());
 /* Return a list */
 List outlist;
 outlist["weight"] = weight;
 outlist["pseudoweight"] = pseudoweight;
 outlist["binrank"] = binrank;
 outlist["modrank"] = modrank;

  return(outlist);
} // end PLANORweightorder


/* ++++++++++++++++++++++++++++++++++++++
PLANORprintLib
Function to display
Interface from R
FUNCTION
 Replace the time consuming R commands:
    for(j in 1:ncol(H)){
      select.kj <- H[,j] != 0
      coeff.char <- paste("^",as.character(H[select.kj,j]), sep="")
      coeff.char[H[select.kj,j] == 1] <- ""
      LIBS.k <- paste(LIBtpf.k[select.kj], coeff.char, sep="")
      cat("1 = ") ; cat(LIBS.k,sep=" ") ; cat("\n")
    }

INPUT
  H: a matrix
  nrow, ncol: dimension of H
  addressofH: address of H: a  matrix
  LIBtpf.k: character vector
  MAXPRINT: maximum number of rows and columns to print
OUTPUT
  Nothing, but print on standard output
CALLED BY
  The R function:  summary.keymatrix
+++++++++++++++++++++++++++++++++++++++++ */

 RcppExport void PLANORprintLib(SEXP SH,
		   SEXP gnrow, SEXP gncol,
		   SEXP gLIBtpf,
		 SEXP gMAXPRINT)
{

  IntegerVector matrice(SH);
  int *nrow = INTEGER_POINTER(gnrow);
  int *ncol = INTEGER_POINTER(gncol);
  IntegerMatrix H = fillMat(matrice, *nrow, *ncol);

  CharacterVector LIBtpf(gLIBtpf);

  int *maxprint = INTEGER_POINTER(gMAXPRINT);
  // Maximum number of columns to display
  // (not all, to limit the amount to display)
  // imin2 is a function of R.h
  int br = MIN(*nrow, *maxprint);
  int bc = MIN(*ncol, *maxprint);
  int zero = (int)0;

  int i,j;

R_CheckUserInterrupt(); // check User interrupt

  for (j=0; j< bc; j++) {
    Rprintf("1 =") ; 
    for (i=0; i < br; i++) {

      //      if (H[i,j] ==0) continue;
      // (NB: H contains positive or null integers)
      if (H(i,j) == zero) {
	continue;
      }


      Rprintf(" %s", CHAR(LIBtpf[i]));

      //      if (H[i,j] !=1) 
      if (H(i,j) != (int)1) {
	Rprintf("^%d ", H[(*nrow)*j+i]);
      }
    } // end i

    Rprintf("\n");
  } // end j

  /*  Warn if all is not printed */
  if ((br <*nrow) || (bc <*ncol)) {
    Rprintf("\n");
  }
  if (br <*nrow) {
    Rprintf("The first %d rows on a total of %d\n", br, *nrow);
  }

  if (bc <*ncol) {
    Rprintf("The first %d columns on a total of %d\n", bc, *ncol);
  }


} // end PLANORprintLib


/* ++++++++++++++++++++++++++++++++++++++
PLANORkernelcheck
Interface from R
FUNCTION
Checks whether any of the N elements in an ineligible set belongs to
the kernel of cbind(PhiStar, admissible[,k]), with PhiStar a p-morphism
matrix and 'admissible' a set of n candidates for making the next
column of PhiStar.
INPUT 
 PhiStar: a (r x s) matrix of integers modulo p
 nrowPhi = r
 ncolPhi = s
 admissible: a (r x n) matrix of integers modulo p
 IneligibleSet: a (s x N) matrix of integers modulo p
 p: a prime number
DETAILS
 For each column A_k of the matrix 'admissible', the function evaluates
 the (r x (s+1)) matrix (PhiStar|A_k). It returns TRUE if no vector
 (I_j'|1)' belongs to the kernel of (PhiStar|A_k), where I_j is the jth column
 of 'IneligibleSet'
CALLED BY
the R function  planor.designkey.recursive
++++++++++++++++++++++++++++++++++++++++ */
RcppExport  SEXP PLANORkernelcheck(SEXP SPhiStar,
				  SEXP nrowPhi,
				  SEXP ncolPhi,
				  SEXP Sadmissible,
				  SEXP nrowadm,
				  SEXP ncoladm,
				  SEXP SIneligibleSet, 
				  SEXP nrowinel,
				  SEXP ncolinel,
				  SEXP Sp) {

  IntegerVector PhiStar(SPhiStar);
  IntegerVector admissible(Sadmissible);
  IntegerVector IneligibleSet(SIneligibleSet);
  int * nrPhi= INTEGER_POINTER(nrowPhi);
  int * ncPhi = INTEGER_POINTER(ncolPhi);
  int * nradm = INTEGER_POINTER(nrowadm);
  int * ncadm = INTEGER_POINTER(ncoladm);
  int * nrinel = INTEGER_POINTER(nrowinel);
  int * ncinel = INTEGER_POINTER(ncolinel);
  int * p = INTEGER_POINTER(Sp);



  IntegerMatrix maPhiStar = fillMat(PhiStar, *nrPhi, *ncPhi);
  IntegerMatrix maadmissible = fillMat(admissible, *nradm, *ncadm);
  IntegerMatrix maIneligibleSet = fillMat(IneligibleSet, *nrinel, *ncinel);


LogicalVector outvec = kernelcheck_basep(maPhiStar, 
				maadmissible, 
					 maIneligibleSet, *p);
 return(wrap(outvec));
} // end PLANORkernelcheck

/* ++++++++++++++++++++++++++++++++++++++
PLANORinverses_basep
Interface from R
Raw calculation of the inverses modulo p
INPUT
p: a prime
CALLED BY
 kernelmatrix.basep, alias.keymatrix
++++++++++++++++++++++++++++++++++++++++ */

RcppExport SEXP PLANORinversesbasep( SEXP Sp) {
  int *pp = INTEGER_POINTER(Sp);
  IntegerVector outvec= inverses_basep(*pp);
  return(wrap(outvec));
} // end PLANORinversesbasep

/* ++++++++++++++++++++++++++++++++++++++
PLANORdesignkeysym
Return a list of design key matrices of size r x s 
with elements in Zp
CALLED BY
 The R function planor.designkey.basep in case of symmetries
and no recursion
++++++++++++++++++++++++++++++++++++++++ */
/* ** HM 2016-09-07 Acceleration de PLANORdesignkeynonsym 
     par les symetries 
     Chercher: *** FAST
** */

RcppExport SEXP PLANORdesignkeysym(SEXP SPhiStar, 
		      SEXP SUstar,
		      SEXP Sineligible_lnz,
		      SEXP Sineligible,
		      SEXP SHmaxj,
		      SEXP SHset,
		      SEXP Ss,
		      SEXP Sf,
		      SEXP Sp,
		      SEXP Smax_sol,
		      SEXP Srandomsearch,
		      SEXP Sverbose) {

  /* Arguments */
  IntegerMatrix PhiStar(SPhiStar);
  IntegerMatrix Ustar(SUstar);
  IntegerMatrix ineligible_lnz(Sineligible_lnz);
  IntegerMatrix ineligible(Sineligible);
  IntegerVector Hmaxj(SHmaxj);
  const List Hset(SHset);
  int *s=INTEGER_POINTER(Ss);
  int *f=INTEGER_POINTER(Sf);
  int smoinsf = *s - *f;
  int *p=INTEGER_POINTER(Sp);
  int * max_sol = INTEGER_POINTER(Smax_sol);
  int *randomsearch= INTEGER_POINTER(Srandomsearch);
  int *verbose= INTEGER_POINTER(Sverbose);

/* Careful: The following structures contain column numbers from 1:
   f,s,ktr,Hmaxj, k siaA, ineligible_lnz,  admiss_j,  Hset */


  /* Return a list of matrices */
  List PhiStar_solution(*max_sol);

  /* *** FAST iaA is no more a list:
  iiaA: initially admissible elements for j [(r x n_j) p-matrices] 
  List iiaA(smoinsf);
  *** */

  /* siaA: status of initially admissible elements for j 
     [n_j-length vectors] */
   List siaA(smoinsf);

   /* wiaA: currently admissible elements for j 
      [vectors of selected indices between 1 and n_j] */
   List wiaA(smoinsf);


   /* ktr: smallest visited index k since the previous visit in j */
  IntegerVector ktr( smoinsf, 0);

  IntegerVector firstvisit (  smoinsf, 1);

  /* *** FAST liaA no more useful:
     liaA:  nbs of initially admissible elements for j 
     IntegerVector liaA(smoinsf, 0);
     *** */

/* Working */
 int i, k,l;
 IntegerMatrix candidates,  iaA, inelig_kj,  newcolj;
 IntegerVector admiss_j,  admiss_j1, admiss_index_j;
 LogicalVector iae_j, select_kj;

/* ****************************************
  BACKTRACK SEARCH
  Remarks:
  1. When j is reached forward, the admissible elements of U*, for
     column f+j  of the key matrix, must be identified. We assume that the ineligible
     elements are all in canonical form, that is, their last non-zero coordinate is
     equal to 1. So we do not need it hence the selection of their first (f+j-1)
     coordinates at step j (forward) below (see Kobilinsky, 2000, end of page 16)
  2. The management of initially admissible elements and hierarchies is different
     from that proposed in Kobilinsky, 2000. Here, new 'initially admissible'
     elements are calculated as soon as hierarchy constraints have changed for
     a given factor A_j.
     *************************************** */


 
 int j = (int)0; // indices begin from 0
 int jprev = (int)-1;
 int nsolution = (int)0;
 int zero = (int)0;

 while(j >= zero){
   int fplusj = *f +j  ; 
   int fplusjmoins1 = *f +j -1; 


  PhiStar = PhiStar(_, Range(0, fplusjmoins1));

/* ** STEP 1: IF column j is reached forward 
      -> identification and updating of the admissible elements ** */

  if(jprev < j){
 /* **
      STEP 1.A: generation of the INITIALLY admissible elements (iae) for j
       -> performed ONLY IF this is the first visit to j with the present
          set of Hset[[j]] columns
       -> this part will remain valid during all the backtrack search except
          for factors A_j subject to changing hierarchy constraints
          (Hmaxj[j] > f)
      ---
      The condition just below is satisfied if this is the first visit to j
      or if A_j is subject to hierarchy constraints that have changed since
      the previous visit to j
** */

  	   if(ktr[j] <= Hmaxj[j]){
	     if(*verbose && (ktr[j] == zero)){ Rprintf("      first visit to column %d\n", fplusj+1);
		   }
	     /* ** 1.A.a: elements of U* that satisfy the hierarchy constraints
        ---
        The IF condition just below is satisfied only if A_j is subject to
        hierarchy constraints and they have changed since the previous visit
** */

		   if(zero < Hmaxj[j]){
		     //  b.aux = b.PhiStar[,Hset[[j]],drop=FALSE]

		     IntegerVector Hsetj=Hset[j];

		     IntegerVector indices(Hsetj.size());
		     // Copy into indices the content of Hset
		     std::copy(Hsetj.begin(), Hsetj.end(),
			      indices.begin());

		     indices = indices - 1; // index form 0
		     IntegerMatrix aux= takeCol(PhiStar, indices);
		     candidates= subgroup(aux, *p, true);
	     	} // end (0 < Hmaxj[j])
		   else {
		     /*  ELSE -> no hierarchy constraint */
		     /* *** FASTAND FURIOUS : 
			same b.iaA for all j, constructed at j=1; 
			the rest is managed by siaA
			Instead of: 
			candidates= Ustar; 
			*** */
		     if (j == zero) {
		       candidates= Ustar;
		     } else {
		       candidates=iaA;
		     }
		   } // end else  -> no hierarchy constraint

		   /* ** 1.A.b: selection based on the relevant 
		      ineligible subset * **/
		     k=MAX(*f, Hmaxj[j]);
		     //  select.kj <- (ineligible.lnz[1,] <= k) & (ineligible.lnz[2,] == (f+j))
		     select_kj = (ineligible_lnz(0,_) <= k) & (ineligible_lnz(1,_)== (fplusj+1)); 


		     if (is_true(all(select_kj== false))) {
		       /* *** TODO: to test this case
			  HM, 14 May 2011: tested and changed: 
			  length -> ncol *** */
		       iae_j = rep(true,candidates.ncol());
		     } else {

		       // b.inelig.kj <- as.matrix(b.ineligible[ seq_len(f+j-1),
		       //                 select.kj, drop=FALSE ])
		       inelig_kj = ineligible( Range(0, fplusjmoins1), _);
		       inelig_kj = takeColLog(inelig_kj, select_kj);

		       iae_j =
			 kernelcheck_basep(PhiStar,  candidates,
					   inelig_kj, *p);

		     } // end else (all(select_kj== false))

		     if (is_true(all(iae_j == false))) {
		       Rprintf("No solution for column %d of the design key\n", fplusj+1);

		       if (nsolution > zero) { 
			 return(wrap(PhiStar_solution[Range(0, nsolution-1)])); 
		       } else {
			 return(wrap(IntegerVector(0)));
				}
		     }
		       
		     /* ** 1.A.c: Information storage after the visit to j ** */
		     /* *** FAST: add the test on j
			Instead of 
			iiaA[j] = takeColLog(candidates, iae_j);
			*** */
		     if (j == zero) {
		       //  b.iaA <- b.candidates[, iae.j, drop=FALSE]
		       iaA = takeColLog(candidates, iae_j);
		       siaA[j] = IntegerVector(rep(*s, iaA.ncol()));
		     } else {
		       //  siaA[[j]] <- apply(rbind(siaA[[j-1]], s*iae.j), 2, min)

		       IntegerVector siaAj1 = IntegerVector(siaA[j-1]); 

		       // pmin is in R sugar
		       siaA[j] = IntegerVector(pmin(siaAj1, IntegerVector(iae_j)* *s));

		     } // end j!=0


		     ktr[j] = k+1;

		     admiss_j =  IntegerVector(siaA[j]);

		     /* *** FAST *** */
		     if (*verbose) {
		       LogicalVector vaux= (admiss_j  == *s);
       		       int auxx=std::accumulate(vaux.begin(), vaux.end(), 0);
		        Rprintf("    ---    col. %d ( j = %d) %d initial candidates\n",
		       			       fplusj+1, j+1, auxx);
		     } 

	   } // end (ktr[j] <= Hmaxj[j])

	   /* ** 
STEP 1.B: UPDATING of the STATUS of the initially admissible elements
                for factor A_j
      -> based on the ineligible elements with last-but-one non-zero element
         between ktr[j] and j-1 and with last non-zero element equal to j
      ---
      The IF condition just below is satisfied when kernel checks are necessary
      ** */

	   if (ktr[j] < fplusj+1) {

	     /* Make elements eligible if their status may have changed */
	     admiss_j =  IntegerVector(siaA[j]);

	     int maxktr = MAX(1, ktr[j]); /* *** FAST *** */
	     admiss_j1 = IntegerVector(siaA[j-1]);

	     for (i=0; i< admiss_j.size(); i++) {
	     /* *** FAST *** */
	     // admiss.j[ admiss.j >= max(1,ktr[j]) ] <- s
	       if (admiss_j[i] >= maxktr) {
		 admiss_j[i] = *s;
	       }
	       /* *** FAST *** */
	       // admiss.j <- apply(rbind(siaA[[j-1]], admiss.j), 2, min) 
	       if (admiss_j[i] > admiss_j1[i]) {
		 admiss_j[i] =admiss_j1[i];
	       }
	     } // end i

	     /*  Loop over the last-but-one non-zero indices */

	     for (int kk= ktr[j]; kk<= fplusj; kk++) {
	       /* relevant ineligible elements */
	       select_kj =  (ineligible_lnz(0,_) == kk) &(ineligible_lnz(1,_)== fplusj+1) ;

	       /* relevant initially admissible elements */
	       // tocheck <- seq_len(liaA)[admiss.j == s]
	       /* *** FAST instead of:
		  tocheck <- seq_len(liaA[[j]])[admiss.j == s] 
		  *** */

	       IntegerVector tocheck = findIndexVal(admiss_j, *s);

	       /* *** HM, 30/06/10 : condition sur length(tocheck) 
		  et sur inelig.kj *** */

	       if ((tocheck.size()>0) & is_true(any(select_kj))) {
		 // b.inelig.kj <- as.matrix(b.ineligible[ seq_len(f+j-1),                                                      select.kj, drop=FALSE ])
		 inelig_kj = ineligible(Range(0, fplusjmoins1), _);
		 inelig_kj = takeColLog(inelig_kj, select_kj);
		 candidates= takeCol(iaA, tocheck);
		 iae_j =
			 kernelcheck_basep(PhiStar,  candidates,
					   inelig_kj, *p);


		 // admiss.j[ tocheck[!ae.j] ] <- k
		 IntegerVector iaefalse = findIndexLog(iae_j, false);
		 for (i=0; i< iaefalse.size(); i++) {
		   l = tocheck[iaefalse[i]];
		   admiss_j[l] =kk;
		 }

	       } // end (any(tocheck == true)

	       /* *** HM, 30/06/10 : fin condition sur length(tocheck) et sur inelig.kj *** */
	     } // end loop kk

	     siaA[j] = IntegerVector(admiss_j);

	     if (firstvisit[j] <= 2) {
	       if (*verbose) {
		 //	cat("    ---    col.", f+j, " ( j =",j, ") ", sum(admiss.j == s), " selected candidates\n")
		 LogicalVector vaux= (admiss_j  == *s);
		 int auxx=std::accumulate(vaux.begin(), vaux.end(), 0);
		 Rprintf("    ---    col. %d ( j = %d) %d selected candidates\n",
		     fplusj+1, j+1, auxx);

	       } 
	       firstvisit[j] = firstvisit[j]+1;
	     } // end (firstvisit[j] <= 2)

	   } // end (ktr[j] < *f+j)

	   /* ** STEP 1.C: UPDATED admissible elements for factor A_j ** */

	   // admiss.index.j <- seq_len(liaA)[siaA[[j]]==s] 
	   /* *** FAST instead of
	      admiss.index.j <- seq_len(liaA[j])[siaA[[j]]==s] 
	      ** */
	   IntegerVector siaAj=IntegerVector(siaA[j]);
	   LogicalVector ad = ( siaAj == *s);
	   admiss_index_j = findIndexLog(ad, true);


	   if (*randomsearch == 1) {
	     // Call the Rcpp function "runif"
	     int nindexsample = (int)admiss_index_j.size() - 1;
	     NumericVector index_sample= runif(nindexsample);
	     admiss_index_j = takeElt(admiss_index_j, index_sample);
	   } 
	   admiss_index_j = admiss_index_j + 1; // index from 1
	   wiaA[j] = admiss_index_j;
	   ktr[j] = fplusj+1;

   } // end  (jprev < j)

/* **
  END OF STEP 1

  STEP 2: INCLUSION OF THE NEXT ADMISSIBLE ELEMENT IN PhiStar
  ** */
  IntegerVector wj= IntegerVector(wiaA[j]);


 /*  CASE 2.A: If there is a next admissible element, add it to PhiStar ** */

  if (wj.size() > zero ) {
    // b.newcolj <- (b.iaA)[,wj[1]]
    /* *** FAST instead of
       b.newcolj <- (b.iaA[[j]])[,wj[1]]
       *** */
    int auxw = wj[0] - (int)1; // index from 0
    // Note : Range is required to not build a column-matrix
    newcolj = iaA(_, Range(auxw, auxw));
    // wiaA[[j]] <- wj[-1]
    if (wj.size() >= (int)2) {
      wiaA[j]=wj[Range(1, wj.size()-1)];
    } else {
      wiaA[j]=IntegerVector(0);
    }

    //   b.PhiStar <- cbind(b.PhiStar, b.newcolj, deparse.level=0)
     PhiStar=cbind(PhiStar, newcolj);

     /* **  CASE 2.A.a: If this was the last column to fill in PhiStar ** */
     if ( (j+1) == smoinsf) {
       /* Save the solution */
      PhiStar_solution[nsolution++] = PhiStar;
      /* Return if enough solutions have been found */
      if (nsolution == *max_sol) {
	Rprintf("The search is closed: max.sol = %d solution(s) found \n", nsolution);
	return wrap(PhiStar_solution[Range(0, nsolution-1)]);
      }
      jprev = j ;
    } // end  if(j == (s-f))
    else {
      /* ** CASE 2.A.b: If this was not the last column, go forward to next j ** */
      jprev = j ; j=j+1;
    }
  } // end if(length(wj) > 0)

  else {
    /* **  CASE 2.B: If there is no next admissible element ** */
    /* Update information */
    /* *** HM: modif 2015
       ktr[ (seq_len(s-f) >= j) & (ktr > 0) ] <- j-1
       changed into :
    //  ktr.change <- (seq(s-f) >= j) & (ktr > j-1)
    //  ktr[ ktr.change ] <- j-1 *** */
    for (int i=j; i< smoinsf; i++) {
      if (ktr[i] > j) {
	ktr[i] = j; // ktr contains indices from 1
      }
    }
    /* Then go backward */
    jprev = j;
    j =  j-1;
  }
 } // end while(j > 0)

 Rprintf("The search is closed: %d solution(s) found \n",
	 nsolution);

 if (nsolution == zero) {
   return(wrap(IntegerVector(0)));
 } else {
   return wrap(PhiStar_solution[Range(0, nsolution-1)]);
 }

} // end PLANORdesignkeysym



/* ++++++++++++++++++++++++++++++++++++++
PLANORdesignkeynonsym
Interface with R
Return a list of design key matrices of size r x s 
with elements in Zp
CALLED BY
 The R function planor.designkey.basep when there are no symmetries
nor recursion
++++++++++++++++++++++++++++++++++++++++ */

RcppExport SEXP PLANORdesignkeynonsym(SEXP SPhiStar, 
		      SEXP SUstar,
		      SEXP Sineligible_lnz,
		      SEXP Sineligible,
		      SEXP SHmaxj,
		      SEXP SHset,
		      SEXP Ss,
		      SEXP Sf,
		      SEXP Sp,
		      SEXP Smax_sol,
		      SEXP Srandomsearch,
		      SEXP Sverbose) {

  /* Arguments */
  IntegerMatrix PhiStar(SPhiStar);
  IntegerMatrix Ustar(SUstar);
  IntegerMatrix ineligible_lnz(Sineligible_lnz);
  IntegerMatrix ineligible(Sineligible);
  IntegerVector Hmaxj(SHmaxj);
  const List Hset(SHset);
  int *s=INTEGER_POINTER(Ss);
  int *f=INTEGER_POINTER(Sf);
  int smoinsf = *s - *f;
  int *p=INTEGER_POINTER(Sp);
  int * max_sol = INTEGER_POINTER(Smax_sol);
  int *randomsearch= INTEGER_POINTER(Srandomsearch);
  int *verbose= INTEGER_POINTER(Sverbose);


  /* Return a list of matrices */
  List PhiStar_solution(*max_sol);
  List nosolution(0);

  /* iiaA: initially admissible elements for j 
     [(r x n_j) p-matrices] */
  List iiaA(smoinsf);

  /* siaA: status of initially admissible elements for j 
     [n_j-length vectors] */
   List siaA(smoinsf);

   /* wiaA: currently admissible elements for j 
      [vectors of selected indices between 1 and n_j] */
   List wiaA(smoinsf);

   /* ktr: smallest visited index k since the previous visit in j */
  IntegerVector ktr(smoinsf, 0);

  /* liaA:  nbs of initially admissible elements for j */
  IntegerVector liaA(smoinsf, 0);

  /* firstvisit: for print */
  IntegerVector firstvisit (smoinsf, 1);


  /* Working */
 int i, l;
 int k= (int)-1;
 IntegerMatrix candidates,  inelig_kj,  newcolj;
 IntegerVector admiss_j,  admiss_j1, admiss_index_j;
 LogicalVector iae_j, select_kj;

/* Careful: The following structures contain column numbers from 1:
   f,s,ktr,Hmaxj, k siaA, ineligible_lnz,  admiss_j,  Hset */


/* ****************************************
  BACKTRACK SEARCH
  Remarks:
  1. When j is reached forward, the admissible elements of U*, for
     column f+j  of the key matrix, must be identified. We assume that the ineligible
     elements are all in canonical form, that is, their last non-zero coordinate is
     equal to 1. So we do not need it hence the selection of their first (f+j-1)
     coordinates at step j (forward) below (see Kobilinsky, 2000, end of page 16)
  2. The management of initially admissible elements and hierarchies is different
     from that proposed in Kobilinsky, 2000. Here, new 'initially admissible'
     elements are calculated as soon as hierarchy constraints have changed for
     a given factor A_j.
     *************************************** */


 int j=0; // indices begin from 0
 int jprev= (int)-1;
 int nsolution =(int)0;
 int zero = (int)0;

 while(j >= zero){
   int fplusj = *f +j  ; 
   int fplusjmoins1 = *f +j -1; 
  PhiStar = PhiStar(_, Range(0, fplusjmoins1));
  /* ** STEP 1: IF column j is reached forward 
      -> identification and updating of the admissible elements ** */

  if(jprev < j){

    /* **
      STEP 1.A: generation of the INITIALLY admissible elements (iae) for j
       -> performed ONLY IF this is the first visit to j with the present
          set of Hset[[j]] columns
       -> this part will remain valid during all the backtrack search except
          for factors A_j subject to changing hierarchy constraints
          (Hmaxj[j] > f)
      ---
      The condition just below is satisfied if this is the first visit to j
      or if A_j is subject to hierarchy constraints that have changed since
      the previous visit to j
** */
      
  	   if(ktr[j] <= Hmaxj[j]){
	     if(*verbose && (ktr[j] == zero)){ Rprintf("      first visit to column %d\n", fplusj+1);
		   }
	     /* ** 1.A.a: elements of U* that satisfy the hierarchy constraints
        ---
        The IF condition just below is satisfied only if A_j is subject to
        hierarchy constraints and they have changed since the previous visit
** */

		   if(zero < Hmaxj[j]){
		     //  b.aux = b.PhiStar[,Hset[[j]],drop=FALSE]
		     IntegerVector Hsetj=Hset[j];
		     IntegerVector indices(Hsetj.size());

		    std::copy(Hsetj.begin(), Hsetj.end(),
			      indices.begin());
		     indices = indices -1;// indices from 0

		     IntegerMatrix aux= takeCol(PhiStar, indices);
		     candidates= subgroup(aux, *p, true);
		   } // end (0 < Hmaxj[j])
		   else {
		     /*  ELSE -> no hierarchy constraint */
		       candidates= Ustar;
		   } // end else  -> no hierarchy constraint

		   /* ** 1.A.b: selection based on the relevant 
		      ineligible subset * **/
		     k=MAX(*f, Hmaxj[j]);
		     //  select.kj <- (ineligible.lnz[1,] <= k) & (ineligible.lnz[2,] == (f+j))
		     select_kj = (ineligible_lnz(0,_) <= k) & (ineligible_lnz(1,_)== (fplusj+1)); 


		     if (is_true(all(select_kj== false))) {
		       /* *** TODO: to test this case
			  HM, 14 May 2011: tested and changed: 
			  length -> ncol *** */
		       iae_j = rep(true,candidates.ncol());
		     } else {
		       // b.inelig.kj <- as.matrix(b.ineligible[ seq_len(f+j-1),
		       //                 select.kj, drop=FALSE ])
		       inelig_kj = ineligible( Range(0, fplusjmoins1), _);
		       inelig_kj = takeColLog(inelig_kj, select_kj);

		       iae_j =
			 kernelcheck_basep(PhiStar,  candidates,
					   inelig_kj, *p);


		     } // end else (all(select_kj== false))

		     if (is_true(all(iae_j == false))) {
		       Rprintf("No solution for column %d of the design key\n", fplusj+1);
		       if (nsolution > zero) { 
			 return(wrap(PhiStar_solution[Range(0, nsolution-1)])); 
		       } else {
			 return(wrap(nosolution));
				}
		     }
		     // b.iaA[[j]] <- b.candidates[, iae.j, drop=FALSE]		       
		     iiaA[j] = takeColLog(candidates, iae_j);

		     /* ** 1.A.c: Information storage after the visit to j ** */



		     // liaA[j] <- sum(iae.j)
		     liaA[j]= std::accumulate(iae_j.begin(), iae_j.end(), 0);
		     // siaA[[j]] <- rep(s, liaA[j])
		     siaA[j] = IntegerVector(rep(*s, liaA[j]));
		     ktr[j] = k+1;
	   } // end (ktr[j] <= Hmaxj[j])

	   /* ** 
STEP 1.B: UPDATING of the STATUS of the initially admissible elements
                for factor A_j
      -> based on the ineligible elements with last-but-one non-zero element
         between ktr[j] and j-1 and with last non-zero element equal to j
      ---
      The IF condition just below is satisfied when kernel checks are necessary
      ** */

	   if (ktr[j] < fplusj+1) {
	     /* Make elements eligible if their status may have changed */

	     admiss_j =  IntegerVector(siaA[j]);

	     // admiss.j[ admiss.j >= ktr[j] ] <- s
	     for (i=0; i< admiss_j.size(); i++) {
	       if (admiss_j[i] >= ktr[j]) {
		 admiss_j[i] = *s;
	       }
	     } // end i


	     /*  Loop over the last-but-one non-zero indices */
	     for (int kk= ktr[j]; kk<= fplusj; kk++) {
	       /* relevant ineligible elements */
	       select_kj =  (ineligible_lnz(0,_) == kk) &(ineligible_lnz(1,_)== fplusj+1) ;

	       /* relevant initially admissible elements */
	       // tocheck <- seq_len(liaA[j])[admiss.j == s]
	       IntegerVector tocheck = findIndexVal(admiss_j, *s);


	       /* *** HM, 30/06/10 : condition sur length(tocheck) 
		  et sur inelig.kj *** */
	       if ((tocheck.size()>0) & is_true(any(select_kj))) {
		 // b.inelig.kj <- as.matrix(b.ineligible[ seq_len(f+j-1),                                                      select.kj, drop=FALSE ])
		 inelig_kj = ineligible(Range(0, fplusjmoins1), _);
		 inelig_kj = takeColLog(inelig_kj, select_kj);
		 candidates= takeCol(iiaA[j], tocheck);
		 iae_j =
			 kernelcheck_basep(PhiStar,  candidates,
					   inelig_kj, *p);

		 // admiss.j[ tocheck[!ae.j] ] <- k
		 IntegerVector iaefalse = findIndexLog(iae_j, false);
		 for (i=0; i< iaefalse.size(); i++) {
		   l = tocheck[iaefalse[i]];
		   admiss_j[l] =kk;
		 }

	       } // end (any(tocheck == true)

	       /* *** HM, 30/06/10 : fin condition sur length(tocheck) et sur inelig.kj *** */
	     } // end loop kk

	     siaA[j] = IntegerVector(admiss_j);

	     if (firstvisit[j] <= 2) {
	       if (*verbose) {
		 //  cat("    ---    col.", f+j, " ( j =",j, ") ", sum(admiss.j == s), " selected candidates\n")
		 LogicalVector vaux= (admiss_j  == *s);
		 int auxx=std::accumulate(vaux.begin(), vaux.end(), 0);
		 Rprintf("    ---    col. %d ( j = %d) %d selected candidates\n",
		     fplusj+1, j+1, auxx);

	       } 
	       firstvisit[j] = firstvisit[j]+1;
	     } // end (firstvisit[j] <= 2)

	   } // end (ktr[j] < *f+j)



	   /* ** STEP 1.C: UPDATED admissible elements for factor A_j ** */
	   // admiss.index.j <- seq_len(liaA[j])[siaA[[j]]==s]

	   IntegerVector siaAj=IntegerVector(siaA[j]);
	   LogicalVector ad = ( siaAj == *s);
	   admiss_index_j = findIndexLog(ad, true);

	   if (*randomsearch == 1) {
	     // Call the Rcpp function "runif"
	     int nindexsample = (int)admiss_index_j.size() - 1;
	     NumericVector index_sample= runif(nindexsample);
	     admiss_index_j = takeElt(admiss_index_j, index_sample);
	   } 
	   admiss_index_j = admiss_index_j + 1; // index from 1
	   wiaA[j] = admiss_index_j;
	   ktr[j] = fplusj+1;
   } // end (jprev < j)

  /* **
  END OF STEP 1

  STEP 2: INCLUSION OF THE NEXT ADMISSIBLE ELEMENT IN PhiStar
  ** */
  IntegerVector wj= IntegerVector(wiaA[j]);


  /*  CASE 2.A: If there is a next admissible element, add it to PhiStar ** */

  if (wj.size() > zero ) {
    // b.newcolj <- (b.iaA[[j]])[,wj[1]]
    int auxw = wj[0] -1; // index from 0
    IntegerMatrix iiaA1 = iiaA[j];
    // Note : Range is required to not build a column-matrix
    newcolj = iiaA1(_, Range(auxw, auxw));
    // wiaA[[j]] <- wj[-1]
    if (wj.size() > 1) {
      wiaA[j]=wj[Range(1, wj.size()-1)];
    } else {
      wiaA[j]=IntegerVector(0);
    }

    //   b.PhiStar <- cbind(b.PhiStar, b.newcolj, deparse.level=0)
     PhiStar=cbind(PhiStar, newcolj);

     /* **  CASE 2.A.a: If this was the last column to fill in PhiStar ** */
     if ( (j+1) == smoinsf) {
       /* Save the solution */
      PhiStar_solution[nsolution++] = PhiStar;
      /* Return if enough solutions have been found */
      if (nsolution == *max_sol) {
	Rprintf("The search is closed: max.sol = %d solution(s) found \n", nsolution);
	return wrap(PhiStar_solution[Range(0, nsolution-1)]);
      }
      jprev = j ;
    } // end  if(j == (s-f))


    else {
      /* ** CASE 2.A.b: If this was not the last column, go forward to next j ** */
      jprev = j ; j=j+1;
    }
  } // end if(length(wj) > 0)

  else {
     /* **  CASE 2.B: If there is no next admissible element ** */

    /* Update information */
    /* *** HM: modif 2015
       ktr[ (seq_len(s-f) >= j) & (ktr > 0) ] <- j-1
       changed into :
    //  ktr.change <- (seq(s-f) >= j) & (ktr > j-1)
    //  ktr[ ktr.change ] <- j-1 *** */
    for (int i=j; i< smoinsf; i++) {
      if (ktr[i] > j) {
	ktr[i] = j; // ktr contains indices from 1
      }
    }
    /* Then go backward */
    jprev = j;
    j =  j-1;
  }
 } // end while(j > 0)

 Rprintf("The search is closed: %d solution(s) found \n",
	 nsolution);

 if (nsolution == zero) {
  return(wrap(nosolution));
 } else {
   return wrap(PhiStar_solution[Range(0, nsolution-1)]);
 }

} // end PLANORdesignkeynonsym

/* ++++++++++++++++++++++++++++++++++++++
PLANORsubgroup
Interface from R
FUNCTION
calculates the non null elements of the subgroup H generated by the columns
of mat, considered as vectors in (Zp)^s
ARGUMENTS:
 matrice: a matrix of integers modulo p whose columns are assumed to
      be independent vectors in (Zp)^s
 p: a prime
 all: if TRUE all elements in H are given, if FALSE only elements up
      to multiplication by an integer modulo p are given
OUTPUT:
 a matrix of integers modulo p whose columns are the subgroup elements
DETAILS:
 it is not checked whether the column vectors of mat are independent, so there
 may be several times 
CALLED BY
 summary.keymatrix
++++++++++++++++++++++++++++++++++++++++ */
RcppExport SEXP PLANORsubgroup(SEXP Smatrice,
				  SEXP nrow,
				  SEXP ncol,
				  SEXP Sp,
				  SEXP Sall) {


  IntegerVector matrice(Smatrice);
  int * nr = INTEGER_POINTER(nrow);
  int * nc = INTEGER_POINTER(ncol);

 IntegerMatrix mamatrice = fillMat(matrice, *nr, *nc);


  int *pp = INTEGER_POINTER(Sp);
  int p=*pp;

  int *pall = INTEGER_POINTER(Sall);
  bool all = bool(*pall);


  IntegerMatrix outmat=subgroup( mamatrice, p, all);

  return(wrap(outmat));

} // end PLANORsubgroup


/* ++++++++++++++ INIT +++++++++++++++++++ */

static const R_CallMethodDef callMethods[] = {
  {"PLANORsubgroup", (DL_FUNC) &PLANORsubgroup, 5},
  {"PLANORinversesbasep", (DL_FUNC) &PLANORinversesbasep, 1},
  {"PLANORprintLib", (DL_FUNC) &PLANORprintLib, 5},
  {"PLANORdesignkeynonsym",  (DL_FUNC) &PLANORdesignkeynonsym, 12},
  {"PLANORdesignkeysym",  (DL_FUNC) &PLANORdesignkeysym, 12},
  {"PLANORkernelcheck", (DL_FUNC) &PLANORkernelcheck, 10},
  {"PLANORweightorder", (DL_FUNC) &PLANORweightorder, 5},
   {NULL, NULL, 0}
};

void   R_init_planor(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);

}
