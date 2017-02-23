IntegerMatrix fillMat(IntegerVector mat, int nr, int nc);
int findFnzIndex (IntegerVector vec);
LogicalVector findFnzOne(IntegerMatrix mat);
int findFOneIndex (IntegerVector vec);
IntegerVector findIndexLog(LogicalVector vec, bool val);
IntegerVector findIndexVal(IntegerVector vec, int val);

IntegerMatrix moduloMat (IntegerMatrix mat, const int p);

NumericMatrix moduloMatNum (NumericMatrix mat, const int p);
// DOES NOT WORK template <typename T> T moduloMat( T mat, const int p);

IntegerVector moduloVec (IntegerVector x, const int p);
IntegerMatrix multMatMod(IntegerMatrix matrice,
			 IntegerMatrix coeffs,
			 int p);
IntegerVector myiota(IntegerVector vec, int p);
IntegerVector sequence( int first, int last);
IntegerMatrix takeCol(IntegerMatrix mat,
		      IntegerVector indices);
IntegerMatrix takeColLog(IntegerMatrix mat,
			 LogicalVector indices);
IntegerVector takeElt(IntegerVector vec,
		      NumericVector indices);
IntegerMatrix takeRowLog(IntegerMatrix mat,
			 LogicalVector indices);
/* ++++++++++++++++++++++++++++++++++++++++++++++++
MACRO imax2, imin2 pas trouvé dans R.h (?)
++++++++++++++++++++++++++++++++++++++++++++++++ */
#define MAX(a, b) ( a > b ?a :b)
#define MIN(a, b) ( a < b ?a :b)

