library("planor")
cat("\n")
cat("***************** PLAN DE L'ARTICLE JTB09 *****************\n")
cat("********** Lurette, Touzeau, Lamboni, Monod ***********\n")
cat("\n")
cat("Eighteen 4-level treatment factors\n")
cat("N=2^12\n")
# 
cat("*** RUN ***\n")
cat("\n")

cat("resolution 4\n")
cat("\n")

ST.K <- planor.designkey(factors=LETTERS[1:18], nlevels=4,
                         model=~(A+B+C+D+E+F+G+H+I+J+K+L+M+N+O+P+Q+R),
                         nunits=2^12,
                       base=~A+B+C+D, max.sol=1)
ST.P <- planor.design(key=ST.K, select=1)
st4<- pick(ST.K,1)
# summary ne marche pas, car subgroup detecte une trop grosse
# big.matrix
#summary.designkey(st4)


