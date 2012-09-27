library("planor")
source("compar.R")

cat("\n")
cat("***************** PLAN DE L'ARTICLE JTB09 *****************\n")
cat("********** Lurette, Touzeau, Lamboni, Monod ***********\n")
cat("\n")
cat("Eighteen 4-level treatment factors\n")
cat("N=2^12\n")
# 
cat("*** RUN ***\n")
# ptmtotal <- proc.time()
# Rprof("ST2.prof")
cat("\n")


cat("resolution 2\n")
cat("\n")
#  ptm <- proc.time()
ST.K <- planor.designkey(factors=LETTERS[1:18], nlevels=2,
                         model=~(A+B+C+D+E+F+G+H+I+J+K+L+M+N+O+P+Q+R)^2 ,
                         nunits=2^12,
                       base=~A+B+C+D+E+F+G+H+I+J+K+L, max.sol=2)
# cat("*** TEMPS planor.designkey", proc.time()-ptm,"\n")
cat("\n")
ST.P <- planor.design(key=ST.K, select=2)
# cat("*** TEMPS design", proc.time()-ptm,"\n")
st2<- pick(ST.K,1)
# Increase print amount:
options(planor.max.print=100)
# ptm <- proc.time()
summary(st2)
# cat("*** TEMPS summary(pick)", proc.time()-ptm,"\n")
# cat("*** TEMPS total", proc.time()-ptmtotal,"\n")
#Rprof(NULL)

cat(" Comparer aux sorties initiales\n")
load("hmtestsSTouzeau2")
# Pbe pour appeler compar sur des big.matrix?
###  
###  slots <- slotNames(OR$ST.K)
###  for (s in slots) {
###      leslot = slot(ST.K,s)
###      ORleslot = slot(OR$ST.K,s)
###    if (s ==".Data") {
###      for (l in 1:length(leslot)) {
###        ll=leslot[[l]]
###        ORll=ORleslot[[l]]
###        for (k in 1:length(ll))
###        print(all.equal(ll[[k]][,], ORll[[k]][,]))
###      }
###    }
###    else {
###      print(all.equal(slotNames(leslot), slotNames(ORleslot)))
###      for (ss in slotNames(leslot)) {
###        s1 <- slot(leslot, ss)
###        s2 <- slot(ORleslot, ss)
###        print(all.equal(s1,s2))
###      }
###    }
###  }

print(compar(ST.P, OR$ST.P))

