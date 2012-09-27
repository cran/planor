pkgname <- "planor"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('planor')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("alias-methods")
### * alias-methods

flush(stderr()); flush(stdout())

### Name: alias-methods
### Title: Methods for Function 'alias': summarize the design properties
### Aliases: alias-method alias alias,designkey-method
###   alias,keymatrix-method alias,listofdesignkeys-method
###   alias,listofkeyrings-method
### Keywords: methods design

### ** Examples

### Creation of an object of class "listofkeyrings"
K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
nunits=3^3, base=~A+B+C, max.sol=2)
### alias on an object of class "keymatrix"
alias(K0[[1]][[1]])
### alias on an object of class "designkey"
alias(K0[1])
### alias on an object of class "listofkeyrings"
alias(K0)



cleanEx()
nameEx("bind-methods")
### * bind-methods

flush(stderr()); flush(stdout())

### Name: bind-methods
### Title: Methods for function 'bind': bind two objects
### Aliases: bind-method bind bind,designfactors,designfactors-method
### Keywords: methods design

### ** Examples

F1 <- planor.factors(factors=c("block",LETTERS[1:4]), nlevels=c(6,6,4,2,6))
F2 <- planor.factors(factors=c("block",LETTERS[11:12]), nlevels=c(6,6,4))
### Method bind on 'designfactors' objects
F3 <- bind(F1,F2)
names(F3)



cleanEx()
nameEx("designfactors-class")
### * designfactors-class

flush(stderr()); flush(stdout())

### Name: designfactors-class
### Title: Class '"designfactors"' and methods of the class
### Aliases: designfactors-class [,designfactors,ANY,ANY,ANY-method
###   length,designfactors-method names,designfactors-method
### Keywords: classes design

### ** Examples

F1 <- planor.factors(factors=c("block",LETTERS[1:4]), nlevels=c(6,6,4,2,6))
F2 <- planor.factors(factors=c("block",LETTERS[11:12]), nlevels=c(4,6,6))
## Method bind - see the warning because two factors in F1 and F2 have
## the same name
F3 <- bind(F1,F2) 
names(F3)
length(F3)
F3@levels
F3.trt <- F3[c(2:5,7,8)]
names(F3.trt)



cleanEx()
nameEx("designkey-class")
### * designkey-class

flush(stderr()); flush(stdout())

### Name: designkey-class
### Title: Class "designkey" and methods of the class
### Aliases: designkey-class
### Keywords: classes design

### ** Examples

### Creation of a 'designkey' object
K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
nunits=3^3, base=~A+B+C, max.sol=2)
print(K0[1])



cleanEx()
nameEx("getDesign-methods")
### * getDesign-methods

flush(stderr()); flush(stdout())

### Name: getDesign-methods
### Title: Methods for function 'getDesign': extract a design
### Aliases: getDesign-method getDesign getDesign,planordesign-method
### Keywords: methods design

### ** Examples

### Creation of a 'planordesign' object
K0 <- planor.designkey(factors=c("R","C","U","A","B1","B2"),
 nlevels=c(3,2,2,3,2,2), model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2,
 nunits=12, base=~R+C+U, max.sol=2)
P0 <- planor.design(key=K0, select=1)
## Method getDesign on the 'planordesign' object
 show(getDesign(P0))



cleanEx()
nameEx("keymatrix-class")
### * keymatrix-class

flush(stderr()); flush(stdout())

### Name: keymatrix-class
### Title: Class "keymatrix" and methods of the class
### Aliases: keymatrix-class
### Keywords: classes design

### ** Examples

showClass("keymatrix")
### Creation of a 'listofkeyrings' object
K0 <- planor.designkey(factors=c("block", LETTERS[1:4]), nlevels=rep(3,5),
  model=~block + (A+B+C+D)^2, estimate=~A+B+C+D,
  nunits=3^3, base=~A+B+C, max.sol=2)
# Method show on a 'keymatrix' of K0
show(K0[[1]][[1]])



cleanEx()
nameEx("keyring-class")
### * keyring-class

flush(stderr()); flush(stdout())

### Name: keyring-class
### Title: Class "keyring" and methods of the class
### Aliases: keyring-class
### Keywords: classes design

### ** Examples

showClass("keyring")
### Creation of a 'listofkeyrings' object
K0 <- planor.designkey(factors=c("block", LETTERS[1:4]), nlevels=rep(3,5),
  model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
  nunits=3^3, base=~A+B+C, max.sol=2)
## Method show applied on a 'keyring' component of K0
show(K0[[1]]) 



cleanEx()
nameEx("listofdesignkeys-class")
### * listofdesignkeys-class

flush(stderr()); flush(stdout())

### Name: listofdesignkeys-class
### Title: Class "listofdesignkeys" and methods of the class
### Aliases: listofdesignkeys-class [,listofdesignkeys,ANY,ANY,ANY-method
### Keywords: classes design

### ** Examples

showClass("listofdesignkeys")
### Creation of a "listofdesignkeys" object
K0 <- planor.designkey(factors=c("R","C","U","A","B1","B2"),
 nlevels=c(3,2,2,3,2,2), model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2,
 nunits=12, base=~R+C+U, max.sol=2)
# Show the object
show(K0)
## Method length
length(K0)
## Extraction: the following two commands are equivalent 
K <- K0[2]
K <- pick(K0,2)



cleanEx()
nameEx("listofkeyrings-class")
### * listofkeyrings-class

flush(stderr()); flush(stdout())

### Name: listofkeyrings-class
### Title: Class "listofkeyrings" and methods of the class
### Aliases: listofkeyrings-class [,listofkeyrings,ANY,ANY,ANY-method
### Keywords: classes design

### ** Examples

showClass("listofkeyrings")
### Creation of a 'listofkeyrings' objct
K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
   model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
   nunits=3^3, base=~A+B+C, max.sol=2)
show(K0)



cleanEx()
nameEx("pick-methods")
### * pick-methods

flush(stderr()); flush(stdout())

### Name: pick-methods
### Title: Methods for function 'pick': extract a single 'designkey' from a
###   list
### Aliases: pick-method pick pick,listofkeyrings-method
###   pick,listofdesignkeys-method
### Keywords: methods design

### ** Examples

# Creation of an object of class "listofdesignkeys"
K2 <- planor.designkey(factors=c("R","C","U","A","B1","B2"),
nlevels=c(3,2,2,3,2,2),  model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2 , nunits=12,
base=~R+C+U, max.sol=2)
# Method 'pick' applied on the "listofdesignkeys" object
K2.1 <- pick(K2,1)
K2.1 <- K2[1] ## Another way of extracting ([ is synonym of pick)

# Creation of an object of class "listofkeyrings"
K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"),
nlevels=rep(3,5), model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
nunits=3^3, base=~A+B+C, max.sol=2)
# Method 'pick' applied on the "listofkeyrings" object
K0.1 <- pick(K0,1)
K0.1 <- K0[1] ## the same



cleanEx()
nameEx("planor-package")
### * planor-package

flush(stderr()); flush(stdout())

### Name: planor-package
### Title: Generation of regular factorial designs
### Aliases: planor-package planor
### Keywords: package design

### ** Examples

# DESIGN SPECIFICATIONS
# Treatments: four 3-level factors A, B, C, D
# Units: 27 in 3 blocks of size 9
# Non-negligible factorial terms:
#   block + A + B + C + D + A:B + A:C + A:D + B:C + B:D + C:D
# Factorial terms to estimate:
#   A + B + C + D
# 1. DIRECT GENERATION, USING 'regular.design'
mydesign <- regular.design(factors=c("block", LETTERS[1:4]),
  nlevels=rep(3,5), model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
  nunits=3^3, randomize=~block/UNITS)
print(mydesign)
# DUMMY ANALYSIS
# Here we omit two-factor interactions from the model, so they are 
# confounded with the residuals (but not with ABCD main effects)
set.seed(123)
mydesign$Y <- runif(27)
mydesign.aov <- aov(Y ~ block + A + B + C + D, data=mydesign)
summary(mydesign.aov)
# 2. STEP-BY-STEP GENERATION, USING 'planor.designkey'
F0 <- planor.factors(factors=c( "block", LETTERS[1:4]), nlevels=rep(3,5),
  block=~block)
M0 <- planor.model(model=~block+(A+B+C+D)^2, estimate=~A+B+C+D) 
K0 <- planor.designkey(factors=F0, model=M0, nunits=3^3, max.sol=2)
summary(K0)
mydesign.S4 <- planor.design(key=K0, select=2)



cleanEx()
nameEx("planor.design-methods")
### * planor.design-methods

flush(stderr()); flush(stdout())

### Name: planor.design-methods
### Title: Methods for function 'planor.design': build a design from a
###   design-key solution
### Aliases: planor.design-method planor.design
###   planor.design,designkey-method planor.design,listofdesignkeys-method
###   planor.design,listofkeyrings-method planor.design,numeric-method
### Keywords: methods design

### ** Examples

### Creation of a 'listofdesignkeys' object
K0 <- planor.designkey(factors=c("R","C","U","A","B1","B2"),
 nlevels=c(3,2,2,3,2,2), model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2,
 nunits=12, base=~R+C+U, max.sol=2)
## Method planor.design applied on the 'listofdesignkeys' object
P0 <- planor.design(key=K0, select=1)
## Method planor.design applied on a designkey' object
P0 <- planor.design(K0[1])


### Creation of a 'listofkeyrings' object
K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
nunits=3^3, base=~A+B+C, max.sol=2, verbose=TRUE)
## Method planor.design applied on a designkey' object
P0 <- planor.design(K0[1])
P0.R <- planor.design(K0[1], randomize=~A+B+C+D) ## randomize the final design


cleanEx()
nameEx("planor.designkey")
### * planor.designkey

flush(stderr()); flush(stdout())

### Name: planor.designkey
### Title: Search for a design key or a collection of design keys
### Aliases: planor.designkey
### Keywords: design

### ** Examples

K0 <- planor.designkey(factors=c("block", LETTERS[1:4]),
  nlevels=rep(3,5), model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
  nunits=3^3, base=~A+B+C, max.sol=2)
## With automatic model generation
Km <- planor.designkey(factors=c("block", LETTERS[1:4]),
  nlevels=rep(2,5), resolution=3, nunits=2^4)



cleanEx()
nameEx("planor.factors")
### * planor.factors

flush(stderr()); flush(stdout())

### Name: planor.factors
### Title: Create an object of class 'designfactors'
### Aliases: planor.factors
### Keywords: design

### ** Examples
planor.factors(c("A","B","C","P"),c(2,3,6,3))
planor.factors(LETTERS[1:12],2)
planor.factors(12,2)
planor.factors( c("A","B","Block"), 3, block=~Block )
zz <- planor.factors( c("A","B","Block"), c(2,3,5))
zz@levels$A <- c("plus","moins")
planor.factors(factors=list(A=c("plus","moins"), B=1:3, Block=1:5))
AB <- data.frame( A=c(rep(c("a","b"),3)), B=rep(c("z","zz","zzz"),rep(2,3)), C=1:6  )
planor.factors(factors=AB)


cleanEx()
nameEx("planor.harmonize")
### * planor.harmonize

flush(stderr()); flush(stdout())

### Name: planor.harmonize
### Title: Harmonize the factors
### Aliases: planor.harmonize
### Keywords: design design

### ** Examples

F2 <- planor.factors(factors=c("block",LETTERS[1:4]), nlevels=c(6,6,6,4,2))
M2 <- planor.model( model=~block+(A+B+C)^2, estimate=~A+B+C )
F2.h <- planor.harmonize(factors=F2, model=M2, base=~A+B)
names(F2)
names(F2.h)



cleanEx()
nameEx("planor.model")
### * planor.model

flush(stderr()); flush(stdout())

### Name: planor.model
### Title: Model and estimate specifications for a design search
### Aliases: planor.model
### Keywords: design

### ** Examples

# Basic example
planor.model(model=~block + (A+B+C)^2, estimate=~(A+B+C)^2)
# Resolution: both calls to 'planor.model' below are equivalent
planor.model(model=~(A+B+C+D)^2, estimate=~A+B+C+D)
myfactors <- planor.factors(factors=c(LETTERS[1:4]), nlevels=rep(2,4))
planor.model(resolution=4, factors=myfactors)
# Complicated examples
planor.model(~A+B+C+D+A:B, ~A+B+C+D, listofmodels=list(c(~E+F,~E)))
planor.model(~A+B+C+D+A:B,~A+B+C+D, listofmodels=
                              list(c(~E+F,~E), ~G, ~H, c(~M+N,~N)))



cleanEx()
nameEx("planor.randomize")
### * planor.randomize

flush(stderr()); flush(stdout())

### Name: planor.randomize
### Title: A function to randomize a factorial design according to an
###   orthogonal block structure
### Aliases: planor.randomize
### Keywords: design

### ** Examples
## Block design
Design <- data.frame(block=rep(1:4,rep(2,4)),
treatment=c("A1","B1","A2","B2","A3","B3","A4","B4"))
planor.randomize(~block, data=Design)       ##  no within-block randomization
planor.randomize(~block/UNITS, data=Design) ##  blocks and units within blocks randomization
## Row-Column design
RowColDes <- data.frame(row=rep(1:3,rep(3,3)),col=rep(1:3,3),
treatment=LETTERS[c(1:3,2,3,1,3,1,2)],
oldRow=rep(1:3,rep(3,3)),oldCol=rep(1:3,3))
planor.randomize(~row*col, data=RowColDes)


cleanEx()
nameEx("planordesign-class")
### * planordesign-class

flush(stderr()); flush(stdout())

### Name: planordesign-class
### Title: Class "planordesign" and methods of the class
### Aliases: planordesign-class [,planordesign,ANY,ANY,ANY-method
### Keywords: classes design

### ** Examples

showClass("planordesign")
### Creation of a 'listofdesignkeys' object
K0 <- planor.designkey(factors=c("R","C","U","A","B1","B2"),
 nlevels=c(3,2,2,3,2,2), model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2,
 nunits=12, base=~R+C+U, max.sol=2)
## Creation of a 'planordesign' object from K0
P0 <- planor.design(key=K0, select=1)
show(P0)



cleanEx()
nameEx("regular.design")
### * regular.design

flush(stderr()); flush(stdout())

### Name: regular.design
### Title: Construct and randomize a regular factorial design
### Aliases: regular.design
### Keywords: design

### ** Examples
mydesign <- regular.design(factors=c("block", LETTERS[1:4]),
  nlevels=rep(3,5), model=~block + (A+B+C+D)^2, estimate=~A+B+C+D,
  nunits=3^3, randomize=~block/UNITS)
print(mydesign)



cleanEx()
nameEx("show-methods")
### * show-methods

flush(stderr()); flush(stdout())

### Name: show-methods
### Title: Methods for function 'show' in package 'planor'
### Aliases: show-method show,designkey-method show,keymatrix-method
###   show,keyring-method show,listofdesignkeys-method
###   show,listofkeyrings-method
### Keywords: methods design

### ** Examples

# Creation of a "listofdesignkeys" object
K0 <- planor.designkey(factors=c("R","C","U","A","B1","B2"),
 nlevels=c(3,2,2,3,2,2), model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2,
 nunits=12, base=~R+C+U, max.sol=2)
## Method show applied on a "keymatrix" object
show(K0[[1]][[1]])
## Method show applied on a "designkey"  object
show(K0[1])
## Method show applied on the "listofdesignkeys" object
show(K0)
K0 # the same

### Creation of a "listofkeyrings" object
K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
nunits=3^3, base=~A+B+C, max.sol=2)
## Method show applied on a "keyring" object
show(K0[[1]]) 
print(K0[[1]]) # the same
K0[[1]] # the same
## Method show applied on the "listofkeyrings" object
show(K0)



cleanEx()
nameEx("summary-methods")
### * summary-methods

flush(stderr()); flush(stdout())

### Name: summary-methods
### Title: Methods for function 'summary' in package 'planor'
### Aliases: summary-method summary,designkey-method
###   summary,keymatrix-method summary,keyring-method
###   summary,listofdesignkeys-method summary,listofkeyrings-method
### Keywords: methods

### ** Examples

### Creation of a "listofdesignkeys" object
K0 <- planor.designkey(factors=c("R","C","U","A","B1","B2"),
  nlevels=c(3,2,2,3,2,2), model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2,
  nunits=12, base=~R+C+U, max.sol=2)
## Method summary applied on a "keymatrix" object
r <- summary(K0[[1]][[1]])
## Method summary applied on a "designkey"  object
summary(K0[1], save=NULL)
# Method summary applied on the "listofdesignkeys" object
r <-summary(K0, show="dt")

### Creation of a "listofkeyrings" object
K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
   model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
   nunits=3^3, base=~A+B+C, max.sol=2)
# Method summary applied on the "keymatrix" object
r <-summary(K0[[1]][[1]])
# Method summary applied on the "keyring" object
r <-summary(K0[[1]])
# Method summary applied on the "listofkeyrings" object
r <- summary(K0, show="dtb", save ="k")
print(r)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
