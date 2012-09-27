\name{planor-package}
\docType{package}
\alias{planor-package}
\alias{planor}
\title{Generation of regular factorial designs}
\description{A package dedicated to the automatic generation of regular
  factorial designs, including fractional designs, orthogonal block designs,
  row-column designs and split-plots.}
\details{The user describes the factors to be controlled in the
  experiment and the anova model to be used when the results will be
  analysed. He or she also specifies the size of the design, that is,
  the number of available experimental units. Then \pkg{planor} looks for a
  design satisfying these specifications and possibly randomizes it. The
  core of the algorithm is the search for the key matrix, an integer
  matrix which determines the aliasing in the resulting factorial
  design.

  The user may use the function \code{\link{regular.design}} where all
  these steps are integrated, and transparent by default. Alternatively,
  the steps can be decomposed by using successively the functions
  \code{\link{planor.factors}}, \code{\link{planor.model}},
  \code{\link{planor.designkey}} and
  \code{\link[=planor.design-method]{planor.design}}. For the expert
  user, the function \code{\link{planor.designkey}} can give several key
  matrix solutions. Alias and summary methods allow to study and compare
  these solutions, in order to select the most appropriate one for the
  final design.}

\note{\itemize{
 \item{The \bold{bigmemory} package should have been installed
  previously to the \bold{planor} installation ; if installed in a non
  standard library, \bold{planor} must be installed in the same
  directory}
 \item{Option \code{planor.max.print} can be set to limit the amount of
  the matrices that are printed, to \code{planor.max.print} rows and
  columns. Default: 20.}}}
  
\author{Hervé Monod, Annie Bouvier, André Kobilinsky
  (Applied Mathematics and Informatics Unit, INRA UR 341 -
  Jouy-en-Josas, France. URL: \url{http://www.jouy.inra.fr/mia_eng/})}

\references{ see \code{citation(planor)} }

\keyword{package}
\keyword{design}

\examples{
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
}
