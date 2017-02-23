###################################################################
# planor R package
# Copyright INRA 2017
# INRA, UR1404, Research Unit MaIAGE
# F78350 Jouy-en-Josas, France.
#
# URL: http://www3.jouy.inra.fr/miaj/public/logiciels/planor/
#
# This file is part of planor R package.
# planor is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# See the GNU General Public License at:
# http://www.gnu.org/licenses/
#
###################################################################
makedesignkey <- function(keys, primes){
## A function to create an object of class designkey directly
## from a list of simple matrices
  ## ARGUMENTS
  ## - keys : a list of n (integer) matrices
  ## - primes : a vector of n prime numbers
  ## OUTPUT
  ## an object of class designkey
  ## DETAILS
  ## the names of the factors are extracted from the matrix column names
  ## the associated model formula is the additive model with all factors
  ## the estimate formula is the constant
  ## EXAMPLE
  ## mat1 <- cbind(diag(3),1)
  ## colnames(mat1) <- c("A","B","C","D")
  ## mat2 <- cbind(diag(2),c(1,2))
  ## colnames(mat2) <- c("E","F","G")
  ## mat.dk <- makedesignkey(list(mat1), primes=c(2))
  ## print(mat.dk)
  ## summary(mat.dk)
  ## alias(mat.dk)
  ## mat.plan <- planor.design(mat.dk)
  ## -----------------------------------------------------
  ## information extraction
  fact.names <- unlist( lapply(keys, colnames) )
  fact.levels <- rep( primes, unlist(lapply(keys, ncol)) )
  tp.fact <- planor.factors(factors=fact.names, nlevels=fact.levels)
  tp.fact@fact.info$model <- TRUE
  tp.fact@pseudo.info$model <- TRUE
  ##
  nunits <- prod( primes^unlist(lapply(keys, nrow)) )
  
  if (nunits > .Machine$integer.max) {
        stop(paste("makedesignkey. Overflow.", nunits, "greater than the maximum integer", .Machine$integer.max))
      }
  

  ##
  names(keys) <- as.character(primes)
  for(p in seq(keys)){
    keys[[p]] <- new("keymatrix", keys[[p]], p=primes[p])
  }
  ##
  tp.model <- as.formula(paste("~", paste(fact.names,collapse="+"), sep=""))
  tp.estim <- as.formula(paste("~", paste(fact.names,collapse="+"), sep=""))
  tp.mod <- planor.model(model=tp.model, estimate=tp.estim)

  ## KEY GENERATION
  newkey <- new("designkey",
                .Data=keys,
                factors=tp.fact,
                nunits=nunits,
                model=tp.mod)

  return(newkey)
}
