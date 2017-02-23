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
#---------------------------------------------------------------------------
#   UTILITIES FOR PLANOR, in particular matrix algebra modulo p,
#   for p a prime
#---------------------------------------------------------------------------
# 1. Miscellaneous
# crossing <- function(n,start=1)
# symmdiff <- function(x,y)
# cross.designs <- function(designs)
# 2. Matrix algebra modulo p, with p a prime
# convertinto.basep <- function (x, p)
# convertfrom.basep <- function (x, p)
# representative.basep <- function(mat,p)
# kernelmatrix.basep <- function(mat,p)
#---------------------------------------------------------------------------
#  1. MISCELLANEOUS
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
crossing <- function(n,start=1){
  # Generates all n1 x n2 x ... x ns combinations of size s with n1,...,ns integers
  # ARGUMENTS
  # - n: a vector of integers of length s
  # - start: integer from where to start the series of integers
  # RETURN
  #  an integer matrix with prod(n) rows and s columns giving all
  #  combinations along the rows, in lexicographic order
  # EXAMPLE
  #  crossing(rep(2,3)) # Internal
  # ----------------------------------------------------------
  N <- prod(n)
  s <- length(n)
  
  if (N > .Machine$integer.max) {
    stop(paste("crossing. Overflow.", N, "greater than the maximum integer", .Machine$integer.max))
  }

  n <- c(n,1)
  crosses <- matrix(NA, N, s)
  for(i in seq_len(s))
    {
      motif <- start + seq_len(n[s+1-i])-1
      repet1 <- rep( prod(n[s+1-i+seq_len(i)]), n[s+1-i] )
      if(i==s){ repet2 <- 1 }
      else{ repet2 <- prod(n[seq_len(s-i)]) }
      crosses[,s-i+1] <- rep( rep( motif, repet1 ), repet2 )
    }
  return(crosses)
}   # end crossing
#---------------------------------------------------------------------------
symmdiff <- function(x,y){
  # Calculates the symmetric differences between two collections Cx and Cy of
  # subsets of S, where S is a set of size n
  # cf. the Delta operation in Kobilinsky 2000, p.5
  # ARGUMENTS
  #  - x: a n-row 0-1 matrix, each column of which defines a subset in Cx
  #  - y: a n-row 0-1 matrix, each column of which defines a subset in Cy
  # RETURN
  #  a n-row 0-1 matrix with one column per distinct symmetric difference
  #  between a subset in Cx and a subset in Cy
  # EXAMPLE: (Internal)
  #  M1 <- diag(3)
  #  M2 <- cbind( c(1,1,0), c(1,0,1) )
  #  a <-symmdiff(M1,M2)
  # print(a)
  # -----------------------------------------------------------------

  Nf <- nrow(x)
  Nx <- ncol(x)
  Ny <- ncol(y)
  IDeltaJ <- array(NA, dim=c(Ny, Nx, Nf))
  for (i in seq_len(Ny)){
    for(j in seq_len(Nx)){
      IDeltaJ[i,j,] <- y[,i] != x[,j]
    }
  }
  dim(IDeltaJ) <- c( prod(dim(IDeltaJ)[c(1,2)]), Nf)
  z.codes <- unique( convertfrom.basep(1*IDeltaJ,2) )
  z.codes <- z.codes[z.codes!=0]
  b.z <- t(convertinto.basep(z.codes,2))
  rownames(b.z) <- rownames(x)
  return(b.z)
} # fin symmdiff

#---------------------------------------------------------------------------
cross.designs <- function(designs){
  # generates a design by crossing one, two or more subdesigns
  # ARGUMENTS
  #  designs: a list of design matrices
  # RETURN
  #  the new design  matrix;
  # EXAMPLE
  #  pl1 <- crossing(c(2,3))
  #  colnames(pl1) <- c("A","B")
  #  pl2 <- crossing(5)
  #  colnames(pl2) <- c("C")
  #  a <- cross.designs(list(pl1,pl2))
  # print(a)
# -------------------------------------------------------

  nrows <- unlist(lapply(designs,nrow))
  ncols <- unlist(lapply(designs,ncol))
  prod.nrows <- prod(nrows)
  sum.ncols <- sum(ncols)
  b.crossdesign <- matrix(0, nrow= prod.nrows, ncol=sum.ncols)


  
  
  colnames(b.crossdesign) <- as.character(rep("NA", sum.ncols))
  i.col <- 0
  for(i in seq_along(designs) ){
    il <- seq_len(i)
    rowindices <- rep( seq_len(nrows[i]), rep(prod.nrows/prod(nrows[il]),nrows[i]) )
    rowindices <- rep( rowindices, prod(nrows[il])/nrows[i] )
    colindices <- sum(ncols[il]) - ncols[i] + seq_len(ncols[i])
    b.crossdesign[, colindices] <- designs[[i]][rowindices,]
    if(!is.null(colnames(designs[[i]]))){
      colnames(b.crossdesign)[colindices] <- colnames(designs[[i]])
    }
  }

  storage.mode(b.crossdesign)  <- "integer"
  

  
  return(b.crossdesign)
}  # end cross.designs
#---------------------------------------------------------------------------
# 2. MATRIX ALGEBRA MODULO p, with p a prime
#---------------------------------------------------------------------------


#---------------------------------------------------------------------------
go1convertinto.basep <- function (x, p) {
 # Conversion of an integer  into base p
 # Called by convertinto.basep and goconvertinto.basep
 # RETURN
 # c(n1, n2, .. nx) such as x=n1+n2*p+n3*p**2+n4*p**3
  # -------------------------------------------------

if (x != round(x) || x < 0)
    return(x)
  val <- x%%p
  while ( (x <- x%/%p) > 0 ) {
    newval <- x%%p
    val <- c(val,newval)
  }
  
  return(val)
} # end go1convertinto.basep



#---------------------------------------------------------------------------
goconvertinto.basep <- function (x, p) {
  # called by convertinto.basep: loop when its argument is vector
  # Call go1convertinto.basep on  max(x) to know
  # the maximal number of columns of the result

  
  
  
  
  
  # ----------------------------------------------------------------

    l <- matrix(0, length(x), length(go1convertinto.basep(max(x),p)))
    for(i in seq_along( x)){
      dec.i <- go1convertinto.basep(x[i],p)
      ## dec.i is an integer64. Conversion into integer
      ## to make 0 and 1 exactly 0 and 1
      u <- as.integer(dec.i)
      l[i, seq_along(dec.i) ] <- u
    }


    return(l)
  } # end goconvertinto.basep


#---------------------------------------------------------------------------

convertinto.basep <- function (x, p) {
  # Conversion of an integer or integer vector x into base p
  # The coefficients are ordered by increasing powers of p
  # ARGUMENTS
  #  - x: an integer or an integer variate
  #  - p: a prime
  # RETURN
  #  if x is an integer, a vector of elements in Zp
  #  if x is a vector, a matrix with length(x) rows
  #       and the adequate number of columns
  # EXAMPLE
  #  convertinto.basep(1:10, 3)
  # -------------------------------------------------------

  if (length(x) > 1)
    return(goconvertinto.basep(x, p))
  else
     return(go1convertinto.basep(x, p))
} # end convertinto.basep

#---------------------------------------------------------------------------

#---------------------------------------------------------------------------

goconvertfrom.basep <- function (x, p) {
  # Function called by convertfrom.basep when x is a vector
  # See convertfrom.basep
  
  ret <- sum( as.integer64(x * p^(seq_along(x)-1) ))

  if (is.na(ret)) {
    stop(paste("convertfrom.basep. Overflow. ", ret,
               " is greater than the maximum integer ",
               ". Generation of a number from a sequence including ", p,
               "^", length(x)-1, " not possible.", sep=""))
  }

  return(ret)
}# end goconvertfrom.basep





#---------------------------------------------------------------------------
convertfrom.basep <- function (x, p) {
  # Conversion of integers x coded as vectors of coefficients in base p
  # to classical integers in base 10
  # ARGUMENTS
  #  - x: a vector of coefficients ordered by increasing p powers
  #     or a n-row matrix with each row interpreted as such a vector;
  #     the elements of x must be between 0 and p-1
  #  - p: a prime
  # RETURN
  #  a scalar integer if x is a vector, or an integer n-vector if x is a matrix
  # EXAMPLE:
  #  vec3 <- convertinto.basep(1:10, 3)
  #  convertfrom.basep( vec3, 3 )
  
  
  
  
  
  
  




# ------------------------------------------------------------
  if (is.matrix(x)) {
    
    
    l <- integer64( nrow(x))
    
    for(i in seq_along( l)){
      l[i] <- goconvertfrom.basep(x[i,],p)
    }
  }
  else {
    l <- goconvertfrom.basep(x,p)
  }
 
    return(l)

} # fin convertfrom.basep


#---------------------------------------------------------------------

#---------------------------------------------------------------------



representative.basep <- function(b.mat,p){
  # generates the minimal set of representatives in base p
  # of the columns x of matrix mat
  # ARGUMENTS
  #  - b.mat : a 0-1  matrix
  #  - p : a prime
  # RETURN
  #  a matrix whose columns are the minimal representatives, with 1
  #  in the last non-zero position of x and all possible combinations
  #  of 1,...,p-1 in the other non-zero positions of x
  # EXAMPLE
  #  a <-representative.basep( as.matrix(c(1,0,1,1,0)), 3 )
  # print(a)
  #  a <- representative.basep( as.matrix(cbind( c(1,1,0),c(1,1,1) )), 3 )
  # ------------------------------------------------------

  if(p==2) return(b.mat %% 2)

  b.representative <- NULL

  for(j in seq_len(ncol(b.mat))){
    x <- b.mat[,j]
    select <- seq_along(x)[x != 0]
    nbtocross <- length(select)-1
    if( nbtocross <= 0 ) mat.j <- as.matrix(x)
    else{
      select <- select[seq_len(nbtocross)]
      N <- (p-1)^nbtocross
      
      if (N > .Machine$integer.max) {
        stop(paste("representative.basep. Overflow.", N, "greater than the maximum integer", .Machine$integer.max))
      }

      mat.j <- matrix(x, nrow(b.mat), N)
      mat.j[select,] <- t( crossing(rep(p-1,nbtocross),start=1) )
    }
    b.representative <- cbind(b.representative, mat.j)
  } # end j
  ret <- b.representative %%p

  storage.mode(ret) <- "integer"

   return(ret)
} # end representative.basep
#---------------------------------------------------------------------------
kernelmatrix.basep <- function(mat,p){
  # Calculates a kernel basis of a p-morphism
  # ARGUMENTS
  #  - mat : the n x s matrix associated with the p-morphism, with elements in Zp;
  #  - p : a prime
  # RETURN
  #  an s x q  matrix of q independent kernel generators, with elements in Zp;
  #  if mat is a regular matrix, then the output is an s x 0 simple matrix
  # DETAILS
  #  done by an adaptation of the Gauss-Jordan elimination algorithm
  # EXAMPLE
  #  mat[,] <- matrix(
  #  c(1,1,0,0, 2,0,1,0,  0,1,2,4, 0,1,2,3, 0,0,0,1) , nrow=4)
  #  matker <- kernelmatrix.basep(mat,5)
  #   Note: mat %*% matker[,] %%5 is a null matrix
# -------------------------------------------------------------------

  nr <- nrow(mat)
  nc <- ncol(mat)
  pseudonames <- colnames(mat)
  mat <- mat %%p
  # special cases when mat=Id or mat=(Id|B): direct solutions
  if( nr == nc ){
    if( all( mat[,seq_len(nr)]==diag(nr) ) ){
      b.mat.kernel <- matrix(0,nr,0)
      storage.mode(b.mat.kernel) <- "integer"      
      # matrix with 0 column
      rownames(b.mat.kernel) <- pseudonames
      return(b.mat.kernel) }}
  if( nr < nc ){
    if( all( mat[,seq_len(nr)]==diag(nr) ) ){
  b.mat.kernel <- rbind( -mat[,seq(nr+1,nc),drop=FALSE], diag(nc-nr) )%%p
      storage.mode(b.mat.kernel) <- "integer"
      rownames(b.mat.kernel) <- pseudonames
      return(b.mat.kernel) }}
  # general case: elimination algorithm
  # - initialisation
  
  # Raw calculation of the inverses modulo p
  invp <- .Call("PLANORinversesbasep", as.integer(p))
  
  
  firstgoodrow <- function(column, rowindex){
    # specific function to indicate the position of the first non-0 element
    # above position "rowindex" in a vector "column"
    # equal to NA if there is no non-0 element above "rowindex" in "column"
    rows <- seq_along(column)
    goodrows <- ( rows >= rowindex) & (column != 0)
    if(sum(goodrows) == 0) firstgood <- NA
    else firstgood <- rows[goodrows][1]
    return(firstgood)
  }
  rowindices <- seq_len(nr)
  colindices <- seq_len(nc)
  i <- 1  ;  continue <- TRUE
  # - main "while" loop
  while(continue){
    # zz : first next row with a non-zero value in the current column
    zz <- firstgoodrow(column=mat[,colindices[i]], rowindex=i)
    # if zz is NA: try next column ...
    if(is.na(zz)){
      if(colindices[i] < nc){
        colindices[i:nc] <- colindices[c((i+1):nc, i)] }
      # ... or terminate main loop if there is no new next column
      else {
        ifinal <- i-1
        continue <- FALSE }
    }
    # if zz is NOT NA: elimination step
    else{
      # row exchange if needed
      if(i != zz)  mat[c(i,zz),] <- mat[c(zz,i),]
      # normalisation
# Raw calculation of the inverses modulo p
      inv <- invp[mat[i,colindices[i]]]
      mat[i,] <- (inv * mat[i,]) %% p
      # orthogonalisation
      for(ii in rowindices[rowindices != i]){
        aik <- mat[ii,colindices[i]]
        mat[ii,] <- (mat[ii,] - aik * mat[i,]) %% p
      }
      # next row if possible, otherwise terminate main loop
      if(i < min(nr,nc)){
        i <- i+1 }
      else{
        ifinal <- nr
        if(nr < nc){ colindices[(nr+1):nc] <- colindices[sort((nr+1):nc)] }
        continue <- FALSE }
    }
  }
  # end of the main "while" loop
  #
  # if there remain columns, construct the kernel generators

  if(ifinal < nc){
    b.mat.kernel <- matrix(NA, nrow=nc, ncol=nc-ifinal)

    genpart <- colindices[seq_len(ifinal)]
    dpdtpart <- colindices[ifinal+seq_len(nc-ifinal)]
    b.mat.kernel[genpart,] <- -mat[seq_len(ifinal),dpdtpart,drop=FALSE]
    b.mat.kernel[dpdtpart,] <- diag(nc-ifinal)
    # Calcul du modulo p
    b.mat.kernel <- b.mat.kernel %% p
  }
  # otherwise there is no kernel generator
  else{
    b.mat.kernel <- matrix(0, nc, 0)
  }
  rownames(b.mat.kernel) <- pseudonames
  storage.mode(b.mat.kernel) <- "integer"
  return(b.mat.kernel)
} # end kernelmatrix.basep


