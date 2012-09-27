#    options(bigmemory.typecast.warning=FALSE)
#---------------------------------------------------------------------------
#   UTILITIES FOR PLANOR, in particular matrix algebra modulo p,
#   for p a prime
#---------------------------------------------------------------------------
# 1. Miscellaneous
# crossing <- function(n,start=1)
# crossingBig <- function(n,start=1) : idem, return a big.matrix
# symmdiff <- function(x,y)
# cross.designs <- function(designs)
# 2. Matrix algebra modulo p, with p a prime
# convertinto.basep <- function (x, p)
# convertfrom.basep <- function (x, p)
# converttfromBig.basep <- function (x, p)
# representative.basep <- function(mat,p)
# kernelmatrix.basep <- function(mat,p)
# subgroup.basep <- function(mat,p,all=FALSE)
#---------------------------------------------------------------------------
#  1. MISCELLANEOUS
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
crossing <- function(n,start=1){
  # Generates all n1 x n2 x ... x ns combinations of size s with n1,...,ns integers
  # ARGUMENTS:
  #  n: a vector of integers of length s
  #  start: integer from where to start the series of integers
  # OUTPUT:
  #  an integer matrix with prod(n) rows and s columns giving all
  #  combinations along the rows, in lexicographic order
  # EXAMPLE:
  #  crossing(rep(2,3))
  N <- prod(n)
  s <- length(n)
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
}   # fin crossing
#---------------------------------------------------------------------------
crossingBig <- function(n,start=1){
  # Version big.matrix of crossing
  # Generates all n1 x n2 x ... x ns combinations of size s with n1,...,ns integers
  # ARGUMENTS:
  #  n: a vector of positive integers of length s
  #  start: integer from where to start the series of integers
  # OUTPUT:
  #  an integer matrix with prod(n) rows and s columns giving all
  #  combinations along the rows, in lexicographic order;
  # big.matrix
  # EXAMPLE:
  #  a <- crossingBig(rep(2,3))
  # print(a)
  
  N <- prod(n)
  s <- length(n)
  n <- c(n,1)
  b.crosses <- big.matrix( N, s, type="short")
  for(i in seq_len(s))
    {
      motif <- start + seq_len(n[s+1-i])-1
      repet1 <- rep( prod(n[s+1-i+seq_len(i)]), n[s+1-i] )
      if(i==s){ repet2 <- 1 }
      else{ repet2 <- prod(n[seq_len(s-i)]) }
      b.crosses[,s-i+1] <- rep( rep( motif, repet1 ), repet2 )
    }
  return(b.crosses)
} # fin crossingBig
#---------------------------------------------------------------------------
symmdiff <- function(x,y){
  # Calculates the symmetric differences between two collections Cx and Cy of
  # subsets of S, where S is a set of size n
  # cf. the Delta operation in Kobilinsky 2000, p.5
  # ARGUMENTS:
  #  x: a n-row 0-1 matrix, each column of which defines a subset in Cx
  #  y: a n-row 0-1 matrix, each column of which defines a subset in Cy
  # OUTPUT:
  #  a n-row 0-1 matrix with one column per distinct symmetric difference
  #  between a subset in Cx and a subset in Cy
  # a big.matrix
  # EXAMPLE:
  #  M1 <- diag(3)
  #  M2 <- cbind( c(1,1,0), c(1,0,1) )
  #  a <-symmdiff(M1,M2)
  # print(a)

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
  b.z <- tconvertintoBig.basep(z.codes,2)
  rownames(b.z) <- rownames(x)
  return(b.z)
}
#---------------------------------------------------------------------------
cross.designs <- function(designs){
  # generates a design by crossing one, two or more subdesigns
  # ARGUMENTS:
  #  designs: a list of design matrices
  # OUTPUT:
  #  the new design matrix;
  # a big.matrix
  # EXAMPLE:
  #  pl1 <- crossing(c(2,3))
  #  colnames(pl1) <- c("A","B")
  #  pl2 <- crossing(5)
  #  colnames(pl2) <- c("C")
  #  a <- cross.designs(list(pl1,pl2))
  # print(a)


  nrows <- unlist(lapply(designs,nrow))
  ncols <- unlist(lapply(designs,ncol))
  prod.nrows <- prod(nrows)
  sum.ncols <- sum(ncols)
  b.crossdesign <- big.matrix( prod.nrows, sum.ncols, type="short")
  
  
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
  return(b.crossdesign)
}  # fin cross.designs
#---------------------------------------------------------------------------
# 2. MATRIX ALGEBRA MODULO p, with p a prime
#---------------------------------------------------------------------------


#---------------------------------------------------------------------------
go1convertinto.basep <- function (x, p) {
   # Conversion of an integer  into base p
  # Called by convertinto.basep and goconvertinto.basep
 # return c(n1, n2, .. nx) such as x=n1+n2*p+n3*p**2+n4*p**3
  
  if (x != round(x) || x < 0) 
    return(x)
  val <- x%%p
  while ( (x <- x%/%p) > 0 ) {
    newval <- x%%p
    val <- c(val,newval)
  }
    return(val)
} #fin go1convertinto.basep

#---------------------------------------------------------------------------

goconvertinto.basep <- function (x, p) {
  # called by convertinto.basep: loop when its argument is vector

  # Call go1convertinto.basep on  max(x) to know
  # the maximal number of columns of the result

  
  
  
  
  
  
    l <- matrix(0, length(x), length(go1convertinto.basep(max(x),p)))
    for(i in seq_along( x)){
      dec.i <- go1convertinto.basep(x[i],p)
      l[i, seq_along(dec.i) ] <- dec.i
    }
    return(l)
  } # fin goconvertinto.basep


#---------------------------------------------------------------------------
convertinto.basep <- function (x, p) {
  # Conversion of an integer or integer vector x into base p
  # The coefficients are ordered by increasing powers of p
  # ARGUMENTS:
  #  x an integer or an integer variate
  #  p a prime
  # OUTPUT:
  #  if x is an integer, a vector of elements in Zp
  #  if x is a vector, a matrix with length(x) rows 
  #       and the adequate number of columns
  # EXAMPLE:
  #  convertinto.basep(1:10, 3)

  
  if (length(x) > 1) 
    return(goconvertinto.basep(x, p))
  else 
     return(go1convertinto.basep(x, p))
} # fin convertinto.basep

#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
goconvertfrom.basep <- function (x, p) {
  # Function called by convertfrom.basep when x is a vector
  # See convertfrom.basep
  return( sum( x * p^(seq_along(x)-1) ))
}



#---------------------------------------------------------------------------
convertfrom.basep <- function (x, p) {
  # Conversion of integers x coded as vectors of coefficients in base p
  # to classical integers in base 10
  # ARGUMENTS:
  #  x: a vector of coefficients ordered by increasing p powers
  #     or a n-row matrix with each row interpreted as such a vector;
  #     the elements of x must be between 0 and p-1
  #  p: a prime
  # Output:
  #  a scalar integer if x is a vector, or an integer n-vector if x is a matrix
  # EXAMPLE:
  #  vec3 <- convertinto.basep(1:10, 3)
  #  convertfrom.basep( vec3, 3 )
  
  
  
  
  
  
  
  




  
  if (is.matrix(x)) {
    
    
    l <- rep(NA, nrow(x))
    for(i in seq_along( l)){
      l[i] <- goconvertfrom.basep(x[i,],p)
    }
    return(l)
  }
  else
    return(goconvertfrom.basep(x,p))
}

#---------------------------------------------------------------------------
converttfromBig.basep <- function (x, p) {
  # Version of convertfrom where the argument x, when it is a matrix,
  # is the transposed of the one of convertfrom
  # (To avoid the transposition in the call)
  # EXAMPLE
  #  vec3 <- convertinto.basep(1:10, 3)
  #  aref <- convertfrom.basep( vec3, 3 )
  # vv=as.big.matrix(t(vec3), type="short")  # no matter the type
  # a<- converttfromBig.basep(vv,3)
  # all.equal(aref,a)
  
  
  
  




  
  if (is.big.matrix(x)) {
    
    
    l <- rep(NA, ncol(x))
    for(i in seq_along( l)){
      l[i] <- goconvertfrom.basep(x[,i],p)
    }
    return(l)
  }
  else
    return(convertfrom.basep(t(x),p))
}


#---------------------------------------------------------------------------
representative.basep <- function(b.mat,p){
  # generates the minimal set of representatives in base p
  # of the columns x of matrix mat 
  # ARGUMENTS:
  #  mat : a 0-1 matrix ; big.matrix
  #  p : a prime
  # OUTPUT:
  #  a matrix whose columns are the minimal representatives, with 1
  #  in the last non-zero position of x and all possible combinations
  #  of 1,...,p-1 in the other non-zero positions of x; big.matrix
  # EXAMPLE:
  #  a <-representative.basep( as.big.matrix(c(1,0,1,1,0)), 3 )
  # print(a)
  #  a <- representative.basep( as.big.matrix(cbind( c(1,1,0),c(1,1,1) )), 3 )

  if(p==2) return(modBig(b.mat,2))
  #
  b.representative <- NULL

  for(j in seq_len(ncol(b.mat))){
    x <- b.mat[,j]
    select <- seq_along(x)[x != 0]
    nbtocross <- length(select)-1
    if( nbtocross <= 0 ) mat.j <- as.matrix(x)
    else{
      select <- select[seq_len(nbtocross)]
      N <- (p-1)^nbtocross
      mat.j <- matrix(x, nrow(b.mat), N)
#       mat.j[select,] <- t( crossing(rep(p-1,nbtocross),start=1) )
      retcr <- crossingBig(rep(p-1,nbtocross),start=1)
      mat.j[select,] <- applyBig(retcr, 1, function(X) {return(X)})
    }
    b.aux <- big.matrix(nrow(mat.j), ncol(mat.j) , type="short",
                        dimnames=dimnames(mat.j))
    b.aux[,] <- mat.j
    b.representative <- cbindBig(b.representative, b.aux)
  } # fin j
#  return(representative %%p)
  return(modBig(b.representative,p))
} # fin representative.basep
#---------------------------------------------------------------------------
kernelmatrix.basep <- function(mat,p){
  # Calculates a kernel basis of a p-morphism
  # ARGUMENTS:
  #  mat : the n x s matrix associated with the p-morphism, with elements in Zp; 
  #  p : a prime
  # OUTPUT:
  #  an s x q matrix of q independent kernel generators, with elements in Zp; a big matrix
  #  if mat is a regular matrix, then the output is an s x 0 matrix
  # (in this case, it is not a big.matrix)
  # DETAILS:
  #  done by an adaptation of the Gauss-Jordan elimination algorithm
  # EXAMPLE: 
  #  mat[,] <- matrix(
  #  c(1,1,0,0, 2,0,1,0,  0,1,2,4, 0,1,2,3, 0,0,0,1) , nrow=4)
  #  matker <- kernelmatrix.basep(mat,5)
  #  a <- multBigSmod(mat, matker,5) # mat %*% matker %%5
  # print(a) # null matrix
  nr <- nrow(mat)
  nc <- ncol(mat)
  pseudonames <- colnames(mat)
  mat <- mat %%p
  # special cases when mat=Id or mat=(Id|B): direct solutions
  if( nr == nc ){
    if( all( mat[,seq_len(nr)]==diag(nr) ) ){
      b.mat.kernel <- matrix(0,nr,0)
    # In this case, the returned structure is not a big.matrix
    # because:
    # A big.matrix instance must have at least one row and one column
      rownames(b.mat.kernel) <- pseudonames
      return(b.mat.kernel) }}
  if( nr < nc ){
    if( all( mat[,seq_len(nr)]==diag(nr) ) ){
 mat.kernel <- rbind( -mat[,seq(nr+1,nc),drop=FALSE], diag(nc-nr) )%%p








 b.mat.kernel <- as.big.matrix(mat.kernel, type="short")
      rownames(b.mat.kernel) <- pseudonames
      return(b.mat.kernel) }}
  # general case: elimination algorithm
  # - initialisation
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
      else{
        ifinal <- i-1
        continue <- FALSE }
    }
    # if zz is NOT NA: elimination step
    else{
      # row exchange if needed
      if(i != zz)  mat[c(i,zz),] <- mat[c(zz,i),]
      # normalisation
# Raw calculation of the inverses modulo p
      inv <- 0
      inv <- .C("PLANORinv", as.integer(p),
                   as.integer(mat[i,colindices[i]]),
                   inv= as.integer(inv))$inv
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
    b.mat.kernel <- big.matrix(nc,nc-ifinal, type="short", init=NA)
    genpart <- colindices[seq_len(ifinal)]
    dpdtpart <- colindices[ifinal+seq_len(nc-ifinal)]
    b.mat.kernel[genpart,] <- -mat[seq_len(ifinal),dpdtpart,drop=FALSE]
    b.mat.kernel[dpdtpart,] <- diag(nc-ifinal)
    # Calcul du modulo p
    #mat.kernel <- mat.kernel %% p
    b.mat.kernel <- modBig(b.mat.kernel, p)
  }
  # otherwise there is no kernel generator
  else{
    b.mat.kernel <- matrix(0, nc, 0)
    # In this case, the returned structure is not a big.matrix
    # because:
    # A big.matrix instance must have at least one row and one column

  }
  rownames(b.mat.kernel) <- pseudonames

  return(b.mat.kernel)
} # fin kernelmatrix.basep


#---------------------------------------------------------------------------
subgroup.basep <- function(b.mat,p,all=FALSE){
  # FUNCTION
  # calculates the non null elements of the subgroup H generated by the columns
  # of b.mat, considered as vectors in (Zp)^s
  # ARGUMENTS:
  #  b.mat: a matrix of integers modulo p whose columns are assumed to
  #       be independent vectors in (Zp)^s;
  #       a big.matrix
  # Ne peut etre appelé sur une matrice à 0 colonne
  #  p: a prime
  #  all: if TRUE all elements in H are given, if FALSE only elements up
  #       to multiplication by an integer modulo p are given
  # OUTPUT:
  #  a matrix of integers modulo p whose columns are the subgroup elements; a big.matrix
  # DETAILS:
  #  it is not checked whether the column vectors of b.mat are independent, so there
  #  may be several times
  
  
  
  
  
  
  
  # EXAMPLES:
  #  a <- subgroup.basep( as.big.matrix(cbind(c(1,2,0),c(2,0,1)), type="short"), 3)
  # print(a)
  # a <- subgroup.basep( as.big.matrix(cbind(c(1,2,0),c(2,0,1)), type="short"), 3, all=TRUE)
  # print(a)

  nbg <- ncol(b.mat)
  # coeffs = matrix to contain the coefficients of the linear combinations of generators
  N <- p^nbg
  
  if (!is.finite(N))
    stop(paste("Matrix too big in subgroup.basep\n",
               "Need", p, "^", nbg,"\n"))
  # ---------------------------------------------
  # On vérifie que l'on peut indicer la matrice résultat,
  # i.e que l'indice max est inf. à l'entier max:
  if ((nrow(b.mat) * (N-1)) > (2**31-1))
    stop(paste("Matrix too big in subgroup.basep\n",
               "Need", (nrow(b.mat) * (N-1)),
               "(nbg=", nbg, "p=", p, "nrow(b.mat)=", nrow(b.mat), ")",
               "max", (2**31-1),"\n"))
  # ---------------------------------------------
#  b.outmat a au plus N-1 colonnes (N>=2)
    b.outmat <- big.matrix( nrow=nrow(b.mat), ncol=(N-1),
                       type="short", init=0,
                       dimnames=list(rownames(b.mat), NULL))

  # Calcul de b.outmat, colonne par colonne
  # En sortie de PLANORsubgroup, b.outmat est modifiée
    .Call("PLANORsubgroup",
                  as.integer(nrow(b.mat)),
                  as.integer(nbg),
                  as.integer(N),
                  b.mat@address,
                  as.integer(p),
                  as.integer(all),
                  b.outmat@address)
  # Eliminer les -1 tout en gardant les dim
  # S'il n'y a que des -1, c'est que la 1ière valeur
  # non nulle, sur chaque ligne de coeffs, est tjrs !=1
  if (all(b.outmat[1,]<0))
    stop("Internal error: subgroup.basep\n")
  
  b.outmat <- b.outmat[,b.outmat[1,]>=0, drop=FALSE]
  b.outmat <- as.big.matrix(b.outmat,   type="short")
  return(b.outmat)
} # fin subgroup.basep
#---------------------------------------------------------------------------
