
# +++++++++++++++++++++++++++++++++++++++++++++++++++
# FUNCTIONS to manipulate big.matrices
# modBig(big, p) = big %%p
# cbindBig(big1, big2) = cbind(big1, big2)
# rbindBig(big1, big2) = rbind(big1, big2)
# multBigmod(big1, big2, p, signe=+1) = (signe*big1  %*% big2) %%p
# multBig(big1, big2) = (big1  %*% big2)
# applyBig(big, dimension, fun)= apply(big, dimension, fun)
# col0Big(big) = 1 if one column at least is entirely equal to zero
# exchangeColBig(big, indT, indS, signe=+1) = res[,indT] <- signe*big[,indS] and return res
# crossingBig(n,start=1) = crossing(n,start=1) return a big.matrix
# converttfromBig.basep= converttfrom.basep
# NOTE: the dimensions of the matrices in argument
#   are not verified (these functions are internal functions)
#   De-comment the sentences with occurence of VERIF to change this.
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# VERIF <- FALSE # To verify the arguments
modBig <- function(big, p)
  {
    # Return a big.matrix= big%%p
    res <- big.matrix(nrow(big),ncol(big), type="short")
    .Call("PLANORmodBig",
          big@address,
          as.integer(p),
          res@address)
    return(res)
  } # end modBig
# ----------------------------------------------------
cbindBig <- function(big1, big2, deparse.level=1)
  {
    # Return a big.matrix=cbind(big1, big2)
    # big1 may be NULL
    if (is.null(big1)) {
      ncol1 <- 0
      address1 <- big2@address # any address: will be ignored
    }    else {
##     if (VERIF) {
##        if (nrow(big1) != nrow(big2))
##          stop("Error in cbindBig: the matrices should have the same number of rows")
##      } # end VERIF
      ncol1 <- ncol(big1)
      address1 <- big1@address
    }
     ret <- big.matrix(nrow(big2), (ncol1+ncol(big2)), type="short",
                      init=111) # init with anything
    .Call("PLANORcbindBig", address1, big2@address,
          ret@address,
          as.integer(nrow(big2)), as.integer(ncol1),
          as.integer(ncol(big2)))
    if (deparse.level==1) {
      if (!is.null(colnames(big1)) && !is.null(colnames(big2)))
        colnames(ret) <-  c(dimnames(big1)[[2]], dimnames(big2)[[2]])
      if (!is.null(rownames(big1)))
          rownames(ret) <- rownames(big1)
      else
        rownames(ret) <- rownames(big2)
    } # end deparse
     return(ret)
  } # end cbindBig

#---------------------------------------------------------------------------
rbindBig <- function(big1, big2, deparse.level=1)
  {
    # Return a big.matrix=rbind(big1, big2)
##      if (VERIF) {
##      if (ncol(big1) != ncol(big2))
##        stop("Error in rbindBig: the matrices should have the same number of cols")
##      } # end VERIF

      ncolb <- ncol(big1)
      nrow1 <- nrow(big1)
      nrow2 <- nrow(big2)
     ret <- big.matrix(nrow1+nrow2, ncolb, type="short",
                      init=111)
    .Call("PLANORrbindBig", big1@address, big2@address,
          ret@address,
          as.integer(nrow1), as.integer(nrow2),
          as.integer(ncolb))
    if (deparse.level==1) {
      if (!is.null(rownames(big1)) && !is.null(rownames(big2)))
        rownames(ret) <- c(dimnames(big1)[[1]],dimnames(big2)[[1]])
      if (!is.null(colnames(big1)))
        colnames(ret) <- colnames(big1)
      else
        colnames(ret) <- colnames(big2)
    } # end deparse
     return(ret)
  } # end rbindBig


#---------------------------------------------------------------------------
rbind1Big <- function(big1, mat2, deparse.level=1)
  {
    # Return a big.matrix=rbind(big1, mat2)
##      if (VERIF) {
##     if (ncol(big1) != ncol(mat2))
##        stop("Error in rbindBig: the matrices should have the same number of cols")
##      } # end VERIF

      ncolb <- ncol(big1)
      nrow1 <- nrow(big1)
      nrow2 <- nrow(mat2)
     ret <- big.matrix(nrow1+nrow2, ncolb, type="short",
                      init=111)
    .Call("PLANORrbind1Big", big1@address, mat2,
          ret@address,
          as.integer(nrow1), as.integer(nrow2),
          as.integer(ncolb))
    if (deparse.level==1) {
      if (!is.null(rownames(big1)) && !is.null(rownames(mat2)))
        rownames(ret) <- c(dimnames(big1)[[1]],dimnames(mat2)[[1]])
      if (!is.null(colnames(big1)))
        colnames(ret) <- colnames(big1)
      else
        colnames(ret) <- colnames(mat2)
    } # end deparse
     return(ret)
  } # end rbind1Big

#---------------------------------------------------------------------------

gotconvertintoBig.basep <- function (x, p) {
  # Conversion of an integer vector x into base p
  # The coefficients are ordered by increasing powers of p
  # Called by tconvertintoBig.basep
  # ARGUMENTS
  #  - x: an integer variate
  #  - p: a prime
  # RETURN
  #   a big.matrix with length(x) columns
  #       and the adequate number of rows
    l <- big.matrix(length(go1convertinto.basep(max(x),p)),
                    length(x),
                    type="short", init=0)
    for(i in seq_along( x)){
      dec.i <- go1convertinto.basep(x[i],p)
      l[ seq_along(dec.i), i ] <- dec.i
    }
    return(l)
  } # end gotconvertintoBig.basep


#---------------------------------------------------------------------------
tconvertintoBig.basep <- function (x, p) {
  # Conversion of an integer or integer vector x into base p
  # The coefficients are ordered by increasing powers of p
  # ARGUMENTS
  #  - x: an integer or an integer variate
  #  - p: a prime
  # RETURN
  #  if x is an integer, a big.matrix with 1 column of elements in Zp
  #  if x is a vector, a big.matrix matrix with length(x) columns
  #       and the adequate number of rows
  # NOTE
  #  return
  # t(convertinto.basep(x,p))  into a big.matrix
  # ----------------------------------------------

  if (length(x) > 1)
    return(gotconvertintoBig.basep(x, p))
  else {
    val <- go1convertinto.basep(x, p)
    b.ret <- big.matrix(length(val), 1, type="short")
    b.ret[,1] <- val
    return(b.ret)
  }

} # end tconvertintoBig.basep

#---------------------------------------------------------------------------
multBigmod <- function(big1, big2, p, signe=+1)
  {
    res <- big.matrix(nrow(big1), ncol(big2), type="short")
    # big2 may be a no big matrix
    if (is.big.matrix(big2)) {
    .Call("PLANORmultBigmod",
       big1@address,    big2@address,
          as.integer(nrow(big1)),
          as.integer(ncol(big1)),
          as.integer(ncol(big2)),
          as.integer(p),
          as.integer(signe),
          res@address)
  }
    else {
          .Call("PLANORmultBigSmod",
       big1@address,    as.double(big2),
          as.integer(nrow(big1)),
          as.integer(ncol(big1)),
          as.integer(ncol(big2)),
          as.integer(p),
          as.integer(signe),
          res@address)
  }
    dimnames(res) <- list(rownames(big1), colnames(big2))
    return(res)
  } # end multBigmod
# --------------------------------------------------------
#---------------------------------------------------------------------------
multBig <- function(big1, big2)
  {
    res <- big.matrix(nrow(big1), ncol(big2), type="short")

    .Call("PLANORmultBig",
       big1@address,    big2@address,
          as.integer(nrow(big1)),
          as.integer(ncol(big1)),
          as.integer(ncol(big2)),
          res@address)
    dimnames(res) <- list(rownames(big1), colnames(big2))
    return(res)
  } # end multBig

# --------------------------------------------------------
applyBig <- function(big, dimension, fun)
  {
    
    
    z <- matrix(big[,], ncol=ncol(big))
    return(apply(z, dimension, fun))
  } # end applyBig
# --------------------------------------------------------
col0Big <- function(big)
  {
# Return 1 if one column at least is entirely equal to zero
    res <-0
    res <-.Call("PLANORcol0Big",
          as.integer(nrow(big)),
          as.integer(ncol(big)),
          big@address,
          res=as.integer(res))
    return(res)
  } #fin col0Big

# --------------------------------------------------------
exchangeColBig <- function(big, indT, indS, signe=+1)
  {
    # res[,indT] <- signe*big[,indS] and return res, a big.matrix
    # This function is motivated by the fact that
    #  big[,indices] does not return a big.matrix
    # When indT is NULL, the first columns of the result matrix
    # are filled in
    if (is.null(indT)) {
      if (is.logical(indS)) {
        z <- indS[indS==TRUE]
        indT <- seq_along(z)
      }
      else
        indT <- seq_along(indS)
    } # end null

    res <- big.matrix(nrow(big), length(indT), type="short")
    res[,indT] <- signe*big[,indS]
    if (!is.null(colnames(big)))
      colnames(res)[indT] <- colnames(big)[indS]
    if (!is.null(rownames(big)))
      rownames(res) <- rownames(big)

    return(res)
  } # end exchangeColBig
# ----------------------------------------------------------
crossingBig <- function(n,start=1){
  # Version big.matrix of crossing
  # Generates all n1 x n2 x ... x ns combinations of size s with n1,...,ns integers
  # ARGUMENTS
  #  - n: a vector of positive integers of length s
  #  - start: integer from where to start the series of integers
  # RETURN
  #  an integer matrix with prod(n) rows and s columns giving all
  #  combinations along the rows, in lexicographic order;
  # big.matrix
  # EXAMPLE
  #  a <- crossingBig(rep(2,3))
  # print(a)
  # -----------------------------------------------

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
} # end crossingBig
#---------------------------------------------------------------------------
converttfromBig.basep <- function (x, p) {
  # Version of convertfrom where the argument x, when it is a matrix,
  # is the transposed of the one of convertfrom
  # (To avoid the transposition in the call)
  # EXAMPLE
  #  vec3 <- convertinto.basep(1:10, 3)
  #  aref <- convertfrom.basep( vec3, 3 )
  #  vv=as.big.matrix(t(vec3), type="short")  # no matter the type
  #  a<- converttfromBig.basep(vv,3)
  #  all.equal(aref,a)
  
  
  
  





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
