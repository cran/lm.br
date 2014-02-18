


.onAttach <- function(...)  {
  library(help=lm.br)$info[[1]] -> version
  version <- version[pmatch("Version",version)]
  um <- strsplit(version," ")[[1]]
  version <- um[nchar(um)>0][2]
  hello <- paste( " lm.br  version ", version, 
    ",  '?lm.br' starts help", sep="" )
  packageStartupMessage( hello )
}





lm.br  <- function( formula, type ="LL", data, subset, na.action,
  weights, inverse =FALSE, var.known =FALSE, offset,
  contrasts =NULL, ... )  {

  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match( c( "formula", "data", "subset", "weights",
    "na.action", "offset" ), names(mf), 0L )
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  offset <- as.vector(model.offset(mf))

  type <- toupper(type)

  x <- model.matrix( mt, mf, contrasts )

  y <- model.response( mf, "numeric")
  if( NCOL(y) > 1 )  
    stop( "mutiple response vectors not supported" )
  ynm <- names(mf)[1]

  w <- model.weights(mf)
  w_ <- w
  if( is.vector(w) )  { 
    if(inverse) { 
      if(any(w==0)) 
        stop( "zero variances not allowed" )  
      else
        w <- 1/w 
    }
  }  else
  if( is.matrix(w) )  {
    n <- NROW(y)
    if( any( dim(w) != c(n, n) ) )  
      stop( "dim('weights') invalid" )
    eW <- eigen(w, TRUE)
    d <- eW$values
    if( any(d <= 0) ) stop("'weights' not positive-definite")
    A <- diag(d^ifelse(inverse, -0.5, 0.5)) %*% t(eW$vector)
    Ainv <- eW$vector %*% diag(d^ifelse(inverse, 0.5, -0.5))
    if(!is.null(offset)) 
      offset <- as.vector( A %*% as.matrix(offset) )
  }

# use 'lm.fit' or 'lm.wfit' to check input, check x-rank per 'tolerance',
# and return some of the output list
  z  <-  if( is.null(w) )
      lm.fit( x, y, offset=offset, ... )
    else  {
      if( is.vector(w) )
        lm.wfit( x, y, w, offset=offset,  ... )
      else
        lm.fit( A %*% x, A %*% y, offset=offset, ... )
    }


  if( length(z$coef) > 0  &&  !is.na(z$coef[1]) )  {

# 'x1' is the term with a coefficient changepoint
    dn <- colnames(x)
    xint <- if( dn[1] == "(Intercept)" )  TRUE  else  FALSE
    if(xint) x1c=2 else x1c=1
    x1 <- as.vector( x[,x1c] )
    x1nm <- dn[x1c]


# 'xb' is the 'x' input matrix, but with two vectors for 'x1'
#  before and after changepoint,  and zeroes for columns
#  that were linearly dependent
    nx <- ncol(x)
    nxb <- ncol(x)+1
    xb <- matrix( 0, nrow(x), nxb )
    xb[,1] <- x[,1]
    for(i in x1c:nx) 
      xb[,i+1] <- if( is.na(z$coef[i]) )  0  else  x[,i]

    bnm <- paste( " ", x1nm, "< theta" )
    bpnm <- paste( " ", x1nm, "> theta" )
    rownames(xb) <- rownames(x)
    colnames(xb) <-  if( xint )
      { if( x1c < nx )  
          c("alpha", bnm, bpnm, dn[(x1c+1):nx])  
        else  
          c("alpha", bnm, bpnm) 
      }
      else
      { if( x1c < nx )  
          c(bnm, bpnm, dn[(x1c+1):nx])  
        else
          c(bnm, bpnm) 
      }


#  construct C++ object
#  re-order for 'x1' non-decreasing,  drop rows with  w = 0,  
#  drop columns that are not independent
    x_ <- x[ order(x1), , drop = FALSE]
    y_ <- y
    if(!is.null(offset)) y_ <- y_ - as.vector(model.offset(mf))
    y_ <- y_[ order(x1) ]
    if( !is.null(w_) )  
      w_ <- if(is.matrix(w_))  w_[ order(x1), order(x1) ]  
            else  w_[ order(x1) ]

    if( is.vector(w_) )  if( any(w_==0) )  {
      ok <- w_!=0
      w_ <- w_[ok]
      y_ <- y_[ok]
      x_ <- x_[ok, , drop = FALSE]
    }
    if( is.vector(w_) )  w_ <- as.matrix(w_)
    if( is.null(w_) )  w_ <- as.matrix( -1 )

    x_fullcol <- x_
    x_ <- x_[ , !is.na(z$coef), drop = FALSE ]

    if( type=="LL" ) {
      if( xint )  model_num <- 1  else
        stop( "'alpha'=0 not supported for type \"LL\"" )
    } else {
      if( type=="TL" ) {
        if( xint )  model_num <- 2  else  model_num <- 3
      }  else  {
        if( type=="LT" ) {
          if( xint )  model_num <- -2  else  model_num <- -3 
        }  else  
          stop( "'type' must be \"LL\", \"LT\" or \"TL\"" )
      }
    }


#  loop to drop columns that are dependent at some changepoint value
#  if x-matrix is dependent at  changepoint = 'th'  then  
#  Q*f(th)=0  on lower rows,  where  f(th) = max( x1-th, 0 )
    x_dep <- TRUE

    while( x_dep )  {

      obj <- new( Cpp_Clmbr, y_, x_, w_, model_num, 
            as.integer(inverse), as.integer(var.known) )

      thQfmin <- obj$param()[5]
      xb[ ,x1c] <- if( type=='TL' )  0  else
          { if( is.infinite(thQfmin) )  -1  else  pmin(x1 - thQfmin, 0 ) }
      xb[ ,x1c+1] <- if( type=='LT' )  0  else
          { if ( is.infinite(thQfmin) )  1  else  pmax(x1 - thQfmin, 0 ) }

      z  <-  if( is.null(w) )
          lm.fit( xb, y, offset=offset, ... )
        else  {
          if( is.vector(w) )
            lm.wfit( xb, y, w, offset=offset,  ... )
          else
            lm.fit( A %*% xb, A %*% y, offset=offset, ... )
        }

      xbr <- z$rank
      nb <- if(type=='LL')  ncol(x_) + 1  else  ncol(x_)
      if( xbr < nb )  {
        xb[ , is.na(z$coef) ]  <- 0
        ok <- vector( 'logical', nx )
        ok[] <- TRUE
        for(i in (x1c+1):nx)  if( is.na(z$coef[i+1]) )  ok[i]<- FALSE
        x_ <- x_fullcol[ , ok, drop = FALSE ]
      } else
        x_dep <- FALSE

    }


# output
    par <- obj$param()
    if( is.nan( par[1] ) )  {   # case of perfect line
      xb <- x
      if(xint) colnames(xb)[1] <- "  alpha"
      for(i in 1:nx) if( is.na(z$coef[i+1]) ) xb[,i] <- 0
    }  
    else  {
      xb[ , x1c ]  <-  if( type=='TL' )  0  else
          { if( is.infinite(par[1]) )  1  else  pmin( x1-par[1], 0 ) }
      xb[ , x1c+1 ]  <-  if( type=='LT' )  0  else
          { if ( is.infinite(par[1]) )  1  else  pmax( x1-par[1], 0 ) }
    }

    if( is.infinite(par[1]) )  {  # can occur for models with alpha=0
      if( type=='LT' )
        { colnames(xb)[1] <- "  (Intercept)";  colnames(xb)[2] <- bnm }
      if( type=='TL' )
        { colnames(xb)[1] <- bpnm;  colnames(xb)[2] <- "  (Intercept)" }
    }


    z  <-  if( is.null(w) )
        lm.fit( xb, y, offset=offset,  ... )
      else  {
        if( is.vector(w) )
          lm.wfit( xb, y, w, offset=offset,  ... )
        else
          lm.fit( A %*% xb, A %*% y, offset=offset, ... )
      }


    if( type=='TL' )  z$coefficients[x1c] <- 0
    if( type=='LT' )  z$coefficients[x1c+1] <- 0

#  add 'theta' to the coefficients
    for(i in ncol(xb):1 )  {
      z$coef[i+1] <- z$coef[i]
      names(z$coef)[i+1] <- names(z$coef)[i]
    }
    z$coef[1] <- par[1]
    names(z$coef)[1] <- "theta"


#  accessor functions
    z$CppObj  <-  obj


    z$ci <- function( CL =0.95, method ="clr" )  {
      method <- toupper(method)
      met <- integer(1)
      if( method=="CLR" )  met <- 1  else  {
        if( method=="AF" )  met <- 2  else
          stop( "'method' must be \"CLR\" or \"AF\"" )
      }
      (z$CppObj)$ci( CL, met )
    }


    z$cr <- function( CL =0.95, method ="clr", incr =NULL, output ="G" )  {
      if( is.null(incr) )  incr<- -1  else 
        if( incr <= 0 )  stop("'incr' must be positive")
      method <- toupper(method)
      met <- integer(1)
      if( method=="CLR" )  met <- 1  else  {
        if( method=="AF" )  met <- 2  else
          stop( "'method' must be \"CLR\" or \"AF\"" )
      }
      if(missing(output) && .Device=="null device")  output <- "T"
      output <- toupper(output)
      if( output=="T" )
        (z$CppObj)$cr3( CL, met, incr )
      else  {
        bounds <- (z$CppObj)$cr4( CL, met, incr, as.integer(FALSE) )
        if( output=="V" )
          return( bounds )
        else  {
          if( output=="G" )  {
            nbd <- nrow(bounds)
            cl <- as( round(100*CL,0), "character" )
            title <- paste( cl, 
              "% conf. region for changepoint by ", method, sep="")
            if( (z$xint && length(z$coef)==4) ||
                  (!z$xint && length(z$coef)==3) )  {
#  univariate model
              n <- length(z$x1)
              x <- y <- matrix( NA, max(n,nbd), 3 )
              x[1:n,1] <- z$x1
              y[1:n,1] <- z$y
              x[1:nbd,2:3] <- bounds[,1]
              y[1:nbd,2:3] <- bounds[,2:3]
              matplot( x, y,
                type=c('p','l','l'), pch=4, lty='solid', col='black',
                main=title, xlab=z$x1nm, ylab=z$ynm )
            } else {
#  multivariate
              x <- y <- matrix( NA, nbd, 2 )
              x[,1:2] <- bounds[,1]
              y[,1:2] <- bounds[,2:3]
              xnm <- paste( "theta (", z$x1nm, ")", sep="")
              matplot( x, y,
                type=c('l','l'), lty='solid', col='black',
                main=title, xlab=xnm, ylab="alpha" )
            }
          } else
            stop("'output' must be \"G\", \"T\" or \"V\"")
        }
      }
    }


    z$sl <- function( theta0, alpha0 =NULL,  method ="clr",
              accuracy =0.001, output ="T" )  {
# overload by if-else statements
      if( !is.null(alpha0) && !is.numeric(alpha0) )  {
        if( is.numeric(method) ) {
          if( !missing(accuracy) )  {
            if( is.character(accuracy) )  output <- accuracy
              else  stop( "'output' must be \"T\", \"V\" or \"B\"" )
          }
          accuracy <- method
        }
        method <- alpha0
        alpha0 <- NULL
      }
      method <- toupper(method)
      output <- toupper(output)
      met <- integer(1)
      if( method=="CLR" )  met <- 1  else  {
        if( method=="AF" )  met <- 2  else  {
          if( method=="MC" )  met <- 3  else
            stop( "'method' must be \"CLR\", \"AF\" or \"MC\"" )
        }
      }
      value <- verbose <- logical(1)
      if( output=="T" )  {
        value <- FALSE
        verbose <- TRUE  
      } else {
        if( output=="V" )  {
          value <- TRUE
          verbose <- FALSE
        } else {
          if( output=="B" )
            value <- verbose <- TRUE
          else  stop( "'output' must be \"T\", \"V\" or \"B\"" )
        }
      }
      if(value) {
        result <- double(1)
        result <- if( is.null(alpha0) )
            (z$CppObj)$sl5( met, as.integer(verbose),
              as.integer(value), accuracy, theta0 )
          else
            (z$CppObj)$sl6( met, as.integer(verbose),
              as.integer(value), accuracy, theta0, alpha0 )
        return( result )
      }  else  {
        if( is.null(alpha0) | !is.numeric(alpha0) )
          (z$CppObj)$sl3( met, accuracy, theta0 )
        else
          (z$CppObj)$sl4( met, accuracy, theta0, alpha0 )
      }
    }


    z$mle <- function( )  (z$CppObj)$mle( )


    z$sety <- function( rWy )  {
      if( NCOL(rWy) > 1 | !is.numeric(rWy) )  
        stop( "'rWy' must be a numeric vector" )
      rWy <- drop(rWy)
      if( !is.null(z$offset) )  rWy <- rWy - z$rWoffset
      rWysorted <- rWy[ order(z$x1) ]
      if( is.vector(z$weights) )  if( any(z$weights==0) )  {
        wsorted <- z$weights[ order(z$x1) ]
        rWysorted <- rWysorted[ wsorted!=0 ]
      }
      (z$CppObj)$sety( rWysorted )
    }


#  output list
    p <-  if( is.nan( par[1] ) )  z$rank  else  z$rank + 1
    z$no_of_parameters <- p
    z$df.residual <- nrow(x_) - p
    z$xint <- xint
    z$x1 <- x1
    z$x1nm <- x1nm
  }

  if( is.vector(w) )  {
    if(!is.null(offset)) z$rWoffset <- sqrt(w) * offset
  }  else
  if( is.matrix(w) )  {
     z$fitted.values <- drop(Ainv %*% z$fitted.values)
     z$residuals <- drop(Ainv %*% z$residuals)
    if(!is.null(offset)) z$rWoffset <- offset
  }

  class(z) <- "lm.br"
  z$call <- call
  z$na.action <- attr(mf, "na.action")
  z$contrasts <- attr(x, "contrasts")
  z$terms <- mt
  z$xlevels <- .getXlevels(mt, mf)
  z$offset <- as.vector(model.offset(mf))
  z$y <- y
  z$ynm <- ynm
  z$weights <- model.weights(mf)
  z$type <- type
  z
}




print.lm.br  <-  function ( x, digits = max(3L, getOption("digits") - 3L), ... )  {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse ="\n"),
    "\n\n", sep = "")
  type <- x$type
  cat( "Broken-line type:  ", type, "\n\n" )
  if (length(coef(x))>1  && !is.na(x$coef[2]))  {
    cat( "Significance Level of H0:\"no changepoint\" vs",
      "H1:\"one changepoint\"\n" )
    cat("  ")
    mx1 <- max( x$x1 )
    mn1 <- min( x$x1 )
    ai <- (mx1-mn1)/( length( x$x1 ) - 1 )
    if( type=='LT' && x$xint )
      x$sl( round( min( mx1 + ai*1.5 ), 2) )
    else
      if( type=='TL' && !x$xint )
        x$sl( -Inf )
      else
        if( type=='LT' && !x$xint )
          x$sl( Inf )
        else
          x$sl( round( max( mn1 - ai*1.5 ), 2) )
    cat("\n")
    x$ci()

# print coefficients unless 'sety' has been called
    par <- x$CppObj$param()
    if( !par[6] )  {
      cat( "Changepoint and coefficients:\n" )
      print.default( round(x$coef, 5) )
    }
    else  cat( "Use 'mle()' for parameter estimates after a call to 'sety'.\n" )
  }
  else  cat( "No coefficients\n" )
  cat("\n")
  invisible(x)
}


