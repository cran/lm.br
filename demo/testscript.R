#
#  This R-script tests inferences by 'lm.br' to confirm its exact accuracy.  
#
#  'testrun' generates random models of all types covered by 'lm.br'.  Then
#  'simtest' generates random observations according to each model.  The
#  inferences that 'lm.br' makes from these observations are compared with 
#  the "unknown", underlying, true model.  
#
#  A complete run takes 10-20 hours, depending on your computer system.
#


library( lm.br )


testrun  <-  function( )  {
# loop to generate random models of different types
# and then call function 'simtest'

  cat("\n\nStarting tests of the R package 'lm.br':\n")
  cat("This script generates random models of all supported types.\n")
  cat("For each random model, it checks coverage frequencies and tests\n")
  cat("functions 'sl', 'ci', 'cr'.  A complete run takes 10-20 hours.\n\n\n")

# if 'track' = TRUE then keep an ongoing log of input values
# in local disk files
  track <- FALSE

  if( track )  {
    scname <- tempfile( pat="simcall", fileext=".txt" )
    rWyname <- tempfile( pat="rWy", fileext=".txt" )
    rWy2name <- tempfile( pat="rWy2", fileext=".txt" )
  }

  set.seed(1234)

## test 2 random models for each type of model
  testspertype <- 2
  totaltests <- 5*2*3*2*testspertype
  tstart <- Sys.time()
  ntests <- 0

  for(nmodel in c(-3,-2,1,2,3) )  { for(vk in 0:1)
      { for(wtype in 1:3)  { for(mvx in 0:1)  
        { for(ntest in 1:testspertype) {

    ntests <- ntests + 1
    tcurrent <- Sys.time()
    dtime <- round( 
      as.numeric( difftime(tcurrent,tstart,units="mins")), 0)
    cat("test", ntests, "of", totaltests, ",  elapsed ", 
      dtime, "min\n")

## generate model

    if(nmodel==1)  {
      nmin<-4
      ntmin<-4
      if(vk==0) ntmin<-5
    }
    if(nmodel== 2 || nmodel== -2 )  {
      nmin<-3
      ntmin<-3
      if(vk==0) ntmin<-4
    }
    if( nmodel==3 || nmodel== -3 )  {
      nmin<-2
      ntmin<-2
      if(vk==0) ntmin<-3
    }

    n <- ntmin + round(10*runif(1),0)

    if(mvx) n <- n+3

    x <- 1:n
    ns<-1
    while(ns<nmin) {
      ns<-1
      for (i in 2:n) {
        rinc<- (runif(1)>0.15)
        x[i] <- x[i-1] + ( 2*runif(1) + 0.05 )*rinc;
        x[i] <- round(x[i],4)
        ns <- ns+rinc
      }
    }

    theta <- x[1]-1 + (x[n]-x[1]+2)*runif(1);
    theta <- round(theta,4)

    x <- x[ sample.int(n) ]

    x2coef <- x3coef <- 0
    if(mvx) {
      x2 <- x3 <- 1:n

      ns <- 1
      while( ns < nmin ) {
        ns<-1
        for (i in 2:n) {
          rinc <- (runif(1)>0.15)
          x2[i] <- x2[i-1] + ( 2*runif(1) + 0.05 )*rinc;
          x2[i] <- round(x2[i],4)
          ns <- ns+rinc
        }
      }

      ns <- 1
      while( ns < nmin ) {
        ns<-1
        for (i in 2:n) {
          rinc <- (runif(1)>0.15)
          x3[i] <- x3[i-1] + ( 2*runif(1) + 0.05 )*rinc;
          x3[i] <- round(x3[i],4)
          ns <- ns+rinc
        }
      }
 
      x2 <- x2[ sample.int(n) ]
      x3 <- x3[ sample.int(n) ]
      xm <- matrix(0,n,3)
      xm[,1] <- x
      xm[,2] <- x2
      xm[,3] <- x3
      x <- xm
      x2coef <- round( (runif(1)-0.5)*5, 4 )
      x3coef <- round( (runif(1)-0.5)*5, 4 )
    }

    alpha <- 0.;

    if (nmodel==1) { beta <- -1. }
    if (nmodel==2 || nmodel==3) { beta<-0 }

    betap <- 2. - runif(1)*2.5;
    betap <- round(betap,4)

    if( nmodel== -2 || nmodel== -3 ) { beta<-betap; betap<-0 }

    if (vk==0) {
      var <- 0.1 + 2*runif(1)
      var <- round(var,4)
    } else var <- 1.;

    wgts <- matrix(0,n,n)
    diag( wgts ) <- 1.

    if (wtype==2)  {
      wgts <- as.vector( diag(wgts) )          
      for(i in 1:n) 
        wgts[i] <- round( 1./( 0.04 + 2*runif(1) ), 4 )
      if(n > ntmin) if(runif(1)<0.2) for(i in 1:(n-ntmin)) {
        ir <- floor((n+1)*runif(1))
        if(runif(1)<0.15) wgts[ir] <- 0
      } 
    }

    if (wtype==3) {  
      notpdf<- TRUE
      while( notpdf ) {  
        wgts[,] <- 0
        for(i in 1:n) 
          wgts[i,i] <- round( 1./( 0.04 + 2*runif(1) ), 4 )
        for(j in 1:8) {  ##3
          ir <- floor((n+1)*runif(1))
          jr <- floor((n+1)*runif(1))
          if(ir==jr) {if(ir==n) {ir<-n-1} else {ir<-ir+1}}
          wgti <- wgts[ir,ir]
          wgtj <- wgts[jr,jr]
          wgts[ir,jr] <- (runif(1) - 0.5)*1.9*sqrt(wgti*wgtj)
          wgts[ir,jr] <- round( wgts[ir,jr] , 4 )
          wgts[jr,ir] <- wgts[ir,jr]
        }  ##3
        eW <- eigen(wgts)
        d <- eW$values
        if( any(d < 1.e-5) ) notpdf<- TRUE  else  notpdf<- FALSE
      }  
    }  

    if(nmodel==1)  model <- 'LL'
    if(nmodel==2 || nmodel==3)  model <- 'TL'
    if(nmodel==-2 || nmodel==-3)  model <- 'LT'
    if(abs(nmodel)==3)  xint <- FALSE  else  xint <- TRUE

#  echo call to 'simtest'
    cat("\n  model ", model, ",")
    if(!xint)  cat("  'alpha' known =0,")
    if(vk)  cat("  'var' known =1,")
    cat("\n  ")
    if(mvx) cat("multivariate, ")
    if(wtype==1) cat("no weights")
    if(wtype==2) cat("vector weights")
    if(wtype==3) cat("matrix weights")
    cat("\n  test", ntest, "of", testspertype, "\n\n")
    cat("  random model:\n")
    cat("  theta=", theta, " alpha=", alpha,
      " B=", beta, " Bp=", betap, " var=", var, "\n")
    if(mvx) {
      cat("  x2coef=", x2coef, " x3coef=", x3coef, "\n")
      cat("  'x' matrix\n")
      print(x)
    } else
      cat("  observations at  x=", x, "\n\n")
    if(wtype==2) cat( "  'weights' = ", wgts, "\n" )
    if(wtype==3) {
      cat( "  'weights' = \n" )
      print(wgts, zero.print = ".")
    }
    cat("\n")

    if( track )  {
#  echo to a file
      simcall <- file( scname, "wt")
      cat("\n  'lm.br' tests for model ", model, file= simcall)
      if(!xint)  cat(",  'alpha' known =0", file= simcall)
      if(vk)  cat(",  'var' known =1", file= simcall)
      cat("\n\n", file= simcall)
      cat("  random model:\n", file= simcall)
      cat("  theta=", theta, " alpha=", alpha,
        " B=", beta, " Bp=", betap, " var=", var, "\n", file= simcall)
      if(mvx) {
        cat("  x2coef=", x2coef, " x3coef=", x3coef, "\n", file= simcall)
        cat("  'x' matrix\n", file= simcall)
        for(i in 1:NROW(x)) cat(x[i,], "\n", file= simcall)
      } else
        cat("  observations at  x=", x, "\n\n", file= simcall)
      if(wtype==2) cat( "  'weights' = ", wgts, "\n", file= simcall )
      if(wtype==3) {
        cat( "  'weights' = \n", file= simcall )
        for(i in 1:NROW(x)) cat(wgts[i,], "\n", file= simcall)
      }
      close( simcall )
    }

    simtest( as.matrix(x), wgts, model, xint, x2coef, x3coef,
             vk, theta, alpha, beta, betap, var, 10000,
               track, rWyname, rWy2name  )


  } } } } }

  if( track )  {
    unlink( scname )
    unlink( rWyname )
    unlink( rWy2name )
  }

  cat( "\ntests of 'lm.br' completed successfully\n\n" )
}






simtest <- function( x, W, model, xint, x2coef, x3coef, vk, 
  theta, alpha, B, Bp, var, N =10000, track, rWyname, rWy2name  )
#
#  This function generates sets of random observations according 
#  to the input model.  The inferences that 'lm.br' makes from these 
#  observations are compared with the true model. 
#  Tests coverage frequency of 95% confidence intervals and regions. 
#  Also, on one set of random data, compares 'sl' by CLR-MC and CLR 
#  and tests 'ci' and 'cr'.
#
{
  if( NCOL(x) > 1 )  mvx <- TRUE  else  mvx <- FALSE

# setup generating model from input
  x1 <- x[,1]

  m <- alpha + B*pmin(x1-theta, 0) + Bp*pmax(x1-theta, 0)
  if(mvx)  m <- m + x2coef*x[,2] + x3coef*x[,3]

  if( is.vector(W) )
    m <- m * sqrt(W)
  else {
    eW <- eigen(W, TRUE)
    D <- eW$values
    Q <- eW$vectors
    rD <- diag( sqrt(D) )
    rW <- Q %*% rD %*% t(Q)
    m <- rW %*% m
  }

  n <- NROW(x)
  sigma <- sqrt(var)


# construct an 'lm.br' object with arbitrary 'y' values
  y <- x1
  mod  <-  if( xint ) 
      lm.br( y ~ x, model, w = W, var.k = vk )
    else
      lm.br( y ~ x+0, model, w = W, var.k = vk )

  cat("\n Monte Carlo simulation test of conf. interval",
       "for 'theta':\n" )
  cat("     no. of     coverage frequency of the", 
    "0.95-confidence interval by\n" )
  cat("   iterations                CLR                    AF",
    "\n" )
  flush.console()
  cicov <- 0.95
  passed <- FALSE
  ntrials <- 0
  while( !passed && ntrials < 3 ) {
    countCLR <- countAF <- 0
    for( i in 1:N ) {

      rWy <- m + rnorm(n,0,sigma)
      rWy <- round(rWy,4)
      if(track) {
        yset <- file( rWyname, "wt")
        cat(rWy,file=yset)
        close(yset)
      }
      mod$sety(rWy)
      stest <- mod$sl(theta,"clr",.0001,FALSE)
      if(stest>0.05) countCLR <- countCLR + 1
      stest <- mod$sl(theta,"af",.0001,FALSE)
      if(stest>0.05) countAF <- countAF + 1

      if(i/1000 - floor(i/1000) == 0) {
        cat( format(i,width=10),
          format(countCLR/i, digits=4, nsmall=4, width=22),
          format( countAF/i, digits=4, nsmall=4, width=22),
          "\n" )
        flush.console()
      }
    }

    cat("\n")
    cicov <- countCLR/N
    passed <- if( cicov < .945 )  FALSE  else  TRUE
    ntrials <- ntrials + 1
  }
  if( !passed && ntrials==3 )
    stop("fails interval coverage test")
  if( cicov > .97) 
    warning( 
      gettextf("conf. interval coverage frequency  %d%%", 
      as.integer(100*cicov)) 
    )

  cat("\n Monte Carlo simulation test of conf. region",
        "for '(theta,alpha)':\n")
  cat("     no. of     coverage frequency of the", 
    "0.95-confidence region by\n")
  cat("   iterations                CLR                    AF",
    "\n")
  if(xint)  passed <- FALSE  else  {
    cat(" not applicable for this model\n\n")
    passed <- TRUE
  }
  flush.console()
  crcov <- .95
  ntrials <- 0
  while( !passed && ntrials < 3 ) {
    countCLR <- countAF <- 0
    for( i in 1:N ) {

      rWy <- m + rnorm(n,0,sigma)
      rWy <- round(rWy,4)
      if(track) {
        yset2 <- file( rWy2name, "wt")
        cat(rWy,file=yset2)
        close(yset2)
      }
      mod$sety(rWy)
      stest <- mod$sl(theta,alpha,"clr",.001,FALSE)
      if(stest>0.05) countCLR <- countCLR + 1
      stest <- mod$sl(theta,alpha,"af",.001,FALSE)
      if(stest>0.05) countAF <- countAF + 1

      if(i/1000 - floor(i/1000) == 0) {
        cat( format(i,width=10),
          format(countCLR/i, digits=4, nsmall=4, width=22),
          format( countAF/i, digits=4, nsmall=4, width=22),
          "\n" )
        flush.console()
      }
    }

    cat("\n")
    crcov <- countCLR/N
    passed <- if( crcov < .945 )  FALSE  else  TRUE
    ntrials <- ntrials + 1
  }
  if( !passed && ntrials==3 ) 
    stop("fails region coverage test")
  if( crcov > .97) 
    warning( 
      gettextf("conf. region coverage frequency  %d%%", 
      as.integer(100*crcov)) 
    )



# use final random 'rWy' data 
  cat("\ntesting 'sl', 'ci' and 'cr' with  'rWy' = \n")
  cat( "  ", rWy )
  cat("\n\n")

  cat("\ncompare 'sl' by methods CLR-MC and CLR:\n\n")
  cat("for sl(theta):\n")
  slmc <- mod$sl(theta,'mc',.001,TRUE)
  cat("\n\n")
  slclr <- mod$sl(theta,'clr',.001,TRUE)
  dif <- slclr - slmc
  if( dif/slmc < -0.1  ||  ( slmc < .1  &&  dif > 0.01) )
  {
    cat(" re-check CLR-MC:\n")
    slmc <- mod$sl(theta,'mc',.0001,TRUE)
    dif <- slclr - slmc
    if( dif/slmc < -0.1  ||  ( slmc < .1  &&  dif > 0.01) )
      warning(
      gettextf("SL by mc %d%% different from clr %d%%",
        as.integer(100*slmc), as.integer(100*slclr)) )
    cat("\n")
  }


  if( xint )  {
    cat("\n\nfor sl(theta,alpha):\n")
    slmc <- mod$sl(theta,alpha,'mc',.001,TRUE)
    cat("\n\n")
    slclr <- mod$sl(theta,alpha,'clr',.001,TRUE)
    dif <- slclr - slmc
    if( dif/slmc < -0.1  ||  ( slmc < .1  &&  dif > 0.01) )
    {
      cat(" re-check CLR-MC:\n")
      slmc <- mod$sl(theta,alpha,'mc',.0001,TRUE)
      dif <- slclr - slmc
      if( dif/slmc < -0.1  || ( slmc < .1  &&  dif > 0.01) )
        warning(
          gettextf("SL by mc %d%% different from clr %d%%",
            as.integer(100*slmc), as.integer(100*slclr)) 
        )
      cat("\n")
    }
  }

  cat("\n\ntest 'ci'\n")
  mod$ci()

  cat("test 'cr'\n")
  mod$cr(output="T")
  cat("\n")
}


testrun()




