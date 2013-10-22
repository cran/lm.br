//

#include "lmbr.h"




RCPP_MODULE(Clmbr){

	using namespace Rcpp ;


	void (Clmbr::*sl1) (int, double, double) = &Clmbr::slR ;
	void (Clmbr::*sl2) (int, double, double, double) = &Clmbr::slR ;
	double (Clmbr::*sl3) (int, int, int, double, double) = &Clmbr::slR ;
	double (Clmbr::*sl4) (int, int, int, double, double, double) = &Clmbr::slR ;

	void (Clmbr::*ci) (double, int) = &Clmbr::ciR ;

	void (Clmbr::*cr1) (double, int, double) = &Clmbr::crR ;
	NumericMatrix (Clmbr::*cr2) (double, int, double, int) = &Clmbr::crR ;

	void (Clmbr::*mle) (void)  const  = &Clmbr::MLE ;

	NumericVector (Clmbr::*param) (void)  const  = &Clmbr::PARAM ;

	void (Clmbr::*sety) (NumericVector) = &Clmbr::SET_rWy ;


	class_<Clmbr>( "Cpp_Clmbr" )
	
	    // expose the constructor
	    .constructor< NumericVector, NumericMatrix, NumericMatrix, int, int, int >()    

		// expose methods
		.method( "sl", sl1 , "SL for th0 by specified method, accuracy" )
		.method( "sl", sl2 , "SL for (th0,a0) by specified method, accuracy" )
		.method( "sl", sl3 , "SL for th0 by specified method, accuracy and output flag" )
		.method( "sl", sl4 , "SL for (th0,a0) by specified method, accuracy and output flag" )
		.method( "ci", ci , "confidence interval for theta by specified method" )
		.method( "cr", cr1 , 
			"confidence region for (theta,alpha) by specified method and increment" )
		.method( "cr", cr2 , 
			"confidence region for (theta,alpha) by specified method and increment, return matrix" )
		.method( "mle", mle , "printout maximum likelihood estimates of parameters" )
		.method( "param", param , "return maximum likelihood estimates of parameters" )
		.method( "sety", sety , 
			"reset values for square-root of weights times y-vector" )
	;

}                     


