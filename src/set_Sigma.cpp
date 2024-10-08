//


#include "lmbr.h"




void Clmbr::set_Sigma( void ) 
// 'Sigma' = covariate or covariance vector or matrix,  errors ~ N( 0, var * Sigma )
// check the input 'weights' vector or matrix,  make symmetric if almost symmetric, 
// then get 'rS' = square root of Sigma, and 'irS' = inverse square root of Sigma 
{
	int i,j;

	if( vectorS )  {
		for (i=0;i<n;i++) {
			const double  wi = *(w_in + i);
			if ( !R_FINITE(wi) )  stop( _("'weights' has invalid entries") );
			if ( wi <= 0 )  stop( _("'weights' has invalid entries") );
		}

	}  else  {
		for (i=0;i<n;i++) for (j=0;j<n;j++) {
			double  wij = *(w_in+j*n+i),  wji = *(w_in+i*n+j);
			if ( !R_FINITE(wij) )  stop( _("'weights' has invalid entries") );
			if( fabs(wij - wji) < zero_eq )  if( i < j)  *(w_in+j*n+i) = wji;
		}
	}


	double  maxD =0.,  minD =Inf;

	if ( vectorS )  {

		for (i=0;i<n;i++) {
			double  Di;
			if( model_in > 0 )  Di= *( w_in + i );  else  Di= *( w_in + (n-1-i) );
			if( Di > maxD )  maxD= Di;
			if( Di < minD )  minD= Di;
			const double  rDi = sqrt( Di );
			if( inverse )
				{  *( rS + i ) = rDi;  *( irS + i ) = 1./rDi; }
			else
				{ *( irS + i ) = rDi;   *( rS + i ) = 1./rDi; }
		}

		if( minD <= 0. )  stop( _("zero or negative 'weights' not allowed") );
		if( minD/maxD < 1.e-7 )  Rf_warning( "%s", _("weights vector might be ill-conditioned for 'clr' method") );

	}  else  {

// use LAPACK routine DSYEVR to get eigenvalues and eigenvectors of Sigma 
		double *  D= R_Calloc( n, double ),  * Q_= R_Calloc( n*n, double );
		
		{
			const char  job ='V',  range ='A',  uplo ='L';
			const int  id =0;	
			const double  tol = 0, dd = 0;
			int  ne, lwork= -1, itmp[1], liwork =-1, info =0;
			int*  isuppZ= R_Calloc(2*n,int);
			double  tmp[1];

//  use 'irS' as the working input matrix
			for(i=0;i<n;i++) for(j=0;j<=i;j++)  
				if( model_in > 0 )  
					*(irS + j*n + i) = *(w_in + j*n + i);  
				else
					*(irS + j*n + i) = *(w_in + (n-1-j)*n + n-1-i);

			F77_CALL(dsyevr)( &job, &range, &uplo, &n_int, irS, &n_int, &dd, &dd, &id, &id, &tol,
								&ne, D, Q_, &n_int, isuppZ, tmp, &lwork, itmp, &liwork, &info FCONE FCONE FCONE);

			if( info )  stop( _("LAPACK routine 'dsyevr' failed") );  else  { lwork= *tmp; liwork= *itmp; }
			double *  work= R_Calloc( lwork, double );
			int*  iwork= R_Calloc( liwork, int );

			F77_CALL(dsyevr)( &job, &range, &uplo, &n_int, irS, &n_int, &dd, &dd, &id, &id, &tol,
								&ne, D, Q_, &n_int, isuppZ, work, &lwork, iwork, &liwork, &info FCONE FCONE FCONE);

			if( info || ne < n )  stop( _("LAPACK routine 'dsyevr' failed") );
			R_Free( isuppZ );  R_Free( work );  R_Free( iwork );
		}

		double *  rD= R_Calloc(n,double);

		for (i=0;i<n;i++)  {
			if( D[i] <= 0. )  stop( _("'weights' matrix not positive-definite") );
			if( D[i] > maxD )  maxD= D[i];
			if( D[i] < minD )  minD= D[i];
			rD[i] =  sqrt( D[i] );
		}

// 'rS' = Q*sqrt(D)*t(Q)  and  'irS' = Q*1/sqrt(D)*t(Q)
		for (i=0;i<n;i++)  for(j=0;j<n;j++)  {
			*(rS + j*n + i) = 0.;
			*(irS + j*n + i) = 0.;
			for(int k=0;k<n;k++)  
				if( inverse )  {
					 *(rS + j*n + i) +=  *(Q_+k*n+i) * rD[k] * ( *(Q_+k*n+j) );
					*(irS + j*n + i) +=  *(Q_+k*n+i) / rD[k] * ( *(Q_+k*n+j) );
				}  else  {
					 *(rS + j*n + i) +=  *(Q_+k*n+i) / rD[k] * ( *(Q_+k*n+j) );
					*(irS + j*n + i) +=  *(Q_+k*n+i) * rD[k] * ( *(Q_+k*n+j) );
				}
		}

		if( minD/maxD < 1.e-7 )  Rf_warning( "%s", _("weights matrix might be ill-conditioned for 'clr' method") );
		R_Free( D );  R_Free( Q_ );  R_Free( rD );
	}


	return;
}

