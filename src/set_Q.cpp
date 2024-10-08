//

#include "lmbr.h"



void  Clmbr::set_Q( void )
//  Canonically reduce the model, using an orthogonal matrix Q for change-of-basis.
{
	int i,j;

	double *  X= R_Calloc( n*xrank, double );
	for (i=0;i<n;i++)  for (j=0;j<xrank;j++)  
		if( model_in > 0 )  *(X+j*n+i) = *(x_in+j*n+i);  else  *(X+j*n+i) = *(x_in+j*n+(n-1-i));

	int xcol;
	if(Model==M3) xcol=0; else xcol=1;
	if( model_in < 0 )  for (i=0;i<n;i++)  *(X+xcol*n+i) *= -1.;


//	compute  X =  irS * X :
	double *  temp= R_Calloc( n, double );
	if( vectorS )  for(j=0;j<xrank;j++)  for(i=0;i<n;i++)  *(X+j*n+i) *=  *( irS + i );
	if( matrixS )  for(j=0;j<xrank;j++)  {
		for(i=0;i<n;i++)  temp[i] = *(X+j*n+i);
		for(i=0;i<n;i++) {
			*(X+j*n+i) = 0.;
			for( int k=0; k<n; k++) *(X+j*n+i) += *( irS + k*n + i ) * temp[k];
		}
	} 
	R_Free( temp );

// 'Q' and 'tau' are global arrays, declared in "lmbr.h"
// set the pre-orthogonalized columns of 'Q' to be the columns of  irS * X ,  
// but with the  1- and x1- vectors as the final two
	for(i=0;i<n;i++)  for(j=xcol+1;j<xrank;j++)  *(Q + (j-xcol-1)*n + i ) =  *(X+j*n+i);
	for(i=0;i<n;i++)  for(j=0;j<xcol+1;j++)  *(Q + (j+xrank-xcol-1)*n + i ) =  *(X+j*n+i);
	R_Free( X );

//  Use LAPACK routines DGEQRF to generate Q and DORMQR to multiply by Q .
	int  lwork,  info;
	double  tmp[1];
	{
		lwork = -1;
		int  xrank_int= static_cast<int>( xrank );

		F77_CALL(dgeqrf)( &n_int, &xrank_int, Q, &n_int, tau, tmp, &lwork, &info );

		if( info )  stop( _("LAPACK routine 'dgeqrf' failed") );  else  lwork=  static_cast<int>( *tmp ); 
		double *  work= R_Calloc( lwork, double );

		F77_CALL(dgeqrf)( &n_int, &xrank_int, Q, &n_int, tau, work, &lwork, &info );

		if( info )  stop( _("LAPACK routine 'dgeqrf' failed") );
		R_Free( work );
	}


	return;
}


