//
// constructor and destructor for class Clmbr
//


#include "lmbr.h"




Clmbr::Clmbr(  NumericVector  yR,  NumericMatrix  xR,  NumericMatrix  wR,  int model_num,
						int  inv,  int  var_k )
// constructor
{
	w_in = x_in = y_in = xs = NULL;
	is = NULL;
	rS = irS = Q = tau = NULL;
	q11 = qx1 = qxx = ck = qff = NULL;
	q10 = qx0 = a0 = b0 = NULL;
    f01 = f0x = NULL;
    B = C = NULL;
	px = psig1 = psigx = pv1h = pxh = NULL;
	nan_m1 = pnse1 = pnuse1 = pusen = NULL;
	nan_m = puqe1 = puqen = puqx = NULL;
	ps1 = psx = NULL;
	pq1 = pqx = NULL;
	pmq1 = pm1h = NULL; 
	py = psy = NULL;
	pqy = NULL;


	R_xlen_t  i, j;
	R_xlen_t  zero= static_cast<R_xlen_t>(0),  one= static_cast<R_xlen_t>(1) ;

	n =  yR.size() ;
	n_int =  static_cast<int>( n ) ;

	xrank = xR.ncol();

	model_in = model_num;

	variance_unknown= !( static_cast<bool>(var_k) );

	inverse= static_cast<bool>(inv);

	bool  cov_matrix_I = true;
	cov_matrix_diagonal = true;

	if( wR(zero,zero) > -0.5 )  {		// flag for null 'weights'
		if( wR.ncol()==one )  {
			for (i=zero;i<n;i++)  if( fabs( wR(i,zero) - 1. ) > zero_eq )  cov_matrix_I = false;
		}  else  {
			for (i=zero;i<n;i++) for (j=zero;j<n;j++) {
				if (i==j &&  fabs( wR(i,j) - 1. )>zero_eq) cov_matrix_I = false;
				if (i!=j && fabs( wR(i,j) )>zero_eq) {cov_matrix_I = false; cov_matrix_diagonal = false;}
			}
		}
	}

	vectorS = false;
	matrixS = false;
	if( !cov_matrix_I )  {
		if( cov_matrix_diagonal )  vectorS = true;  else  matrixS = true;
	}


	y_in = R_Calloc( n, double );  
	x_in = R_Calloc( n*xrank, double ); 
	if( vectorS )  w_in = R_Calloc( n, double );  
	if( matrixS )  w_in = R_Calloc( n*n, double );  
 

// store input values

	for (i=zero;i<n;i++) {
		y_in[i] =  yR[i];
		for(j=zero; j< xrank; j++)  *(x_in+j*n+i) = xR(i,j);
		if( vectorS )  { if( wR.ncol()==one )  *(w_in+i) = wR(i,zero);  else  *(w_in+i) = wR(i,i); }
		if( matrixS )  for(j=zero; j< n; j++)  *(w_in+j*n+i) = wR(i,j);
	}



	initialize();
}





Clmbr::Clmbr( const Clmbr  &initM )
//copy constructor
{
	w_in = x_in = y_in = xs = NULL;
	is = NULL;
	rS = irS = Q = tau = NULL;
	q11 = qx1 = qxx = ck = qff = NULL;
	q10 = qx0 = a0 = b0 = NULL;
    f01 = f0x = NULL;
    B = C = NULL;
	px = psig1 = psigx = pv1h = pxh = NULL;
	nan_m1 = pnse1 = pnuse1 = pusen = NULL;
	nan_m = puqe1 = puqen = puqx = NULL;
	ps1 = psx = NULL;
	pq1 = pqx = NULL;
	pmq1 = pm1h = NULL; 
	py = psy = NULL;
	pqy = NULL;

	
	model_in = initM.model_in;

	variance_unknown= initM.variance_unknown;

	inverse= initM.inverse;

	n = initM.n;
	n_int = initM.n_int;

	xrank = initM.xrank;

	cov_matrix_diagonal = initM.cov_matrix_diagonal;
	vectorS = initM.vectorS;
	matrixS = initM.matrixS;

	y_in = R_Calloc( n, double );  
	x_in = R_Calloc( n*xrank, double ); 
	if( vectorS )  w_in = R_Calloc( n, double );  
	if( matrixS )  w_in = R_Calloc( n*n, double );  


// store input values

	for (int i=0;i<n;i++) {
		y_in[i] = (*initM.py)[i];
		for(int j=0; j< xrank; j++)  *(x_in+j*n+i) = *(initM.x_in+j*n+i);
		if( vectorS )  w_in[i] = initM.w_in[i];
		if( matrixS )  for(int j=0; j< n; j++)  *(w_in+j*n+i) = *(initM.w_in+j*n+i);
	}


	initialize();
}




Clmbr::~Clmbr()
// destructor
{
	const Vector<double>  zeroVec(0);
	
	R_xlen_t  zero= static_cast<R_xlen_t>(0),  one= static_cast<R_xlen_t>(1) ;

	*px = *pv1h = *pxh = *psig1 = *psigx = *nan_m1 = *nan_m = zeroVec;
	*pnse1 = *pnuse1 = *pusen = *puqe1 = *puqen = *puqx = zeroVec;
	*py = *psy = *pqy = zeroVec;
	for(R_xlen_t i=zero; i<ns+one; i++) {
		ps1[i]= zeroVec;  psx[i]= zeroVec;  pq1[i]= zeroVec;  pqx[i]= zeroVec;  pmq1[i]= zeroVec;
	}
	if(Model==M3)  *pm1h = zeroVec;

	R_Free( w_in );  R_Free( x_in );  R_Free( y_in );  R_Free( xs );
	R_Free( px );
	R_Free( rS );  R_Free( irS );  R_Free( Q ); R_Free( tau );
	R_Free( is );
	R_Free( q11 );  R_Free( qx1 );  R_Free( qxx );  R_Free( ck );  R_Free( qff );
	R_Free( q10 );  R_Free( qx0 );  R_Free( a0 );  R_Free( b0 );
	R_Free( f01 );  R_Free( f0x );
	R_Free( B );  R_Free( C );
	R_Free( psig1 );  R_Free( psigx );  R_Free( pv1h );  R_Free( pxh );
	R_Free( nan_m1 );  R_Free( pnse1 );  R_Free( pnuse1 );  R_Free( pusen );
	R_Free( nan_m );  R_Free( puqe1 );  R_Free( puqen );  R_Free( puqx );
	R_Free( ps1 );  R_Free( psx );
	R_Free( pq1 );  R_Free( pqx );
	R_Free( pmq1 ); R_Free( pm1h );
	R_Free( py );  R_Free( psy );
	R_Free( pqy );

}


