#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>
#include <stdarg.h>
#include <stddef.h>

const double NEARZERO = 1.0e-10;       // interpretation of "zero"

// Prototypes
//void print_vec( char title[], double *v , int n);
//void print_mat( char title[], double **A, int m, int n );
//void matrixTimesVector( double **A, double *v, double *C, int n );
void vectorCombination( double a, double *U, double b, double *V , double *W, int n);
//double dotProduct( double *U, double *V, int n );
double vectorNorm( double *V , int n );
//void copyVector( double *V, double *W , int n);
//void subVector( double *x, double *V, double *W , int n);
//void addVector( double *x, double *V, double *W , int n);
//void saxpy(double a, double * x, double * y, int n);

//void conjugateGradientSolver( double **A, double *B, double *X, int m, int n );


//======================================================================


int main(){

	double ** A;
	double * b;
	double * x;
	double * x0;
	double * x_check;

	int i,n,m;

	n = 2;
	m = n;

	A = (double**) malloc(n* sizeof(double*));
	b = (double*) malloc(n* sizeof(double));
	x = (double*) malloc(n* sizeof(double));
	x0 = (double*) malloc(n* sizeof(double));
	x_check = (double*) malloc(n* sizeof(double));
	for (i = 0; i<n; i++){A[i] = (double*) malloc(n*sizeof(double));x[i]=0.;x_check[i]=0.;}

	//   {{ 4, 1},{ 1, 3}}

	A[0][0] = 4.;
	A[0][1] = 1.;
	A[1][0] = 1.;
	A[1][1] = 3.;



	// vec B = { 1, 2 };
	b[0] = 1.;
	b[1] = 2.;

	x0[0] = 2;
	x0[1] = 1;

	print_mat( "matrice A:\n", A,m,n );
	print_vec( "vecteur b: \n", b,n );

	conjugateGradientSolver( A, b, x0, m, n );

	printf("Solves AX = B\n");
	print_vec( "\nX:", x,n );
	matrixTimesVector( A, x, x_check, n ) ;
	print_vec( "\nCheck AX:\n", x_check, n);

	for (i = 0; i<n; i++){free(A[i]);}
	free(A);
	free(b);
	free(x);
	free(x0);
	free(x_check);

}


//======================================================================


void print_vec( char title[] , double *v, int n )
{
	double x;
	printf("%s",title);

	for ( int i = 0; i < n; i++ ){
		x = v[i];   if ( abs( x ) < NEARZERO ) x = 0.0;
		printf("%f\t",x);
	}
	printf("\n");
}


//======================================================================


void print_mat( char title[], double **A, int m, int n )
{
	double x;
	printf("%s",title);

	for ( int i = 0; i < m; i++ ){
		for ( int j = 0; j < n; j++ ){
			x = A[i][j];   if ( abs( x ) < NEARZERO ) x = 0.0;
			printf("%f\t",x);
		}
		printf("\n");
	}
}


//======================================================================


void matrixTimesVector( double **A, double *v, double *C, int n )     // Matrix times vector
{
	int i,j;
	double tmp;
	for ( i = 0; i < n; i++ ){
		tmp = 0.;
		for ( j = 0; j < n; j++ ){
			tmp += A[i][j] * v[j];
		}
		C[i] = tmp;
	}
}


//======================================================================


void vectorCombination( double a, double *U, double b, double *V , double *W, int n)        // Linear combination of vectors
{
	int j;
	double tmp1, tmp2;
	for (j = 0; j < n; j++ ){
		tmp1 = U[j]; tmp2 = V[j];
		W[j] = a * tmp1 + b * tmp2;
	}
}


//======================================================================


double dotProduct( double *U, double *V, int n )  
{
	double ret = 0.;
	int i;
	for (i = 0;i < n ; i++){ret += U[i] * V[i];}
	return ret;
}


//======================================================================


double vectorNorm( double *V , int n)
{
   return sqrt( dotProduct( V, V , n) );
}

//======================================================================

void copyVector( double *V, double *W , int n)
{
	int i;
	for (i = 0;i < n ; i++){V[i] = W[i];}
}

//======================================================================

void subVector( double *x, double *V, double *W , int n)
{
	// x = V + W
	int i;
	for (i = 0;i < n ; i++){x[i] = V[i] - W[i];}
}

//======================================================================

void addVector( double *x, double *V, double *W , int n)
{
	// x = V + W
	int i;
	for (i = 0;i < n ; i++){x[i] = V[i] + W[i];}
}




//======================================================================


void conjugateGradientSolver( double **A, double *B, double *X, int m, int n )
{
	double TOLERANCE = 1.0e-10;

	double * R;
	double * P;
	double * Rold;
	double * AP;
	double alpha;
	double beta;

	int k = 0;

	R = (double*) malloc(n* sizeof(double));
	P = (double*) malloc(n* sizeof(double));
	Rold = (double*) malloc(n* sizeof(double));
	AP = (double*) malloc(n* sizeof(double));

//	copyVector(R,B, n);

	matrixTimesVector( A, X, AP, n );
	subVector(R, B, AP, n); 
	copyVector(P,R, n);
	while ( k < n ){
		copyVector(Rold,R, n);
		matrixTimesVector( A, P, AP, n );

		alpha = dotProduct( R, R, n ) / fmax( dotProduct( P, AP, n ), NEARZERO );

		vectorCombination( 1.0, X, alpha, P, X, n );            // Next estimate of solution
		vectorCombination( 1.0, R, -alpha, AP, Rold, n );          // Residual

		print_vec("Residuel :\n",Rold,n);
		print_vec("X :\n",X,n);

		printf("\t[Iteration %d] : vectorNorm(R) = %f and dotProduct = %f\n",k,vectorNorm( R, n ), dotProduct( R, R, n ) );
		if ( vectorNorm( R, n ) < TOLERANCE ) break;             // Convergence test

		beta = dotProduct( R, R, n ) / fmax( dotProduct( Rold, Rold, n ), NEARZERO );
		vectorCombination( 1.0, R, beta, P, P, n );             // Next gradient
		printf("\t[Iteration %d] : alpha = %f, beta = %f\n",k,alpha, beta);
		k++;
	}
/*
	free(R);
	free(P);
	free(Rold);
	free(AP);
*/
//return X;
}
