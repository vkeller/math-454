#include <math.h>

static void simple_dgemm(int n, double alpha, const double *A, const double *B,
                         double beta, double *C)
{
    int i;
    int j;
    int k;
    for (i = 0; i < n; ++i)
      {
        for (j = 0; j < n; ++j)
          {
            double prod = 0;
            for (k = 0; k < n; ++k)
              {
                prod += A[k * n + i] * B[j * n + k];
              }
            //std::cout << prod << std::endl;
            C[j * n + i] = alpha*prod + beta*C[j*n + i];
          }
      }
}



double verifyResult(const double *mat, const double *mat_ref, int M, int N)
{
  	double norm = 0.0;
	int i;
	int j;

    for (i = 0; i < M; i++)
    	{
        for (j = 0; j < N; j++)
          {
            //mat[i+j*M] << " " << mat_ref[i + j*M] << std::endl;
            if (abs((double)mat[i + j*M ] - (double)mat_ref[i + j*M ]) > norm)
              {
                norm = abs((double)mat[i + j*M] - (double)mat_ref[i + j*M]);
		if (norm != 0.)
			printf("%d %d = %f\n", i, j, norm); 
              }
          }
//        std::cout << "----" << std::endl;
      }
    return norm;
}


static void fill(double *A, int size, double v)
{
	int ii;
  	for (ii = 0; ii < size; ++ii)
    	A[ii] = (double) v;
}

static void eye(double *A, int lda, int n)
{
   	fill(A, lda*n, 0.);
	int ii;

   	for (ii = 0; ii < lda; ++ii)
    	A[ii*lda + ii] = 1.;
}
