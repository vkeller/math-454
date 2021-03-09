/*
 * mxv.hpp
 *
 *  Created on: Mar 21, 2011
 *      Author: vkeller
 */

#ifndef MXV_HPP_
#define MXV_HPP_

#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <limits>

#include <stdio.h>
#include <stdlib.h>

#ifdef WRITEINFILE
#include <fstream>
#include <string>
#include <sstream>
#endif

#ifdef INTEL
	#include <smmintrin.h>
	#include "mkl_cblas.h"
	#define _aligned_free(a) _mm_free(a)
	#define _aligned_malloc(a, b) _mm_malloc(a, b)
	#define ALIGN_16 __attribute__((aligned(16)))
#elif BLUEGENE
	#include <builtins.h>
	#include "essl.h"
	#define _aligned_free(a) free(a)
	#define _aligned_malloc(a, b) malloc(a)
	#define ALIGN_16

#else
	#define _aligned_free(a) free(a)
	#define _aligned_malloc(a, b) malloc(a)
	#define ALIGN_16
#endif



//#include <mxm.hpp>


using namespace std;

/*
void cblas_dgemv(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TA,
                 const int M, const int N, const double alpha, const double *A,
                 const int lda, const double *X, const int incX,
                 const double beta, double *Y, const int incY)

////////////////// 2D //////////////////

template<typename T> void cblas_call(const int M, const int N, const T alpha, const T *A,
        const int lda, const T *X, const int incX,
        const T beta, T *Y, const int incY){}

template <> void cblas_call<double>(const int M, const int N, const double alpha, const double *A,
        const int lda, const double *X, const int incX,
        const double beta, double *Y, const int incY){
	cblas_dgemv(CblasRowMajor,CblasNoTrans, M, N, alpha, A, lda, X, incX, beta, Y, incY);
}
template <> void cblas_call<float>(const int M, const int N, const float alpha, const float *A,
        const int lda, const float *X, const int incX,
        const float beta, float *Y, const int incY){
	cblas_sgemv(CblasRowMajor,CblasNoTrans, M, N, alpha, A, lda, X, incX, beta, Y, incY);
}

////////////////// 1D //////////////////

template<typename T> void cblas_call2(const int M, const int N, const T alpha, const T *A,
        const int lda, const T *X, const int incX,
        const T beta, T *Y, const int incY){}

template <> void cblas_call2<double>(const int M, const int N, const double alpha, const double *A,
        const int lda, const double *X, const int incX,
        const double beta, double *Y, const int incY){
	cblas_dgemv(CblasRowMajor,CblasNoTrans, M, N, alpha, A, lda, X, incX, beta, Y, incY);
}
template <> void cblas_call2<float>(const int M, const int N, const float alpha, const float *A,
        const int lda, const float *X, const int incX,
        const float beta, float *Y, const int incY){
	cblas_sgemv(CblasRowMajor,CblasNoTrans, M, N, alpha, A, lda, X, incX, beta, Y, incY);
}

*/

class MXV{

public:
	void compute();

private:
	template<typename T> void initMat(T **mat, T value, int x, int y){
		int i, j;
		for (i=0;i<x;i++){
			for (j=0;j<y;j++){
				mat[i][j] = value;
			}
		}
	}

	template<typename T> T verifyMat(T **mat, int x, int y){
		T ret=(T)0;
		int i, j;
		for (i=0;i<x;i++){
			for (j=0;j<y;j++){
				ret = ret+mat[i][j];
			}
		}
		return ret;
	}

	template<typename T> void initVec(T *vec, T value, int x){
		int i;
		for (i=0;i<x;i++){
			vec[i] = value;
		}
	}
	template<typename T> T verifyVec(T *vec, int x){
		T ret=(T)0;
		int i;
		for (i=0;i<x;i++){
			ret = ret+vec[i];
		}
		return ret;
	}

/*
	template<typename T> T matmult(int nn, int ndiag, int * jcoef, T *DUa, T** DAMAT, T *DU){
	    struct timeval tv1,tv2;
	    double elapsedtime;
	    int it, i,j, ni, jc,jt;
	    double flops, gflops;
	    T sum=(T)0;

	    int nnm1=nn-1;
	    int nnm2=nn-2;
	    int nnjc;
	    int jcm2;
	    int jcm1;

		cout << "Start computation (matmult1) " << endl;
	    gettimeofday(&tv1, (struct timezone*)0);
	    ni=100000000/nn;
		for (it=0;it<ni;it++){

			MXV::initVec(DUa,(T)0,nn);

			// First three diagonals (5,0,1)
			DUa[0] = DUa[0] + DAMAT[0][0]*DU[0]+DAMAT[0][1]*DU[1];
			for (j=1;j<nnm1;j++){
				DUa[j] = DUa[j] + DAMAT[j][0]*DU[0]+DAMAT[j][1]*DU[1] + DAMAT[j][5]*DU[j-1];
			}
			DUa[nnm1]=DUa[nnm1]+DAMAT[nnm1][0]*DU[nnm1]+ DAMAT[nnm1][5]*DU[nnm2];



			// Second three diagonals (2,3,4)
			jc = jcoef[4] +1;
			nnjc=nn-jc;
			jcm2 = jc-2;
			for (j=0;j<nnjc;j++){
				jt = jcm2+j;
				DUa[jt] = DUa[jt] + DAMAT[jt][2]*DU[jt-2]+DAMAT[jt][3]*DU[jt-1] + DAMAT[jt][4]*DU[jt];
			}
			jt=nn-jc+2;
			DUa[jt]=DUa[jt]+DAMAT[jt][2]*DU[nn-2]+DAMAT[jt][3]*DU[nn-1];
			DUa[jt+1]=DUa[jt+1]+DAMAT[jt+1][2]*DU[nn-1];

			//Next three columns (7,8,9)
			jc=-jcoef[8]+1;
			jcm2 = jc-2;
			jcm1 = jc-1;
			DUa[jcm2]=DUa[jcm2]+DAMAT[jcm2][6]*DU[0];
			DUa[jcm1]=DUa[jcm1]+DAMAT[jcm1][6]*DU[1]+DAMAT[jcm1][7]*DU[0];
			for (i=jc;i<nnm1;i++){
				j=i-jc+1;
				DUa[i]=DUa[i]+DAMAT[i][6]*DU[j+2]+DAMAT[i][7]*DU[j+1]+DAMAT[i][8]*DU[j];
			}
		}
		gettimeofday(&tv2, (struct timezone*)0);
		sum=MXV::verifyVec(DUa, nn);

		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
	    flops=(double)18*(double)ni*(double)nn;
	    gflops=flops/elapsedtime/1.e9;
		cout << "End computation" << endl;
		printRes(flops,elapsedtime,gflops);
		return sum;
	}

*/

	template<typename T> T matmult2D(int nn, int ndiag, T *DUa, T** DAMAT, T *DU, int ni){
	    struct timeval tv1,tv2;
	    double elapsedtime = (double)0;
	    int it, i,j, jc,jt;
	    double flops, gflops;
	    T sum=(T)0;
		for (it=0;it<ni;it++){

			MXV::initVec(DUa,(T)0,nn);

		    gettimeofday(&tv1, (struct timezone*)0);
			for (i=0;i<nn;i++){
				for (j=0;j<ndiag;j++){
					DUa[i] = DUa[i] + DAMAT[i][j]*DU[j];
				}
			}
			gettimeofday(&tv2, (struct timezone*)0);
			elapsedtime+=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;

		}
		sum=MXV::verifyVec(DUa, nn);

		flops = (double)2*(double)nn*(double)ndiag*(double)ni;
	    gflops=flops/elapsedtime/1.e9;
		printRes(flops,elapsedtime,gflops);
		return sum;
	}

	template<typename T> double matmult2D_GF(int nn, int ndiag, T *DUa, T** DAMAT, T *DU, int ni){
	    struct timeval tv1,tv2;
	    double elapsedtime = (double)0;
	    int it, i,j;
//	    int jc,jt;
	    double flops, gflops;
//	    T sum=(T)0;
		for (it=0;it<ni;it++){

			MXV::initVec(DUa,(T)0,nn);

		    gettimeofday(&tv1, (struct timezone*)0);
			for (i=0;i<nn;i++){
				for (j=0;j<ndiag;j++){
					DUa[i] = DUa[i] + DAMAT[i][j]*DU[j];
				}
			}
			gettimeofday(&tv2, (struct timezone*)0);
			elapsedtime+=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;

		}

		flops = (double)2*(double)nn*(double)ndiag*(double)ni;
	    gflops=flops/elapsedtime/1.e9;
		return gflops;
	}



	template<typename T> T matmult1D(int nn, int ndiag, T *Ua, T* AMAT, T *U, int ni){
	    struct timeval tv1,tv2;
	    double elapsedtime = (double)0;
	    int it, i,j, jc,jt;
	    double flops, gflops;
	    T sum=(T)0;
		for (it=0;it<ni;it++){

			MXV::initVec(Ua,(T)0,nn);

		    gettimeofday(&tv1, (struct timezone*)0);
			for (j=0;j<ndiag;j++){
				for (i=0;i<nn;i++){
					Ua[i] = Ua[i] + AMAT[i+j*nn]*U[i];
				}
			}
			gettimeofday(&tv2, (struct timezone*)0);
			elapsedtime+=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;

		}
		sum=MXV::verifyVec(Ua, nn);

		flops = (double)2*(double)nn*(double)ndiag*(double)ni;
	    gflops=flops/elapsedtime/1.e9;
		printRes(flops,elapsedtime,gflops);
		return sum;
	}

	template<typename T> double matmult1D_GF(int nn, int ndiag, T *Ua, T* AMAT, T *U, int ni){
	    struct timeval tv1,tv2;
	    double elapsedtime = (double)0;
	    int it, i,j;
	    double flops, gflops;
//	    T sum=(T)0;
		for (it=0;it<ni;it++){

			MXV::initVec(Ua,(T)0,nn);

		    gettimeofday(&tv1, (struct timezone*)0);
			for (j=0;j<ndiag;j++){
				for (i=0;i<nn;i++){
					Ua[i] = Ua[i] + AMAT[i+j*nn]*U[i];
				}
			}
			gettimeofday(&tv2, (struct timezone*)0);
			elapsedtime+=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;

		}
		flops = (double)2*(double)nn*(double)ndiag*(double)ni;
	    gflops=flops/elapsedtime/1.e9;
		return gflops;
	}



	/*

	y := alpha*A*x + beta*y

	void cblas_dgemv(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TA,
	                 const int M, const int N, const double alpha, const double *A,
	                 const int lda, const double *X, const int incX,
	                 const double beta, double *Y, const int incY)
	template<typename T> T matmult_gemv(int nn, int ndiag, T *Ua, T** AMAT, T *U, int ni){
		T sum = (T)0;
	    struct timeval tv1,tv2;
	    double elapsedtime;
	    int it, i,j, jc,jt;
	    double flops, gflops;

	    double some_double = (double)0;
	    float some_float = (float)0;

//		cout << "Start computation (gemv version) " << endl;

	    for (it=0;it<ni;it++){
			MXV::initVec(Ua,(T)0,nn);
		    gettimeofday(&tv1, (struct timezone*)0);
			cblas_call(ndiag,nn,(T)1.0,&AMAT[0][0],nn,U,1,(T)0.0,Ua,1);
			gettimeofday(&tv2, (struct timezone*)0);
			elapsedtime+=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		}

		sum=MXV::verifyVec(Ua, nn);

//		cblas_call(some_double);
//		cblas_call(some_float);

//	    flops=(double)18*(double)ni*(double)nn;
		flops = (double)2*(double)nn*(double)ndiag*(double)ni;
	    gflops=flops/elapsedtime/1.e9;
//		cout << "End computation" << endl;
		printRes(flops,elapsedtime,gflops);
		return sum;
	}


	template<typename T> T matmult2_gemv(int nn, int ndiag, T *Ua, T* AMAT, T *U, int ni){
		T sum = (T)0;
	    struct timeval tv1,tv2;
	    double elapsedtime;
	    int it, i,j, jc,jt;
	    double flops, gflops;

	    double some_double = (double)0;
	    float some_float = (float)0;

	    gettimeofday(&tv1, (struct timezone*)0);

	    for (it=0;it<ni;it++){
			MXV::initVec(Ua,(T)0,nn);
			cblas_call2(ndiag,nn,(T)1.0,AMAT,nn,U,1,(T)0.0,Ua,1);
		}

		gettimeofday(&tv2, (struct timezone*)0);
		sum=MXV::verifyVec(Ua, nn);

		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
	    flops=(double)18*(double)ni*(double)nn;
	    gflops=flops/elapsedtime/1.e9;
		printRes(flops,elapsedtime,gflops);
		return sum;
	}

*/

	void printRes(double flops, double elapsedtime, double gflops);


	template<typename T> void printSum(T res, T sum){
		cout << "Exact res       = " << (int)res << endl;
		cout << "Computed res    = " << (int)sum << endl;
		cout << "Diff (r_e-r_c)  = " << (int)(res-sum) << endl;
	}

};



#endif /* MXV_HPP_ */
