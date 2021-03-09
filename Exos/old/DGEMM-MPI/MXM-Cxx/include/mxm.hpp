/*
 * mxm.hpp
 *
 *  Created on: Mar 17, 2011
 *      Author: vkeller
 */

#ifndef MXM_HPP_
#define MXM_HPP_


#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <iomanip>

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

//#define sX 256
//#define sY 256


using namespace std;

class MXM{

public:
	void compute();

private:

	template<typename T> void initMat1D(T *mat, T value, int x, int y){
		int i;
		int t=x*y;
		for (i=0;i<t;i++){
			mat[i] = value;
		}
	}


	template<typename T> T verifyMat1D(T *mat, int x, int y){
		T ret=(T)0;
		int i;
		int t=x*y;
		for (i=0;i<t;i++){
			ret = ret+mat[i];
		}
		return ret;
	}


	template<typename T> T verifyMat2D(T *mat[], int x, int y){
		T ret=(T)0;
		int i, j;
		for (i=0;i<x;i++){
			for (j=0;j<y;j++){
				ret = ret+mat[i][j];
			}
		}
		return ret;
	}
	template<typename T> void initMat2D(T *mat[], T value, int x, int y){
		int i, j;
		for (i=0;i<x;i++){
			for (j=0;j<x;j++){
				mat[i][j] = (T)value;
			}
		}
	}

	void isAligned(void *p, int line);

	template<typename T> void matmult2D(T *SAp[],T *SBp[],T *SCp[], int sizeX, int sizeY){

		double elapsedtime;
		double gflops;
		double flops;
		int i, j, k;
	    struct timeval tv1,tv2;
	    T ver;

		cout << "Order\t   GF/s\t\t time\t\t FLOPS\t\t verif" << endl;

	    cout << "I-J-K\t";
		initMat2D(SAp,(T)4.0, sizeX, sizeY);
		initMat2D(SBp,(T)7.0, sizeX, sizeY);
		initMat2D(SCp,(T)0.0, sizeX, sizeY);
		gettimeofday(&tv1, (struct timezone*)0);
		for (i=0;i<sizeX;i++){
			for (j=0;j<sizeY;j++){
	//						#pragma ivdep
	//						#pragma vector aligned
				for (k=0;k<sizeY;k++){
					SCp[i][j]=SCp[i][j]+SAp[i][k]*SBp[k][j];
				}
			}
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		flops=(double)sizeX*(double)sizeY*(double)sizeY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));
		ver=(T)0;
		ver = verifyMat2D(SCp, sizeX, sizeY);
		cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  ver << endl;


	    cout << "I-K-J\t";
		initMat2D(SCp,(T)0.0, sizeX, sizeY);
	    gettimeofday(&tv1, (struct timezone*)0);
		for (i=0;i<sizeX;i++){
			for (k=0;k<sizeY;k++){
	//#pragma ivdep
	//#pragma vector aligned
				for (j=0;j<sizeY;j++){
					SCp[i][j]=SCp[i][j]+SAp[i][k]*SBp[k][j];
				}
			}
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		flops=(double)sizeX*(double)sizeY*(double)sizeY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));
		ver=(T)0;
		ver = verifyMat2D(SCp, sizeX, sizeY);
		cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  ver << endl;


	    cout << "J-I-K\t";
		initMat2D(SCp,(T)0.0, sizeX, sizeY);
	    gettimeofday(&tv1, (struct timezone*)0);
		for (j=0;j<sizeY;j++){
			for (i=0;i<sizeX;i++){
	//#pragma ivdep
				//#pragma vector aligned
				for (k=0;k<sizeY;k++){
					SCp[i][j]=SCp[i][j]+SAp[i][k]*SBp[k][j];
				}
			}
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		flops=(double)sizeX*(double)sizeY*(double)sizeY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));
		ver=(T)0;
		ver = verifyMat2D(SCp, sizeX, sizeY);
		cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  ver << endl;


	    cout << "J-K-I\t";
		initMat2D(SCp,(T)0.0, sizeX, sizeY);
	    gettimeofday(&tv1, (struct timezone*)0);
		for (j=0;j<sizeY;j++){
			for (k=0;k<sizeY;k++){
	//#pragma ivdep
	//#pragma vector aligned
				for (i=0;i<sizeX;i++){
					SCp[i][j]=SCp[i][j]+SAp[i][k]*SBp[k][j];
				}
			}
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		flops=(double)sizeX*(double)sizeY*(double)sizeY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));
		ver=(T)0;
		ver = verifyMat2D(SCp, sizeX, sizeY);
		cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  ver << endl;


	    cout << "K-I-J\t";
		initMat2D(SCp,(T)0.0, sizeX, sizeY);
	    gettimeofday(&tv1, (struct timezone*)0);
		for (k=0;k<sizeY;k++){
			for (i=0;i<sizeX;i++){
				//#pragma ivdep
				//#pragma vector aligned
				for (j=0;j<sizeY;j++){
					SCp[i][j]=SCp[i][j]+SAp[i][k]*SBp[k][j];
				}
			}
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		flops=(double)sizeX*(double)sizeY*(double)sizeY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));
		ver=(T)0;
		ver = verifyMat2D(SCp, sizeX, sizeY);
		cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  ver << endl;

	    cout << "K-J-I\t";
		initMat2D(SCp,(T)0.0, sizeX, sizeY);
	    gettimeofday(&tv1, (struct timezone*)0);
		for (k=0;k<sizeY;k++){
			for (j=0;j<sizeY;j++){
	//#pragma ivdep
	//#pragma vector aligned
				for (i=0;i<sizeX;i++){
					SCp[i][j]=SCp[i][j]+SAp[i][k]*SBp[k][j];
				}
			}
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		flops=(double)sizeX*(double)sizeY*(double)sizeY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));
		ver=(T)0;
		ver = verifyMat2D(SCp, sizeX, sizeY);
		cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  ver << endl;

	}



	template<typename T> void matmult1D(T *SAp1D,T *SBp1D,T *SCp1D, int sizeX, int sizeY){

		double elapsedtime;
		double gflops;
		double flops;
		int i, j, k;
	    struct timeval tv1,tv2;
	    T ver=0.0;
		cout << "Order\t   GF/s\t\t time\t\t FLOPS\t\t verif " << endl;


	    cout << "I-J-K\t";
		initMat1D(SAp1D,(T)4.0, sizeX, sizeY);
		initMat1D(SBp1D,(T)7.0, sizeX, sizeY);
		initMat1D(SCp1D,(T)0.0, sizeX, sizeY);
		gettimeofday(&tv1, (struct timezone*)0);
		for (i=0;i<sizeX;i++){
			for (j=0;j<sizeY;j++){
				for (k=0;k<sizeY;k++){
					SCp1D[i*sizeY+j]=SCp1D[i*sizeY+j]+SAp1D[i*sizeY+k]*SBp1D[k*sizeX+j];
				}
			}
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		flops=(double)sizeX*(double)sizeY*(double)sizeY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));
		ver=(T)0;
		ver = verifyMat1D(SCp1D, sizeX, sizeY);
		cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  ver << endl;


		cout << "I-K-J\t";
		initMat1D(SCp1D,(T)0.0, sizeX, sizeY);
		gettimeofday(&tv1, (struct timezone*)0);
		for (i=0;i<sizeX;i++){
			for (k=0;k<sizeY;k++){
				for (j=0;j<sizeY;j++){
					SCp1D[i*sizeY+j]=SCp1D[i*sizeY+j]+SAp1D[i*sizeY+k]*SBp1D[k*sizeX+j];
				}
			}
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		flops=(double)sizeX*(double)sizeY*(double)sizeY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));
		ver=(T)0;
		ver = verifyMat1D(SCp1D, sizeX, sizeY);
		cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  ver << endl;


		cout << "J-I-K\t";
		initMat1D(SCp1D,(T)0.0, sizeX, sizeY);
		gettimeofday(&tv1, (struct timezone*)0);
		for (j=0;j<sizeY;j++){
			for (i=0;i<sizeX;i++){
				for (k=0;k<sizeY;k++){
					SCp1D[i*sizeY+j]=SCp1D[i*sizeY+j]+SAp1D[i*sizeY+k]*SBp1D[k*sizeX+j];
				}
			}
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		flops=(double)sizeX*(double)sizeY*(double)sizeY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));
		ver=(T)0;
		ver = verifyMat1D(SCp1D, sizeX, sizeY);
		cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  ver << endl;


		cout << "J-K-I\t";
		initMat1D(SCp1D,(T)0.0, sizeX, sizeY);
		gettimeofday(&tv1, (struct timezone*)0);
		for (j=0;j<sizeY;j++){
			for (k=0;k<sizeY;k++){
				for (i=0;i<sizeX;i++){
					SCp1D[i*sizeY+j]=SCp1D[i*sizeY+j]+SAp1D[i*sizeY+k]*SBp1D[k*sizeX+j];
				}
			}
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		flops=(double)sizeX*(double)sizeY*(double)sizeY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));
		ver=(T)0;
		ver = verifyMat1D(SCp1D, sizeX, sizeY);
		cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  ver << endl;




		cout << "K-I-J\t";
		initMat1D(SCp1D,(T)0.0, sizeX, sizeY);
		gettimeofday(&tv1, (struct timezone*)0);
		for (k=0;k<sizeY;k++){
			for (i=0;i<sizeX;i++){
				for (j=0;j<sizeY;j++){
					SCp1D[i*sizeY+j]=SCp1D[i*sizeY+j]+SAp1D[i*sizeY+k]*SBp1D[k*sizeX+j];
				}
			}
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		flops=(double)sizeX*(double)sizeY*(double)sizeY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));
		ver=(T)0;
		ver = verifyMat1D(SCp1D, sizeX, sizeY);
		cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  ver << endl;


		cout << "K-J-I\t";
		initMat1D(SCp1D,(T)0.0, sizeX, sizeY);
		gettimeofday(&tv1, (struct timezone*)0);
		for (k=0;k<sizeY;k++){
			for (j=0;j<sizeY;j++){
				for (i=0;i<sizeX;i++){
					SCp1D[i*sizeY+j]=SCp1D[i*sizeY+j]+SAp1D[i*sizeY+k]*SBp1D[k*sizeX+j];
				}
			}
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		flops=(double)sizeX*(double)sizeY*(double)sizeY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));
		ver=(T)0;
		ver = verifyMat1D(SCp1D, sizeX, sizeY);
		cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  ver << endl;


	}































/*
 * ----------------------------------------------------
 * Versions able to write in a file
 * ----------------------------------------------------
 */


	template<typename T> void matmult2D_GF(T *SAp[],T *SBp[],T *SCp[], int sizeX, int sizeY, ofstream * file){

		double elapsedtime;
		double gflops;
		double flops;
		int i, j, k;
	    struct timeval tv1,tv2;
	    T ver;


//	    cout << "I-J-K\t";
		initMat2D(SAp,(T)4.0, sizeX, sizeY);
		initMat2D(SBp,(T)7.0, sizeX, sizeY);
		initMat2D(SCp,(T)0.0, sizeX, sizeY);
		gettimeofday(&tv1, (struct timezone*)0);
		for (i=0;i<sizeX;i++){
			for (j=0;j<sizeY;j++){
	//						#pragma ivdep
	//						#pragma vector aligned
				for (k=0;k<sizeY;k++){
					SCp[i][j]=SCp[i][j]+SAp[i][k]*SBp[k][j];
				}
			}
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		flops=(double)sizeX*(double)sizeY*(double)sizeY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));

		*file << gflops << "\t" ;

		ver=(T)0;
		ver = verifyMat2D(SCp, sizeX, sizeY);
//		cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  ver << endl;


//	    cout << "I-K-J\t";
		initMat2D(SCp,(T)0.0, sizeX, sizeY);
	    gettimeofday(&tv1, (struct timezone*)0);
		for (i=0;i<sizeX;i++){
			for (k=0;k<sizeY;k++){
	//#pragma ivdep
	//#pragma vector aligned
				for (j=0;j<sizeY;j++){
					SCp[i][j]=SCp[i][j]+SAp[i][k]*SBp[k][j];
				}
			}
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		flops=(double)sizeX*(double)sizeY*(double)sizeY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));
		*file << gflops << "\t" ;
		ver=(T)0;
		ver = verifyMat2D(SCp, sizeX, sizeY);
//		cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  ver << endl;


//	    cout << "J-I-K\t";
		initMat2D(SCp,(T)0.0, sizeX, sizeY);
	    gettimeofday(&tv1, (struct timezone*)0);
		for (j=0;j<sizeY;j++){
			for (i=0;i<sizeX;i++){
	//#pragma ivdep
				//#pragma vector aligned
				for (k=0;k<sizeY;k++){
					SCp[i][j]=SCp[i][j]+SAp[i][k]*SBp[k][j];
				}
			}
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		flops=(double)sizeX*(double)sizeY*(double)sizeY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));
		*file << gflops << "\t" ;
		ver=(T)0;
		ver = verifyMat2D(SCp, sizeX, sizeY);
//		cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  ver << endl;


//	    cout << "J-K-I\t";
		initMat2D(SCp,(T)0.0, sizeX, sizeY);
	    gettimeofday(&tv1, (struct timezone*)0);
		for (j=0;j<sizeY;j++){
			for (k=0;k<sizeY;k++){
	//#pragma ivdep
	//#pragma vector aligned
				for (i=0;i<sizeX;i++){
					SCp[i][j]=SCp[i][j]+SAp[i][k]*SBp[k][j];
				}
			}
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		flops=(double)sizeX*(double)sizeY*(double)sizeY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));
		*file << gflops << "\t" ;
		ver=(T)0;
		ver = verifyMat2D(SCp, sizeX, sizeY);
//		cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  ver << endl;


//	    cout << "K-I-J\t";
		initMat2D(SCp,(T)0.0, sizeX, sizeY);
	    gettimeofday(&tv1, (struct timezone*)0);
		for (k=0;k<sizeY;k++){
			for (i=0;i<sizeX;i++){
				//#pragma ivdep
				//#pragma vector aligned
				for (j=0;j<sizeY;j++){
					SCp[i][j]=SCp[i][j]+SAp[i][k]*SBp[k][j];
				}
			}
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		flops=(double)sizeX*(double)sizeY*(double)sizeY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));
		*file << gflops << "\t" ;
		ver=(T)0;
		ver = verifyMat2D(SCp, sizeX, sizeY);
//		cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  ver << endl;

//	    cout << "K-J-I\t";
		initMat2D(SCp,(T)0.0, sizeX, sizeY);
	    gettimeofday(&tv1, (struct timezone*)0);
		for (k=0;k<sizeY;k++){
			for (j=0;j<sizeY;j++){
	//#pragma ivdep
	//#pragma vector aligned
				for (i=0;i<sizeX;i++){
					SCp[i][j]=SCp[i][j]+SAp[i][k]*SBp[k][j];
				}
			}
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		flops=(double)sizeX*(double)sizeY*(double)sizeY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));
		*file << gflops << "\t" << endl ;
		ver=(T)0;
		ver = verifyMat2D(SCp, sizeX, sizeY);
//		cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  ver << endl;

	}



	template<typename T> void matmult1D_GF(T *SAp1D,T *SBp1D,T *SCp1D, int sizeX, int sizeY, ofstream * file){

		double elapsedtime;
		double gflops;
		double flops;
		int i, j, k;
	    struct timeval tv1,tv2;
	    T ver=0.0;


//	    cout << "I-J-K\t";
		initMat1D(SAp1D,(T)4.0, sizeX, sizeY);
		initMat1D(SBp1D,(T)7.0, sizeX, sizeY);
		initMat1D(SCp1D,(T)0.0, sizeX, sizeY);
		gettimeofday(&tv1, (struct timezone*)0);
		for (i=0;i<sizeX;i++){
			for (j=0;j<sizeY;j++){
				for (k=0;k<sizeY;k++){
					SCp1D[i*sizeY+j]=SCp1D[i*sizeY+j]+SAp1D[i*sizeY+k]*SBp1D[k*sizeX+j];
				}
			}
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		flops=(double)sizeX*(double)sizeY*(double)sizeY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));
		*file << gflops << "\t" ;
		ver=(T)0;
		ver = verifyMat1D(SCp1D, sizeX, sizeY);
//		cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  ver << endl;


//		cout << "I-K-J\t";
		initMat1D(SCp1D,(T)0.0, sizeX, sizeY);
		gettimeofday(&tv1, (struct timezone*)0);
		for (i=0;i<sizeX;i++){
			for (k=0;k<sizeY;k++){
				for (j=0;j<sizeY;j++){
					SCp1D[i*sizeY+j]=SCp1D[i*sizeY+j]+SAp1D[i*sizeY+k]*SBp1D[k*sizeX+j];
				}
			}
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		flops=(double)sizeX*(double)sizeY*(double)sizeY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));
		*file << gflops << "\t" ;
		ver=(T)0;
		ver = verifyMat1D(SCp1D, sizeX, sizeY);
//		cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  ver << endl;


//		cout << "J-I-K\t";
		initMat1D(SCp1D,(T)0.0, sizeX, sizeY);
		gettimeofday(&tv1, (struct timezone*)0);
		for (j=0;j<sizeY;j++){
			for (i=0;i<sizeX;i++){
				for (k=0;k<sizeY;k++){
					SCp1D[i*sizeY+j]=SCp1D[i*sizeY+j]+SAp1D[i*sizeY+k]*SBp1D[k*sizeX+j];
				}
			}
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		flops=(double)sizeX*(double)sizeY*(double)sizeY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));
		*file << gflops << "\t" ;
		ver=(T)0;
		ver = verifyMat1D(SCp1D, sizeX, sizeY);
//		cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  ver << endl;


//		cout << "J-K-I\t";
		initMat1D(SCp1D,(T)0.0, sizeX, sizeY);
		gettimeofday(&tv1, (struct timezone*)0);
		for (j=0;j<sizeY;j++){
			for (k=0;k<sizeY;k++){
				for (i=0;i<sizeX;i++){
					SCp1D[i*sizeY+j]=SCp1D[i*sizeY+j]+SAp1D[i*sizeY+k]*SBp1D[k*sizeX+j];
				}
			}
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		flops=(double)sizeX*(double)sizeY*(double)sizeY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));
		*file << gflops << "\t" ;
		ver=(T)0;
		ver = verifyMat1D(SCp1D, sizeX, sizeY);
//		cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  ver << endl;




//		cout << "K-I-J\t";
		initMat1D(SCp1D,(T)0.0, sizeX, sizeY);
		gettimeofday(&tv1, (struct timezone*)0);
		for (k=0;k<sizeY;k++){
			for (i=0;i<sizeX;i++){
				for (j=0;j<sizeY;j++){
					SCp1D[i*sizeY+j]=SCp1D[i*sizeY+j]+SAp1D[i*sizeY+k]*SBp1D[k*sizeX+j];
				}
			}
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		flops=(double)sizeX*(double)sizeY*(double)sizeY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));
		*file << gflops << "\t" ;
		ver=(T)0;
		ver = verifyMat1D(SCp1D, sizeX, sizeY);
//		cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  ver << endl;


//		cout << "K-J-I\t";
		initMat1D(SCp1D,(T)0.0, sizeX, sizeY);
		gettimeofday(&tv1, (struct timezone*)0);
		for (k=0;k<sizeY;k++){
			for (j=0;j<sizeY;j++){
				for (i=0;i<sizeX;i++){
					SCp1D[i*sizeY+j]=SCp1D[i*sizeY+j]+SAp1D[i*sizeY+k]*SBp1D[k*sizeX+j];
				}
			}
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		flops=(double)sizeX*(double)sizeY*(double)sizeY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));
		*file << gflops << endl ;
		ver=(T)0;
		ver = verifyMat1D(SCp1D, sizeX, sizeY);
//		cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  ver << endl;


	}






};



#endif /* MXM_HPP_ */
