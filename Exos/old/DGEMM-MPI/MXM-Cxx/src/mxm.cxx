/*
 * mxm.cxx
 *
 *  Created on: Mar 17, 2011
 *      Author: vkeller
 */

#include <mxm.hpp>

/*
void MXM::isAligned(void *p, int line){
	if (((int)p % 16 != 0)) {
		cout << "Alignement error: pointer to " << p << " at source code line " << line << endl;
		exit(1);
	}
}
*/

#define SIZE_MAX 512
#define STEP 1

void MXM::compute(){


	double elapsedtime;
	double gflops;
	double flops;
	float Sver;
	double Dver;
	long double DDver;

	int i, j, k, iterator;

	int sX = 256;
	int sY = 256;



	float **SAp ALIGN_16, **SBp ALIGN_16, **SCp ALIGN_16;
	float *SAp1D ALIGN_16, *SBp1D ALIGN_16, *SCp1D ALIGN_16;

	double **DAp ALIGN_16, **DBp ALIGN_16, **DCp ALIGN_16;
	double *DAp1D ALIGN_16, *DBp1D ALIGN_16, *DCp1D ALIGN_16;

	long double **DDAp ALIGN_16, **DDBp ALIGN_16, **DDCp ALIGN_16;
	long double *DDAp1D ALIGN_16, *DDBp1D ALIGN_16, *DDCp1D ALIGN_16;


	float **SAp_na, **SBp_na, **SCp_na;
	float *SAp1D_na, *SBp1D_na, *SCp1D_na;

	double **DAp_na, **DBp_na, **DCp_na;
	double *DAp1D_na, *DBp1D_na, *DCp1D_na;

//	long double **DDAp_na, **DDBp_na, **DDCp_na;
//	long double *DDAp1D_na, *DDBp1D_na, *DDCp1D_na;



    struct timeval tv1,tv2;

//    float st1,st2,st3;
//    float st4,st5,st6;
//    double dt1,dt2,dt3;

//    int kt;



#ifdef WRITEINFILE
	bool printres=false;
#else
	bool printres=true;
#endif

#ifdef WRITEINFILE

	cout << "Open file" << endl;
    ofstream myfile_32_static("dataMXM-32-stat.dat",ios::out);
    myfile_32_static << "Size\tstatic" << endl;

    ofstream myfile_32_2D_align("dataMXM-32-2D-align.dat",ios::out);
    myfile_32_2D_align << "Size\tI-J-K\tI-K-J\tJ-I-K\tJ-K-I\tK-I-J\tK-J-I" << endl;

    ofstream myfile_32_2D_non_align("dataMXM-32-2D-non-align.dat",ios::out);
    myfile_32_2D_non_align << "Size\tI-J-K\tI-K-J\tJ-I-K\tJ-K-I\tK-I-J\tK-J-I" << endl;

    ofstream myfile_32_1D_align("dataMXM-32-1D-align.dat",ios::out);
    myfile_32_1D_align << "Size\tI-J-K\tI-K-J\tJ-I-K\tJ-K-I\tK-I-J\tK-J-I" << endl;

    ofstream myfile_32_1D_non_align("dataMXM-32-1D-non-align.dat",ios::out);
    myfile_32_1D_non_align << "Size\tI-J-K\tI-K-J\tJ-I-K\tJ-K-I\tK-I-J\tK-J-I" << endl;

    ofstream myfile_32_dgemm("dataMXM-32-dgemm.dat",ios::out);
    myfile_32_dgemm << "Size\tDGEMM" << endl;

    ofstream myfile_64_static("dataMXM-64-stat.dat",ios::out);
    myfile_64_static << "Size\tstatic" << endl;

    ofstream myfile_64_2D_align("dataMXM-64-2D-align.dat",ios::out);
    myfile_64_2D_align << "Size\tI-J-K\tI-K-J\tJ-I-K\tJ-K-I\tK-I-J\tK-J-I" << endl;

    ofstream myfile_64_2D_non_align("dataMXM-64-2D-non-align.dat",ios::out);
    myfile_64_2D_non_align << "Size\tI-J-K\tI-K-J\tJ-I-K\tJ-K-I\tK-I-J\tK-J-I" << endl;

    ofstream myfile_64_1D_align("dataMXM-64-1D-align.dat",ios::out);
    myfile_64_1D_align << "Size\tI-J-K\tI-K-J\tJ-I-K\tJ-K-I\tK-I-J\tK-J-I" << endl;

    ofstream myfile_64_1D_non_align("dataMXM-64-1D-non-align.dat",ios::out);
    myfile_64_1D_non_align << "Size\tI-J-K\tI-K-J\tJ-I-K\tJ-K-I\tK-I-J\tK-J-I" << endl;

    ofstream myfile_64_dgemm("dataMXM-64-dgemm.dat",ios::out);
    myfile_64_dgemm << "Size\tDGEMM" << endl;


#endif


    for (iterator=10; iterator<SIZE_MAX;iterator=iterator+STEP){

    	cout << "size : " << iterator << endl;

		// ------------------------------------------------------------------------------------------------------------------------------------------
		// ------------------------------------------------------------------------------------------------------------------------------------------
		// FLOAT VERSION
		// ------------------------------------------------------------------------------------------------------------------------------------------
		// ------------------------------------------------------------------------------------------------------------------------------------------

    	sX = iterator;
    	sY = iterator;

    	float SA[sX][sY], SB[sX][sY], SC[sX][sY];
    	double DA[sX][sY], DB[sX][sY], DC[sX][sY];
    	long double DDA[sX][sY], DDB[sX][sY], DDC[sX][sY];


		if (printres){
			cout << "Simple precision (32 bits) MXM (STATIC)" << endl;
			cout << "Order\t   GF/s\t\t time\t\t FLOPS\t\t verif" <<endl;
		}

		for (i=0;i<sX;i++){
			for (j=0;j<sY;j++){
				SA[i][j] = (float)4.0;
				SB[i][j] = (float)7.0;
				SC[i][j] = (float)0.0;
			}
		}

		gettimeofday(&tv1, (struct timezone*)0);

		for (i=0;i<sX;i++){
			for (j=0;j<sY;j++){
				for (k=0;k<sY;k++){
					SC[i][j]=SC[i][j]+SA[i][k]*SB[k][j];
				}
			}
		}

		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;

		flops=(double)sX*(double)sY*(double)sY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));

	#ifdef WRITEINFILE
		myfile_32_static << sX << "\t" << gflops << endl;
	#endif
		if (printres)
			cout << "Static\t";

		Sver = (float)0;
		for (i=0;i<sX;i++){
			for (j=0;j<sY;j++){
				Sver = Sver+SC[i][j];
			}
		}
		if (printres)
			cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  Sver << endl;









		// ---------------------------------------------------
		// Allocate memory for matrices

		if (printres)
			cout << "Simple precision (32 bits) MXM (ALIGNED PTR in 2D)" << endl;

		SAp = (float**) _aligned_malloc(sX * sizeof(float*), 16);
		SBp = (float**) _aligned_malloc(sX * sizeof(float*), 16);
		SCp = (float**) _aligned_malloc(sX * sizeof(float*), 16);

		for(i=0; i<sX; i++){
			SAp[i] = (float*) _aligned_malloc(sY * sizeof(float), 16);
			SBp[i] = (float*) _aligned_malloc(sY * sizeof(float), 16);
			SCp[i] = (float*) _aligned_malloc(sY * sizeof(float), 16);
		}

	#ifdef BLUEGENE
		#pragma disjoint (**SAp, **SBp, **SCp)
		__alignx(16,SAp);
		__alignx(16,SBp);
		__alignx(16,SCp);
	#endif

	#ifdef WRITEINFILE
		myfile_32_2D_align << sX << "\t";
		matmult2D_GF(SAp,SBp,SCp,sX,sY, &myfile_32_2D_align);
	#else
		matmult2D(SAp,SBp,SCp,sX,sY);
	#endif
		// Deallocate memory ;-)
		for(i=0; i<sX; i++){
			_aligned_free(SAp[i]);
			_aligned_free(SBp[i]);
			_aligned_free(SCp[i]);
		}
		_aligned_free(SAp);
		_aligned_free(SBp);
		_aligned_free(SCp);
	// ---------------------------------------------------





		// ---------------------------------------------------
		// Allocate memory for matrices

		if (printres)
			cout << "Simple precision (32 bits) MXM (NON ALIGNED PTR in 2D)" << endl;

		SAp_na = (float**) malloc(sX * sizeof(float*));
		SBp_na = (float**) malloc(sX * sizeof(float*));
		SCp_na = (float**) malloc(sX * sizeof(float*));

		for(i=0; i<sX; i++){
			SAp_na[i] = (float*) malloc(sY * sizeof(float));
			SBp_na[i] = (float*) malloc(sY * sizeof(float));
			SCp_na[i] = (float*) malloc(sY * sizeof(float));
		}

#ifdef WRITEINFILE
	myfile_32_2D_non_align << sX << "\t";
	matmult2D_GF(SAp_na,SBp_na,SCp_na,sX,sY, &myfile_32_2D_non_align);
#else
	matmult2D(SAp_na,SBp_na,SCp_na,sX,sY);
#endif

		// Deallocate memory ;-)
		for(i=0; i<sX; i++){
			free(SAp_na[i]);
			free(SBp_na[i]);
			free(SCp_na[i]);
		}
		free(SAp_na);
		free(SBp_na);
		free(SCp_na);
	// ---------------------------------------------------






		// ---------------------------------------------------
		// Allocate memory for matrices

		if (printres)
			cout << "Simple precision (32 bits) MXM (ALIGNED PTR in 1D)" << endl;

		SAp1D = (float*) _aligned_malloc(sX * sY * sizeof(float), 16);
		SBp1D = (float*) _aligned_malloc(sX * sY * sizeof(float), 16);
		SCp1D = (float*) _aligned_malloc(sX * sY * sizeof(float), 16);

	#ifdef BLUEGENE
		#pragma disjoint (*SAp1D, *SBp1D, *SCp1D)
		__alignx(16,SAp1D);
		__alignx(16,SBp1D);
		__alignx(16,SCp1D);
	#endif

#ifdef WRITEINFILE
	myfile_32_1D_align << sX << "\t";
	matmult1D_GF(SAp1D,SBp1D,SCp1D, sX,sY, &myfile_32_1D_align);
#else
	matmult1D(SAp1D,SBp1D,SCp1D, sX,sY);
#endif

		// Deallocate memory ;-)
			_aligned_free(SAp1D);
			_aligned_free(SBp1D);
			_aligned_free(SCp1D);
		// ---------------------------------------------------


			// ---------------------------------------------------
			// Allocate memory for matrices

		if (printres)
			cout << "Simple precision (32 bits) MXM (NON ALIGNED PTR in 1D)" << endl;

		SAp1D_na = (float*) malloc(sX * sY * sizeof(float));
		SBp1D_na = (float*) malloc(sX * sY * sizeof(float));
		SCp1D_na = (float*) malloc(sX * sY * sizeof(float));

#ifdef WRITEINFILE
	myfile_32_1D_non_align << sX << "\t";
	matmult1D_GF(SAp1D_na,SBp1D_na,SCp1D_na, sX,sY, &myfile_32_1D_non_align);
#else
	matmult1D(SAp1D_na,SBp1D_na,SCp1D_na, sX,sY);
#endif

		// Deallocate memory ;-)
		free(SAp1D_na);
		free(SBp1D_na);
		free(SCp1D_na);
		// ---------------------------------------------------


		// ---------------------------------------------------
		// Allocate memory for matrices

		if (printres)
			cout << "Simple precision (32 bits) MXM (NON ALIGNED PTR in 1D) GEMM version" << endl;

		SAp1D_na = (float*) malloc(sX * sY * sizeof(float));
		SBp1D_na = (float*) malloc(sX * sY * sizeof(float));
		SCp1D_na = (float*) malloc(sX * sY * sizeof(float));

		MXM::initMat1D(SAp1D_na,(float)4.0, sX, sY);
		MXM::initMat1D(SBp1D_na,(float)7.0, sX, sY);
		MXM::initMat1D(SCp1D_na,(float)0.0, sX, sY);


		gettimeofday(&tv1, (struct timezone*)0);
	#ifdef INTEL
		cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,sY,sX,sY,(float)1, SAp1D_na, max(sX, sY),SBp1D_na, max(sX, sY), (float)0, SCp1D_na, max(sX, sY));
	#elif BLUEGENE
		sgemm("n", "n", sY,sX,sY, (float)1,SAp1D_na,max(sX, sY),SBp1D_na,max(sX, sY), (float)0, SCp1D_na,max(sX, sY));
	#endif
		gettimeofday(&tv2, (struct timezone*)0);

		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;

		flops=(double)sX*(double)sY*(double)sY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));
#ifdef WRITEINFILE
	myfile_32_dgemm << sX << "\t" << gflops << endl;
#endif

		if (printres)
			cout << "GEMM  \t";

		Sver = (float)0;
		for (i=0;i<sX;i++){
			for (j=0;j<sY;j++){
				Sver = Sver+SCp1D_na[i*sY+j];
			}
		}
		if (printres)
			cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  Sver << endl;


		// Deallocate memory ;-)
		free(SAp1D_na);
		free(SBp1D_na);
		free(SCp1D_na);
		// ---------------------------------------------------







			// ------------------------------------------------------------------------------------------------------------------------------------------
			// ------------------------------------------------------------------------------------------------------------------------------------------
			// DOUBLE VERSION
			// ------------------------------------------------------------------------------------------------------------------------------------------
			// ------------------------------------------------------------------------------------------------------------------------------------------



		if (printres){
			cout << "Double precision (64 bits) MXM (STATIC)" << endl;
			cout << "Order\t   GF/s\t\t time\t\t FLOPS\t\t verif " << endl;
		}

		for (i=0;i<sX;i++){
			for (j=0;j<sY;j++){
				DA[i][j] = (double)4.0;
				DB[i][j] = (double)7.0;
				DC[i][j] = (double)0.0;
			}
		}

		gettimeofday(&tv1, (struct timezone*)0);

		for (i=0;i<sX;i++){
			for (j=0;j<sY;j++){
				for (k=0;k<sY;k++){
					DC[i][j]=DC[i][j]+DA[i][k]*DB[k][j];
				}
			}
		}

		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;

		flops=(double)sX*(double)sY*(double)sY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));

#ifdef WRITEINFILE
	myfile_64_static << sX << "\t" << gflops << endl;
#endif

		if (printres)
			cout << "Static\t";

		Dver = (double)0;
		for (i=0;i<sX;i++){
			for (j=0;j<sY;j++){
				Dver = Dver+DC[i][j];
			}
		}

		if (printres)
			cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  Dver << endl;






		// ---------------------------------------------------
		// Allocate memory for matrices

		if (printres)
			cout << "Double precision (64 bits) MXM (ALIGNED PTR 2D)" << endl;


		DAp = (double**) _aligned_malloc(sX * sizeof(double*), 16);
		DBp = (double**) _aligned_malloc(sX * sizeof(double*), 16);
		DCp = (double**) _aligned_malloc(sX * sizeof(double*), 16);

		for(i=0; i<sX; i++){
			DAp[i] = (double*) _aligned_malloc(sY * sizeof(double), 16);
			DBp[i] = (double*) _aligned_malloc(sY * sizeof(double), 16);
			DCp[i] = (double*) _aligned_malloc(sY * sizeof(double), 16);
		}
	#ifdef BLUEGENE
		#pragma disjoint (**DAp, **DBp, **DCp)
		__alignx(16,DAp);
		__alignx(16,DBp);
		__alignx(16,DCp);
	#endif

#ifdef WRITEINFILE
		myfile_64_2D_align << sX << "\t";
		matmult2D_GF(DAp,DBp,DCp,sX,sY, &myfile_64_2D_align);
#else
		matmult2D(DAp,DBp,DCp,sX,sY);
#endif

		// Deallocate memory ;-)
		for(i=0; i<sX; i++){
			_aligned_free(DAp[i]);
			_aligned_free(DBp[i]);
			_aligned_free(DCp[i]);
		}
		_aligned_free(DAp);
		_aligned_free(DBp);
		_aligned_free(DCp);
	// ---------------------------------------------------



		// ---------------------------------------------------
		// Allocate memory for matrices

		if (printres)
			cout << "Double precision (64 bits) MXM (NON ALIGNED PTR 2D)" << endl;


		DAp_na = (double**) malloc(sX * sizeof(double*));
		DBp_na = (double**) malloc(sX * sizeof(double*));
		DCp_na = (double**) malloc(sX * sizeof(double*));

		for(i=0; i<sX; i++){
			DAp_na[i] = (double*) malloc(sY * sizeof(double));
			DBp_na[i] = (double*) malloc(sY * sizeof(double));
			DCp_na[i] = (double*) malloc(sY * sizeof(double));
		}

#ifdef WRITEINFILE
		myfile_64_2D_non_align << sX << "\t";
		matmult2D_GF(DAp_na,DBp_na,DCp_na,sX,sY, &myfile_64_2D_non_align);
#else
		matmult2D(DAp_na,DBp_na,DCp_na,sX,sY);
#endif

		// Deallocate memory ;-)
		for(i=0; i<sX; i++){
			free(DAp_na[i]);
			free(DBp_na[i]);
			free(DCp_na[i]);
		}
		free(DAp_na);
		free(DBp_na);
		free(DCp_na);
	// ---------------------------------------------------



		// ---------------------------------------------------
		// Allocate memory for matrices

		if (printres)
			cout << "Double precision (64 bits) MXM (ALIGNED PTR in 1D)" << endl;

		DAp1D = (double*) _aligned_malloc(sX * sY * sizeof(double), 16);
		DBp1D = (double*) _aligned_malloc(sX * sY * sizeof(double), 16);
		DCp1D = (double*) _aligned_malloc(sX * sY * sizeof(double), 16);

	#ifdef BLUEGENE
		#pragma disjoint (*DAp1D, *DBp1D, *DCp1D)
		__alignx(16,DAp1D);
		__alignx(16,DBp1D);
		__alignx(16,DCp1D);
	#endif

#ifdef WRITEINFILE
		myfile_64_1D_align << sX << "\t";
		matmult1D_GF(DAp1D,DBp1D,DCp1D, sX,sY, &myfile_64_1D_align);
#else
		matmult1D(DAp1D,DBp1D,DCp1D, sX,sY);
#endif

		// Deallocate memory ;-)
			_aligned_free(DAp1D);
			_aligned_free(DBp1D);
			_aligned_free(DCp1D);
		// ---------------------------------------------------


		// ---------------------------------------------------
		// Allocate memory for matrices

		if (printres)
			cout << "Double precision (64 bits) MXM (NON ALIGNED PTR in 1D)" << endl;

		DAp1D_na = (double*) malloc(sX * sY * sizeof(double));
		DBp1D_na = (double*) malloc(sX * sY * sizeof(double));
		DCp1D_na = (double*) malloc(sX * sY * sizeof(double));

#ifdef WRITEINFILE
		myfile_64_1D_non_align << sX << "\t";
	    matmult1D_GF(DAp1D_na,DBp1D_na,DCp1D_na, sX,sY, &myfile_64_1D_non_align);
#else
	    matmult1D(DAp1D_na,DBp1D_na,DCp1D_na, sX,sY);
#endif

		// Deallocate memory ;-)
		free(DAp1D_na);
		free(DBp1D_na);
		free(DCp1D_na);
		// ---------------------------------------------------


		// ---------------------------------------------------
		// Allocate memory for matrices


		// ---------------------------------------------------

		if (printres)
			cout << "Double precision (64 bits) MXM (NON ALIGNED PTR in 1D) GEMM version" << endl;

		DAp1D_na = (double*) malloc(sX * sY * sizeof(double));
		DBp1D_na = (double*) malloc(sX * sY * sizeof(double));
		DCp1D_na = (double*) malloc(sX * sY * sizeof(double));

		MXM::initMat1D(DAp1D_na,(double)4.0, sX, sY);
		MXM::initMat1D(DBp1D_na,(double)7.0, sX, sY);
		MXM::initMat1D(DCp1D_na,(double)0.0, sX, sY);

		gettimeofday(&tv1, (struct timezone*)0);
	#ifdef INTEL
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,sY,sX,sY,(double)1, DAp1D_na,max(sX, sY),DBp1D_na,max(sX, sY),(double)0, DCp1D_na,max(sX, sY));
	#elif BLUEGENE
		dgemm("n", "n", sY,sX,sY, (double)1,DAp1D_na,max(sX, sY),DBp1D_na,max(sX, sY), (double)0, DCp1D_na,max(sX, sY));
	#endif
		gettimeofday(&tv2, (struct timezone*)0);

		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;

		flops=(double)sX*(double)sY*(double)sY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));

		if (printres)
			cout << "GEMM  \t";

		Dver = (double)0;
		for (i=0;i<sX;i++){
				for (j=0;j<sY;j++){
						Dver = Dver+DCp1D_na[i*sY+j];
				}
		}

		if (printres)
			cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  Dver << endl;


		// Deallocate memory ;-)
		free(DAp1D_na);
		free(DBp1D_na);
		free(DCp1D_na);
		// ---------------------------------------------------









/*


		// ------------------------------------------------------------------------------------------------------------------------------------------
		// ------------------------------------------------------------------------------------------------------------------------------------------
		// LONG DOUBLE VERSION
		// ------------------------------------------------------------------------------------------------------------------------------------------
		// ------------------------------------------------------------------------------------------------------------------------------------------

	#ifdef INTEL

		if (printres){
			cout << "Long double precision (128 bits) MXM (STATIC)" << endl;
			cout << "Order\t   GF/s\t\t time\t\t FLOPS\t\t verif " << endl;
		}

		for (i=0;i<sX;i++){
			for (j=0;j<sY;j++){
				DDA[i][j] = (long double)4.0;
				DDB[i][j] = (long double)7.0;
				DDC[i][j] = (long double)0.0;
			}
		}

		gettimeofday(&tv1, (struct timezone*)0);

		for (i=0;i<sX;i++){
			for (j=0;j<sY;j++){
				for (k=0;k<sY;k++){
					DDC[i][j]=DDC[i][j]+DDA[i][k]*DDB[k][j];
				}
			}
		}

		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;

		flops=(double)sX*(double)sY*(double)sY*2.0;
		gflops=flops/(elapsedtime*pow((double)10,9));

		if (printres)
			cout << "Static\t";

		DDver = (long double)0;
		for (i=0;i<sX;i++){
			for (j=0;j<sY;j++){
				DDver = DDver+DDC[i][j];
			}
		}

		if (printres)
			cout << "   "<< gflops <<"\t"  << setprecision(6) << elapsedtime << "\t" << setprecision(6) << flops << "\t" <<setprecision(6) <<  DDver << endl;



		if (printres)
			cout << "Long double precision (128 bits) MXM (ALIGNED PTR 2D)" << endl;


		DDAp = (long double**) _aligned_malloc(sX * sizeof(long double*), 16);
		DDBp = (long double**) _aligned_malloc(sX * sizeof(long double*), 16);
		DDCp = (long double**) _aligned_malloc(sX * sizeof(long double*), 16);

		for(i=0; i<sX; i++){
			DDAp[i] = (long double*) _aligned_malloc(sY * sizeof(long double), 16);
			DDBp[i] = (long double*) _aligned_malloc(sY * sizeof(long double), 16);
			DDCp[i] = (long double*) _aligned_malloc(sY * sizeof(long double), 16);
		}

		matmult2D(DDAp,DDBp,DDCp,sX,sY);


		// Deallocate memory ;-)
		for(i=0; i<sX; i++){
			_aligned_free(DDAp[i]);
			_aligned_free(DDBp[i]);
			_aligned_free(DDCp[i]);
		}
		_aligned_free(DDAp);
		_aligned_free(DDBp);
		_aligned_free(DDCp);
	// ---------------------------------------------------


		// ---------------------------------------------------
		// Allocate memory for matrices

		if (printres)
			cout << "Long double precision (128 bits) MXM (ALIGNED PTR in 1D)" << endl;

		DDAp1D = (long double*) _aligned_malloc(sX * sY * sizeof(long double), 16);
		DDBp1D = (long double*) _aligned_malloc(sX * sY * sizeof(long double), 16);
		DDCp1D = (long double*) _aligned_malloc(sX * sY * sizeof(long double), 16);

		matmult1D(DDAp1D,DDBp1D,DDCp1D, sX,sY);

		// Deallocate memory ;-)
			_aligned_free(DDAp1D);
			_aligned_free(DDBp1D);
			_aligned_free(DDCp1D);
		// ---------------------------------------------------


	#elif BLUEGENE
		if (printres)
			cout << "Long double version is not supported on the BlueGene/P machine" << endl;
	#endif


*/


    }



#ifdef WRITEINFILE
    myfile_32_static.close();
    myfile_32_2D_align.close();
    myfile_32_2D_non_align.close();
    myfile_32_1D_align.close();
    myfile_32_1D_non_align.close();
    myfile_32_dgemm.close();

    myfile_64_static.close();
    myfile_64_2D_align.close();
    myfile_64_2D_non_align.close();
    myfile_64_1D_align.close();
    myfile_64_1D_non_align.close();
    myfile_64_dgemm.close();
	cout << "Close file" << endl;
#endif




}

