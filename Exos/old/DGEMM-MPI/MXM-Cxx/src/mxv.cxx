/*
 * mxv.cxx
 *
 *  Created on: Mar 21, 2011
 *      Author: vkeller
 *
 *
 *      n1 = sizeX
 *      n2 = sizeY
 */

#include <mxv.hpp>

#define SIZE_MAX 100
#define STEP 2

void MXV::compute(){

//	int n1=3999, n2=3999, ndiag=9, nn=(n1+1)*(n2+1), it, ni;
//	int n1=1999, n2=1999, ndiag=9, nn=(n1+1)*(n2+1), it, ni;
	int n1, n2, ndiag=9, nn, it, ni;

#ifdef WRITEINFILE
	bool printres=false;
#else
	bool printres=true;
#endif

	printres=false;

	int i, j, size;

	double **DAMATp ALIGN_16;
	double *DUp ALIGN_16, *DUap  ALIGN_16, *DAMAT1Dp ALIGN_16;
	double Dsum=(double)0, Dres=(double)0;

	float **SAMATp ALIGN_16;
	float *SUp ALIGN_16, *SUap  ALIGN_16, *SAMAT1Dp ALIGN_16;
	float Ssum=(float)0, Sres=(float)0;

	double **DAMATp_na;
	double *DUp_na, *DUap_na, *DAMAT1Dp_na;

	float **SAMATp_na;
	float *SUp_na, *SUap_na, *SAMAT1Dp_na;


//	double s, cst ;

    struct timeval tv1,tv2;
    double elapsedtime;

    double flops=0, gflops=0;

//    int * jcoef;

#ifdef WRITEINFILE

	cout << "I write in a file" << endl;

    ofstream myfile;
    ostringstream outfn;
    outfn << "dataMXV-" << SIZE_MAX << "-" << STEP << ".dat";
	string filename = outfn.str();
	myfile.open (filename.c_str());

	// write header
	myfile << "size\t32_S\t32_2D_A\t32_2D_NA\t32_1D_A\t32_1D_NA\t32_GEMV\t64_S\t64_2D_A\t64_2D_NA\t64_1D_A\t64_1D_NA\t64_GEMV" << endl;
#endif
    for (size=2;size<SIZE_MAX;size++){

    	n1=size*STEP;
    	n2=size*STEP;
    	nn=(n1+1)*(n2+1);


		double DAMAT[nn][ndiag];
		double DU[nn];
		double DUa[nn];

		float SAMAT[nn][ndiag];
		float SU[nn];
		float SUa[nn];

#ifdef WRITEINFILE
		myfile << nn<<"\t";
#endif
	// =================================================================
	// INITIALIZATION

		if (printres)
			cout << "Init matrices vectors and other stuff" << endl;

		for (i=0;i<nn;i++){
			for (j=0;j<ndiag;j++){
				DAMAT[i][j] = (double)1.0;
				SAMAT[i][j] = (float)1.0;
			}
			DU[i]=(double)1.0;
			DUa[i]=(double)0.0;
			SU[i]=(float)1.0;
			SUa[i]=(float)0.0;
		}


		Dres = (double)nn*(double)ndiag;
		Sres = (float)nn*(float)ndiag;

		ni = 10;
		flops = (double)2*(double)nn*(double)ndiag*(double)ni;



		if (printres){
			cout << "=============================================" << endl;
			cout << "================== 32 bits ==================" << endl;
			cout << "=============================================" << endl;
			cout << "Size = " << nn << " (n1 = " << n1 << " ; n2 = " << n2 << ")" << endl;
		}
	// =================================================================

		if (printres)
			cout << "MXV (STATIC) in 32 bits" << endl;

		elapsedtime = (double)0;

		for (it=0;it<ni;it++){
			MXV::initVec(SUa,(float)0,nn);

			gettimeofday(&tv1, (struct timezone*)0);
			for (i=0;i<nn;i++){
				for (j=0;j<ndiag;j++){
					SUa[i] = SUa[i] + SAMAT[i][j]*SU[i];
				}
			}
			gettimeofday(&tv2, (struct timezone*)0);
			elapsedtime+=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		}

		Ssum=MXV::verifyVec(SUa, nn);

		gflops=flops/elapsedtime/1.e9;
		if (printres){
			cout << "End computation" << endl;
			printRes(flops,elapsedtime,gflops);
		}

#ifdef WRITEINFILE
		myfile <<gflops<<"\t";
#endif

		//		if (printres)
		//			MXV::printSum(Sres, Ssum);





	// =================================================================

		if (printres)
			cout << "MXV (PTR ALIGNED 2D _aligned_malloc) in 32 bits" << endl;

		SAMATp = (float**) _aligned_malloc(nn * sizeof(float*), 16);
		SUp = (float*) _aligned_malloc(nn * sizeof(float), 16);
		SUap = (float*) _aligned_malloc(nn * sizeof(float), 16);


		for(i=0; i<nn; i++){
			SAMATp[i] = (float*) _aligned_malloc(ndiag * sizeof(float), 16);
		}
	#ifdef BLUEGENE
		#pragma disjoint (*SUap, *SUp, **SAMATp)
		__alignx(16,SUp);
		__alignx(16,SUap);
		__alignx(16,SAMATp);
	#endif

		MXV::initMat(SAMATp,(float)1.0,nn, ndiag);
		MXV::initVec(SUp,(float)1.0,nn);
		MXV::initVec(SUap,(float)0.0,nn);

#ifdef WRITEINFILE
		Ssum = MXV::matmult2D_GF(nn,ndiag,SUap, SAMATp,SUp, ni);
		myfile <<Ssum<<"\t";
#else
		Ssum = MXV::matmult2D(nn,ndiag,SUap, SAMATp,SUp, ni);
#endif
//		if (printres)
//			MXV::printSum(Sres, Ssum);

		for(i=0; i<nn; i++){
			_aligned_free(SAMATp[i]);
		}
		_aligned_free(SAMATp);
		_aligned_free(SUp);
		_aligned_free(SUap);








	// =================================================================

		if (printres)
			cout << "MXV (PTR NON ALIGNED 2D malloc) in 32 bits" << endl;

		SAMATp_na = (float**) malloc(nn * sizeof(float*));
		SUp_na = (float*) malloc(nn * sizeof(float));
		SUap_na = (float*) malloc(nn * sizeof(float));

		for(i=0; i<nn; i++){
			SAMATp_na[i] = (float*) malloc(ndiag * sizeof(float));
		}

		MXV::initMat(SAMATp_na,(float)1.0,nn, ndiag);
		MXV::initVec(SUp_na,(float)1.0,nn);
		MXV::initVec(SUap_na,(float)0.0,nn);

#ifdef WRITEINFILE
		Ssum = MXV::matmult2D_GF(nn,ndiag,SUap_na, SAMATp_na,SUp_na, ni);
		myfile <<Ssum<<"\t";
#else
		Ssum = MXV::matmult2D(nn,ndiag,SUap_na, SAMATp_na,SUp_na, ni);
#endif
//		if (printres)
//			MXV::printSum(Sres, Ssum);

		for(i=0; i<nn; i++){
			free(SAMATp_na[i]);
		}
		free(SAMATp_na);
		free(SUp_na);
		free(SUap_na);


		// =================================================================

		if (printres)
			cout << "MXV (PTR ALIGNED 1D _aligned_malloc) in 32 bits " << endl;

		SAMAT1Dp = (float*) _aligned_malloc(nn * ndiag * sizeof(float), 16);
		SUp = (float*) _aligned_malloc(nn * sizeof(float), 16);
		SUap = (float*) _aligned_malloc(nn * sizeof(float), 16);

		#ifdef BLUEGENE
		#pragma disjoint (*SUap, *SUp, *SAMAT1Dp)
		__alignx(16,SUp);
		__alignx(16,SUap);
		__alignx(16,SAMAT1Dp);
		#endif


		MXV::initVec(SAMAT1Dp,(float)1.0,(nn*ndiag));
		MXV::initVec(SUp,(float)1.0,nn);
		MXV::initVec(SUap,(float)0.0,nn);

#ifdef WRITEINFILE
		Ssum=MXV::matmult1D_GF(nn,ndiag,SUap, SAMAT1Dp,SUp, ni);
		myfile <<Ssum<<"\t";
#else
		Ssum=MXV::matmult1D(nn,ndiag,SUap, SAMAT1Dp,SUp, ni);
#endif
//		if (printres)
//			MXV::printSum(Sres, Ssum);

		_aligned_free(SAMAT1Dp);
		_aligned_free(SUp);
		_aligned_free(SUap);


		// =================================================================

		if (printres)
			cout << "MXV (PTR NON ALIGNED 1D malloc) in 64 bits " << endl;

		SAMAT1Dp_na = (float*) malloc(nn * ndiag * sizeof(float));
		SUp_na = (float*) malloc(nn * sizeof(float));
		SUap_na = (float*) malloc(nn * sizeof(float));

		MXV::initVec(SAMAT1Dp_na,(float)1.0,(nn*ndiag));
		MXV::initVec(SUp_na,(float)1.0,nn);
		MXV::initVec(SUap_na,(float)0.0,nn);

#ifdef WRITEINFILE
		Ssum=MXV::matmult1D_GF(nn,ndiag,SUap_na, SAMAT1Dp_na,SUp_na, ni);
		myfile <<Ssum<<"\t";
#else
		Ssum=MXV::matmult1D(nn,ndiag,SUap_na, SAMAT1Dp_na,SUp_na, ni);
#endif
//		if (printres)
//			MXV::printSum(Sres, Ssum);

		free(SAMAT1Dp_na);
		free(SUp_na);
		free(SUap_na);







	// =================================================================

		if (printres)
			cout << "MXV (PTR ALIGNED 1D malloc) in 32 bits GEMV" << endl;

		SAMAT1Dp = (float*) malloc(nn*ndiag * sizeof(float));
		SUp = (float*) malloc(nn * sizeof(float));
		SUap = (float*) malloc(nn * sizeof(float));

		MXV::initVec(SAMAT1Dp,(float)1.0,(nn*ndiag));
		MXV::initVec(SUp,(float)1.0,nn);
		MXV::initVec(SUap,(float)0.0,nn);

		elapsedtime = (double)0;

		gettimeofday(&tv1, (struct timezone*)0);
		for (it=0;it<ni;it++){
			MXV::initVec(SUap,(float)0,nn);
	#ifdef INTEL
			cblas_sgemv(CblasRowMajor,CblasNoTrans, nn, ndiag, 1, SAMAT1Dp, ndiag, SUp, 1, 0, SUap, 1);
	#elif BLUEGENE
			sgemv("n", nn, ndiag, 1, SAMAT1Dp, max(1,nn), SUp, 1, 0, SUap, 1);
	#endif
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;

		Ssum=MXV::verifyVec(SUa, nn);

		gflops=flops/elapsedtime/1.e9;
		if (printres)
			printRes(flops,elapsedtime,gflops);

#ifdef WRITEINFILE
		myfile <<gflops<<"\t";
#endif
//		if (printres)
//			MXV::printSum(Sres, Ssum);

		free(SAMAT1Dp);
		free(SUp);
		free(SUap);




		if (printres){
			cout << "=============================================" << endl;
			cout << "================== 64 bits ==================" << endl;
			cout << "=============================================" << endl;
			cout << "Size = " << nn << " (n1 = " << n1 << " ; n2 = " << n2 << ")" << endl;
		}

	// =================================================================

		if (printres)
			cout << "MXV (STATIC 2D) in 64 bits" << endl;

		elapsedtime=(double)0;

		for (it=0;it<ni;it++){
			MXV::initVec(DUa,(double)0,nn);

			gettimeofday(&tv1, (struct timezone*)0);
			for (i=0;i<nn;i++){
				for (j=0;j<ndiag;j++){
					DUa[i] = DUa[i] + DAMAT[i][j]*DU[i];
				}
			}
			gettimeofday(&tv2, (struct timezone*)0);
			elapsedtime+=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;
		}

		gflops=flops/elapsedtime/1.e9;
		if (printres)
			printRes(flops,elapsedtime,gflops);

		Dsum=MXV::verifyVec(DUa, nn);
#ifdef WRITEINFILE
		myfile <<gflops<<"\t";
#endif
//		if (printres)
//			MXV::printSum(Dres, Dsum);



	// =================================================================

		if (printres)
			cout << "MXV (PTR ALIGNED 2D _aligned_malloc) in 64 bits " << endl;

		DAMATp = (double**) _aligned_malloc(nn * sizeof(double*), 16);
		DUp = (double*) _aligned_malloc(nn * sizeof(double), 16);
		DUap = (double*) _aligned_malloc(nn * sizeof(double), 16);

		for(i=0; i<nn; i++){
			DAMATp[i] = (double*) _aligned_malloc(ndiag * sizeof(double), 16);
		}


	#ifdef BLUEGENE
		#pragma disjoint (*DUap, *DUp, **DAMATp)
		__alignx(16,DUp);
		__alignx(16,DUap);
		__alignx(16,DAMATp);
	#endif


		MXV::initMat(DAMATp,(double)1.0,nn, ndiag);
		MXV::initVec(DUp,(double)1.0,nn);
		MXV::initVec(DUap,(double)0.0,nn);

#ifdef WRITEINFILE
		Dsum=MXV::matmult2D_GF(nn,ndiag,DUap, DAMATp,DUp, ni);
		myfile <<Dsum<<"\t";
#else
		Dsum=MXV::matmult2D(nn,ndiag,DUap, DAMATp,DUp, ni);
#endif
//		if (printres)
//			MXV::printSum(Dres, Dsum);

		for(i=0; i<nn; i++){
			_aligned_free(DAMATp[i]);
		}
		_aligned_free(DAMATp);
		_aligned_free(DUp);
		_aligned_free(DUap);




	// =================================================================

		if (printres)
			cout << "MXV (PTR NON ALIGNED 2D malloc) in 64 bits" << endl;

		DAMATp_na = (double**) malloc(nn * sizeof(double*));
		DUp_na = (double*) malloc(nn* sizeof(double));
		DUap_na = (double*) malloc(nn * sizeof(double));

		for(i=0; i<nn; i++){
			DAMATp_na[i] = (double*) malloc(ndiag * sizeof(double));
		}
		MXV::initMat(DAMATp_na,(double)1.0,nn, ndiag);
		MXV::initVec(DUp_na,(double)1.0,nn);
		MXV::initVec(DUap_na,(double)0.0,nn);

#ifdef WRITEINFILE
		Dsum=MXV::matmult2D_GF(nn,ndiag,DUap_na, DAMATp_na,DUp_na, ni);
		myfile <<Dsum<<"\t";
#else
		Dsum=MXV::matmult2D(nn,ndiag,DUap_na, DAMATp_na,DUp_na, ni);
#endif
//		if (printres)
//			MXV::printSum(Dres, Dsum);

		for(i=0; i<nn; i++){
			free(DAMATp_na[i]);
		}
		free(DAMATp_na);
		free(DUp_na);
		free(DUap_na);


		// =================================================================

		if (printres)
			cout << "MXV (PTR ALIGNED 1D _aligned_malloc) in 64 bits " << endl;

		DAMAT1Dp = (double*) _aligned_malloc(nn * ndiag * sizeof(double), 16);
		DUp = (double*) _aligned_malloc(nn * sizeof(double), 16);
		DUap = (double*) _aligned_malloc(nn * sizeof(double), 16);

		#ifdef BLUEGENE
		#pragma disjoint (*DUap, *DUp, **DAMATp)
		__alignx(16,DUp);
		__alignx(16,DUap);
		__alignx(16,DAMATp);
		#endif


		MXV::initVec(DAMAT1Dp,(double)1.0,(nn*ndiag));
		MXV::initVec(DUp,(double)1.0,nn);
		MXV::initVec(DUap,(double)0.0,nn);

#ifdef WRITEINFILE
		Dsum=MXV::matmult1D_GF(nn,ndiag,DUap, DAMAT1Dp,DUp, ni);
		myfile <<Dsum<<"\t";
#else
		Dsum=MXV::matmult1D(nn,ndiag,DUap, DAMAT1Dp,DUp, ni);
#endif
//		if (printres)
//			MXV::printSum(Dres, Dsum);

		_aligned_free(DAMAT1Dp);
		_aligned_free(DUp);
		_aligned_free(DUap);


		// =================================================================

		if (printres)
			cout << "MXV (PTR NON ALIGNED 1D malloc) in 64 bits " << endl;

		DAMAT1Dp_na = (double*) malloc(nn * ndiag * sizeof(double));
		DUp_na = (double*) malloc(nn * sizeof(double));
		DUap_na = (double*) malloc(nn * sizeof(double));

		MXV::initVec(DAMAT1Dp_na,(double)1.0,(nn*ndiag));
		MXV::initVec(DUp_na,(double)1.0,nn);
		MXV::initVec(DUap_na,(double)0.0,nn);

#ifdef WRITEINFILE
		Dsum=MXV::matmult1D_GF(nn,ndiag,DUap_na, DAMAT1Dp_na,DUp_na, ni);
		myfile <<Dsum<<"\t";
#else
		Dsum=MXV::matmult1D(nn,ndiag,DUap_na, DAMAT1Dp_na,DUp_na, ni);
#endif
//		if (printres)
//			MXV::printSum(Dres, Dsum);

		free(DAMAT1Dp_na);
		free(DUp_na);
		free(DUap_na);

	// =================================================================

		if (printres)
			cout << "MXV (PTR ALIGNED 1D malloc) in 64 bits GEMV" << endl;

		DAMAT1Dp = (double*) malloc(nn*ndiag * sizeof(double));
		DUp = (double*) malloc(nn * sizeof(double));
		DUap = (double*) malloc(nn * sizeof(double));

		MXV::initVec(DAMAT1Dp,(double)1.0,(nn*ndiag));
		MXV::initVec(DUp,(double)1.0,nn);
		MXV::initVec(DUap,(double)0.0,nn);

		elapsedtime=(double)0;

		gettimeofday(&tv1, (struct timezone*)0);
		for (it=0;it<ni;it++){
			MXV::initVec(DUap,(double)0,nn);
	#ifdef INTEL
			cblas_dgemv(CblasRowMajor,CblasNoTrans, nn, ndiag, 1, DAMAT1Dp, ndiag, DUp, 1, 0, DUap, 1);
	#elif BLUEGENE
			dgemv("n", nn, ndiag, 1, DAMAT1Dp, max(1,nn), DUp, 1, 0, DUap, 1);
	#endif
		}
		gettimeofday(&tv2, (struct timezone*)0);
		elapsedtime=((tv2.tv_sec-tv1.tv_sec)*1000000.0 + (tv2.tv_usec-tv1.tv_usec))/1000000.0;

		Dsum=MXV::verifyVec(DUa, nn);

		gflops=flops/elapsedtime/1.e9;
		if (printres)
			printRes(flops,elapsedtime,gflops);

#ifdef WRITEINFILE
		myfile <<gflops << endl;
#endif
//		if (printres)
//			MXV::printSum(Dres, Dsum);

		free(DAMAT1Dp);
		free(DUp);
		free(DUap);

#ifdef WRITEINFILE
		cout << "Size = " << nn << " (n1 = " << n1 << " ; n2 = " << n2 << ")" << endl;
#endif

    }

#ifdef WRITEINFILE
	myfile.close();
#endif
//	cout << "General memory free" << endl;

}




void MXV::printRes(double flops, double elapsedtime, double gflops){
	cout << "Flops           = " << flops << endl;
	cout << "Time            = " << elapsedtime << endl;
	cout << "GF/s            = " << gflops << endl;
}
