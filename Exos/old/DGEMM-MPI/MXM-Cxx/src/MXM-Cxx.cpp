//============================================================================
// Name        : MXM-Cxx.cpp
// Author      : Vince
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
// Add
//
// -align -axSSE4.2 -vec-report3 -ansi-alias-check -ansi-alias -falias -O3 -ipo
//
// in the Makefile to get all opt
//============================================================================

#include <stdio.h>
#include <stdlib.h>

#include <mxm.hpp>
#include <mxv.hpp>
#include <mpi.h>

int main(int argc, char *argv[]) {

	int choix=1;

	int myrank, mysize;


	MPI_Init (&argc, &argv);	/* starts MPI */
	MPI_Comm_rank (MPI_COMM_WORLD, &myrank);	/* get current process id */
	MPI_Comm_size (MPI_COMM_WORLD, &mysize);	/* get number of processes */

	if (myrank==0) {

		cout << "MYRANK IS " << myrank << " and we are " << mysize << endl;

		cout << "Ciao. Bienvenue sur le ch'tit bench de Vince" << endl;
		cout << " 1 : MXM uniquement" << endl;
		cout << " 2 : MXV uniquement" << endl;
		cout << " 3 : MXM  & MXV" << endl;
		cout << " Faites votre choix: " << endl;
		cin >> choix ;

		switch (choix){

			case 1: {
				MXM *mxm = new MXM();
				mxm->compute();
				free(mxm);
				break;
			}
			case 2: {
				MXV *mxv = new MXV();
				mxv->compute();
				free(mxv);
				break;
			}
			case 3: {
				MXM *mxm = new MXM();
				mxm->compute();
				free(mxm);
				MXV *mxv = new MXV();
				mxv->compute();
				free(mxm);
				break;
			}
			default: {
				cout << "Hé ho, j'ai dit 1, 2 ou 3... Rrogntudjû ! " << endl;
			}

		}
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
