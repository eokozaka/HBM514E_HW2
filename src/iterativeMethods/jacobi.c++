#include "general.h"
#include "stdio.h"
#include "math.h"
#include "mpi.h"
void jacobi(int n, int nprocs, int iam){
//	MPI_Barrier(MPI_COMM_WORLD);
//	printf("I am %d in jacobi\n",iam);

	MPI_Request sreq, rreq;
	MPI_Status status;

//  Calculate the number of rows per processor
	int nRowsLocal = n / nprocs;

//  Calculate the global beginning and ending indices for each process.
	int iBeginGlob = iam * nRowsLocal;
	int iEndGlob   = iBeginGlob + (nRowsLocal - 1);

//	printf("I am %d and my global indices are from %d to %d\n",iam, iBeginGlob 
//		   ,iEndGlob);
	   
/***********************************************
* Allocate the a matrix, b vector and x vector *
***********************************************/

 	double *A = new double[nRowsLocal*n];
	double *b = new double[nRowsLocal];
	double *xnew = new double[nRowsLocal];
	double *xold = new double[nRowsLocal];
	double *xres = new double[nRowsLocal];
	double maxRes=0.0,tolerance=1.0;

	for (int i = 0; i< nRowsLocal; i++){
		xold[i] = 0;
	}

	MPI_Datatype xChunk; // Datatype to pass the chunks of vectors x & b.
    MPI_Type_contiguous(nRowsLocal, MPI_DOUBLE, &xChunk);
    MPI_Type_commit(&xChunk);
/* Find the upper and lower neighbours. */
//
	int myup, mydown;

    mydown = (iam + 1) % nprocs;
	if( iam == 0)
    myup = nprocs - 1;
	else 
    myup = (iam - 1) % nprocs;

//	printf("I am %d and my up is %d and my down is %d \n",iam,myup,mydown);

/************************************************/ 

/* Initialize x vector for previous solution locally */
	for (int i = 0; i< nRowsLocal; i++){
		xold[i] = 0; 
	}

/* Initialize A matrix and b vector locally. */
	for (int i = 0; i < nRowsLocal; i++){
		int iGlobRow = iam * nRowsLocal + i;
		b[i] = iGlobRow;
		for (int j = 0; j < n; j++){
			int iGlob = iam * nRowsLocal + i;
			if( iGlob != j ){
				A[ i * n + j ] = 0.500 ;
			}else{
				A[ i * n + j ] = n;

			}	
		}
	}	
	
	int counter = 1;
	int maxIter = 1000;
	double t0 = 0, t1 = 0;
	t0 = MPI_Wtime();
	while(tolerance>1e-04 && counter <= maxIter){
		for (int iproc = 0; iproc < nprocs; iproc++){
			MPI_Isend(xold, 1, xChunk,   myup, 101, MPI_COMM_WORLD, &sreq);
			// Calculation of the sums
			for (int i = 0; i < nRowsLocal;i++){
				int iGlobRow = iam*nRowsLocal + i;
	//			if( iam == 0){
	//				printf("Local Row Number %d, Global Row Number %d\n", i,iGlobRow);
	//			}
				for (int j = 0; j < nRowsLocal; j++){
	
					int iLocalA  = i*n + ((iproc + iam)%nprocs)*nRowsLocal + j;
					int	iGlobA   = iam * nRowsLocal * n + iLocalA;
					int iGlobCol = iLocalA % n;
					
					if( iGlobRow != iGlobCol)		
					xnew[i] -= A[iLocalA]*xold[j];
				}
			}
	//		//Receive a chunk of x vector from bottom.
			MPI_Recv(xold, 1, xChunk, mydown, 101, MPI_COMM_WORLD, &status);
			MPI_Wait(&sreq, &status);
		}//end shift loop

		for (int i = 0; i<nRowsLocal;i++){
			int iLocalDiag = iam*nRowsLocal + i*(n+1);

			xnew[i] = (b[i] + xnew[i]) / A[iLocalDiag];
			xres[i] = abs(xnew[i] - xold[i]);

			if (xres[i]>maxRes) maxRes = xres[i];
			xold[i] = xnew[i];
		}//end calculate X, local max residual and swap

		//find the maximum tolerance in between processors.
		MPI_Allreduce(&maxRes, &tolerance, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
//		if(iam == 0) printf("Iteration %d Residual: %f \n",counter,tolerance);
		counter++;
		maxRes = 0;
	}// end while loop

	t1 = MPI_Wtime();

	printVector(xnew, nRowsLocal);
    MPI_Barrier(MPI_COMM_WORLD); 
	MPI_Type_free(&xChunk); 
	if(iam == 0){
		printf("Parallel Time: %f\n",t1-t0);
	}
 	delete [] A; 
	delete [] b; 
	delete [] xnew;
	delete [] xold;
	delete [] xres;
}
