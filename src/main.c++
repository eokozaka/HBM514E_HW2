#include"stdio.h"
#include"math.h"
#include"iostream"
#include"mainConfig.h"
#include <random>
#include "mpi.h"
#ifdef USE_MYSORT
#include "sorting.h"
#endif
#include "general.h"
#include "directMethods.h"
double* matMul(double* a,double* b, double* c, int n){
	for (int i=0; i<n;i++){
		for (int j=0; j<n;j++){
			for (int k=0; k<n;k++){
				c[ i * n + j] += a[ i * n + k] * b[ k*n + j ];  
			}
		}
	}
	return c;
}
double* matVec(double* a,double* b, double* c, int n){
	for (int i=0; i<n;i++){
		for (int j=0; j<n;j++){
				c[i] += a[ i * n + j] * b[j];  
		}
	}
	return c;
}



void rowPointers(double* a, int n){

//	for(int i=0;i<n;i++){
//		std::cout << "Row number " << i << " Address: "<< &a[i*n] << " Value "
//			<< a[i*n]<<"\n";
//	}

}

int* transpose(int* a, int* at, int n){
	return at;
}

int main(int argc, char* argv[]){
	int n;
	n = atoi(argv[1]);
	int nprocs, iam;
	MPI_Request sreq, rreq;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &iam);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	if(iam == 0) std::cout << " Version " << VERSION_MAJOR << "." << VERSION_MINOR << "\n";

//  Calculate the number of rows per processor
	int nRowsLocal = n / nprocs;
//  Calculate the local beginning and ending indices for each process.
	int iBeginGlob = iam * nRowsLocal;
	int iEndGlob = iBeginGlob + (nRowsLocal - 1);
//	printf("I am %d and my global indices are from %d to %d\n",iam, iBeginGlob 
//		   ,iEndGlob);
/*	***********************************************
 *	Allocate the a matrix, b vector and x vector  *
**************************************************/
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

	int myup, mydown;

	
    mydown = (iam + 1) % nprocs;
	if( iam == 0)
    myup = nprocs - 1;
	else 
    myup = (iam - 1) % nprocs;
//	printf("I am %d and my up is %d and my down is %d \n",iam,myup,mydown);
//
	for (int i = 0; i< nRowsLocal; i++){
		xold[i] = 0; //(iam + 1) * i;
		std::cout << "\r" << iam << "** " << i << " " << xold[i] << "\n";
	}

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

//	printNSMatrix(A,nRowsLocal,n);


//	for (int i = 0; i< nRowsLocal; i++){
//		printf("I am %d and my xold[%d] = %f\n",iam, i, xold[i]);
//	}


//	for (int i = 0; i< nRowsLocal; i++){
//		printf("I am %d and after comms my xold[%d] = %f\n",iam, i, xold[i]);
//	}
	
	int counter = 1;
	int maxIter = 1000;
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
//					int iLocalA  = i*n + iproc*nRowsLocal + j;
					int	iGlobA   = iam * nRowsLocal * n + iLocalA;
					int iGlobCol = iLocalA % n;
	
	//				if( iam == 0 && iGlobRow != iGlobCol){
	//					printf("%d,%d,%d,%d ",iLocalA, iGlobA,iGlobRow,iGlobCol); // Local index of A
	//				}
					
					if( iGlobRow != iGlobCol)		
			//		xnew[i] -= A[i*n+iproc*nRowsLocal+j]*xold[j];
					xnew[i] -= A[iLocalA]*xold[j];
				}
			//		if( iam == 0) printf("\n");
			}
	//		//Receive a chunk of x vector from bottom.
			MPI_Recv(xold, 1, xChunk, mydown, 101, MPI_COMM_WORLD, &status);
//			MPI_Irecv(xold, 1, xChunk, mydown, 101, MPI_COMM_WORLD, &rreq);
	//	
			MPI_Wait(&sreq, &status);
//			MPI_Wait(&rreq, &status);
		}//end shift loop
//		MPI_Isend(xold, 1, xChunk,   myup, 101, MPI_COMM_WORLD, &sreq);
//		MPI_Irecv(xold, 1, xChunk, mydown, 101, MPI_COMM_WORLD, &rreq);
//	//	
//		MPI_Wait(&sreq, &status);
//		MPI_Wait(&rreq, &status);


		for (int i = 0; i<nRowsLocal;i++){
			int iLocalDiag = iam*nRowsLocal + i*(n+1);
	//		printf("I am %d, local row is %d and the diag is %f\n",iam,i,A[iLocalDiag]);
			xnew[i] = (b[i] + xnew[i]) / A[iLocalDiag];
			xres[i] = abs(xnew[i] - xold[i]);
			std::cout << "\r" << iam << 
				" index: " << i << 
				" counter: " << counter << 
				" b: "    << b[i] << 
				" diag: " <<  A[iLocalDiag]	<< 
				" xold[i] " << xold[i] << 
				" xnew[i] " << xnew[i] << 
				" tol " << tolerance << 
				" xres[i] " << xres[i] <<  "\n";
			if (xres[i]>maxRes) maxRes = xres[i];
			xold[i] = xnew[i];
		}//end calculate X, local max residual and swap

		//find the maximum tolerance in between processors.
		MPI_Allreduce(&maxRes, &tolerance, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
//		if(iam == 0) printf("Iteration %d Residual: %f \n",counter,tolerance);
		counter++;
		maxRes = 0;
	}// end while loop

	printVector(xnew, nRowsLocal);

 	delete [] A; 
	delete [] b; 
	delete [] xnew;
	delete [] xold;
	delete [] xres;
	MPI_Finalize();
}
