#include "general.h"
#include <stdio.h>
#include"mpi.h"
void gaussElimination2DBlockCyclic(int nprocs, int size){
	int dims[2] = {0,0};
	int iam;
	int blockSize = 10; // This is the block size of a single block.
	MPI_Comm_size(MPI_COMM_WORLD, &iam);
	MPI_Dims_create(nprocs, 2,dims); // Let the MPI decide number of virtual blocks.
	printf("I am %d and dimensions are: row %d and col %d", iam, dims[0], dims[1]);
	int rowBlockSize = blockSize * dims[0]; // The number of rows on a single processor
	int colBlockSize = blockSize * dims[1]; // The number of cols on a single processor
	int nRowBlocks = 0, nColBlocks = 0;

	if (size % rowBlockSize == 0){
		nRowBlocks = size / rowBlockSize; // Number of "ROW" blocks 
		printf("Number of row blocks: %d\n", nRowBlocks);                       
	}else{
		printf("Row size is not divisible by blocksize\n");                     
	} 

	if (size % colBlockSize == 0){
		nColBlocks = size / colBlockSize; // Number of "COL" blocks 
		printf("Number of col blocks: %d\n", nRowBlocks);                       
	}else{
		printf("Row size is not divisible by blocksize\n");                     
	} 

    int nLocalSize = blockSize * nRowBlocks * blockSize * nColBlocks;
	// ELIMINATION STAGE.
	// LOOP OVER GLOBAL INDEX.
	for(int iGlob = 0; iGlob < size; iGlob++){
		int rootProc;// Need a method to find root processor.
		int rowCoord;// Row coordinate of the root processor.
//		if (coord[0] == rowCoord && ){
			// broadcast
			// Find my global minimum col index.
			// broadcast my diagonal element with its local row index
			// to the row comm &
//		}
		// - FIND THE GLOBAL I & J
		// - CHECK IF THE ELEMENT IS ON DIAGONAL. 
		// - WHICH PROC DOES THE DIAG BELONG??
		// - THEN DIVIDE THE LOCAL ROW WITH THAT ELEMENT
		// - ROW BROADCAST THAT ELEMENT AND DIVIDE THE SAME ROW
		//   ON ALL ROW PROCS.
		// - COLUMN BROADCAST THE DIVIDED ROW AND SUBTRACT FROM 
		//   ALL ROWS !!BELOW!!.
	}
}
double partialPivoting(double *a, int nrow, int ncol, int icol){
// find the maximum element and its index in icol
	int imax,amax;
// How to replace rows without cache misses.
	for (int irow=icol;irow<nrow;irow++){
		if(irow==icol){
			amax = a[irow*ncol+icol];
			imax = 0;
		}
		int k0 =irow*ncol+icol;
//		printf("Relevant elements %d %f %d\n",k0,a[k0],amax);
		if(a[k0]>amax){
			amax = a[irow];
			imax = irow;
		}
//		printf("Find max %d %d %d \n",irow,imax,amax);
	}
	return imax;
}

double* gaussElimination(double* A, double* b, double* x, int n){
	int ncol = n+1,nrow = n;
	double *aug = new double[ncol*nrow];
	int kdiag = 0, k1 = 0, k2 = 0, k3 = 0;
	double fact = 0.0;
// Create the augmented matrix
	for (int i=0; i<n; i++){
		for (int j=0; j<n+1;j++){
			if(j<n){
				aug[i*ncol + j] = A[i*n + j];
			}else{	
				aug[i*ncol + j] = b[i];
			}	
		}

	}
//

// Forward elimination
//	printNSMatrix(aug,nrow,ncol);
	for (int i=0;i<nrow-1;i++){
		int rowMax;
		rowMax = partialPivoting(aug,nrow,ncol,i);
//		printf("****** %d\n",rowMax);
		fact = 1;
		kdiag = i * ncol + i;
		if(aug[kdiag] == 0){
//			printf("Invalid Matrix!\n");
			return 0;}
		for(int j=i+1;j<nrow;j++){
			k1 = j * ncol + i;
			fact =double(aug[k1])/double(aug[kdiag]);
//			printf("%d %d %d %f %f %f\n",j,i,k1,fact,aug[k1],aug[kdiag]);
			for(int k = 0; k<ncol;k++){
				k2 = j*ncol + k;
				k3 = i*ncol + k;
				aug[k2] -= fact*aug[k3];
			}
//			printf("\n");
//			printNSMatrix(aug,n,n+1);
		}
	}
////////////////////////////////////


// Backward substitution.

// ////////////////////
	delete [] aug;
	return x;
}

void gaussEliminationRowBlockCyclic(int iam,int nprocs, int size){
	double t0,t1;

	// Number of rows per processor
	int nRows = size / nprocs;
	double *ALocal = new double[nRows*size]; 
	double *bLocal = new double[nRows]; 
	double *xLocal = new double[nRows]; 
	double *sendBuf = new double[size + 1]; 

	// Initialize ALocal and bLocal
	//
	//
	// Create data type for sending rows below
	// +1 for the augmentation
	MPI_Datatype rowType;
    MPI_Type_contiguous(size+1, MPI_DOUBLE, &rowType);
    MPI_Type_commit(&rowType);

	//INIT MATRIX
	for (int i=0;i<nRows;++i){
		int globalI = iam + i * nprocs;
	//	if(iam == 0) printf("I am %d and my globalI is %d\n",iam,globalI);
		for (int j=0;j<size;++j){
			int localIndex = i*size + j;
			if(j == globalI){
				ALocal[localIndex] = size;
			}else{
				ALocal[localIndex] = 0.1;
			}
		}
		bLocal[i] = i+1;
	}
	for(int j=0;j<size+1;j++) sendBuf[j] = 0.0;

	for (int i = 0; i<size;i++){ //Global loop
		int localProc = i % nprocs;
		int localRowIndex = i / nprocs;
		int pivotIndex = localRowIndex * size + i;
		if( iam == localProc){
			for (int j = i+1; j<size; j++){
				int local1Dindex = localRowIndex * size + j;
				ALocal[local1Dindex] = ALocal[local1Dindex] / ALocal[pivotIndex]; // Row normalization
			}
			bLocal[localRowIndex] = bLocal[localRowIndex] / ALocal[pivotIndex];
			ALocal[pivotIndex] = 1.0;
			for (int j = 0;j<size;j++){
				int local1Dindex = localRowIndex * size + j;
				sendBuf[j] = ALocal[local1Dindex];
			}
			sendBuf[size] = bLocal[localRowIndex];
		}
//			printf("Iam %d size %d, localRowIndex %d global row index %d\n",iam,size,localRowIndex,i);
		MPI_Bcast(sendBuf, 1, rowType, localProc, MPI_COMM_WORLD);
//		}else{
		printf("***** IAM ******* %d\n ",iam);
		printVector(sendBuf, size+1);
		if(iam<=localProc ){
			for(int j = localRowIndex + 1; j<nRows; j++){
				for(int k=i+1;k<size;k++){
				int local1Dindex = localRowIndex * size + k;
				int pivotIndex = localRowIndex * size + i;
				ALocal[local1Dindex] -= ALocal[pivotIndex] * sendBuf[k];
				// Elimination for the local proc and procs before the local.
				}
			}
		}else if(iam>localProc){
			for(int j = localRowIndex; j<nRows; j++){
				for(int k=i+1;k<size;k++){
				// Elimination for the procs after the local.
				int local1Dindex = localRowIndex * size + k;
				int pivotIndex = localRowIndex * size + i;
				ALocal[local1Dindex] -= ALocal[pivotIndex] * sendBuf[k];
				}
			}
		}
//		if(iam==0)printf("Global index is %d and belongs to proc %d %d\n",i,localProc,localRowIndex);
	}
//	printNSMatrix(ALocal,nRows,size);
//	Back substitution.
}
