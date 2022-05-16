#include "general.h"
#include <stdio.h>

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
		printf("Relevant elements %d %f %d\n",k0,a[k0],amax);
		if(a[k0]>amax){
			amax = a[irow];
			imax = irow;
		}
		printf("Find max %d %d %d \n",irow,imax,amax);
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
	printNSMatrix(aug,nrow,ncol);
	for (int i=0;i<nrow-1;i++){
		int rowMax;
		rowMax = partialPivoting(aug,nrow,ncol,i);
		printf("****** %d\n",rowMax);
		fact = 1;
		kdiag = i * ncol + i;
		if(aug[kdiag] == 0){
			printf("Invalid Matrix!\n");
			return 0;}
		for(int j=i+1;j<nrow;j++){
			k1 = j * ncol + i;
			fact =double(aug[k1])/double(aug[kdiag]);
			printf("%d %d %d %f %f %f\n",j,i,k1,fact,aug[k1],aug[kdiag]);
			for(int k = 0; k<ncol;k++){
				k2 = j*ncol + k;
				k3 = i*ncol + k;
				aug[k2] -= fact*aug[k3];
			}
			printf("\n");
			printNSMatrix(aug,n,n+1);
		}
	}
////////////////////////////////////


// Backward substitution.

// ////////////////////
	delete [] aug;
	return x;
}
