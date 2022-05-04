#include"stdio.h"
#include"math.h"
#include"iostream"
#include"mainConfig.h"
#include <random>
#ifdef USE_MYSORT
#include "sorting.h"
#endif
void printMatrix(int* a,int n){
	int k;
	for (int i=0; i<n;i++){
		for (int j=0; j<n;j++){
			k = i * n + j;
			printf("%d\t ",a[k]);
		}
		printf("\n");
	}
	printf("\n");
}
void printNSMatrix(int *a,int m, int n){
	int k;
	for (int i=0; i<m;i++){
		for (int j=0; j<n;j++){
			k = i * n + j;
			printf("%d\t ",a[k]);
		}
		printf("\n");
	}
	printf("\n");
}
void printVector(int* a,int n){
	for (int i=0; i<n;i++){
			printf("%d\t ",a[i]);
		}
		printf("\n");
		printf("\n");
}
int* init3BSymMatrix(int* a,int n){
	int k;
	for (int i=0; i<n;i++){
		for (int j=0; j<n;j++){
			k = i * n + j;
			if ( i == j ){
				a[k] = i + 1;
			}else if( i == j-1){
				a[k] = -2;
			}else if( i == j+1){
				a[k] = -2;
			}else{
				a[k] = 0;
			}
		}
	}
	return a;
}
int* initEyeMatrix(int* a,int n){
	int k;
	for (int i=0; i<n;i++){
		for (int j=0; j<n;j++){
			k = i * n + j;
			if ( i == j ){
				a[k] = 1;
			}else{
				a[k] = 0;
			}
		}
	}
	return a;
}
int* initRandomFullMatrix(int* a,int n){
	int k;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<>distr(1,3);


	for (int i=0; i<n;i++){
		for (int j=0; j<n;j++){
			k = i * n + j;
			a[k] = distr(gen);
		}
	}
	return a;
}
int* initRandomVector(int* a,int n){
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<>distr(1,3);


	for (int i=0; i<n;i++){
			a[i] = distr(gen);
	}
	return a;
}
int* matMul(int* a,int* b, int* c, int n){
	for (int i=0; i<n;i++){
		for (int j=0; j<n;j++){
			for (int k=0; k<n;k++){
				c[ i * n + j] += a[ i * n + k] * b[ k*n + j ];  
			}
		}
	}
	return c;
}
int* matVec(int* a,int* b, int* c, int n){
	for (int i=0; i<n;i++){
		for (int j=0; j<n;j++){
				c[i] += a[ i * n + j] * b[j];  
		}
	}
	return c;
}

int* gaussElimination(int* A, int* b, int* x, int n){
	int ncol = n+1,nrow = n;
	int *aug = new int[ncol*nrow];
	int  kdiag = 0, k1 = 0, k2 = 0, k3 = 0;
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
		fact = 1;
		kdiag = i * ncol + i;
		if(aug[kdiag] == 0){
			printf("Invalid Matrix!\n");
			return 0;}
		for(int j=i+1;j<nrow;j++){
			k1 = j * ncol + i;
			fact =double(aug[k1])/double(aug[kdiag]);
			printf("%d %d %d %f %d %d\n",j,i,k1,fact,aug[k1],aug[kdiag]);
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

int pivoting(int *a, int nrow, int ncol, int icol){
// find the maximum element and its index in icol
	int imax,amax;
// How to replace rows without cache misses.
	amax = a[icol];
	imax = 0;
	for (int irow=1;irow<nrow;irow++){
		int k0 =nrow*irow+icol;
		printf("Relevant elements %d %d %d\n",k0,a[k0],amax);
		if(a[k0]>amax){
			amax = a[irow];
			imax = irow;
		}
		printf("Find max %d %d %d \n",irow,imax,amax);
	}
	return imax;
}

void rowPointers(int* a, int n){

	for(int i=0;i<n;i++){
		std::cout << "Row number " << i << " Address: "<< &a[i*n] << " Value "
			<< a[i*n]<<"\n";
	}

}

int* transpose(int* a, int* at, int n){
	return at;
}

int main(int argc, char* argv[]){
	int n = 5;
	int *a = new int[n*n];
	int *b = new int[n*1];
	int *c = new int[n*1];

	std::cout << " Matrix Multiplication Version " << VERSION_MAJOR << "." << VERSION_MINOR << "\n";

	initRandomFullMatrix(a,n);
	printMatrix(a,n);
	rowPointers(a,n);
#ifdef USE_MYSORT
	sort(b, n);
#endif

//	max = pivoting(a,n,n-2);
//	printf("*** %d ***\n",max);
//	initRandomVector(b,n);
//
//	printMatrix(a,n);
//	printVector(b,n);
//
//	matVec(a, b, c, n);
//
//	printVector(c,n);
//
//	gaussElimination(a, b, c,n);

	delete [] a;
	delete [] b;
	delete [] c;
}
