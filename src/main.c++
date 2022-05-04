#include"stdio.h"
#include"math.h"
#include"iostream"
#include"mainConfig.h"
#include <random>
#ifdef USE_MYSORT
#include "sorting.h"
#endif
void printMatrix(double* a,int n){
	int k;
	for (int i=0; i<n;i++){
		for (int j=0; j<n;j++){
			k = i * n + j;
			printf("%f\t ",a[k]);
		}
		printf("\n");
	}
	printf("\n");
}
void printNSMatrix(double *a,int m, int n){
	int k;
	for (int i=0; i<m;i++){
		for (int j=0; j<n;j++){
			k = i * n + j;
			printf("%f\t ",a[k]);
		}
		printf("\n");
	}
	printf("\n");
}
void printVector(double * a,int n){
	for (int i=0; i<n;i++){
			printf("%f\t ",a[i]);
		}
		printf("\n");
		printf("\n");
}
double* init3BSymMatrix(double * a,int n){
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
double* initEyeMatrix(double* a,int n){
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
double* initRandomFullMatrix(double* a,int n){
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
double* initRandomVector(double* a,int n){
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<>distr(1,3);


	for (int i=0; i<n;i++){
			a[i] = distr(gen);
	}
	return a;
}
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


void rowPointers(double* a, int n){

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
	double *a = new double[n*n];
	double *b = new double[n*1];
	double *c = new double[n*1];

	std::cout << " Matrix Multiplication Version " << VERSION_MAJOR << "." << VERSION_MINOR << "\n";

	initRandomFullMatrix(a,n);
	printMatrix(a,n);
	rowPointers(a,n);
	initRandomVector(b, n);
	initRandomVector(c, n);
#ifdef USE_MYSORT
	sort(b, n);
#endif

//	max = pivoting(a,n,n-2);
//	printf("*** %d ***\n",max);
//	initRandomVector(b,n);
//
//	printMatrix(a,n);
	printVector(b,n);
//
//	matVec(a, b, c, n);
//
	printVector(c,n);
//
	gaussElimination(a, b, c,n);

	delete [] a;
	delete [] b;
	delete [] c;
}
