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
	int nproc, iam;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &iam);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

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
	int kk = index1d(2,3,5);
	printf("index for 2 3 5 is: %d\n",kk);
		

	delete [] a;
	delete [] b;
	delete [] c;
	MPI_Finalize();
}
