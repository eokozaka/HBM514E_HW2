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
#include "iterativeMethods.h"
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
	n = atoi(argv[1]); //Matrix size
	int nprocs, iam;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &iam);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	if(iam == 0) std::cout << " Version " << VERSION_MAJOR << "." << VERSION_MINOR << "\n";
	jacobi(n,nprocs,iam);
	MPI_Finalize();
}
