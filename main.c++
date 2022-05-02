#include"stdio.h"
#include "math.h"
//#include "H5Cpp.h"
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
int* initMatrix(int* a,int n){
	int k;
	for (int i=0; i<n;i++){
		for (int j=0; j<n;j++){
			k = i * n + j;
			a[k] = i;
//			printf("%d ",k);
		}
	}
	return a;
}
int* matmul(int* a,int* b, int* c, int n){
	for (int i=0; i<n;i++){
		for (int j=0; j<n;j++){
			for (int k=0; k<n;k++){
				c[ i * n + j] += a[ i * n + k] * b[ k*n + j ];  
			}
//			printf("%d ",k);
		}
	}
	return c;
}

int main(int argc, char* argv[]){
	int n = 4;
	int k = 0;
	int l = 0;
	int m = 0;
	int* a = new int[n*n];
	int* b = new int[n*n];
	int* c = new int[n*n];

	initMatrix(a,n);
	initMatrix(b,n);

	matmul(a, b, c, n);

	printMatrix(a,n);
	printMatrix(b,n);
	printMatrix(c,n);


}
