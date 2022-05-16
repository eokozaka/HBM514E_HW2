#include "stdio.h"
#include "math.h"
#include <random>

int index(int i, int j, int n){
	int k = i*n + j;
	return k;
}
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
