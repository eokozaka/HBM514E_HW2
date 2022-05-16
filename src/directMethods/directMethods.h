#include "general.h"
#include <stdio.h>
double partialPivoting(double *a, int nrow, int ncol, int icol);
void lu(double *a, int n, double *l, double *u);
void gaussElimination(double *a, double *b, double* x,int n);
