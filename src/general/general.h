inline int index1d(int i, int j, int n){
	int k = i * n + j;
	return k;
}

void printMatrix(double* a, int n);
void printNSMatrix(double *a, int m, int n);
void printVector(double *a, int n);
double* init3BSymMatrix(double *a, int n);
double* initEyeMatrix(double* a,int n);
double* initRandomFullMatrix(double* a,int n);
double* initRandomVector(double* a,int n);
