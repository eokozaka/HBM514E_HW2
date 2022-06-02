double partialPivoting(double *a, int nrow, int ncol, int icol);
void lu(double *a, int n, double *l, double *u);
void gaussElimination2DBlockCyclic(int nprocs ,int size);
void gaussEliminationRowBlockCyclic(int nprocs ,int size);
void gaussEliminationSerial(int size);
