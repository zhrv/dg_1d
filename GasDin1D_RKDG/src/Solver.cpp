#include "Solver.h"

void Solver::init(int cellsCount, int bSize) {
	blockSize = bSize;
	int n = cellsCount*blockSize;
	a = new CSRMatrix(n);
	b = new double[n];
	x = new double[n];
}

void Solver::zero() {
	memset(x, 0, sizeof(double)*a->n);
	memset(b, 0, sizeof(double)*a->n);
	a->zero();
}


Solver::~Solver() {
	delete a;
	delete[] b;
	delete[] x;
}

void Solver::setMatrElement(int i, int j, double** matr9) {
	for (int ii = 0; ii < blockSize; ii++) {
		for (int jj = 0; jj < blockSize; jj++) {
			a->set(ii + i*blockSize, ii + j * blockSize, matr9[ii][jj]);
		}
	}
}

void Solver::setRightElement(int i, double* vect9) {
	for (int ii = 0; ii < blockSize; ii++) {
		b[ii + i*blockSize] = vect9[ii];
	}
}


void Solver::addMatrElement(int i, int j, double** matr9) {
	for (int ii = 0; ii < blockSize; ii++) {
		for (int jj = 0; jj < blockSize; jj++) {
			a->add(ii + i*blockSize, ii + j * blockSize, matr9[ii][jj]);
		}
	}
}

void Solver::addRightElement(int i, double* vect9) {
	for (int ii = 0; ii < blockSize; ii++) {
		b[ii + i*blockSize] += vect9[ii];
	}
}


void SolverZeidel::solve(double eps, int& maxIter) {
	double aii;
	double err = 1.0;
	int step = 0;
	double tmp;
	//memset(x, 0, sizeof(double)*a->n);
	while(err > eps && step < maxIter)
	{
		step++;
		for (int i = 0; i < a->n; i++)
		{
			tmp = 0.0;
			aii = 0;
			for (int k = a->ia[i]; k < a->ia[i+1]; k++)
			{
				if (i == a->ja[k])
				{
					aii = a->a[k];
				} else {
					tmp += a->a[k]*x[a->ja[k]];
				}
			}
			if (aii == 0) 
			{
				printf("ZEIDEL_SOLVER: error: a[%d, %d] = 0\n", i, i);

			}
			x[i] = (-tmp+b[i])/aii;
		}
		err = 0.0;
		for (int i = 0; i < a->n; i++)
		{
			tmp = 0.0;
			for (int k = a->ia[i]; k < a->ia[i+1]; k++)
			{
				tmp += a->a[k]*x[a->ja[k]];
			}
			err += fabs(tmp-b[i]);
		}
		//int qqqqq = 0; // ZHRV_WARN
	}
	if (step >= maxIter)
	{
		printf("ZEIDEL_SOLVER: (warning) maximum iterations done (%d); error: %e\n", step, err);
	}
	maxIter = step;
}

void SolverJacobi::solve(double eps, int& maxIter)
{
	printf("JACOBI_SOLVER: not released...\n");
}

