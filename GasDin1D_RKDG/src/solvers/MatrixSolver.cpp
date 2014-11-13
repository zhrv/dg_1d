#include "MatrixSolver.h"
#include "global.h"
#include "SolverHypreBoomerAMG.h"
#include "SolverHyprePcg.h"
#include "SolverHypreGmres.h"
#include "SolverHypreFlexGmres.h"
#include "SolverHypreFlexGmresPrecAMG.h"
#include "SolverZeidel.h"
#include "SolverHypreCustomSeidel.h"



MatrixSolver* MatrixSolver::create(const char* solverName)
{
	if (strcmp(solverName, "HYPRE_BoomerAMG") == 0) {
		return new SolverHypreBoomerAmg();
	}
	else
	if (strcmp(solverName, "HYPRE_PCG") == 0) {
		return new SolverHyprePcg();
	}
	else
	if (strcmp(solverName, "HYPRE_FlexGMRES") == 0) {
		return new SolverHypreFlexGmres();
	}
	else
	if (strcmp(solverName, "HYPRE_FlexGMRESPrecAMG") == 0) {
		return new SolverHypreFlexGmresPrecAMG();
	}
	else
	if (strcmp(solverName, "HYPRE_GMRES") == 0) {
		return new SolverHypreGmres();
	}
	else
	if (strcmp(solverName, "CUSTOM_HYPRE_SEIDEL") == 0) {
		return new SolverHypreCustomSeidel();
	}
	else
	if (strcmp(solverName, "ZEIDEL") == 0) {
		return new SolverZeidel();
	}
	else {
		log("ERROR (SolverFactory): wrong solver name, used HYPRE Flexible GMRES solver...\n");
		return new SolverHypreFlexGmres();
	}
}





void MatrixSolver::init(int cellsCount, int blockDimension)
{
	blockDim = blockDimension;
	int n = cellsCount*blockDim;
	a = new CSRMatrix(n);
	b = new double[n];
	x = new double[n];
}

void MatrixSolver::zero() {
	memset(x, 0, sizeof(double)*a->n);
	memset(b, 0, sizeof(double)*a->n);
	a->zero();
}

MatrixSolver::~MatrixSolver()
{
	delete a;
	delete[] b;
	delete[] x;
}

void MatrixSolver::setMatrElement(int i, int j, double** matrDim)
{
	for (int ii = 0; ii < blockDim; ++ii)
	{
		for (int jj = 0; jj < blockDim; ++jj)
		{
			a->set(ii+i*blockDim, jj+j*blockDim, matrDim[ii][jj]);
		}
	}
}

void MatrixSolver::setRightElement(int i, double* vectDim)
{
	for (int ii = 0; ii < blockDim; ++ii)
	{
		b[ii+i*blockDim] = vectDim[ii];
	}
}


void MatrixSolver::addMatrElement(int i, int j, double** matrDim)
{
	for (int ii = 0; ii < blockDim; ++ii)
	{
		for (int jj = 0; jj < blockDim; ++jj)
		{
			a->add(ii+i*blockDim, jj+j*blockDim, matrDim[ii][jj]);
		}
	}
}

void MatrixSolver::createMatrElement(int i, int j) {
	for (int ii = 0; ii < blockDim; ++ii)
	{
		for (int jj = 0; jj < blockDim; ++jj)
		{
			a->init(ii + i*blockDim, jj + j*blockDim);
		}
	}

}

void MatrixSolver::initCSR()
{
	a->assemble();
}

void MatrixSolver::addRightElement(int i, double* vectDim)
{
	for (int ii = 0; ii < blockDim; ii++)
	{
		b[ii+i*blockDim] += vectDim[ii];
	}
}

void MatrixSolver::printToFile(const char* fileName)
{
	a->printToFile(fileName);
}


int SolverJacobi::solve(double eps, int& maxIter)
{
	if (!tempXAlloc) {
		tempX = new double [a->n];
		tempXAlloc = true;
	}
	double	aii;
	double	err = 1.0;
	int		step = 0;
	double	tmp;
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
				if (i == a->ja[k])	// i == j
				{
					aii = a->a[k];
				} else {
					tmp += a->a[k]*x[a->ja[k]];
				}
			}
			if (fabs(aii) <= eps*eps) 
			{
				log("JACOBI_SOLVER: error: a[%d, %d] = 0\n", i, i);
				return MatrixSolver::RESULT_ERR_ZERO_DIAG;
			}
			//x[i] = (-tmp+b[i])/aii;
			tempX[i] = (-tmp+b[i])/aii;
		}
		err = 0.0;
		for (int i = 0; i < a->n; i++)
		{
			tmp = 0.0;
			for (int k = a->ia[i]; k < a->ia[i+1]; k++)
			{
				x[a->ja[k]] = tempX[a->ja[k]];
				tmp += a->a[k]*x[a->ja[k]];
			}
			err += fabs(tmp-b[i]);
		}
		int qqqqq = 0; // ZHRV_WARN
	}
	if (step >= maxIter)
	{
		log("JACOBI_SOLVER: (warning) maximum iterations done (%d); error: %e\n", step, err);
		maxIter = step;
		return MatrixSolver::RESULT_ERR_MAX_ITER;
	}
	maxIter = step;
	return MatrixSolver::RESULT_OK;
}


