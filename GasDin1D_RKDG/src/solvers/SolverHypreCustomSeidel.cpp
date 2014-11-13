#include "SolverHypreCustomSeidel.h"
#include "global.h"
#include "_hypre_utilities.h"

SolverHypreCustomSeidel::SolverHypreCustomSeidel()
{
}


SolverHypreCustomSeidel::~SolverHypreCustomSeidel()
{
	HYPRE_IJMatrixDestroy(A);
	delete[] x;
	delete[] b;
	delete[] cols;
	delete[] values;
}


int SolverHypreCustomSeidel::solve(double eps, int& maxIter)
{
	int result = MatrixSolver::RESULT_OK;

	/* Assemble after setting the coefficients */
	HYPRE_IJMatrixAssemble(A);

	/* Get the parcsr matrix object to use */
	HYPRE_IJMatrixGetObject(A, (void**)&parcsr_A);


	double	aii;
	double	err = 1.0;
	int		step = 0;
	double	tmp;
	double	*values;
	int		*cols;
	int		size = 0;


	while (err > eps && step < maxIter)
	{
		step++;
		for (int i = ilower; i <= iupper; i++)
		{
			HYPRE_ParCSRMatrixGetRow(parcsr_A, i, &size, &cols, &values);
			HYPRE_ParCSRMatrixRestoreRow(parcsr_A, i, &size, &cols, &values);
			tmp = 0.0;
			aii = 0;
			for (int k = 0; k < size; k++) {
				if (cols[k] == i) {
					aii = values[k];
				}
				else {
					tmp += values[k] * x[cols[k]];
				}
			}
			if (fabs(aii) <= eps*eps)
			{
				log("ZEIDEL_SOLVER: error: a[%d, %d] = 0\n", i, i);
				return MatrixSolver::RESULT_ERR_ZERO_DIAG;
			}
			x[i] = (-tmp + b[i]) / aii;
		}
		err = 0.0;
		for (int i = ilower; i <= iupper; i++)
		{
			HYPRE_ParCSRMatrixGetRow(parcsr_A, i, &size, &cols, &values);
			HYPRE_ParCSRMatrixRestoreRow(parcsr_A, i, &size, &cols, &values);
			tmp = 0.0;
			for (int k = 0; k < size; k++) {
				tmp += values[k] * x[cols[k]];
			}
			err += fabs(tmp - b[i]);
		}
		if (PRINT_LEVEL > 0) {
			printf("SEIDEL SOLVER: step = %5d\terr = %16.8e\n", step, err);
		}
	}
	if (step >= maxIter)
	{
		log("ZEIDEL_SOLVER: (warning) maximum iterations done (%d); error: %e\n", step, err);
		maxIter = step;
		return MatrixSolver::RESULT_ERR_MAX_ITER;
	}
	maxIter = step;
	return MatrixSolver::RESULT_OK;
}



void SolverHypreCustomSeidel::addRightElement(int i, double* vectDim)
{
	MatrixSolver::addRightElement(i, vectDim);
}

void SolverHypreCustomSeidel::setRightElement(int i, double* vectDim)
{
	MatrixSolver::setRightElement(i, vectDim);
}


void SolverHypreCustomSeidel::init(int cellsCount, int blockDimension)
{
	blockDim = blockDimension;
	int n = cellsCount*blockDim;
	ilower = 0;
	iupper = n - 1;
	local_size = iupper - ilower + 1;

	cols	= new int[blockDim];
	values	= new double[local_size];
	x		= new double[local_size];
	b		= new double[local_size];
	
	memset(x, 0, local_size*sizeof(double));
	memset(b, 0, local_size*sizeof(double));

	HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);
	HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
	HYPRE_IJMatrixInitialize(A);
}

void SolverHypreCustomSeidel::zero() {
	HYPRE_IJMatrixDestroy(A);

	HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);
	HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
	HYPRE_IJMatrixInitialize(A);

	memset(x, 0, local_size*sizeof(double));
	memset(b, 0, local_size*sizeof(double));
}



