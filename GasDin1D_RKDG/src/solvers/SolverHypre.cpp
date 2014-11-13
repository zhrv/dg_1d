#include "SolverHypre.h"


void SolverHypre::initMatrVectors()
{
	HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);
	HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
	HYPRE_IJMatrixInitialize(A);

	/* Create the rhs and solution */
	HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &bb);
	HYPRE_IJVectorSetObjectType(bb, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(bb);

	HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &xx);
	HYPRE_IJVectorSetObjectType(xx, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(xx);
}

void SolverHypre::init(int cellsCount, int blockDimension)
{
	blockDim = blockDimension;
	int n = cellsCount*blockDim;
	ilower = 0;
	iupper = n - 1;
	local_size = iupper - ilower + 1;

	cols	= new int[blockDim];
	values	= new double[n];
	x		= new double[n];

	initMatrVectors();

}

void SolverHypre::zero() {
	HYPRE_IJMatrixDestroy(A);
	HYPRE_IJVectorDestroy(bb);
	HYPRE_IJVectorDestroy(xx);

	HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);
	HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
	HYPRE_IJMatrixInitialize(A);

	/* Create the rhs and solution */
	HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &bb);
	HYPRE_IJVectorSetObjectType(bb, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(bb);

	HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &xx);
	HYPRE_IJVectorSetObjectType(xx, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(xx);

	memset(x, 0, local_size*sizeof(double));
}

SolverHypre::~SolverHypre()
{
	HYPRE_IJMatrixDestroy(A);
	HYPRE_IJVectorDestroy(bb);
	HYPRE_IJVectorDestroy(xx);
	delete[] x;
}

void SolverHypre::setMatrElement(int i, int j, double** matrDim)
{
	int row;

	for (int ii = 0; ii < blockDim; ++ii)
	{
		row = ii + i*blockDim;
		for (int jj = 0; jj < blockDim; ++jj)
		{
			values[jj] = matrDim[ii][jj];
			cols[jj] = jj + j*blockDim;
		}
		HYPRE_IJMatrixSetValues(A, 1, &blockDim, &row, cols, values);
	}


}

void SolverHypre::addMatrElement(int i, int j, double** matrDim)
{
	int row;

	for (int ii = 0; ii < blockDim; ++ii)
	{
		row = ii + i*blockDim;
		for (int jj = 0; jj < blockDim; ++jj)
		{
			values[jj] = matrDim[ii][jj];
			cols[jj] = jj + j*blockDim;
		}
		HYPRE_IJMatrixAddToValues(A, 1, &blockDim, &row, cols, values);
	}

}


void SolverHypre::setRightElement(int i, double* vectDim)
{
	for (int ii = 0; ii < blockDim; ++ii)
	{
		cols[ii] = ii + i*blockDim;
	}
	HYPRE_IJVectorSetValues(bb, blockDim, cols, vectDim);
}


void SolverHypre::addRightElement(int i, double* vectDim)
{
	for (int ii = 0; ii < blockDim; ++ii)
	{
		cols[ii] = ii + i*blockDim;
	}
	HYPRE_IJVectorAddToValues(bb, blockDim, cols, vectDim);
}

void SolverHypre::setParameter(const char* name, int val)
{
	if (strcmp(name, "PRINT_LEVEL") == 0) {
		PRINT_LEVEL = val;
	}
}


void SolverHypre::printToFile(const char* fileName)
{
	int n = local_size;
	int nc = 1;
	double * x = new double[n];
	int * cols = new int[n];
	FILE * fp = fopen(fileName, "w");
	for (int i = 0; i < n; i++) {
		cols[i] = ilower + i;
	}

	for (int row = 0; row < n; row++) {
		HYPRE_IJMatrixGetValues(A, 1, &nc, &row, &row, x);
		for (int i = 0; i < nc; i++) {
			fprintf(fp, "%25.16e  ", x[i]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n\n=============================================================================================\n\n\n");
	HYPRE_IJVectorGetValues(xx, local_size, cols, x);
	for (int i = 0; i < n; i++) {
		fprintf(fp, "%25.16e  ", x[i]);
	}

	fclose(fp);
	delete[] x;
	delete[] cols;
}

