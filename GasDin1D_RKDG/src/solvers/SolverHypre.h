#pragma once
#include "MatrixSolver.h"
class SolverHypre :
	public MatrixSolver
{
public:
	~SolverHypre();
	
	virtual void init(int cellsCount, int blockDimension);

	virtual void zero();

	virtual void setMatrElement(int i, int j, double** matrDim);
	virtual void setRightElement(int i, double* vectDim);
	virtual void addMatrElement(int i, int j, double** matrDim);
	virtual void addRightElement(int i, double* vectDim);
	virtual void createMatrElement(int i, int j) {}
	virtual void setParameter(const char* name, int val);
	virtual void printToFile(const char* fileName);

	virtual void initCSR() {};

protected:
	void initMatrVectors();

protected:
	int* cols;
	double* values;
	HYPRE_Solver solver, precond;
	HYPRE_Int ilower, iupper, local_size;
	HYPRE_IJMatrix A;
	HYPRE_ParCSRMatrix parcsr_A;
	HYPRE_IJVector bb;
	HYPRE_ParVector par_bb;
	HYPRE_IJVector xx;
	HYPRE_ParVector par_xx;

	int PRINT_LEVEL = 0;
};

