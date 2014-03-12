#pragma once;
#include <cmath> 
#include "CSR.h"

const int SOLVER_TYPE_ZEIDEL = 1;
const int SOLVER_TYPE_JACOBI = 2;


class Solver 
{
public:
	CSRMatrix* a;
	double* b;
	double* x;
	virtual void solve(double eps, int& maxIter) = 0;
	void init(int cellsCount);
	void zero();
	void setMatrElement(int i, int j, double** matr3);
	void setRightElement(int i, double* vect3);
	void addMatrElement(int i, int j, double** matr3);
	void addRightElement(int i, double* vect3);
	~Solver();

	void printMatr() { a->print(); }
};

class SolverZeidel : public Solver 
{
	virtual void solve(double eps, int& maxIter);
};

class SolverJacobi : public Solver  
{
	virtual void solve(double eps, int& maxIter);
};






class SolverFactory 
{
public:
	static Solver* create(int type) {
		switch (type) {
		case SOLVER_TYPE_ZEIDEL:
			return new SolverZeidel();
			break;
		case SOLVER_TYPE_JACOBI:
			return new SolverJacobi();
			break;
		}
	}
};