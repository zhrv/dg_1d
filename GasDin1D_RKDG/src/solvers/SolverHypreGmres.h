#pragma once
#include "SolverHypre.h"
class SolverHypreGmres : public SolverHypre
{
public:
	virtual int solve(double eps, int& maxIter);
	virtual char* getName() { return "HYPRE GMRES"; }
	void setParameter(const char* name, int val);
private:
	int KRYLOV_DIM = 30;
};

