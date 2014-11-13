#pragma once
#include "SolverHypre.h"

class SolverHypreFlexGmres : public SolverHypre
{
public:
	virtual int solve(double eps, int& maxIter);
	virtual char* getName() { return "HYPRE Flexible GMRES"; }
	void setParameter(const char* name, int val);
private:
	int KRYLOV_DIM = 30;
};

