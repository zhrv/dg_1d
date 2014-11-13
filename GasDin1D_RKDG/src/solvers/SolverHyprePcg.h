#pragma once
#include "SolverHypre.h"
class SolverHyprePcg : public SolverHypre
{
public:
	virtual int solve(double eps, int& maxIter);
	virtual char* getName() { return "HYPRE PCG"; }
};

