#include "SolverZeidel.h"
#include "global.h"
#include <algorithm>
#include <ctime>

int SolverZeidel::solve(double eps, int& maxIter)
{
	double	aii;
	double	err = 1.0;
	int		step = 0;
	double	tmp;
	
	memset(x, 0, sizeof(double)*a->n);
	while (err > eps && step < maxIter)
	{
		step++;
		for (int i = 0; i < a->n; i++)
		{
			tmp = 0.0;
			aii = 0;
			for (int k = a->ia[i]; k < a->ia[i+1]; k++) {
				if (a->ja[k] == i) {
					aii = a->a[k];
				}
				else {
					tmp += a->a[k] * x[a->ja[k]];
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
		for (int i = 0; i < a->n; i++)
		{
			tmp = 0.0;
			for (int k = a->ia[i]; k < a->ia[i + 1]; k++) {
				tmp += a->a[k] * x[a->ja[k]];
			}
			err += fabs(tmp - b[i]);
		}
		//int qqqqq = 0; // ZHRV_WARN
		//printf("SEIDEL SOLVER: step = %5d\terr = %16.8e\n", step, err);
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
