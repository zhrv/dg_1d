#include "SolverHypreGmres.h"
#include "global.h"
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"

int SolverHypreGmres::solve(double eps, int& maxIter)
{
	int result = MatrixSolver::RESULT_OK;
	//static bool firstInit = true;

	/* Set the solution to zero */
	//if (firstInit)
	{
		int    *rows;

		rows = (int*)calloc(local_size, sizeof(int));

		for (int i = 0; i < local_size; i++)
		{
			rows[i] = ilower + i;
			x[i] = 0.0;
		}
		free(rows);

		HYPRE_IJVectorSetValues(xx, local_size, rows, x);
		
		//firstInit = false;
	}


	/* Assemble after setting the coefficients */
	HYPRE_IJMatrixAssemble(A);
	HYPRE_IJVectorAssemble(bb);
	HYPRE_IJVectorAssemble(xx);



	/* Get the parcsr matrix object to use */
	HYPRE_IJMatrixGetObject(A, (void**)&parcsr_A);
	HYPRE_IJVectorGetObject(bb, (void **)&par_bb);
	HYPRE_IJVectorGetObject(xx, (void **)&par_xx);


	/* Choose a solver and solve the system */

	/* GMRES */
	{
		double final_res_norm;

		/* Create solver */
		HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &solver);

		/* Set some parameters (See Reference Manual for more parameters) */
		//HYPRE_GMRESSetSkipRealResidualCheck(solver, 1);
		HYPRE_GMRESSetKDim(solver,			KRYLOV_DIM);
		HYPRE_GMRESSetMaxIter(solver,		maxIter);			/* max iterations */
		HYPRE_GMRESSetTol(solver,			eps);				/* conv. tolerance */
		HYPRE_GMRESSetAbsoluteTol(solver,	0.0);
		HYPRE_GMRESSetPrintLevel(solver,	PRINT_LEVEL);		/* print solve info */
		HYPRE_GMRESSetLogging(solver,		1);					/* needed to get run info later */


		/* Now setup and solve! */
		HYPRE_ParCSRGMRESSetup(solver, parcsr_A, par_bb, par_xx);
		HYPRE_ParCSRGMRESSolve(solver, parcsr_A, par_bb, par_xx);

		
		/* Run info - needed logging turned on */
		int initMaxIter = maxIter;
		HYPRE_GMRESGetNumIterations(solver, &maxIter);
		HYPRE_GMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
		if (initMaxIter <= maxIter) {
			result |= MatrixSolver::RESULT_ERR_MAX_ITER;
		}
		if (final_res_norm >= eps || !isfinite(final_res_norm)) {
			result |= MatrixSolver::RESULT_ERR_CONVERG;
		}

		/* Destory solver */
		HYPRE_ParCSRGMRESDestroy(solver);
	}

	if (result == MatrixSolver::RESULT_OK) {
		int *rows = (int*)calloc(local_size, sizeof(int));
		for (int i = 0; i < local_size; i++)
			rows[i] = ilower + i;

		/* get the local solution */
		HYPRE_IJVectorGetValues(xx, local_size, rows, x);

		delete[] rows;
	}


	return result;
}


void SolverHypreGmres::setParameter(const char* name, int val)
{
	if (strcmp(name, "KRYLOV_DIM") == 0) {
		KRYLOV_DIM = val;
	}
	else {
		SolverHypre::setParameter(name, val);
	}
}

