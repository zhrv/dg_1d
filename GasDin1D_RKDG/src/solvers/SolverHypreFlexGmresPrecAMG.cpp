#include "SolverHypreFlexGmresPrecAMG.h"
#include "global.h"
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"

int SolverHypreFlexGmresPrecAMG::solve(double eps, int& maxIter)
{
	int result = MatrixSolver::RESULT_OK;
	/* Set the solution to zero */
	{
		int    *rows;

		rows = (int*)calloc(local_size, sizeof(int));

		for (int i = 0; i < local_size; i++)
		{
			rows[i] = ilower + i;
			x[i] = 0.0;
		}


		HYPRE_IJVectorSetValues(xx, local_size, rows, x);

		free(rows);
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

	/* Flexible GMRES */
	{
		int    num_iterations;
		double final_res_norm;
		int    restart = 30;
		int    modify = 1;


		/* Create solver */
		HYPRE_ParCSRFlexGMRESCreate(MPI_COMM_WORLD, &solver);

		/* Set some parameters (See Reference Manual for more parameters) */
		HYPRE_FlexGMRESSetKDim(solver, restart);
		HYPRE_FlexGMRESSetMaxIter(solver, maxIter);		/* max iterations */
		HYPRE_FlexGMRESSetTol(solver, eps);				/* conv. tolerance */
		HYPRE_FlexGMRESSetPrintLevel(solver, PRINT_LEVEL);		/* print solve info */
		HYPRE_FlexGMRESSetLogging(solver, 1);			/* needed to get run info later */

		/* Now set up the AMG preconditioner and specify any parameters */
		HYPRE_BoomerAMGCreate(&precond);
		HYPRE_BoomerAMGSetPrintLevel(precond, PRINT_LEVEL); /* print amg solution info */
		HYPRE_BoomerAMGSetCoarsenType(precond, 6);
		HYPRE_BoomerAMGSetRelaxType(precond, 6); /* Sym G.S./Jacobi hybrid */
		HYPRE_BoomerAMGSetNumSweeps(precond, 1);
		HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero */
		HYPRE_BoomerAMGSetMaxIter(precond, 1); /* do only one iteration! */
		
		/* Set the FlexGMRES preconditioner */
		HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup, precond);

		//if (modify)
		//	/*	this is an optional call  - if you don't call it, hypre_FlexGMRESModifyPCDefault
		//		is used - which does nothing.  Otherwise, you can define your own, similar to
		//		the one used here */
		//	HYPRE_FlexGMRESSetModifyPC(solver,
		//			                       (HYPRE_PtrToModifyPCFcn)hypre_FlexGMRESModifyPCAMGExample);

		/* Now setup and solve! */
		HYPRE_ParCSRFlexGMRESSetup(solver, parcsr_A, par_bb, par_xx);
		HYPRE_ParCSRFlexGMRESSolve(solver, parcsr_A, par_bb, par_xx);

		/* Run info - needed logging turned on */
		int initMaxIter = maxIter;
		HYPRE_FlexGMRESGetNumIterations(solver, &maxIter);
		HYPRE_FlexGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
		if (initMaxIter <= maxIter) {
			result |= MatrixSolver::RESULT_ERR_MAX_ITER;
		}
		if (final_res_norm >= eps || !isfinite(final_res_norm)) {
			result |= MatrixSolver::RESULT_ERR_CONVERG;
		}

		/* Destory solver and preconditioner */
		HYPRE_ParCSRFlexGMRESDestroy(solver);
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


void SolverHypreFlexGmresPrecAMG::setParameter(const char* name, int val)
{
	if (strcmp(name, "KRYLOV_DIM") == 0) {
		KRYLOV_DIM = val;
	}
	else {
		SolverHypre::setParameter(name, val);
	}
}
