#include "SolverHyprePcg.h"
#include "global.h"
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"

int SolverHyprePcg::solve(double eps, int& maxIter)
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

	//printToFile("matr.txt");
	//{
	//	FILE *fp = fopen("matr.txt", "w");
	//	int *rows = (int*)calloc(local_size, sizeof(int));
	//	for (int i = 0; i < local_size; i++)
	//		rows[i] = ilower + i;

	//	for (int i = 0; i < local_size; i++) {
	//		HYPRE_IJMatrixGetValues(A, 1, &local_size, &i, rows, x);
	//		for (int j = 0; j < local_size; j++) {
	//			fprintf(fp, "%16.8e ", x[j]);
	//		}
	//		fprintf(fp, "\n");
	//	}
	//	
	//	
	//	HYPRE_IJVectorGetValues(bb, local_size, rows, x);
	//	for (int i = 0; i < local_size; i++) {
	//		fprintf(fp, "%16.8e\n", x[i]);
	//	}
	//	
	//	fclose(fp);
	//	delete[] rows;
	//}



	/* PCG */
	{
		double final_res_norm;

		/* Create solver */
		HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);

		/* Set some parameters (See Reference Manual for more parameters) */
		HYPRE_PCGSetMaxIter(solver, maxIter); /* max iterations */
		HYPRE_PCGSetTol(solver, eps); /* conv. tolerance */
		HYPRE_PCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
		HYPRE_PCGSetPrintLevel(solver, PRINT_LEVEL); /* prints out the iteration info */
		HYPRE_PCGSetLogging(solver, 1); /* needed to get run info later */

		/* Now setup and solve! */
		HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_bb, par_xx);
		HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_bb, par_xx);

		/* Run info - needed logging turned on */
		int initMaxIter = maxIter;
		HYPRE_PCGGetNumIterations(solver, &maxIter);
		HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
		if (initMaxIter <= maxIter) {
			result |= MatrixSolver::RESULT_ERR_MAX_ITER;
		}
		if (final_res_norm >= eps || !isfinite(final_res_norm)) {
			result |= MatrixSolver::RESULT_ERR_CONVERG;
			log("HYPRE PCG solver errror: Final Relative Residual Norm = %e; Config Error = %e\n", final_res_norm, eps);
		}
		//if (myid == 0)
		//{
		//	log("\n");
		//	log("Iterations = %d\n", maxIter);
		//	log("Final Relative Residual Norm = %e\n", final_res_norm);
		//	log("\n");
		//}

		/* Destroy solver */
		HYPRE_ParCSRPCGDestroy(solver);
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


