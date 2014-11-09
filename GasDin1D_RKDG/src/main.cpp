#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "Solver.h"
#include "functions.h"
#include <float.h>

const int		FUNC_COUNT = 3;
const int		MATR_BLOCK = 3 * FUNC_COUNT;

const double	CFL		= 2.5e-3;

const double	LIM_ALPHA = 2.0;

const int		N		= 100;
const double	XMIN	= -1.0; 
const double	XMAX	=  1.0;
const double	EPS		= 1.0e-5;
const double	GAM		= 1.4;
const double	AGAM	= GAM-1.0;
const double	TMAX	= 0.2;

const int		MAX_ITER	= 5000;
const int		SAVE_STEP	= 100;
const int		PRINT_STEP	= 1;

double **ro, **ru, **re;
double **ro_, **ru_, **re_;
double *r, *u, *e, *p;
double *dr, *du, *de;


double **cellGP, **cellGW, *cellC;

double *smMu;

double h	=	(XMAX-XMIN)/N;
double tau = CFL*h; // <= 1.e-5

#define FLUX     roe_orig
#define FLUX_RHS flux_rim

#define LIMITER_I  calcLimiterEigenv
#define LIMITER_II calcLimiter_II

const bool USE_LIMITER_I	= true;
const bool USE_LIMITER_II	= false;
const bool USE_SMOOTHER		= false;

Solver *S;

/*  базисные функции и поля  */
double getField(int i, int iCell, double x);
double getF(int i, int iCell, double x);
double getDF(int i, int iCell, double x);
double getFieldsAvg(int idField, int iCell, double x);
/*  функции  */
void init();
void done();

void primToCons(double r, double p, double u, double &ro, double &ru, double &re);
void consToPrim(double &r, double &p, double &u, double ro, double ru, double re);

void calcIntegral();
void calcMatrWithTau();
void calcMatrFlux();
void calcRHS();
void calcLimiterCons();
void calcLimiterEigenv();
void calcLimiter_II();
void calcSmoother();

int main(int argc, char** argv) 
{
#ifdef _DEBUG
	_controlfp(~(_MCW_EM & (~_EM_INEXACT) & (~_EM_UNDERFLOW)), _MCW_EM);
#endif
	init();
	double t = 0.0;
	int step = 0;
	while (t < TMAX) {
		t += tau;
		step++;
		S->zero();
		
		calcIntegral();			// вычисляем интеграл от(dF / dU)*deltaU*dFi / dx
		calcMatrWithTau();		// вычисляем матрицы перед производной по времени
		calcMatrFlux();			// Вычисляем потоковые величины 
		calcRHS();				// Вычисляем столбец правых членов

		int maxIter = MAX_ITER;
		
		S->solve(EPS, maxIter);
		
		for (int iCell = 0, ind = 0; iCell < N; iCell++, ind += MATR_BLOCK) {
			int shift = 0;
			for (int j = 0; j < FUNC_COUNT; j++)
				ro[iCell][j] += S->x[ind + (shift++)];
			for (int j = 0; j < FUNC_COUNT; j++)
				ru[iCell][j] += S->x[ind + (shift++)];
			for (int j = 0; j < FUNC_COUNT; j++)
				re[iCell][j] += S->x[ind + (shift++)];
		}

		if (USE_LIMITER_I)  LIMITER_I();
		if (USE_LIMITER_II) LIMITER_II();
		if (USE_SMOOTHER)   calcSmoother();

		for (int iCell = 0, ind = 0; iCell < N; iCell++, ind += 9) {
			consToPrim(r[iCell], p[iCell], u[iCell], ro[iCell][0], ru[iCell][0], re[iCell][0]);
		}


		if (step % PRINT_STEP == 0) {
			printf("%10d | ITER: %4d | RO: ... | RU: ... | RE ... \n", step, maxIter);
		}
		if (step % SAVE_STEP == 0) {
			char str[50];
			sprintf(str, "res_%010d.csv", step);
			FILE * fp = fopen(str, "w");
			fprintf(fp, "x,r,p,u\n");
			for (int i = 0; i < N; i++) {
				fprintf(fp, "%25.15e, %25.15e, %25.15e, %25.15e\n", XMIN+h*i, r[i], p[i], u[i]);
			}
			fclose(fp);
			printf("           | File '%s' is written...\n", str);
		}
	}

	return 0;
}



void memAlloc() {
	ro = new double*[N];
	ru = new double*[N];
	re = new double*[N];
	ro_ = new double*[N];
	ru_ = new double*[N];
	re_ = new double*[N];
	matrM = new double**[N];
	cellGP = new double*[N];
	cellGW = new double*[N];
	cellC = new double[N];
	for (int i = 0; i < N; i++) {
		ro[i] = new double[FUNC_COUNT];
		ru[i] = new double[FUNC_COUNT];
		re[i] = new double[FUNC_COUNT];
		ro_[i] = new double[FUNC_COUNT];
		ru_[i] = new double[FUNC_COUNT];
		re_[i] = new double[FUNC_COUNT];
		cellGP[i] = new double[2];
		cellGW[i] = new double[2];
		matrM[i] = new double*[FUNC_COUNT];
		for (int j = 0; j < FUNC_COUNT; j++) {
			matrM[i][j] = new double[FUNC_COUNT];
		}
	}
	r = new double[N];
	u = new double[N];
	p = new double[N];
	vect = new double[FUNC_COUNT];
	matr = new double*[FUNC_COUNT];
	matr1 = new double*[FUNC_COUNT];
	for (int i = 0; i < FUNC_COUNT; i++) {
		matr[i] = new double[FUNC_COUNT];
		matr1[i] = new double[FUNC_COUNT];
	}
	A = new double*[3];
	L = new double*[3];
	R = new double*[3];
	Am = new double*[3];
	Ap = new double*[3];
	for (int i = 0; i < 3; i++) {
		A[i] = new double[3];
		Am[i] = new double[3];
		Ap[i] = new double[3];
		L[i] = new double[3];
		R[i] = new double[3];
	}
	matr9 = new double*[MATR_BLOCK];
	for (int i = 0; i < MATR_BLOCK; i++) {
		matr9[i] = new double[MATR_BLOCK];
	}

	mA = new double**[2];
	for (int i = 0; i < 2; i++) {
		mA[i] = new double*[3];
		for (int j = 0; j < 3; j++) {
			mA[i][j] = new double[3];
		}
	}
	smMu = new double[N + 5];
}

void memFree() {
	for (int i = 0; i < N; i++) {
		delete[] ro[i];
		delete[] ru[i];
		delete[] re[i];
		delete[] ro_[i];
		delete[] ru_[i];
		delete[] re_[i];
		delete[] cellGP[i];
		delete[] cellGW[i];
		for (int j = 0; j < FUNC_COUNT; j++) {
			delete[] matrM[i][j];
		}
		delete[] matrM[i];
	}
	delete[] ro;
	delete[] ru;
	delete[] re;
	delete[] ro_;
	delete[] ru_;
	delete[] re_;
	delete[] matrM;
	delete[] cellC;

	delete[] r;
	delete[] u;
	delete[] p;
	for (int i = 0; i < FUNC_COUNT; i++) {
		delete[] A[i];
		delete[] Am[i];
		delete[] Ap[i];
		delete[] L[i];
		delete[] R[i];
		delete[] matr[i];
		delete[] matr1[i];
	}
	delete[] A;
	delete[] Ap;
	delete[] Am;
	delete[] L;
	delete[] R;
	delete[] matr;
	delete[] matr1;
	delete[] vect;
	for (int i = 0; i < 3 * FUNC_COUNT; i++) {
		delete[] matr9[i];
	}
	delete[] matr9;

	for (int i = 0; i < FUNC_COUNT; i++) {
		for (int j = 0; j < FUNC_COUNT; j++) {
			delete[] mA[i][j];
		}
		delete[] mA[i];
	}
	delete[] mA;
	delete[] smMu;
}


void calcMassMatr() {
	for (int iCell = 0; iCell < N; iCell++){

		for (int i = 0; i < FUNC_COUNT; i++) {
			for (int j = 0; j < FUNC_COUNT; j++){
				matrM[iCell][i][j] = 0.0;
				for (int iGP = 0; iGP < 2; iGP++) {
					double x = cellGP[iCell][iGP];
					matrM[iCell][i][j] = matrM[iCell][i][j] + cellGW[iCell][iGP]*getF(i, iCell, x)*getF(j, iCell, x);
				}
			}
		}
	}
}


void init() {
	memAlloc();
	
	// узлы квадратур Гаусса и центры ячеек
	for (int i = 0; i < N; i++){
		double x1 = XMIN + i*h;
		double x2 = x1 + h;

		cellGW[i][0] = h*0.5;
		cellGW[i][1] = h*0.5;
		double tmp = h / sqrt(3.0);
		cellGP[i][0] = 0.5*((x1 + x2) - tmp);
		cellGP[i][1] = 0.5*((x1 + x2) + tmp);
		cellC[i] = 0.50*(x1 + x2);
	}

	S = SolverFactory::create(SOLVER_TYPE_ZEIDEL);
	S->init(N, MATR_BLOCK);
	for (int i = 0; i < N; i++) {
		double x = cellC[i];
		// Sod  !!! XMIN = -1.0; XMAX = 1.0;
		if (x < 0.0) {
			r[i] = 1.0;
			p[i] = 1.0;
			u[i] = 0.0;
		} else {
			r[i] = 0.125;
			p[i] = 0.1;
			u[i] = 0.0;
		}
		// Lax  !!! XMIN = -1.0; XMAX = 1.0;
		//if (x < 0.0) {
		//	r[i] = 0.445;
		//	p[i] = 3.528;
		//	u[i] = 0.698;
		//}
		//else {
		//	r[i] = 0.5;
		//	p[i] = 0.571;
		//	u[i] = 0.0;
		//}
		// Blast waves !!! XMIN = 0.0; XMAX = 1.0;
		//if (x <= 0.1) {
		//	r[i] = 1.0;
		//	p[i] = 2500.0 * r[i] * AGAM;
		//	u[i] = 0.0;
		//}
		//else if (x >= 0.9) {
		//	r[i] = 1.0;
		//	p[i] = 250.0 * r[i] * AGAM;
		//	u[i] = 0.0;
		//}
		//else {
		//	r[i] = 1.0;
		//	p[i] = 0.025 * r[i] * AGAM;
		//	u[i] = 0.0;
		//}
		//r[i] = 0.125;
		//p[i] = 0.1;
		//u[i] = 1.0;

		primToCons(r[i], p[i], u[i], ro[i][0], ru[i][0], re[i][0]);
		for (int j = 1; j < FUNC_COUNT; j++) {
			ro[i][j] = ru[i][j] = re[i][j] = 0.0;
		}
		//double dt = CFL*h / (fabs(u[i]) + sqrt(GAM*p[i] / r[i]));
		//if (dt < tau) tau = dt;
	}
	printf("TAU = %25.16e\n\n", tau);


	calcMassMatr();

}

void done() {
	delete S;
	memFree();
}



void primToCons(double r, double p, double u, double &ro, double &ru, double &re) {
	ro = r;
	ru = r*u;
	re = (p / AGAM + r*u*u*0.5);
}

void consToPrim(double &r, double &p, double &u, double ro, double ru, double re) {
	r = ro;
	u = ru / r;
	p = AGAM*(re - r*u*u*0.5);
}



void calcIntegral() {


	for (int iCell = 0; iCell < N; iCell++) {

		for (int i = 0; i < MATR_BLOCK; i++) {
			for (int j = 0; j < MATR_BLOCK; j++) {
				matr9[i][j] = 0.0;
			}
		}

		for (int k = 0; k < 2; k++) {
			double x = cellGP[iCell][k];
			double fRO = getField(0, iCell, x);
			double fRU = getField(1, iCell, x);
			double fRE = getField(2, iCell, x);
			double fU = fRU / fRO;
			double fP, fE;
			fE = fRE / fRO - 0.50*(fU*fU);
			URS(1, fRO, fP, fE, GAM);
			double c2 = GAM*fP / fRO;
			//calcMatrA(fU, c2, GAM, mA[k]);
			calcMatrA(c2, fU, GAM, mA[k]);
		}

		for (int ii = 0; ii < 3; ii++) {
			for (int jj = 0; jj < 3;jj++){
				for (int i = 0; i < FUNC_COUNT; i++) {
					for (int j = 0; j < FUNC_COUNT; j++) {
						matr[i][j] = 0.0;
						for (int k = 0; k < 2; k++) {
							double x = cellGP[iCell][k];
							matr[i][j] = matr[i][j] - cellGW[iCell][k]*mA[k][ii][jj]*getDF(i, iCell, x)*getF(j, iCell, x);
						}
					}
				}
				//call Solver_AddMatr3(m9, m3, ii, jj)
				addMatr3ToMatr9(matr9, matr, ii, jj, FUNC_COUNT);
			}
		}

		//call Solver_AddMatr9(matr, m9, iCell, iCell)
		S->addMatrElement(iCell, iCell, matr9);
	}

}



void calcMatrWithTau() {
	for (int iCell = 0; iCell < N; iCell++) {

		for (int i = 0; i < MATR_BLOCK; i++) {
			for (int j = 0; j < MATR_BLOCK; j++) {
				matr9[i][j] = 0.0;
			}
		}


		for (int i = 0; i < FUNC_COUNT; i++) {
			for (int j = 0; j < FUNC_COUNT; j++) {
				matr[i][j] = matrM[iCell][i][j] / tau;
			}
		}

		for (int ii = 0; ii < 3; ii++) {
			addMatr3ToMatr9(matr9, matr, ii, ii, FUNC_COUNT);
		}

		S->addMatrElement(iCell, iCell, matr9);

	}

}


void calcMatrFlux() { 
	double ri, ei, pi, ui, vi, wi, rb, pb, ub, vb, wb, re, pe, ue, ve, we;
	vb = 0.0; ve = 0.0; wb = 0.0; we = 0.0;

	for (int iCell = 1; iCell < N - 1; iCell++) {

		{ //!< вычисляем потоки на границе iCell+1/2

			for (int i = 0; i < MATR_BLOCK; i++){
				for (int j = 0; j < MATR_BLOCK; j++){
					matr9[i][j] = 0.0;
				}
			}

			double x1 = cellC[iCell] + 0.50*h;

			double ROl = getField(0, iCell, x1);
			double RUl = getField(1, iCell, x1);
			double REl = getField(2, iCell, x1);
			double ROr = getField(0, iCell + 1, x1);
			double RUr = getField(1, iCell + 1, x1);
			double REr = getField(2, iCell + 1, x1);
			consToPrim(rb, pb, ub, ROl, RUl, REl);
			consToPrim(re, pe, ue, ROr, RUr, REr);
			FLUX(ri, ei, pi, ui, vi, wi,
				rb, pb, ub, vb, wb,
				re, pe, ue, ve, we, GAM);


			calcMatrAPM(GAM*pi / ri, ui, GAM, Am, Ap);



			for (int i = 0; i < MATR_BLOCK; i++) {
				for (int j = 0; j < MATR_BLOCK; j++){
					matr9[i][j] = 0.0;
				}
			}

			for (int ii = 1; ii < 3; ii++){
				for (int jj = 1; jj < 3; jj++){
					for (int i = 0; i < FUNC_COUNT; i++){
						for (int j = 0; j < FUNC_COUNT; j++){
							matr[i][j] = Am[ii][jj] * getF(i, iCell, x1)*getF(j, iCell + 1, x1);
						}
					}
					addMatr3ToMatr9(matr9, matr, ii, jj, FUNC_COUNT);
				}
			}

			//call Solver_AddMatr9(matr, m9, iCell, iCell)
			S->addMatrElement(iCell, iCell+1, matr9);


			for (int i = 0; i < MATR_BLOCK; i++) {
				for (int j = 0; j < MATR_BLOCK; j++) {
					matr9[i][j] = 0.0;
				}
			}
			for (int ii = 0; ii < 3; ii++){
				for (int jj = 0; jj < 3; jj++){
					for (int i = 0; i < FUNC_COUNT; i++){
						for (int j = 0; j < FUNC_COUNT; j++){
							matr[i][j] = Ap[ii][jj] * getF(i, iCell, x1)*getF(j, iCell, x1);
						}
					}
					addMatr3ToMatr9(matr9, matr, ii, jj, FUNC_COUNT);
				}
			}
			S->addMatrElement(iCell, iCell, matr9);
			//call Solver_AddMatr9(matr, m9, iCell, iCell + 1)

		} // iCell+1/2



		{ //!< вычисляем потоки на границе iCell-1/2

			for (int i = 0; i < MATR_BLOCK; i++){
				for (int j = 0; j < MATR_BLOCK; j++){
					matr9[i][j] = 0.0;
				}
			}

			//осреднение по Роу
			double x1 = cellC[iCell] + 0.50*h;

			double ROl = getField(0, iCell-1, x1);
			double RUl = getField(1, iCell-1, x1);
			double REl = getField(2, iCell-1, x1);
			double ROr = getField(0, iCell, x1);
			double RUr = getField(1, iCell, x1);
			double REr = getField(2, iCell, x1);
			consToPrim(rb, pb, ub, ROl, RUl, REl);
			consToPrim(re, pe, ue, ROr, RUr, REr);
			FLUX(ri, ei, pi, ui, vi, wi,
				rb, pb, ub, vb, wb,
				re, pe, ue, ve, we, GAM);


			calcMatrAPM(GAM*pi / ri, ui, GAM, Am, Ap);



			for (int i = 0; i < MATR_BLOCK; i++) {
				for (int j = 0; j < MATR_BLOCK; j++){
					matr9[i][j] = 0.0;
				}
			}

			for (int ii = 1; ii < 3; ii++){
				for (int jj = 1; jj < 3; jj++){
					for (int i = 0; i < FUNC_COUNT; i++){
						for (int j = 0; j < FUNC_COUNT; j++){
							matr[i][j] = -Am[ii][jj] * getF(i, iCell, x1)*getF(j, iCell, x1);
						}
					}
					addMatr3ToMatr9(matr9, matr, ii, jj, FUNC_COUNT);
				}
			}

			//call Solver_AddMatr9(matr, m9, iCell, iCell)
			S->addMatrElement(iCell, iCell, matr9);


			for (int i = 0; i < MATR_BLOCK; i++) {
				for (int j = 0; j < MATR_BLOCK; j++) {
					matr9[i][j] = 0.0;
				}
			}
			for (int ii = 0; ii < 3; ii++){
				for (int jj = 0; jj < 3; jj++){
					for (int i = 0; i < FUNC_COUNT; i++){
						for (int j = 0; j < FUNC_COUNT; j++){
							matr[i][j] = -Ap[ii][jj] * getF(i, iCell, x1)*getF(j, iCell - 1, x1);
						}
					}
					addMatr3ToMatr9(matr9, matr, ii, jj, FUNC_COUNT);
				}
			}
			S->addMatrElement(iCell, iCell-1, matr9);
			//call Solver_AddMatr9(matr, m9, iCell, iCell + 1)


		} // iCell-1/2



	} // for (iCell)



	{	int iCell = 0;

		{ //!< вычисляем потоки на границе iCell+1/2

			for (int i = 0; i < MATR_BLOCK; i++){
				for (int j = 0; j < MATR_BLOCK; j++){
					matr9[i][j] = 0.0;
				}
			}

			//осреднение по Роу
			double x1 = cellC[iCell] + 0.50*h;

			double ROl = getField(0, iCell, x1);
			double RUl = getField(1, iCell, x1);
			double REl = getField(2, iCell, x1);
			double ROr = getField(0, iCell + 1, x1);
			double RUr = getField(1, iCell + 1, x1);
			double REr = getField(2, iCell + 1, x1);
			consToPrim(rb, pb, ub, ROl, RUl, REl);
			consToPrim(re, pe, ue, ROr, RUr, REr);
			FLUX(ri, ei, pi, ui, vi, wi,
				rb, pb, ub, vb, wb,
				re, pe, ue, ve, we, GAM);


			calcMatrAPM(GAM*pi / ri, ui, GAM, Am, Ap);



			for (int i = 0; i < MATR_BLOCK; i++) {
				for (int j = 0; j < MATR_BLOCK; j++){
					matr9[i][j] = 0.0;
				}
			}

			for (int ii = 1; ii < 3; ii++){
				for (int jj = 1; jj < 3; jj++){
					for (int i = 0; i < FUNC_COUNT; i++){
						for (int j = 0; j < FUNC_COUNT; j++){
							matr[i][j] = Am[ii][jj] * getF(i, iCell, x1)*getF(j, iCell + 1, x1);
						}
					}
					addMatr3ToMatr9(matr9, matr, ii, jj, FUNC_COUNT);
				}
			}

			//call Solver_AddMatr9(matr, m9, iCell, iCell)
			S->addMatrElement(iCell, iCell+1, matr9);


			for (int i = 0; i < MATR_BLOCK; i++) {
				for (int j = 0; j < MATR_BLOCK; j++) {
					matr9[i][j] = 0.0;
				}
			}
			for (int ii = 0; ii < 3; ii++){
				for (int jj = 0; jj < 3; jj++){
					for (int i = 0; i < FUNC_COUNT; i++){
						for (int j = 0; j < FUNC_COUNT; j++){
							matr[i][j] = Ap[ii][jj] * getF(i, iCell, x1)*getF(j, iCell, x1);
						}
					}
					addMatr3ToMatr9(matr9, matr, ii, jj, FUNC_COUNT);
				}
			}
			S->addMatrElement(iCell, iCell, matr9);
			//call Solver_AddMatr9(matr, m9, iCell, iCell + 1)

		} // iCell+1/2



		{ //!< вычисляем потоки на границе iCell-1/2

			for (int i = 0; i < MATR_BLOCK; i++){
				for (int j = 0; j < MATR_BLOCK; j++){
					matr9[i][j] = 0.0;
				}
			}

			//осреднение по Роу
			double x1 = cellC[iCell] + 0.50*h;

			double ROr = getField(0, iCell, x1);
			double RUr = getField(1, iCell, x1);
			double REr = getField(2, iCell, x1);
			//consToPrim(rb, pb, ub, ROl, RUl, REl);
			consToPrim(re, pe, ue, ROr, RUr, REr);
			rb = re; pb = pe; ub = ue;
			FLUX(ri, ei, pi, ui, vi, wi,
				rb, pb, ub, vb, wb,
				re, pe, ue, ve, we, GAM);


			calcMatrAPM(GAM*pi / ri, ui, GAM, Am, Ap);



			for (int i = 0; i < MATR_BLOCK; i++) {
				for (int j = 0; j < MATR_BLOCK; j++){
					matr9[i][j] = 0.0;
				}
			}

			for (int ii = 1; ii < 3; ii++){
				for (int jj = 1; jj < 3; jj++){
					for (int i = 0; i < FUNC_COUNT; i++){
						for (int j = 0; j < FUNC_COUNT; j++){
							matr[i][j] = -Am[ii][jj] * getF(i, iCell, x1)*getF(j, iCell, x1);
						}
					}
					addMatr3ToMatr9(matr9, matr, ii, jj, FUNC_COUNT);
				}
			}

			//call Solver_AddMatr9(matr, m9, iCell, iCell)
			S->addMatrElement(iCell, iCell, matr9);


			//for (int i = 0; i < 9; i++) {
			//	for (int j = 0; j < 9; j++) {
			//		matr9[i][j] = 0.0;
			//	}
			//}
			//for (int ii = 0; ii < 3; ii++){
			//	for (int jj = 0; jj < 3; jj++){
			//		for (int i = 0; i < 3; i++){
			//			for (int j = 0; j < 3; j++){
			//				matr[i][j] = -Ap[ii][jj] * getF(i, iCell, x1)*getF(j, iCell - 1, x1);
			//			}
			//		}
			//		addMatr3ToMatr9(matr9, matr, ii, jj);
			//	}
			//}
			//S->addMatrElement(iCell, iCell - 1, matr9);
			////call Solver_AddMatr9(matr, m9, iCell, iCell + 1)


		} // iCell-1/2



	} // (iCell=0)


	{	int iCell = N - 1;

		{ //!< вычисляем потоки на границе iCell+1/2

			for (int i = 0; i < MATR_BLOCK; i++){
				for (int j = 0; j < MATR_BLOCK; j++){
					matr9[i][j] = 0.0;
				}
			}

			//осреднение по Роу
			double x1 = cellC[iCell] + 0.50*h;

			double ROl = getField(0, iCell, x1);
			double RUl = getField(1, iCell, x1);
			double REl = getField(2, iCell, x1);
			//double ROr = getField(0, iCell + 1, x1);
			//double RUr = getField(1, iCell + 1, x1);
			//double REr = getField(2, iCell + 1, x1);
			consToPrim(rb, pb, ub, ROl, RUl, REl);
			//consToPrim(re, pe, ue, ROr, RUr, REr);
			re = rb; pe = pb; ue = ub;
			FLUX(ri, ei, pi, ui, vi, wi,
				rb, pb, ub, vb, wb,
				re, pe, ue, ve, we, GAM);


			calcMatrAPM(GAM*pi / ri, ui, GAM, Am, Ap);



			//for (int i = 0; i < 9; i++) {
			//	for (int j = 0; j < 9; j++){
			//		matr9[i][j] = 0.0;
			//	}
			//}

			//for (int ii = 1; ii < 3; ii++){
			//	for (int jj = 1; jj < 3; jj++){
			//		for (int i = 0; i < 3; i++){
			//			for (int j = 0; j < 3; j++){
			//				matr[i][j] = Am[ii][jj] * getF(i, iCell, x1)*getF(j, iCell + 1, x1);
			//			}
			//		}
			//		addMatr3ToMatr9(matr9, matr, ii, jj);
			//	}
			//}

			////call Solver_AddMatr9(matr, m9, iCell, iCell)
			//S->addMatrElement(iCell, iCell + 1, matr9);


			for (int i = 0; i < MATR_BLOCK; i++) {
				for (int j = 0; j < MATR_BLOCK; j++) {
					matr9[i][j] = 0.0;
				}
			}
			for (int ii = 0; ii < 3; ii++){
				for (int jj = 0; jj < 3; jj++){
					for (int i = 0; i < FUNC_COUNT; i++){
						for (int j = 0; j < FUNC_COUNT; j++){
							matr[i][j] = Ap[ii][jj] * getF(i, iCell, x1)*getF(j, iCell, x1);
						}
					}
					addMatr3ToMatr9(matr9, matr, ii, jj, FUNC_COUNT);
				}
			}
			S->addMatrElement(iCell, iCell, matr9);
			//call Solver_AddMatr9(matr, m9, iCell, iCell + 1)

		} // iCell+1/2



		{ //!< вычисляем потоки на границе iCell-1/2

			for (int i = 0; i < MATR_BLOCK; i++){
				for (int j = 0; j < MATR_BLOCK; j++){
					matr9[i][j] = 0.0;
				}
			}

			//осреднение по Роу
			double x1 = cellC[iCell] + 0.50*h;

			double ROl = getField(0, iCell-1, x1);
			double RUl = getField(1, iCell-1, x1);
			double REl = getField(2, iCell-1, x1);
			double ROr = getField(0, iCell, x1);
			double RUr = getField(1, iCell, x1);
			double REr = getField(2, iCell, x1);
			consToPrim(rb, pb, ub, ROl, RUl, REl);
			consToPrim(re, pe, ue, ROr, RUr, REr);
			FLUX(ri, ei, pi, ui, vi, wi,
				rb, pb, ub, vb, wb,
				re, pe, ue, ve, we, GAM);


			calcMatrAPM(GAM*pi / ri, ui, GAM, Am, Ap);



			for (int i = 0; i < MATR_BLOCK; i++) {
				for (int j = 0; j < MATR_BLOCK; j++){
					matr9[i][j] = 0.0;
				}
			}

			for (int ii = 1; ii < 3; ii++){
				for (int jj = 1; jj < 3; jj++){
					for (int i = 0; i < FUNC_COUNT; i++){
						for (int j = 0; j < FUNC_COUNT; j++){
							matr[i][j] = -Am[ii][jj] * getF(i, iCell, x1)*getF(j, iCell, x1);
						}
					}
					addMatr3ToMatr9(matr9, matr, ii, jj, FUNC_COUNT);
				}
			}

			//call Solver_AddMatr9(matr, m9, iCell, iCell)
			S->addMatrElement(iCell, iCell, matr9);


			for (int i = 0; i < MATR_BLOCK; i++) {
				for (int j = 0; j < MATR_BLOCK; j++) {
					matr9[i][j] = 0.0;
				}
			}

			for (int ii = 0; ii < 3; ii++){
				for (int jj = 0; jj < 3; jj++){
					for (int i = 0; i < FUNC_COUNT; i++){
						for (int j = 0; j < FUNC_COUNT; j++){
							matr[i][j] = -Ap[ii][jj] * getF(i, iCell, x1)*getF(j, iCell - 1, x1);
						}
					}
					addMatr3ToMatr9(matr9, matr, ii, jj, FUNC_COUNT);
				}
			}
			S->addMatrElement(iCell, iCell - 1, matr9);
			//call Solver_AddMatr9(matr, m9, iCell, iCell + 1)


		} // iCell-1/2



	} // (iCell = N - 1)





}

void calcRHS() {
	double  ri, ei, pi, ui, vi, wi,
			rb, pb, ub, vb=0.0, wb=0.0,
			re, pe, ue, ve=0.0, we=0.0,
			alpha, fRO, fRU, fRE;
	double *vect = new double[9];
	for (int iCell = 0; iCell < N; iCell++) {
		memset(vect, 0, MATR_BLOCK*sizeof(double));
		for (int k = 0; k < 2; k++) {
			double x = cellGP[iCell][k];
			double w = cellGW[iCell][k];
			fRO = getField(0, iCell, x);
			fRU = getField(1, iCell, x);
			fRE = getField(2, iCell, x);
			consToPrim(ri, pi, ui, fRO, fRU, fRO);
			fRO = ri*ui;
			fRU = ri*ui*ui + pi;
			fRE = (pi/AGAM+ri*ui*ui*0.5 +pi)*ui;
			int shift = 0;
			for (int j = 0; j < FUNC_COUNT; j++)
				vect[shift++] += w*fRO*getDF(j, iCell, x);
			for (int j = 0; j < FUNC_COUNT; j++)
				vect[shift++] += w*fRU*getDF(j, iCell, x);
			for (int j = 0; j < FUNC_COUNT; j++)
				vect[shift++] += w*fRE*getDF(j, iCell, x);
			//vect[0] += w*fRO*getDF(0, iCell, x);
			//vect[1] += w*fRO*getDF(1, iCell, x); 
			//vect[2] += w*fRO*getDF(2, iCell, x); 
			//vect[3] += w*fRU*getDF(0, iCell, x); 
			//vect[4] += w*fRU*getDF(1, iCell, x); 
			//vect[5] += w*fRU*getDF(2, iCell, x); 
			//vect[6] += w*fRE*getDF(0, iCell, x); 
			//vect[7] += w*fRE*getDF(1, iCell, x); 
			//vect[8] += w*fRE*getDF(2, iCell, x); 
		}
		S->addRightElement(iCell, vect);
	}

	for (int iCell = 1; iCell < N-1; iCell++) {
		{ // iCell+1/2
			double x = cellC[iCell]+0.5*h;
			double ROl = getField(0, iCell, x);
			double RUl = getField(1, iCell, x);
			double REl = getField(2, iCell, x);
			double ROr = getField(0, iCell+1, x);
			double RUr = getField(1, iCell+1, x);
			double REr = getField(2, iCell+1, x);
			consToPrim(rb, pb, ub, ROl, RUl, REl);
			consToPrim(re, pe, ue, ROr, RUr, REr);
			FLUX_RHS(fRO, fRU, fRE,		rb, pb, ub,		re, pe, ue, GAM);
			int shift = 0;
			for (int j = 0; j < FUNC_COUNT; j++)
				vect[shift++] = -fRO*getF(j, iCell, x);
			for (int j = 0; j < FUNC_COUNT; j++)
				vect[shift++] = -fRU*getF(j, iCell, x);
			for (int j = 0; j < FUNC_COUNT; j++)
				vect[shift++] = -fRE*getF(j, iCell, x);
			//vect[0] = -fRO*getF(0, iCell, x);
			//vect[1] = -fRO*getF(1, iCell, x);
			//vect[2] = -fRO*getF(2, iCell, x);
			//vect[3] = -fRU*getF(0, iCell, x);
			//vect[4] = -fRU*getF(1, iCell, x);
			//vect[5] = -fRU*getF(2, iCell, x);
			//vect[6] = -fRE*getF(0, iCell, x);
			//vect[7] = -fRE*getF(1, iCell, x);
			//vect[8] = -fRE*getF(2, iCell, x);
			S->addRightElement(iCell, vect);
		}
		{ // iCell-1/2
			double x = cellC[iCell]-0.5*h;
			double ROl = getField(0, iCell-1, x);
			double RUl = getField(1, iCell-1, x);
			double REl = getField(2, iCell-1, x);
			double ROr = getField(0, iCell, x);
			double RUr = getField(1, iCell, x);
			double REr = getField(2, iCell, x);
			consToPrim(rb, pb, ub, ROl, RUl, REl);
			consToPrim(re, pe, ue, ROr, RUr, REr);
			FLUX_RHS(fRO, fRU, fRE, rb, pb, ub, re, pe, ue, GAM);
			int shift = 0;
			for (int j = 0; j < FUNC_COUNT; j++)
				vect[shift++] = fRO*getF(j, iCell, x);
			for (int j = 0; j < FUNC_COUNT; j++)
				vect[shift++] = fRU*getF(j, iCell, x);
			for (int j = 0; j < FUNC_COUNT; j++)
				vect[shift++] = fRE*getF(j, iCell, x);
			//vect[0] = fRO*getF(0, iCell, x);
			//vect[1] = fRO*getF(1, iCell, x);
			//vect[2] = fRO*getF(2, iCell, x);
			//vect[3] = fRU*getF(0, iCell, x);
			//vect[4] = fRU*getF(1, iCell, x);
			//vect[5] = fRU*getF(2, iCell, x);
			//vect[6] = fRE*getF(0, iCell, x);
			//vect[7] = fRE*getF(1, iCell, x);
			//vect[8] = fRE*getF(2, iCell, x);
			S->addRightElement(iCell, vect);
		}
	}
	
	{ 
		int iCell = 0;
		{ // iCell+1/2
			double x = cellC[iCell]+0.5*h;
			double ROl = getField(0, iCell, x);
			double RUl = getField(1, iCell, x);
			double REl = getField(2, iCell, x);
			double ROr = getField(0, iCell+1, x);
			double RUr = getField(1, iCell+1, x);
			double REr = getField(2, iCell+1, x);
			consToPrim(rb, pb, ub, ROl, RUl, REl);
			consToPrim(re, pe, ue, ROr, RUr, REr);
			FLUX_RHS(fRO, fRU, fRE, rb, pb, ub, re, pe, ue, GAM);
			int shift = 0;
			for (int j = 0; j < FUNC_COUNT; j++)
				vect[shift++] = -fRO*getF(j, iCell, x);
			for (int j = 0; j < FUNC_COUNT; j++)
				vect[shift++] = -fRU*getF(j, iCell, x);
			for (int j = 0; j < FUNC_COUNT; j++)
				vect[shift++] = -fRE*getF(j, iCell, x);
			//vect[0] = -fRO*getF(0, iCell, x);
			//vect[1] = -fRO*getF(1, iCell, x);
			//vect[2] = -fRO*getF(2, iCell, x);
			//vect[3] = -fRU*getF(0, iCell, x);
			//vect[4] = -fRU*getF(1, iCell, x);
			//vect[5] = -fRU*getF(2, iCell, x);
			//vect[6] = -fRE*getF(0, iCell, x);
			//vect[7] = -fRE*getF(1, iCell, x);
			//vect[8] = -fRE*getF(2, iCell, x);
			S->addRightElement(iCell, vect);
		}
		{ // iCell-1/2
			double x = cellC[iCell]-0.5*h;
			//double ROl = getField(0, iCell-1, x);
			//double RUl = getField(1, iCell-1, x);
			//double REl = getField(2, iCell-1, x);
			double ROr = getField(0, iCell, x);
			double RUr = getField(1, iCell, x);
			double REr = getField(2, iCell, x);
			//consToPrim(rb, pb, ub, ROl, RUl, REl);
			consToPrim(re, pe, ue, ROr, RUr, REr);
			rb = re; pb = pe; ub = ue;
			FLUX_RHS(fRO, fRU, fRE, rb, pb, ub, re, pe, ue, GAM);
			int shift = 0;
			for (int j = 0; j < FUNC_COUNT; j++)
				vect[shift++] = fRO*getF(j, iCell, x);
			for (int j = 0; j < FUNC_COUNT; j++)
				vect[shift++] = fRU*getF(j, iCell, x);
			for (int j = 0; j < FUNC_COUNT; j++)
				vect[shift++] = fRE*getF(j, iCell, x);
			//vect[0] = fRO*getF(0, iCell, x);
			//vect[1] = fRO*getF(1, iCell, x);
			//vect[2] = fRO*getF(2, iCell, x);
			//vect[3] = fRU*getF(0, iCell, x);
			//vect[4] = fRU*getF(1, iCell, x);
			//vect[5] = fRU*getF(2, iCell, x);
			//vect[6] = fRE*getF(0, iCell, x);
			//vect[7] = fRE*getF(1, iCell, x);
			//vect[8] = fRE*getF(2, iCell, x);
			S->addRightElement(iCell, vect);
		}
	}

	{ 
		int iCell = N-1;
		{ // iCell+1/2
			double x = cellC[iCell]+0.5*h;
			double ROl = getField(0, iCell, x);
			double RUl = getField(1, iCell, x);
			double REl = getField(2, iCell, x);
			//double ROr = getField(0, iCell+1, x);
			//double RUr = getField(1, iCell+1, x);
			//double REr = getField(2, iCell+1, x);
			consToPrim(rb, pb, ub, ROl, RUl, REl);
			//consToPrim(re, pe, ue, ROr, RUr, REr);
			re = rb; pe = pb; ue = ub;
			FLUX_RHS(fRO, fRU, fRE, rb, pb, ub, re, pe, ue, GAM);
			int shift = 0;
			for (int j = 0; j < FUNC_COUNT; j++)
				vect[shift++] = -fRO*getF(j, iCell, x);
			for (int j = 0; j < FUNC_COUNT; j++)
				vect[shift++] = -fRU*getF(j, iCell, x);
			for (int j = 0; j < FUNC_COUNT; j++)
				vect[shift++] = -fRE*getF(j, iCell, x);
			//vect[0] = -fRO*getF(0, iCell, x);
			//vect[1] = -fRO*getF(1, iCell, x);
			//vect[2] = -fRO*getF(2, iCell, x);
			//vect[3] = -fRU*getF(0, iCell, x);
			//vect[4] = -fRU*getF(1, iCell, x);
			//vect[5] = -fRU*getF(2, iCell, x);
			//vect[6] = -fRE*getF(0, iCell, x);
			//vect[7] = -fRE*getF(1, iCell, x);
			//vect[8] = -fRE*getF(2, iCell, x);
			S->addRightElement(iCell, vect);
		}
		{ // iCell-1/2
			double x = cellC[iCell]-0.5*h;
			double ROl = getField(0, iCell-1, x);
			double RUl = getField(1, iCell-1, x);
			double REl = getField(2, iCell-1, x);
			double ROr = getField(0, iCell, x);
			double RUr = getField(1, iCell, x);
			double REr = getField(2, iCell, x);
			consToPrim(rb, pb, ub, ROl, RUl, REl);
			consToPrim(re, pe, ue, ROr, RUr, REr);
			FLUX_RHS(fRO, fRU, fRE, rb, pb, ub, re, pe, ue, GAM);
			int shift = 0;
			for (int j = 0; j < FUNC_COUNT; j++)
				vect[shift++] = fRO*getF(j, iCell, x);
			for (int j = 0; j < FUNC_COUNT; j++)
				vect[shift++] = fRU*getF(j, iCell, x);
			for (int j = 0; j < FUNC_COUNT; j++)
				vect[shift++] = fRE*getF(j, iCell, x);
			//vect[0] = fRO*getF(0, iCell, x);
			//vect[1] = fRO*getF(1, iCell, x);
			//vect[2] = fRO*getF(2, iCell, x);
			//vect[3] = fRU*getF(0, iCell, x);
			//vect[4] = fRU*getF(1, iCell, x);
			//vect[5] = fRU*getF(2, iCell, x);
			//vect[6] = fRE*getF(0, iCell, x);
			//vect[7] = fRE*getF(1, iCell, x);
			//vect[8] = fRE*getF(2, iCell, x);
			S->addRightElement(iCell, vect);
		}
	}
	delete[] vect;
}

inline double _sign_(double a) {
	if (a < 0.0) {
		return -1.0;
	}
	else if (a > 0.0) {
		return 1.0;
	}
	else {
		return 0.0;
	}
}

inline double _min_(double a, double b) {
	if (a < b) return a;
	return b;
}

double minmod(double a, double b, double c) {
	
	if ((_sign_(a) == _sign_(b)) && (_sign_(c) == _sign_(b))) {
		return _sign_(a)*_min_(_min_(fabs(a), fabs(b)), fabs(c));
	}
	return 0.0;
}

void calcLimiterCons() {
	{
		int iCell = 0;
		// проецируем на линейный базис
		double u0, u1, u1l;
		u0 = ro[iCell][0] + ro[iCell][2] / 12.0;
		u1 = ro[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(ro[iCell + 1][0] + ro[iCell + 1][2] / 12.0 - u0),
			LIM_ALPHA*(u0 - ro[iCell][0] - ro[iCell][2] / 12.0));
		if (u1l != u1) {
			ro[iCell][0] = u0;
			ro[iCell][1] = u1l;
			ro[iCell][2] = 0.0;
		}

		u0 = ru[iCell][0] + ru[iCell][2] / 12.0;
		u1 = ru[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(ru[iCell + 1][0] + ru[iCell + 1][2] / 12.0 - u0),
			LIM_ALPHA*(u0 - ru[iCell][0] - ru[iCell][2] / 12.0));
		if (u1l != u1) {
			ru[iCell][0] = u0;
			ru[iCell][1] = u1l;
			ru[iCell][2] = 0.0;
		}

		u0 = re[iCell][0] + re[iCell][2] / 12.0;
		u1 = re[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(re[iCell + 1][0] + re[iCell + 1][2] / 12.0 - u0),
			LIM_ALPHA*(u0 - re[iCell][0] - re[iCell][2] / 12.0));
		if (u1l != u1) {
			re[iCell][0] = u0;
			re[iCell][1] = u1l;
			re[iCell][2] = 0.0;
		}

	}
	for (int iCell = 1; iCell < N - 1; iCell++) {
		// проецируем на линейный базис
		double u0, u1, u1l;
		u0 = ro[iCell][0] + ro[iCell][2] / 12.0;
		u1 = ro[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(ro[iCell + 1][0] + ro[iCell + 1][2] / 12.0 - u0),
			LIM_ALPHA*(u0 - ro[iCell - 1][0] - ro[iCell - 1][2] / 12.0));
		if (u1l != u1) {
			ro[iCell][0] = u0;
			ro[iCell][1] = u1l;
			ro[iCell][2] = 0.0;
		}

		u0 = ru[iCell][0] + ru[iCell][2] / 12.0;
		u1 = ru[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(ru[iCell + 1][0] + ru[iCell + 1][2] / 12.0 - u0),
			LIM_ALPHA*(u0 - ru[iCell - 1][0] - ru[iCell - 1][2] / 12.0));
		if (u1l != u1) {
			ru[iCell][0] = u0;
			ru[iCell][1] = u1l;
			ru[iCell][2] = 0.0;
		}

		u0 = re[iCell][0] + re[iCell][2] / 12.0;
		u1 = re[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(re[iCell + 1][0] + re[iCell + 1][2] / 12.0 - u0),
			LIM_ALPHA*(u0 - re[iCell - 1][0] - re[iCell - 1][2] / 12.0));
		if (u1l != u1) {
			re[iCell][0] = u0;
			re[iCell][1] = u1l;
			re[iCell][2] = 0.0;
		}

	}
	{
		int iCell = N - 1;
		// проецируем на линейный базис
		double u0, u1, u1l;
		u0 = ro[iCell][0] + ro[iCell][2] / 12.0;
		u1 = ro[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(ro[iCell][0] + ro[iCell][2] / 12.0 - u0),
			LIM_ALPHA*(u0 - ro[iCell - 1][0] - ro[iCell - 1][2] / 12.0));
		if (u1l != u1) {
			ro[iCell][0] = u0;
			ro[iCell][1] = u1l;
			ro[iCell][2] = 0.0;
		}

		u0 = ru[iCell][0] + ru[iCell][2] / 12.0;
		u1 = ru[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(ru[iCell][0] + ru[iCell][2] / 12.0 - u0),
			LIM_ALPHA*(u0 - ru[iCell - 1][0] - ru[iCell - 1][2] / 12.0));
		if (u1l != u1) {
			ru[iCell][0] = u0;
			ru[iCell][1] = u1l;
			ru[iCell][2] = 0.0;
		}

		u0 = re[iCell][0] + re[iCell][2] / 12.0;
		u1 = re[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(re[iCell][0] + re[iCell][2] / 12.0 - u0),
			LIM_ALPHA*(u0 - re[iCell - 1][0] - re[iCell - 1][2] / 12.0));
		if (u1l != u1) {
			re[iCell][0] = u0;
			re[iCell][1] = u1l;
			re[iCell][2] = 0.0;
		}

	}
}

void calcLimiterEigenv() {
	double c2_, u_;
	for (int i = 0; i < N; i++) {
		u_ = (ru[i][0] + ((FUNC_COUNT < 3) ? 0.0 : ru[i][2] / 12.0)) / (ro[i][0] + ((FUNC_COUNT < 3) ? 0.0 : ro[i][2] / 12.0));
		c2_ = ((re[i][0] + ((FUNC_COUNT < 3) ? 0.0 : re[i][2] / 12.0)) / (ro[i][0] + ((FUNC_COUNT < 3) ? 0.0 : ro[i][2] / 12.0)) - u_*u_*0.5)*GAM*AGAM;
		calcMatrL(c2_, u_, GAM, L);
		for (int j = 0; j < FUNC_COUNT; j++) {
			ro_[i][j] = L[0][0] * ro[i][j] + L[0][1] * ru[i][j] + L[0][2] * re[i][j];
			ru_[i][j] = L[1][0] * ro[i][j] + L[1][1] * ru[i][j] + L[1][2] * re[i][j];
			re_[i][j] = L[2][0] * ro[i][j] + L[2][1] * ru[i][j] + L[2][2] * re[i][j];
		}
	}


	{
		int iCell = 0;
		// проецируем на линейный базис
		double u0, u1, u1l;
		u0 = ro_[iCell][0] + ((FUNC_COUNT < 3) ? 0.0 : ro_[iCell][2] / 12.0);
		u1 = ro_[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(ro_[iCell + 1][0] + ((FUNC_COUNT < 3) ? 0.0 : ro_[iCell + 1][2] / 12.0) - u0),
			LIM_ALPHA*(u0 - ro_[iCell][0] - ((FUNC_COUNT < 3) ? 0.0 : ro_[iCell][2] / 12.0)));
		if (u1l != u1) {
			ro_[iCell][0] = u0;
			ro_[iCell][1] = u1l;
			if (FUNC_COUNT > 2) ro_[iCell][2] = 0.0;
		}

		u0 = ru_[iCell][0] + ((FUNC_COUNT < 3) ? 0.0 : ru_[iCell][2] / 12.0);
		u1 = ru_[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(ru_[iCell + 1][0] + ((FUNC_COUNT < 3) ? 0.0 : ru_[iCell + 1][2] / 12.0) - u0),
			LIM_ALPHA*(u0 - ru_[iCell][0] - ((FUNC_COUNT < 3) ? 0.0 : ru_[iCell][2] / 12.0)));
		if (u1l != u1) {
			ru_[iCell][0] = u0;
			ru_[iCell][1] = u1l;
			if (FUNC_COUNT > 2) ru_[iCell][2] = 0.0;
		}

		u0 = re_[iCell][0] + ((FUNC_COUNT < 3) ? 0.0 : re_[iCell][2] / 12.0);
		u1 = re_[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(re_[iCell + 1][0] + ((FUNC_COUNT < 3) ? 0.0 : re_[iCell + 1][2] / 12.0) - u0),
			LIM_ALPHA*(u0 - re_[iCell][0] - ((FUNC_COUNT < 3) ? 0.0 : re_[iCell][2] / 12.0)));
		if (u1l != u1) {
			re_[iCell][0] = u0;
			re_[iCell][1] = u1l;
			if (FUNC_COUNT > 2) re_[iCell][2] = 0.0;
		}

	}
	for (int iCell = 1; iCell < N - 1; iCell++) {
		// проецируем на линейный базис
		double u0, u1, u1l;
		u0 = ro_[iCell][0] + ((FUNC_COUNT < 3) ? 0.0 : ro_[iCell][2] / 12.0);
		u1 = ro_[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(ro_[iCell + 1][0] + ((FUNC_COUNT < 3) ? 0.0 : ro_[iCell + 1][2] / 12.0) - u0),
			LIM_ALPHA*(u0 - ro_[iCell - 1][0] - ((FUNC_COUNT < 3) ? 0.0 : ro_[iCell - 1][2] / 12.0)));
		if (u1l != u1) {
			ro_[iCell][0] = u0;
			ro_[iCell][1] = u1l;
			if (FUNC_COUNT > 2) ro_[iCell][2] = 0.0;
		}

		u0 = ru_[iCell][0] + ((FUNC_COUNT < 3) ? 0.0 : ru_[iCell][2] / 12.0);
		u1 = ru_[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(ru_[iCell + 1][0] + ((FUNC_COUNT < 3) ? 0.0 : ru_[iCell + 1][2] / 12.0) - u0),
			LIM_ALPHA*(u0 - ru_[iCell - 1][0] - ((FUNC_COUNT < 3) ? 0.0 : ru_[iCell - 1][2] / 12.0)));
		if (u1l != u1) {
			ru_[iCell][0] = u0;
			ru_[iCell][1] = u1l;
			if (FUNC_COUNT > 2) ru_[iCell][2] = 0.0;
		}

		u0 = re_[iCell][0] + ((FUNC_COUNT < 3) ? 0.0 : re_[iCell][2] / 12.0);
		u1 = re_[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(re_[iCell + 1][0] + ((FUNC_COUNT < 3) ? 0.0 : re_[iCell + 1][2] / 12.0) - u0),
			LIM_ALPHA*(u0 - re_[iCell - 1][0] - ((FUNC_COUNT < 3) ? 0.0 : re_[iCell - 1][2] / 12.0)));
		if (u1l != u1) {
			re_[iCell][0] = u0;
			re_[iCell][1] = u1l;
			if (FUNC_COUNT > 2) re_[iCell][2] = 0.0;
		}

	}
	{
		int iCell = N - 1;
		// проецируем на линейный базис
		double u0, u1, u1l;
		u0 = ro_[iCell][0] + ((FUNC_COUNT < 3) ? 0.0 : ro_[iCell][2] / 12.0);
		u1 = ro_[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(ro_[iCell][0] + ((FUNC_COUNT < 3) ? 0.0 : ro_[iCell][2] / 12.0) - u0),
			LIM_ALPHA*(u0 - ro_[iCell - 1][0] - ((FUNC_COUNT < 3) ? 0.0 : ro_[iCell - 1][2] / 12.0)));
		if (u1l != u1) {
			ro_[iCell][0] = u0;
			ro_[iCell][1] = u1l;
			if (FUNC_COUNT > 2) ro_[iCell][2] = 0.0;
		}

		u0 = ru_[iCell][0] + ((FUNC_COUNT < 3) ? 0.0 : ru_[iCell][2] / 12.0);
		u1 = ru_[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(ru_[iCell][0] + ((FUNC_COUNT < 3) ? 0.0 : ru_[iCell][2] / 12.0) - u0),
			LIM_ALPHA*(u0 - ru_[iCell - 1][0] - ((FUNC_COUNT < 3) ? 0.0 : ru_[iCell - 1][2] / 12.0)));
		if (u1l != u1) {
			ru_[iCell][0] = u0;
			ru_[iCell][1] = u1l;
			if (FUNC_COUNT > 2) ru_[iCell][2] = 0.0;
		}

		u0 = re_[iCell][0] + ((FUNC_COUNT < 3) ? 0.0 : re_[iCell][2] / 12.0);
		u1 = re_[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(re_[iCell][0] + ((FUNC_COUNT < 3) ? 0.0 : re_[iCell][2] / 12.0) - u0),
			LIM_ALPHA*(u0 - re_[iCell - 1][0] - ((FUNC_COUNT < 3) ? 0.0 : (re_[iCell - 1][2] / 12.0))));
		if (u1l != u1) {
			re_[iCell][0] = u0;
			re_[iCell][1] = u1l;
			if (FUNC_COUNT > 2) re_[iCell][2] = 0.0;
		}

	}

	for (int i = 0; i < N; i++) {
		u_ = (ru[i][0] + ((FUNC_COUNT < 3) ? 0.0 : ru[i][2] / 12.0)) / (ro[i][0] + ((FUNC_COUNT < 3) ? 0.0 : (ro[i][2] / 12.0)));
		c2_ = ((re[i][0] + ((FUNC_COUNT < 3) ? 0.0 : re[i][2] / 12.0)) / (ro[i][0] + ((FUNC_COUNT < 3) ? 0.0 : (ro[i][2] / 12.0))) - u_*u_*0.5)*GAM*AGAM;
		calcMatrR(c2_, u_, GAM, R);
		for (int j = 0; j < FUNC_COUNT; j++) {
			ro[i][j] = R[0][0] * ro_[i][j] + R[0][1] * ru_[i][j] + R[0][2] * re_[i][j];
			ru[i][j] = R[1][0] * ro_[i][j] + R[1][1] * ru_[i][j] + R[1][2] * re_[i][j];
			re[i][j] = R[2][0] * ro_[i][j] + R[2][1] * ru_[i][j] + R[2][2] * re_[i][j];
		}
	}
}

void calcLimiter_II() {
	double x[5];
	for (int i = 0; i < N; i++) {
		x[0] = cellC[i];
		x[1] = cellGP[i][0];
		x[2] = cellGP[i][1];
		x[3] = cellC[i] - 0.5*h;
		x[4] = cellC[i] + 0.5*h;
		for (int k = 0; k < 5; k++) {
			double fRO = getField(0, i, x[k]);
			double fRU = getField(1, i, x[k]);
			double fRE = getField(2, i, x[k]);
			double r = fRO;
			double p = (fRE - 0.5*fRU*fRU / fRO)*AGAM;
			if ((r < EPS) || (p < EPS)) {
				ro[i][0] += ro[i][2] / 12.0;
				ro[i][1] = 0.0;
				if (FUNC_COUNT > 2) ro[i][2] = 0.0;
				ru[i][0] += ru[i][2] / 12.0;
				ru[i][1] = 0.0;
				if (FUNC_COUNT > 2) ru[i][2] = 0.0;
				re[i][0] += re[i][2] / 12.0;
				re[i][1] = 0.0;
				if (FUNC_COUNT > 2) re[i][2] = 0.0;
				printf("WARNING: bad cell #%d\n", i);
				break;
			}
			
		}
	}
}

#define _pow_2_(x) ((x)*(x))

void calcSmoother() {
	//копируем старые значения
	for (int i = 0; i < N; i++) {
		memcpy(ro_[i], ro[i], 3 * sizeof(double));
		memcpy(ru_[i], ru[i], 3 * sizeof(double));
		memcpy(re_[i], re[i], 3 * sizeof(double));
	}
	for (int idField = 0; idField < 3; idField++) {
		double** fld, **fldOld;
		switch (idField) {
		case 0:
			fld		= ro_;
			fldOld	= ro;
			break;
		case 1:
			fld		= ru_;
			fldOld	= ru;
			break;
		case 2:
			fld		= re_;
			fldOld	= re;
			break;
		}

		// i+1/2
		for (int iCell = 0; iCell < N-1; iCell++) {
			int c1 = iCell;
			int c2 = iCell + 1;

			double fIntP1P2 = 0;
			double fIntP1Pavg = 0;
			double fIntP2Pavg = 0;
			double fInt1, fInt2;

			// Integral[ (P1-P2)^2 ]d(V1+V2)
			fInt1 = 0.;
			fInt2 = 0.;
			for (int iGP = 0; iGP < 2; ++iGP) // по точкам Гаусса
			{
				fInt1 += cellGW[c1][iGP] * _pow_2_(getField(idField, c1, cellGP[c1][iGP]) - getField(idField, c2, cellGP[c1][iGP]));
				fInt2 += cellGW[c2][iGP] * _pow_2_(getField(idField, c1, cellGP[c2][iGP]) - getField(idField, c2, cellGP[c2][iGP]));
			}
			fIntP1P2 = fInt1 + fInt2;

			if (fabs(fIntP1P2) < EPS)
			{
				smMu[iCell] = 0.0;
				continue;
			}

			// Integral[ (P1-Pavg)^2 ]d(V1+V2)
			fInt1 = 0.;
			fInt2 = 0.;
			for (int iGP = 0; iGP < 2; ++iGP) // по точкам Гаусса
			{
				fInt1 += cellGW[c1][iGP] * _pow_2_(getField(idField, c1, cellGP[c1][iGP]) - getFieldsAvg(idField, iCell, cellGP[c1][iGP]));
				fInt2 += cellGW[c2][iGP] * _pow_2_(getField(idField, c1, cellGP[c2][iGP]) - getFieldsAvg(idField, iCell, cellGP[c2][iGP]));
			}
			fIntP1Pavg = fInt1 + fInt2;

			// Integral[ (P2-Pavg)^2 ]d(V1+V2)
			fInt1 = 0.;
			fInt2 = 0.;
			for (int iGP = 0; iGP < 2; ++iGP) // по точкам Гаусса
			{
				fInt1 += cellGW[c1][iGP] * _pow_2_(getField(idField, c2, cellGP[c1][iGP]) - getFieldsAvg(idField, iCell, cellGP[c1][iGP]));
				fInt2 += cellGW[c2][iGP] * _pow_2_(getField(idField, c2, cellGP[c2][iGP]) - getFieldsAvg(idField, iCell, cellGP[c2][iGP]));
			}
			fIntP2Pavg = fInt1 + fInt2;
			double ROl = getField(0, iCell, cellC[iCell]);
			double RUl = getField(1, iCell, cellC[iCell]);
			double REl = getField(2, iCell, cellC[iCell]);
			double ROr = getField(0, iCell+1, cellC[iCell]);
			double RUr = getField(1, iCell+1, cellC[iCell]);
			double REr = getField(2, iCell+1, cellC[iCell]);
			double RI, EI, PI, UI, VI, WI;
			roe_orig(RI, EI, PI, UI, VI, WI,
				ROl, (REl - 0.5*RUl*RUl / ROl)*(GAM - 1.0), RUl / ROl, 0, 0,
				ROr, (REl - 0.5*RUl*RUl / ROl)*(GAM - 1.0), RUr / ROr, 0, 0, GAM);
			double CFL = tau*(fabs(UI)+sqrt(GAM*PI/RI)) / h;
			smMu[iCell] = CFL*2.0 * fIntP1P2 / (fIntP1Pavg + fIntP2Pavg);
		}

		for (int i = 1; i < N-1; i++) {
			if (smMu[i] + smMu[i - 1] > 1.0) {
				printf("WARNING: bad cell #%d for smoothing...\n", i);
				continue;
			}
			for (int j = 0; j < 3; j++) {
				fld[i][j] = ((1.0 - smMu[i] - smMu[i - 1])*fldOld[i][j] + smMu[i - 1] * fldOld[i - 1][j] + smMu[i] * fldOld[i + 1][j]);
			}
		}

	}

	for (int i = 0; i < N; i++) {
		memcpy(ro[i], ro_[i], 3 * sizeof(double));
		memcpy(ru[i], ru_[i], 3 * sizeof(double));
		memcpy(re[i], re_[i], 3 * sizeof(double));
	}

}

#undef _pow_2_


double getF(int i, int iCell, double x) {
	if (i == 0) {
		return 1.0;
	}
	else if (i == 1) {
		return (x - cellC[iCell]) / h;
	}
	else if (i == 2) {
		return ((x - cellC[iCell]) / h)*((x - cellC[iCell]) / h);
	}
	printf("ERROR: bad getF() argument!");
	return 0.0;

} // getF

double getDF(int i, int iCell, double x) {
	if (i == 0) {
		return 0.0;
	}
	else if (i == 1) {
		return 1.0 / h;
	}
	else if (i == 2) {
		return 2.0*(x - cellC[iCell]) / h / h;
	}
	printf("ERROR: bad getF() argument!");
	return 0.0;

} // getF


//double getF(int i, int iCell, double x) {
//	if (i == 0) {
//		return 1.0;
//	}
//	else if (i == 1) {
//		return (x - cellC[iCell]);
//	}
//	else if (i == 2) {
//		return ((x - cellC[iCell]))*((x - cellC[iCell]));
//	}
//	printf("ERROR: bad getF() argument!");
//	return 0.0;
//
//} // getF
//
//double getDF(int i, int iCell, double x) {
//	if (i == 0) {
//		return 0.0;
//	}
//	else if (i == 1) {
//		return 1.0;
//	}
//	else if (i == 2) {
//		return 2.0*(x - cellC[iCell]);
//	}
//	printf("ERROR: bad getF() argument!");
//	return 0.0;
//
//} // getF
//



double getField(int i, int iCell, double x){

	double f1 = getF(1, iCell, x);
	double f2 = getF(2, iCell, x);
	switch (i) {
	case 0:
		return ro[iCell][0] + ro[iCell][1] * f1 + ((FUNC_COUNT < 3) ? 0.0 : ro[iCell][2] * f2);
		break;
	case 1:
		return ru[iCell][0] + ru[iCell][1] * f1 +( (FUNC_COUNT < 3) ? 0.0 : ru[iCell][2] * f2);
		break;
	case 2:
		return re[iCell][0] + re[iCell][1] * f1 + ((FUNC_COUNT < 3) ? 0.0 : re[iCell][2] * f2);
		break;
	}
}

/**
 *	Average on iCell+1/2
 */
double getFieldsAvg(int idField, int iCell, double x) {
	if (iCell > -1 && iCell < N - 1) {
		return 0.5*(getField(idField, iCell, x)+getField(idField, iCell + 1, x));
	}
	else {
		if (iCell = -1) {
			return getField(idField, iCell+1, x);
		}
		else {
			return getField(idField, iCell, x);
		}
	}
}




























//////////////////////////////////////////////////////////////////////////////////////////////


/*
for (int iCell = 1; iCell < N; iCell++) {
double ri,ei,pi,ui,vi,wi;
double rb = iCell>1 ? r[iCell-1]+0.5*minmod(r[iCell]-r[iCell-1], r[iCell-1]-r[iCell-2]) : r[iCell-1];
double pb = iCell>1 ? p[iCell-1]+0.5*minmod(p[iCell]-p[iCell-1], p[iCell-1]-p[iCell-2]) : p[iCell-1];
double ub = iCell>1 ? u[iCell-1]+0.5*minmod(u[iCell]-u[iCell-1], u[iCell-1]-u[iCell-2]) : u[iCell-1];
double vb = 0;
double wb = 0;
double re = iCell<N-1 ? r[iCell]-0.5*minmod(r[iCell+1]-r[iCell], r[iCell]-r[iCell-1]) : r[iCell];
double pe = iCell<N-1 ? p[iCell]-0.5*minmod(p[iCell+1]-p[iCell], p[iCell]-p[iCell-1]) : p[iCell];
double ue = iCell<N-1 ? u[iCell]-0.5*minmod(u[iCell+1]-u[iCell], u[iCell]-u[iCell-1]) : u[iCell];
double ve = 0;
double we = 0;
roe_orig(ri,ei,pi,ui,vi,wi,
rb, pb, ub, vb, wb,
re, pe, ue, ve, we, GAM);
double c2 = GAM*pi/ri;
double alpha = _max_(fabs(ub)+sqrt(GAM*pb/rb), fabs(ue)+sqrt(GAM*pe/re));

calcMatrAPM(c2, ui,  Am,  Ap);

// i-1/2
for (int i = 0; i < 3; i++) {
for (int j = 0; j < 3; j++) {
matr[i][j] = -Ap[i][j];
}
}
S->addMatrElement(iCell, iCell-1, matr);

for (int i = 0; i < 3; i++) {
for (int j = 0; j < 3; j++) {
matr[i][j] = -Am[i][j];
}
}
S->addMatrElement(iCell, iCell, matr);

// i+1/2
for (int i = 0; i < 3; i++) {
for (int j = 0; j < 3; j++) {
matr[i][j] = Am[i][j];
}
}
S->addMatrElement(iCell-1, iCell, matr);

for (int i = 0; i < 3; i++) {
for (int j = 0; j < 3; j++) {
matr[i][j] = Ap[i][j];
}
}
S->addMatrElement(iCell-1, iCell-1, matr);

vect[0] = 0.5*(re*ue+rb*ub-alpha*(re-rb));
vect[1] = 0.5*(re*ue*ue+pe+rb*ub*ub+pb-alpha*(re*ue-rb*ub));
vect[2] = 0.5*(((pe/AGAM+re*ue*ue*0.5)+pe)*ue+((pb/AGAM+rb*ub*ub*0.5)+pb)*ub-alpha*((pe/AGAM+re*ue*ue*0.5)-(pb/AGAM+rb*ub*ub*0.5)));
S->addRightElement(iCell, vect);

vect[0] *= -1.0;
vect[1] *= -1.0;
vect[2] *= -1.0;
S->addRightElement(iCell-1, vect);


//rim_orig(ri,ei,pi,ui,vi,wi,
//		 rb, pb, ub, vb, wb,
//		 re, pe, ue, ve, we, GAM);
//
//vect[0] = ri*ui;
//vect[1] = ri*ui*ui+pi;
//vect[2] = ((pi/AGAM+ri*ui*ui*0.5)+pi)*ui;
//S->addRightElement(iCell, vect);

//vect[0] *= -1.0;
//vect[1] *= -1.0;
//vect[2] *= -1.0;
//S->addRightElement(iCell-1, vect);
}

{
int iCell = 0;

double ri,ei,pi,ui,vi,wi;
double rb = r[iCell];
double pb = p[iCell];
double ub = u[iCell];
double vb = 0;
double wb = 0;
double re = r[iCell];
double pe = p[iCell];
double ue = u[iCell];
double ve = 0;
double we = 0;
roe_orig(ri,ei,pi,ui,vi,wi,
rb, pb, ub, vb, wb,
re, pe, ue, ve, we, GAM);
double c2 = GAM*pi/ri;
double alpha = _max_(fabs(ub)+sqrt(GAM*pb/rb), fabs(ue)+sqrt(GAM*pe/re));

calcMatrAPM(c2, ui, Am,  Ap);

for (int i = 0; i < 3; i++) {
for (int j = 0; j < 3; j++) {
matr[i][j] = -Am[i][j];
}
}
S->addMatrElement(iCell, iCell, matr);

vect[0] = 0.5*(re*ue+rb*ub-alpha*(re-rb));
vect[1] = 0.5*(re*ue*ue+pe+rb*ub*ub+pb-alpha*(re*ue-rb*ub));
vect[2] = 0.5*(((pe/AGAM+re*ue*ue*0.5)+pe)*ue+((pb/AGAM+rb*ub*ub*0.5)+pb)*ub-alpha*((pe/AGAM+re*ue*ue*0.5)-(pb/AGAM+rb*ub*ub*0.5)));
S->addRightElement(iCell, vect);

//rim_orig(ri,ei,pi,ui,vi,wi,
//		 rb, pb, ub, vb, wb,
//		 re, pe, ue, ve, we, GAM);
//
//vect[0] = ri*ui;
//vect[1] = ri*ui*ui+pi;
//vect[2] = ((pi/AGAM+ri*ui*ui*0.5)+pi)*ui;
//S->addRightElement(iCell, vect);

}

{
int iCell = N;

double ri,ei,pi,ui,vi,wi;
double rb = r[iCell-1];
double pb = p[iCell-1];
double ub = u[iCell-1];
double vb = 0;
double wb = 0;
double re = r[iCell-1];
double pe = p[iCell-1];
double ue = u[iCell-1];
double ve = 0;
double we = 0;
roe_orig(ri,ei,pi,ui,vi,wi,
rb, pb, ub, vb, wb,
re, pe, ue, ve, we, GAM);
double c2 = GAM*pi/ri;
double alpha = _max_(fabs(ub)+sqrt(GAM*pb/rb), fabs(ue)+sqrt(GAM*pe/re));

calcMatrAPM(c2, ui, Am, Ap);

for (int i = 0; i < 3; i++) {
for (int j = 0; j < 3; j++) {
matr[i][j] = Ap[i][j];
}
}
S->addMatrElement(iCell-1, iCell-1, matr);

vect[0] = -0.5*(re*ue+rb*ub-alpha*(re-rb));
vect[1] = -0.5*(re*ue*ue+pe+rb*ub*ub+pb-alpha*(re*ue-rb*ub));
vect[2] = -0.5*(((pe/AGAM+re*ue*ue*0.5)+pe)*ue+((pb/AGAM+rb*ub*ub*0.5)+pb)*ub-alpha*((pe/AGAM+re*ue*ue*0.5)-(pb/AGAM+rb*ub*ub*0.5)));
S->addRightElement(iCell-1, vect);

//rim_orig(ri,ei,pi,ui,vi,wi,
//		 rb, pb, ub, vb, wb,
//		 re, pe, ue, ve, we, GAM);
//
//vect[0] = -ri*ui;
//vect[1] = -ri*ui*ui-pi;
//vect[2] = -((pi/AGAM+ri*ui*ui*0.5)+pi)*ui;
//S->addRightElement(iCell-1, vect);

}

for (int iCell = 0; iCell < N; iCell++) {
for (int i = 0; i < 3; i++) {
for (int j = 0; j < 3; j++) {
matr[i][j] = (i==j ? h/tau : 0.0);
}
}
S->addMatrElement(iCell, iCell, matr);
}
*/
