#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <ctime>
#include <cstdarg>


const int		FUNC_COUNT	=	2;
const int		MATR_BLOCK	=	3 * FUNC_COUNT;

const double	CFL			=	1.e-2;

const double	LIM_ALPHA	=	2.0;

const int		N			=	100;
const double	XMIN		=	-1.0;
const double	XMAX		=	1.0;
const double	EPS			=	1.0e-5;
const double	GAM			=	5.0 / 3.0;
const double	AGAM		=	GAM - 1.0;
double			TMAX		=	0.07;

const int		MAX_ITER	=	7000;
const int		SAVE_STEP	=	1000;
const int		PRINT_STEP	=	1;

const double	h			=	(XMAX - XMIN) / N;
double			tau			=	1.e-5;//CFL*h; // <= 1.e-5
const double	SIGMA		=	1.e-6;//2.5e-2 * h * h;

double **ro, **ru, **re;
double **ro_, **ru_, **re_;
double *r, *u, *e, *p;
double **dr, **du, **de;
double **matr, *vect, **matr1, **matr9, ***matrM, ***matrInvM;
double **A, **Ap, **Am, **L, **R, ***mA;

double **cellGP, **cellGW, *cellC;

double *smMu;

FILE *hLog;

double L_, R_, P_, U_, TIME_/*, E_, CV_, TIME_, MU_, KP_*/; // параметры обезразмеривания

#define POW_2(X) ((X)*(X))

#define FLUX     rim_orig
#define FLUX_RHS flux_rim

#define LIMITER_I  calcLimiterEigenv
#define LIMITER_II calcLimiter_II

const bool USE_LIMITER_I = true;
const bool USE_LIMITER_II = true;
const bool USE_SMOOTHER = false;

/*  базисные функции и поля  */
double getField(int i, int iCell, double x);
double getF(int i, int iCell, double x);
double getDF(int i, int iCell, double x);
double getFieldsAvg(int idField, int iCell, double x);
/*  функции  */
void init();
void done();
void LOG(char * format, ...);

void primToCons(double r, double p, double u, double &ro, double &ru, double &re);
void consToPrim(double &r, double &p, double &u, double ro, double ru, double re);

void calcIntegral();
void calcFlux();

extern void roe_orig(double& RI, double& EI, double& PI, double& UI, double& VI, double& WI,
	double RB, double PB, double UB, double VB, double WB,
	double RE, double PE, double UE, double VE, double WE, double GAM);

extern void rim_orig(double& RI, double& EI, double& PI, double& UI, double& VI, double& WI,
	double RB, double PB, double UB, double VB, double WB,
	double RE, double PE, double UE, double VE, double WE, double GAM);

void flux_lf(double &fr, double &fu, double &fe,
	double rb, double pb, double ub, double re, double pe, double ue, double GAM);

void flux_cir(double &fr, double &fu, double &fe,
	double rb, double pb, double ub, double re, double pe, double ue, double GAM);

void flux_rim(double &fr, double &fu, double &fe,
	double rb, double pb, double ub, double re, double pe, double ue, double GAM);

void URS(int iCase, double& r, double& p, double& e, double gam);
double _max_(double a, double b);

void calcMatrAPM(double c2, double u, double GAM, double** Am, double** Ap);
void calcMatrA(double c2, double u, double GAM, double** A);
void calcMatrL(double c2, double u, double GAM, double** R);
void calcMatrR(double c2, double u, double GAM, double** R);
void calcLimiterCons();
void calcLimiterEigenv();
void calcLimiter_II();
void calcSmoother();

void multMatrVect(double** M, double* v, double* result, int n);


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
		for (int i = 0; i < N; i++) {
			memcpy(ro_[i], ro[i], FUNC_COUNT*sizeof(double));
			memcpy(ru_[i], ru[i], FUNC_COUNT*sizeof(double));
			memcpy(re_[i], re[i], FUNC_COUNT*sizeof(double));
		}
		int iter = 0;
		double err, maxU;
		do {
			iter++;
			for (int i = 0; i < N; i++) {
				memset(dr[i], 0, FUNC_COUNT*sizeof(double));
				memset(du[i], 0, FUNC_COUNT*sizeof(double));
				memset(de[i], 0, FUNC_COUNT*sizeof(double));
			}
			calcIntegral();
			calcFlux();
			double tmpRO[FUNC_COUNT], tmpRU[FUNC_COUNT], tmpRE[FUNC_COUNT];
			double dRO[FUNC_COUNT], dRU[FUNC_COUNT], dRE[FUNC_COUNT];
			double delta;
			err = 0.0;
			maxU = 0.0;
			for (int i = 0; i < N; i++) {
				//for (int iF = 0; iF < FUNC_COUNT; iF++) {
				//	dRO[iF] = (ro[i][iF] - ro_[i][iF]) / tau;
				//	dRU[iF] = (ru[i][iF] - ru_[i][iF]) / tau;
				//	dRE[iF] = (re[i][iF] - re_[i][iF]) / tau;
				//}
				//multMatrVect(matrM[i], dRO, tmpRO, FUNC_COUNT);
				//multMatrVect(matrM[i], dRU, tmpRU, FUNC_COUNT);
				//multMatrVect(matrM[i], dRE, tmpRE, FUNC_COUNT);
				//for (int j = 0; j < FUNC_COUNT; j++) {
				//	delta = SIGMA*(tmpRO[j] + dr[i][j]);
				//	if (fabs(delta) > err) err = fabs(delta);
				//	ro[i][j] -= delta;

				//	delta = SIGMA*(tmpRU[j] + du[i][j]);
				//	if (fabs(delta) > err) err = fabs(delta);
				//	ru[i][j] -= delta;

				//	delta = SIGMA*(tmpRE[j] + de[i][j]);
				//	if (fabs(delta) > err) err = fabs(delta);
				//	re[i][j] -= delta;
				//}
				multMatrVect(matrInvM[i], dr[i], tmpRO, FUNC_COUNT);
				multMatrVect(matrInvM[i], du[i], tmpRU, FUNC_COUNT);
				multMatrVect(matrInvM[i], de[i], tmpRE, FUNC_COUNT);
				for (int j = 0; j < FUNC_COUNT; j++) {
					delta = SIGMA*(tmpRO[j] + (ro[i][j] - ro_[i][j]) / tau);
					if (fabs(delta) > err) err = fabs(delta);
					if (fabs(ro[i][j]) > maxU) maxU = fabs(ro[i][j]);
					ro[i][j] -= delta;

					delta = SIGMA*(tmpRU[j] + (ru[i][j] - ru_[i][j]) / tau);
					if (fabs(delta) > err) err = fabs(delta);
					if (fabs(ru[i][j]) > maxU) maxU = fabs(ru[i][j]);
					ru[i][j] -= delta;

					delta = SIGMA*(tmpRE[j] + (re[i][j] - re_[i][j]) / tau);
					if (fabs(delta) > err) err = fabs(delta);
					if (fabs(re[i][j]) > maxU) maxU = fabs(re[i][j]);
					re[i][j] -= delta;
				}
			}
		} while (iter < MAX_ITER && err/maxU > EPS);

		//if (LIMITER_I) calcLimiterEigenv();

		if (step%PRINT_STEP == 0) LOG("step = %6d | time = %16.8E | iterations = %6d | err = %16.8E\n", step, t*TIME_, iter, err);
		if (step%SAVE_STEP == 0) {
			char fName[64];
			sprintf(fName, "res_gp_%08d.csv", step);
			FILE * fp = fopen(fName, "w");
			for (int i = 0; i < N; i++) {
				for (int iG = 0; iG < 2; iG++) {
					double x = cellGP[i][iG];
					double fRO = getField(0, i, x);
					double fRU = getField(1, i, x);
					double fRE = getField(2, i, x);
					double r, p, u;
					consToPrim(r, p, u, fRO, fRU, fRE);
					fprintf(fp, "%16.8E, %16.8E, %16.8E, %16.8E, %16.8E\n", x, r, p, u, cellGW[i][iG]);
				}
				fprintf(fp, "\n");
			}
			fclose(fp);
		}
	}

	done();
}

void calcIntegral()
{
	double ri, ui, ei, pi;
	for (int iCell = 0; iCell < N; iCell++) {
		for (int iG = 0; iG < 2; iG++) {
			double x = cellGP[iCell][iG];
			double fRO = getField(0, iCell, x);
			double fRU = getField(1, iCell, x);
			double fRE = getField(2, iCell, x);
			consToPrim(ri, pi, ui, fRO, fRU, fRE);
			URS(3, ri, pi, ei, GAM);
			fRO = ri*ui;
			fRU = fRO*ui+pi;
			fRE = (ri*(ei+0.5*ui*ui)+pi)*ui;
			fRO *= cellGW[iCell][iG];
			fRU *= cellGW[iCell][iG];
			fRE *= cellGW[iCell][iG];
			for (int iF = 0; iF < FUNC_COUNT; iF++) {
				dr[iCell][iF] -= getDF(iF, iCell, x) * fRO;
				du[iCell][iF] -= getDF(iF, iCell, x) * fRU;
				de[iCell][iF] -= getDF(iF, iCell, x) * fRE;
			}
		}
	}
}

void calcFlux()
{
	double ri, ei, pi, ui, vi, wi, rb, pb, ub, vb, wb, re, pe, ue, ve, we;
	vb = 0.0; ve = 0.0; wb = 0.0; we = 0.0;

	for (int iCell = 0; iCell < N - 1; iCell++) {
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

		double FR = ri*ui;
		double FU = FR*ui + pi;
		double FE = FR*(ei + 0.5*ui*ui) + pi*ui;

		for (int iF = 0; iF < FUNC_COUNT; iF++) {
			dr[iCell][iF] += FR*getF(iF, iCell, x1);
			du[iCell][iF] += FU*getF(iF, iCell, x1);
			de[iCell][iF] += FE*getF(iF, iCell, x1);

			dr[iCell + 1][iF] -= FR*getF(iF, iCell + 1, x1);
			du[iCell + 1][iF] -= FU*getF(iF, iCell + 1, x1);
			de[iCell + 1][iF] -= FE*getF(iF, iCell + 1, x1);
		}

	}
	{	int iCell = 0;
		double x1 = cellC[iCell] - 0.50*h;

		double ROr = getField(0, iCell, x1);
		double RUr = getField(1, iCell, x1);
		double REr = getField(2, iCell, x1);
		consToPrim(re, pe, ue, ROr, RUr, REr);
		rb = re; ub = ue; pb = pe; 
		FLUX(ri, ei, pi, ui, vi, wi,
			rb, pb, ub, vb, wb,
			re, pe, ue, ve, we, GAM);

		double FR = ri*ui;
		double FU = FR*ui + pi;
		double FE = FR*(ei + 0.5*ui*ui) + pi*ui;

		for (int iF = 0; iF < FUNC_COUNT; iF++) {
			dr[iCell][iF] -= FR*getF(iF, iCell, x1);
			du[iCell][iF] -= FU*getF(iF, iCell, x1);
			de[iCell][iF] -= FE*getF(iF, iCell, x1);

		}

	}
	{	int iCell = N-1;
		double x1 = cellC[iCell] + 0.50*h;

		double ROl = getField(0, iCell, x1);
		double RUl = getField(1, iCell, x1);
		double REl = getField(2, iCell, x1);
		consToPrim(rb, pb, ub, ROl, RUl, REl);
		re = rb; ue = ub; pe = pb;
		FLUX(ri, ei, pi, ui, vi, wi,
			rb, pb, ub, vb, wb,
			re, pe, ue, ve, we, GAM);

		double FR = ri*ui;
		double FU = FR*ui + pi;
		double FE = FR*(ei + 0.5*ui*ui) + pi*ui;

		for (int iF = 0; iF < FUNC_COUNT; iF++) {
			dr[iCell][iF] += FR*getF(iF, iCell, x1);
			du[iCell][iF] += FU*getF(iF, iCell, x1);
			de[iCell][iF] += FE*getF(iF, iCell, x1);
		}

	}
}

double _max_(double a, double b) {
	if (a > b) return a;
	return b;
}

void inverseMatr(double** a_src, double **am, int N)
{
	int	*	mask;
	double	fmaxval;
	int		maxind;
	int		tmpi;
	double	tmp;
	//double	a[N][N];

	double	**a;

	mask = new int[N];
	a = new double*[N];
	for (int i = 0; i < N; i++)
	{
		a[i] = new double[N];
		for (int j = 0; j < N; j++)
		{
			a[i][j] = a_src[i][j];
		}
	}
	//::memcpy(a, a_src, sizeof(double)*N*N);

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (i == j)
			{
				am[i][j] = 1.0;
			}
			else {
				am[i][j] = 0.0;
			}
		}
	}
	for (int i = 0; i < N; i++)
	{
		mask[i] = i;
	}
	for (int i = 0; i < N; i++)
	{
		maxind = i;
		fmaxval = fabs(a[i][i]);
		for (int ni = i + 1; ni < N; ni++)
		{
			if (fabs(fmaxval) <= fabs(a[ni][i]))
			{
				fmaxval = fabs(a[ni][i]);
				maxind = ni;
			}
		}
		fmaxval = a[maxind][i];
		if (fmaxval == 0)
		{
			printf("ERROR! Determinant of mass matrix is zero...\n");
			printf("-------------------------------------------\n");
			for (int ii = 0; ii < N; ii++) {
				for (int jj = 0; jj < N; jj++) {
					printf("%16.8E ", a_src[ii][jj]);
				}
				printf("\n");
			}
			printf("===========================================\n\n\n");

			return;
		}
		if (i != maxind)
		{
			for (int nj = 0; nj < N; nj++)
			{
				tmp = a[i][nj];
				a[i][nj] = a[maxind][nj];
				a[maxind][nj] = tmp;

				tmp = am[i][nj];
				am[i][nj] = am[maxind][nj];
				am[maxind][nj] = tmp;
			}
			tmpi = mask[i];
			mask[i] = mask[maxind];
			mask[maxind] = tmpi;
		}
		double aii = a[i][i];
		for (int j = 0; j < N; j++)
		{
			a[i][j] = a[i][j] / aii;
			am[i][j] = am[i][j] / aii;
		}
		for (int ni = 0; ni < N; ni++)
		{
			if (ni != i)
			{
				double fconst = a[ni][i];
				for (int nj = 0; nj < N; nj++)
				{
					a[ni][nj] = a[ni][nj] - fconst *  a[i][nj];
					am[ni][nj] = am[ni][nj] - fconst * am[i][nj];
				}
			}
		}
	}
	//for (int i = 0; i < N; i++)
	//{
	//	if (mask[i] != i) 
	//	{
	//		for (int j = 0; j < N; j++) 
	//		{
	//			tmp				= a[i][j];
	//			a[i][j]			= a[mask[i]][j];
	//			a[mask[i]][j]	= tmp;
	//		}
	//	}
	//}
	for (int i = 0; i < N; i++)
	{
		delete[] a[i];
	}
	delete[] a;
	delete[] mask;
	return;
}


void calcMassMatr() {
	for (int iCell = 0; iCell < N; iCell++){
		//for (int i = 0; i < FUNC_COUNT; i++) {
		//	for (int j = 0; j < FUNC_COUNT; j++){
		//		matrM[iCell][i][j] = 0.0;
		//		for (int iGP = 0; iGP < 2; iGP++) {
		//			double x = cellGP[iCell][iGP];
		//			matrM[iCell][i][j] = matrM[iCell][i][j] + cellGW[iCell][iGP] * getF(i, iCell, x)*getF(j, iCell, x);
		//		}
		//	}
		//}
		matrM[iCell][0][0] = h;
		matrM[iCell][0][1] = matrM[iCell][1][0] = 0.0;
		matrM[iCell][1][1] = h/12.0;

		if (FUNC_COUNT == 3) {
			matrM[iCell][0][2] = matrM[iCell][2][0] = h / 12.0;
			matrM[iCell][1][2] = matrM[iCell][2][1] = 0.0;
			matrM[iCell][2][2] = h/80;
		}

		inverseMatr(matrM[iCell], matrInvM[iCell], FUNC_COUNT);
	}
	
}

void memAlloc() {
	ro = new double*[N];
	ru = new double*[N];
	re = new double*[N];
	dr = new double*[N];
	du = new double*[N];
	de = new double*[N];
	ro_ = new double*[N];
	ru_ = new double*[N];
	re_ = new double*[N];
	matrM = new double**[N];
	matrInvM = new double**[N];
	cellGP = new double*[N];
	cellGW = new double*[N];
	cellC = new double[N];
	for (int i = 0; i < N; i++) {
		ro[i] = new double[FUNC_COUNT];
		ru[i] = new double[FUNC_COUNT];
		re[i] = new double[FUNC_COUNT];
		dr[i] = new double[FUNC_COUNT];
		du[i] = new double[FUNC_COUNT];
		de[i] = new double[FUNC_COUNT];
		ro_[i] = new double[FUNC_COUNT];
		ru_[i] = new double[FUNC_COUNT];
		re_[i] = new double[FUNC_COUNT];
		cellGP[i] = new double[2];
		cellGW[i] = new double[2];
		matrM[i] = new double*[FUNC_COUNT];
		matrInvM[i] = new double*[FUNC_COUNT];
		for (int j = 0; j < FUNC_COUNT; j++) {
			matrM[i][j] = new double[FUNC_COUNT];
			matrInvM[i][j] = new double[FUNC_COUNT];
		}
	}
	r = new double[N];
	u = new double[N];
	p = new double[N];
	A = new double*[3];
	L = new double*[3];
	R = new double*[3];
	//Am = new double*[3];
	//Ap = new double*[3];
	for (int i = 0; i < 3; i++) {
		A[i] = new double[3];
		//Am[i] = new double[3];
		//Ap[i] = new double[3];
		L[i] = new double[3];
		R[i] = new double[3];
	}
	//matr9 = new double*[MATR_BLOCK];
	//for (int i = 0; i < MATR_BLOCK; i++) {
	//	matr9[i] = new double[MATR_BLOCK];
	//}

	//mA = new double**[2];
	//for (int i = 0; i < 2; i++) {
	//	mA[i] = new double*[3];
	//	for (int j = 0; j < 3; j++) {
	//		mA[i][j] = new double[3];
	//	}
	//}
	smMu = new double[N + 5];
}

void memFree() {
	for (int i = 0; i < N; i++) {
		delete[] ro[i];
		delete[] ru[i];
		delete[] re[i];
		delete[] dr[i];
		delete[] du[i];
		delete[] de[i];
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
	delete[] dr;
	delete[] du;
	delete[] de;
	delete[] ro_;
	delete[] ru_;
	delete[] re_;
	delete[] matrM;
	delete[] cellC;

	delete[] r;
	delete[] u;
	delete[] p;
	for (int i = 0; i < 3; i++) {
		delete[] A[i];
		//delete[] Am[i];
		//delete[] Ap[i];
		delete[] L[i];
		delete[] R[i];
	}
	delete[] A;
	//delete[] Ap;
	//delete[] Am;
	delete[] L;
	delete[] R;

	delete[] smMu;
}

void init() {
	char *logName = new char[128];
	time_t t = time(0);   // get time now
	struct tm * now = localtime(&t);
	sprintf(logName, "%4d-%02d-%02d_%02d.%02d.%02d", now->tm_year + 1900, now->tm_mon + 1, now->tm_mday, now->tm_hour, now->tm_min, now->tm_sec);
	logName = strcat(logName, ".log");

	hLog = fopen(logName, "w");

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

	for (int i = 0; i < N; i++) {
		double x = cellC[i];
		// Sod  !!! XMIN = -1.0; XMAX = 1.0;
		//if (x < 0.0) {
		//	r[i] = 1.0;
		//	p[i] = 1.0;
		//	u[i] = 0.0;
		//} else {
		//	r[i] = 0.125;
		//	p[i] = 0.1;
		//	u[i] = 0.0;
		//}
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

		// Simple wave
		const double L = 0.2;
		if (fabs(x) < L) {
			r[i] = 1.0 + exp(2.0 - 2.0*((L*L) / ((L*L) - (x*x))));
		}
		else {
			r[i] = 1.0;
		}
		double e = ::pow(r[i], AGAM);
		u[i] = -2.0*sqrt(e*AGAM*GAM) / AGAM;
		p[i] = r[i] * e*AGAM;

		//double dt = CFL*h / (fabs(u[i]) + sqrt(GAM*p[i] / r[i]));
		//if (dt < tau) tau = dt;
	}
	printf("TAU = %25.16e\n\n", tau);

	double maxR = r[0];
	double maxP = p[0];
	for (int i = 1; i < N; i++) {
		if (maxR < r[i]) maxR = r[i];
		if (maxP < r[i]) maxP = p[i];
	}

	// параметры обезразмеривания
	L_ = 1.0;
	R_ = 1.0; //maxR;					// характерная плотность = начальная плотность
	P_ = 1.0; //maxP;					// характерное давление = начальное давление 
	//T_ = maxT;					// характерная температура = начальная температура
	U_ = 1.0; //sqrt(P_ / R_);		// характерная скорость = sqrt( P_ / R_ )
	//E_ = POW_2(U_);			// характерная энергия  = U_**2
	//CV_ = POW_2(U_) / T_;	// характерная теплоёмкость  = U_**2 / T_
	TIME_ = 1.0; //L_ / U_;			// характерное время
	//MU_ = R_ * U_ * L_;		// характерная вязкость = R_ * U_ * L_
	///KP_ = R_ * POW_2(U_) * U_ * L_ / T_;	// коэффициент теплопроводности = R_ * U_**3 * L_ / T_
	//CV_ = POW_2(U_) / T_;	// характерная теплоёмкость  = U_**2 / T_

	for (int i = 0; i < N; i++) {
		r[i] /= R_;
		p[i] /= P_;
		u[i] /= U_;
		primToCons(r[i], p[i], u[i], ro[i][0], ru[i][0], re[i][0]);
		for (int j = 1; j < FUNC_COUNT; j++) {
			ro[i][j] = 0.0;
			ru[i][j] = 0.0;
			re[i][j] = 0.0;
		}
	}

	TMAX /= TIME_;
	tau	 /= TIME_;

	calcMassMatr();

}

void done() {
	fclose(hLog);
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

inline double avg(double *fld)
{
	double avg = fld[0];
	if (FUNC_COUNT == 3) avg += fld[2] / 12.0;
	return avg;
}

void calcLimiterEigenv() {
	double c2_, u_;
	for (int i = 0; i < N; i++) {
		u_ = avg(ru[i]) / avg(ro[i]);
		c2_ = (avg(re[i]) / avg(ro[i]) - u_*u_*0.5)*GAM*AGAM;
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
		u0 = avg(ro_[iCell]);
		u1 = ro_[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(avg(ro_[iCell + 1]) - u0),
			LIM_ALPHA*(u0 - avg(ro_[iCell])));
		if (u1l != u1) {
			ro_[iCell][0] = u0;
			ro_[iCell][1] = u1l;
			if (FUNC_COUNT > 2) ro_[iCell][2] = 0.0;
		}

		u0 = avg(ru_[iCell]);
		u1 = ru_[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(avg(ru_[iCell + 1]) - u0),
			LIM_ALPHA*(u0 - avg(ru_[iCell])));
		if (u1l != u1) {
			ru_[iCell][0] = u0;
			ru_[iCell][1] = u1l;
			if (FUNC_COUNT > 2) ru_[iCell][2] = 0.0;
		}

		u0 = avg(re_[iCell]);
		u1 = re_[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(avg(re_[iCell + 1]) - u0),
			LIM_ALPHA*(u0 - avg(re_[iCell])));
		if (u1l != u1) {
			re_[iCell][0] = u0;
			re_[iCell][1] = u1l;
			if (FUNC_COUNT > 2) re_[iCell][2] = 0.0;
		}

	}
	for (int iCell = 1; iCell < N - 1; iCell++) {
		// проецируем на линейный базис
		double u0, u1, u1l;
		u0 = avg(ro_[iCell]);
		u1 = ro_[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(avg(ro_[iCell + 1]) - u0),
			LIM_ALPHA*(u0 - avg(ro_[iCell - 1])));
		if (u1l != u1) {
			ro_[iCell][0] = u0;
			ro_[iCell][1] = u1l;
			if (FUNC_COUNT > 2) ro_[iCell][2] = 0.0;
		}

		u0 = avg(ru_[iCell]);
		u1 = ru_[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(avg(ru_[iCell + 1]) - u0),
			LIM_ALPHA*(u0 - avg(ru_[iCell - 1])));
		if (u1l != u1) {
			ru_[iCell][0] = u0;
			ru_[iCell][1] = u1l;
			if (FUNC_COUNT > 2) ru_[iCell][2] = 0.0;
		}

		u0 = avg(re_[iCell]);
		u1 = re_[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(avg(re_[iCell + 1]) - u0),
			LIM_ALPHA*(u0 - avg(re_[iCell - 1])));
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
		u0 = avg(ro_[iCell]);
		u1 = ro_[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(avg(ro_[iCell]) - u0),
			LIM_ALPHA*(u0 - avg(ro_[iCell - 1])));
		if (u1l != u1) {
			ro_[iCell][0] = u0;
			ro_[iCell][1] = u1l;
			if (FUNC_COUNT > 2) ro_[iCell][2] = 0.0;
		}

		u0 = avg(ru_[iCell]);
		u1 = ru_[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(avg(ru_[iCell]) - u0),
			LIM_ALPHA*(u0 - avg(ru_[iCell - 1])));
		if (u1l != u1) {
			ru_[iCell][0] = u0;
			ru_[iCell][1] = u1l;
			if (FUNC_COUNT > 2) ru_[iCell][2] = 0.0;
		}

		u0 = avg(re_[iCell]);
		u1 = re_[iCell][1];
		u1l = minmod(u1,
			LIM_ALPHA*(avg(re_[iCell]) - u0),
			LIM_ALPHA*(u0 - avg(re_[iCell - 1])));
		if (u1l != u1) {
			re_[iCell][0] = u0;
			re_[iCell][1] = u1l;
			if (FUNC_COUNT > 2) re_[iCell][2] = 0.0;
		}

	}

	for (int i = 0; i < N; i++) {
		u_ = avg(ru[i]) / avg(ro[i]);
		c2_ = (avg(re[i]) / avg(ro[i]) - u_*u_*0.5)*GAM*AGAM;
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
		memcpy(ro_[i], ro[i], FUNC_COUNT * sizeof(double));
		memcpy(ru_[i], ru[i], FUNC_COUNT * sizeof(double));
		memcpy(re_[i], re[i], FUNC_COUNT * sizeof(double));
	}
	for (int idField = 0; idField < 3; idField++) {
		double** fld, **fldOld;
		switch (idField) {
		case 0:
			fld = ro_;
			fldOld = ro;
			break;
		case 1:
			fld = ru_;
			fldOld = ru;
			break;
		case 2:
			fld = re_;
			fldOld = re;
			break;
		}

		// i+1/2
		for (int iCell = 0; iCell < N - 1; iCell++) {
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
			double ROr = getField(0, iCell + 1, cellC[iCell]);
			double RUr = getField(1, iCell + 1, cellC[iCell]);
			double REr = getField(2, iCell + 1, cellC[iCell]);
			double RI, EI, PI, UI, VI, WI;
			roe_orig(RI, EI, PI, UI, VI, WI,
				ROl, (REl - 0.5*RUl*RUl / ROl)*(GAM - 1.0), RUl / ROl, 0, 0,
				ROr, (REl - 0.5*RUl*RUl / ROl)*(GAM - 1.0), RUr / ROr, 0, 0, GAM);
			double CFL = tau*(fabs(UI) + sqrt(GAM*PI / RI)) / h;
			smMu[iCell] = CFL*2.0 * fIntP1P2 / (fIntP1Pavg + fIntP2Pavg);
		}

		for (int i = 1; i < N - 1; i++) {
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
		memcpy(ro[i], ro_[i], FUNC_COUNT * sizeof(double));
		memcpy(ru[i], ru_[i], FUNC_COUNT * sizeof(double));
		memcpy(re[i], re_[i], FUNC_COUNT * sizeof(double));
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


double getField(int i, int iCell, double x){

	double f1 = getF(1, iCell, x);
	double f2 = getF(2, iCell, x);
	switch (i) {
	case 0:
		return ro[iCell][0] + ro[iCell][1] * f1 + (FUNC_COUNT < 3 ? 0.0 : ro[iCell][2]) * f2;
		break;
	case 1:
		return ru[iCell][0] + ru[iCell][1] * f1 + (FUNC_COUNT < 3 ? 0.0 : ru[iCell][2]) * f2;
		break;
	case 2:
		return re[iCell][0] + re[iCell][1] * f1 + (FUNC_COUNT < 3 ? 0.0 : re[iCell][2]) * f2;
		break;
	}
}

/**
*	Average on iCell+1/2
*/
double getFieldsAvg(int idField, int iCell, double x) {
	if (iCell > -1 && iCell < N - 1) {
		return 0.5*(getField(idField, iCell, x) + getField(idField, iCell + 1, x));
	}
	else {
		if (iCell = -1) {
			return getField(idField, iCell + 1, x);
		}
		else {
			return getField(idField, iCell, x);
		}
	}
}


void calcMatrL(double c2, double u, double GAM, double** L)
{
	double cz, kk, deta;

	cz = sqrt(c2);
	kk = GAM - 1.0;
	deta = 0.5*kk*u*u;

	L[0][0] = -c2 + deta;
	L[0][1] = -kk*u;
	L[0][2] = kk;

	L[1][0] = -cz*u + deta;
	L[1][1] = cz - kk*u;
	L[1][2] = kk;

	L[2][0] = cz*u + deta;
	L[2][1] = -cz - kk*u;
	L[2][2] = kk;

}

void calcMatrR(double c2, double u, double GAM, double** R)
{
	double cz, kk, deta;

	cz = sqrt(c2);
	kk = GAM - 1.0;
	R[0][0] = -1.0 / c2;
	R[0][1] = 0.5 / c2;
	R[0][2] = 0.5 / c2;

	R[1][0] = -u / c2;
	R[1][1] = 0.5*u / c2 + 0.5 / cz;
	R[1][2] = 0.5*u / c2 - 0.5 / cz;

	R[2][0] = -0.5*u*u / c2;
	R[2][1] = 0.25*u*u / c2 + 0.5 / kk + 0.5*u / cz;
	R[2][2] = 0.25*u*u / c2 + 0.5 / kk - 0.5*u / cz;

}

void calcMatrA(double c2, double u, double GAM, double** A) {
	double cz = sqrt(c2);
	double kk = GAM - 1.0;

	A[0][0] = 0.0;
	A[0][1] = 1.0;
	A[0][2] = 0.0;

	A[1][0] = (GAM - 3.0)*u*u / 2.0;
	A[1][1] = (3.0 - GAM)*u;
	A[1][2] = kk;

	A[2][0] = kk*(u*u*u) - u*(c2 / kk + GAM*u*u / 2.0);
	A[2][1] = c2 / kk + GAM*u*u / 2.0 - 1.5*kk*u*u;
	A[2][2] = GAM*u;
}

void calcMatrAPM(double c2, double u, double GAM, double** Am, double** Ap) {
	double matr[3][3], matr1[3][3];

	double cz = sqrt(c2);
	double ll = fabs(u) + cz;

	calcMatrL(c2, u, GAM, L);
	calcMatrR(c2, u, GAM, R);

	//double E1[3][3], E2[3][3];

	//for (int i = 0; i < 3; i++) {
	//	for (int j = 0; j < 3; j++) {
	//		E1[i][j] = 0.0;
	//		E2[i][j] = 0.0;
	//		for (int k = 0; k < 3; k++) {
	//			E1[i][j] += L[i][k] * R[k][j];
	//			E2[i][j] += R[i][k] * L[k][j];
	//		}
	//	}
	//}


	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			if (i != j) matr1[i][j] = 0.0;
		}
	}

	// Am
	matr1[0][0] = 0.5*(u - fabs(u));
	matr1[1][1] = 0.5*(u + cz - fabs(u + cz));
	matr1[2][2] = 0.5*(u - cz - fabs(u - cz));
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++)	{
			matr[i][j] = 0.0;
			for (int k = 0; k < 3; k++) {
				matr[i][j] += R[i][k] * matr1[k][j];
			}
		}
	}
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++)	{
			Am[i][j] = 0.0;
			for (int k = 0; k < 3; k++) {
				Am[i][j] += matr[i][k] * L[k][j];
			}
		}
	}

	// Ap
	matr1[0][0] = 0.5*(u + fabs(u));
	matr1[1][1] = 0.5*(u + cz + fabs(u + cz));
	matr1[2][2] = 0.5*(u - cz + fabs(u - cz));
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++)	{
			matr[i][j] = 0.0;
			for (int k = 0; k < 3; k++) {
				matr[i][j] += R[i][k] * matr1[k][j];
			}
		}
	}
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++)	{
			Ap[i][j] = 0.0;
			for (int k = 0; k < 3; k++) {
				Ap[i][j] += matr[i][k] * L[k][j];
			}
		}
	}

}



void roe_orig(double& RI, double& EI, double& PI, double& UI, double& VI, double& WI,
	double RB, double PB, double UB, double VB, double WB,
	double RE, double PE, double UE, double VE, double WE, double GAM) {

	double AGAM = (GAM - 1.0);

	// Схема ROE

	//double fG = GAM;

	double fSB = sqrt(RB);
	double fSE = sqrt(RE);
	double fS_ = 1.0 / (fSB + fSE);

	RI = fSB * fSE;

	UI = (fSB * UB + fSE * UE) * fS_;
	VI = (fSB * VB + fSE * VE) * fS_;
	WI = (fSB * WB + fSE * WE) * fS_;

	double EB = PB / (RB*AGAM);
	double EE = PE / (RE*AGAM);
	//EI = (fSB * EB + fSE * EE) * fS_;
	//TI = ( fSB * TB + fSE * TE ) * fS_;


	double HB = EB + (UB*UB + VB*VB + WB*WB)*0.5 + PB / RB;
	double HE = EE + (UE*UE + VE*VE + WE*WE)*0.5 + PE / RE;

	double HI = (fSB * HB + fSE * HE) * fS_;

	PI = (HI - (UI*UI + VI*VI + WI*WI)*0.5) * RI * AGAM / GAM;

	// TODO:
	EI = PI / (RI*AGAM);


}

void rim_orig(double& RI, double& EI, double& PI, double& UI, double& VI, double& WI,
	double RB, double PB, double UB, double VB, double WB,
	double RE, double PE, double UE, double VE, double WE, double GAM) {

	double AGAM = (GAM - 1.0);
	double BGAM = (2.0*sqrt(GAM / AGAM));
	double CGAM = (1.0 / GAM);
	double DGAM = (2.0 / AGAM);
	double EGAM = (AGAM / (GAM + 1.0));
	double GGAM = (sqrt(GAM*AGAM));
	double HGAM = (AGAM / 2.0);
	double FGAM = (3.0*GAM - 1.0);
	double OGAM = (AGAM / (2.0*GAM));
	double QGAM = (GAM + 1.0);
	double PGAM = (QGAM / (2.0*GAM));
	double RGAM = (4.0*GAM);
	double SGAM = (GAM*AGAM);
	double TGAM = (QGAM / 2.0);
	double UGAM = (sqrt(AGAM / GAM));

	double RF, RS, EF, ES, SBL, SFL, SSL, SEL, D;
	double    eps = 1.0e-5;
	double    CB = sqrt(GAM*PB / RB);
	double    CE = sqrt(GAM*PE / RE);
	double    EB = CB*CB / SGAM;
	double    EE = CE*CE / SGAM;
	double    RCB = RB*CB;
	double    RCE = RE*CE;
	double    DU = UB - UE;
	if (DU < -2.0*(CB + CE) / AGAM) {
		//printf(" ATTENTION!!!  VACUUM \n");
		RF = 0.0;
		RS = 0.0;
		EF = 0.0;
		ES = 0.0;
		SBL = UB - CB;
		SFL = UB + 2.0*CB / AGAM;
		SSL = UE - 2.0*CE / AGAM;
		SEL = UE + CE;
		goto lbl9;
	}
	double P = (PB*RCE + PE*RCB + DU*RCB*RCE) / (RCB + RCE);
lbl5:
	if (P<eps) P = eps;

	double PPB = P / PB;
	if (PB>P) goto lbl1;
	double PKB = PGAM*PPB + OGAM;
	double ZNB = RCB*sqrt(PKB);
	double F1 = (P - PB) / ZNB;
	double FS1 = (QGAM*PPB + FGAM) / (RGAM*ZNB*PKB);
	goto lbl2;
lbl1:
	double ZFB = CB*exp(log(PPB)*OGAM);
	F1 = DGAM*(ZFB - CB);
	FS1 = ZFB / (GAM*P);
lbl2:
	double PPE = P / PE;
	if (PE>P) goto lbl3;
	double PKE = PGAM*PPE + OGAM;
	double ZNE = RCE*sqrt(PKE);
	double F2 = (P - PE) / ZNE;
	double FS2 = (QGAM*PPE + FGAM) / (RGAM*ZNE*PKE);
	goto lbl4;
lbl3:
	double ZFE = CE*exp(log(PPE)*OGAM);
	F2 = DGAM*(ZFE - CE);
	FS2 = ZFE / (GAM*P);
lbl4:
	double DP = (DU - F1 - F2) / (FS1 + FS2);
	P = P + DP;
	if (fabs(DU - F1 - F2)>eps) goto lbl5;


	PPB = P / PB;
	PPE = P / PE;

	//       ZFB=CB*PPB**OGAM;
	//       ZFE=CE*PPE**OGAM;
	ZFB = CB*exp(log(PPB)*OGAM);
	ZFE = CE*exp(log(PPE)*OGAM);
	if (PB>P) goto lbl6;
	D = UB - sqrt((TGAM*P + HGAM*PB) / RB);
	double UBD = UB - D;
	double RUBD = RB*UBD;
	RF = RUBD*RUBD / (PB - P + RUBD*UBD);
	double UF = D + RUBD / RF;
	EF = P / (AGAM*RF);
	SBL = D;
	SFL = D;
	goto lbl7;
lbl6:
	EF = ZFB*ZFB / SGAM;
	UF = UB + DGAM*(CB - ZFB);
	RF = P / (AGAM*EF);
	SBL = UB - CB;
	SFL = UF - ZFB;
lbl7:
	if (PE>P) goto lbl8;
	D = UE + sqrt((TGAM*P + HGAM*PE) / RE);
	double UED = UE - D;
	double RUED = RE*UED;
	RS = RUED*RUED / (PE - P + RUED*UED);
	double US = D + RUED / RS;
	ES = P / (AGAM*RS);
	SEL = D;
	SSL = D;
	goto lbl9;
lbl8:
	ES = ZFE*ZFE / SGAM;
	US = UE - DGAM*(CE - ZFE);
	RS = P / (AGAM*ES);
	SSL = US + ZFE;
	SEL = UE + CE;
lbl9:
	// 
	// C     compute the interpolation value
	if (SEL <= 0.0) {
		RI = RE;
		EI = EE;
		UI = UE;
		VI = VE;
		WI = WE;
		goto lbl157;
	}

	if (SBL >= 0.0) {
		RI = RB;
		EI = EB;
		UI = UB;
		VI = VB;
		WI = WB;
		goto lbl157;
	}

	if ((SSL >= 0.0) && (SFL <= 0.0)) {
		if (US >= 0.0) {
			RI = RF;
			EI = EF;
			UI = UF;
			VI = VB;
			WI = WB;
		}
		else {
			RI = RS;
			EI = ES;
			UI = US;
			VI = VE;
			WI = WE;
		}
		goto lbl157;
	}

	if (SFL>0.0) {
		UI = (UB + DGAM*GGAM*sqrt(EB)) / (1 + DGAM);
		VI = VB;
		WI = WB;
		EI = (UI*UI) / SGAM;
		RI = RB*exp(log(EI / EB)*(1 / AGAM));
	}
	else {
		UI = (UE - DGAM*GGAM*sqrt(EE)) / (1 + DGAM);
		VI = VE;
		WI = WE;
		EI = (UI*UI) / SGAM;
		RI = RE*exp(log(EI / EE)*(1 / AGAM));
	}

lbl157:
	PI = AGAM*EI*RI;

}

void flux_cir(double &fr, double &fu, double &fe,
	double rb, double pb, double ub, double re, double pe, double ue, double GAM)
{
	double ri, ei, pi, ui, vi, wi, vb = 0, ve = 0, wb = 0, we = 0;
	double AGAM = GAM - 1.0;
	double matr[3][3], matr1[3][3], Ap[3][3], Am[3][3];
	
	roe_orig(ri, ei, pi, ui, vi, wi,
		rb, pb, ub, vb, wb,
		re, pe, ue, ve, we, GAM);

	double c2 = GAM*pi / ri;
	double cz = sqrt(c2);

	calcMatrL(c2, ui, GAM, L);
	calcMatrR(c2, ui, GAM, R);

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			if (i != j) matr1[i][j] = 0.0;
		}
	}
	
	// Am
	matr1[0][0] = fabs(ui);
	matr1[1][1] = fabs(ui + cz);
	matr1[2][2] = fabs(ui - cz);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++)	{
			matr[i][j] = 0.0;
			for (int k = 0; k < 3; k++) {
				matr[i][j] += R[i][k] * matr1[k][j];
			}
		}
	}
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++)	{
			Am[i][j] = 0.0;
			for (int k = 0; k < 3; k++) {
				Am[i][j] += matr[i][k] * L[k][j];
			}
		}
	}

	double U[3];

	ri = re - rb;
	ui = re*ue - rb*ub;
	ei = (pe / AGAM + re*ue*ue*0.5) - (pb / AGAM + rb*ub*ub*0.5);

	U[0] = ri*Am[0][0] + ui*Am[0][1] + ei*Am[0][2];
	U[1] = ri*Am[1][0] + ui*Am[1][1] + ei*Am[1][2];
	U[2] = ri*Am[2][0] + ui*Am[2][1] + ei*Am[2][2];

	fr = 0.5*(rb*ub + re*ue - U[0]);
	fu = 0.5*(rb*ub*ub + pb + re*ue*ue + pe - U[1]);
	fe = 0.5*((pb / AGAM + rb*ub*ub*0.5 + pb)*ub + (pe / AGAM + re*ue*ue*0.5 + pe)*ue - U[2]);
}

void flux_lf(double &fr, double &fu, double &fe,
	double rb, double pb, double ub, double re, double pe, double ue, double GAM)
{
	double AGAM = GAM - 1.0;
	double alpha = _max_(fabs(ub) + sqrt(GAM*pb / rb), fabs(ue) + sqrt(GAM*pe / re));
	fr = 0.5*(rb*ub + re*ue - alpha*(re - rb));
	fu = 0.5*(rb*ub*ub + pb + re*ue*ue + pe - alpha*(re*ue - rb*ub));
	fe = 0.5*((pb / AGAM + rb*ub*ub*0.5 + pb)*ub + (pe / AGAM + re*ue*ue*0.5 + pe)*ue - alpha*((pe / AGAM + re*ue*ue*0.5) - (pb / AGAM + rb*ub*ub*0.5)));
}

void flux_rim(double &fr, double &fu, double &fe,
	double RB, double PB, double UB,
	double RE, double PE, double UE, double GAM)
{
	double ri, ei, pi, ui, vi, wi, VB = 0, VE = 0, WB = 0, WE = 0;
	double AGAM = GAM - 1.0;
	rim_orig(ri, ei, pi, ui, vi, wi,
		RB, PB, UB, VB, WB,
		RE, PE, UE, VE, WE, GAM);
	fr = ri*ui;
	fu = ri*ui*ui + pi;
	fe = (pi / AGAM + ri*ui*ui*0.5 + pi)*ui;

}

void multMatrVect(double** M, double* v, double* result, int n)
{
	memset(result, 0, n*sizeof(double));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			result[i] += M[i][j] * v[j];
		}
	}
}

void URS(int iCase, double& r, double& p, double& e, double gam) {
	switch (iCase) {
	case 1:
		p = r*e*(gam - 1);
		break;
	case 2:
		r = p / (e*(gam - 1));
		break;
	case 3:
		e = p / (r*(gam - 1));
		break;
	default:
		printf("ERROR: bad URS() parameter\n");
		break;
	}
}// URS


void LOG(char * format, ...)
{
	va_list arglist;

	va_start(arglist, format);

	vprintf(format, arglist);
	vfprintf(hLog, format, arglist);

	va_end(arglist);

	fflush(stdout);
	fflush(hLog);
}
