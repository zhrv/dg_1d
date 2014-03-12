#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "Solver.h"
#include "functions.h"

const double	CFL		= 1000.0;

const int		N		= 500;
const double	XMIN	= -1.0; 
const double	XMAX	=  1.0;
const double	EPS		= 1.0e-5;
const double	GAM		= 1.4;
const double	AGAM	= GAM-1.0;
const double	TMAX	= 0.2;

const int		MAX_ITER	= 1000;
const int		SAVE_STEP	= 100;

double **ro, **ru, **re;
double *r, *u, *e, *p; 
double *dr, *du, *de;

double **matr, *vect, **matr1;
double **A, **Ap, **Am, **L, **R;

double h	=	(XMAX-XMIN)/N;
double tau	=	1.0e-4;


void memAlloc() {
	ro = new double[N];
	ru = new double[N];
	re = new double[N];
	r = new double[N];
	u = new double[N];
	p = new double[N];
	vect = new double[3];
	matr = new double*[3];
	matr1 = new double*[3];
	A = new double*[3];
	L = new double*[3];
	R = new double*[3];
	Am = new double*[3];
	Ap = new double*[3];
	for (int i = 0; i < 3; i++) {
		matr[i]  = new double[3];
		matr1[i]  = new double[3];
		A[i]  = new double[3];
		Am[i] = new double[3];
		Ap[i] = new double[3];
		L[i] = new double[3];
		R[i] = new double[3];
	}
}

void memFree() {
	delete[] ro;
	delete[] ru;
	delete[] re;
	delete[] r;
	delete[] u;
	delete[] p;
	for (int i = 0; i < 3; i++) {
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
}



void primToCons(double r, double p, double u, double &ro, double &ru, double &re) {
	ro = r;
	ru = r*u;
	re = (p/AGAM+r*u*u*0.5);
}

void consToPrim(double &r, double &p, double &u, double ro, double ru, double re) {
	r = ro;
	u = ru/r;
	p = AGAM*(re-r*u*u*0.5);
}

void calcMatrL(double c2, double u, double** L)
{	
	double cz, kk, deta;

    cz = sqrt(c2);
    kk = GAM-1.0;
	deta = 0.5*kk*u*u;

	L[0][0] = -c2+deta;
	L[0][1] = -kk*u;
	L[0][2] = kk;

	L[1][0] = -cz*u+deta;
	L[1][1] = cz-kk*u;
	L[1][2] = kk;

	L[2][0] = cz*u+deta;
	L[2][1] = -cz-kk*u;
	L[2][2] = kk;
	
}

void calcMatrR(double c2, double u, double** R)
{	
	double cz, kk, deta;

    cz = sqrt(c2);
    kk = GAM-1.0;
	R[0][0] = -1.0/c2;
	R[0][1] = 0.5/c2;
	R[0][2] = 0.5/c2;
	
	R[1][0] = -u/c2;
	R[1][1] = 0.5*u/c2 + 0.5/cz;
	R[1][2] = 0.5*u/c2 - 0.5/cz;
	
	R[2][0] = -0.5*u*u/c2;
	R[2][1] = 0.25*u*u/c2 + 0.5/kk + 0.5*u/cz;
	R[2][2] = 0.25*u*u/c2 + 0.5/kk - 0.5*u/cz;
	
}

void calcMatrAPM(double c2, double u, double** Am, double** Ap) {
    double cz = sqrt(c2);
	double ll = fabs(u)+cz;
	
	calcMatrL(c2, u, L);
	calcMatrR(c2, u, R);
	
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			if (i != j) matr1[i][j] = 0.0;
	
	// Am
	matr1[0][0] = 0.5*(u-fabs(u));
	matr1[1][1] = 0.5*(u+cz-fabs(u+cz));
	matr1[2][2] = 0.5*(u-cz-fabs(u-cz));
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++)	{
			matr[i][j] = 0.0;
			for (int k = 0; k < 3; k++) {
				matr[i][j] += R[i][k]*matr1[k][j];
			}
		}
	}
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++)	{
			Am[i][j] = 0.0;
			for (int k = 0; k < 3; k++) {
				Am[i][j] += matr[i][k]*L[k][j];
			}
		}
	}

	// Ap
	matr1[0][0] = 0.5*(u+fabs(u));
	matr1[1][1] = 0.5*(u+cz+fabs(u+cz));
	matr1[2][2] = 0.5*(u-cz+fabs(u-cz));
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++)	{
			matr[i][j] = 0.0;
			for (int k = 0; k < 3; k++) {
				matr[i][j] += R[i][k]*matr1[k][j];
			}
		}
	}
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++)	{
			Ap[i][j] = 0.0;
			for (int k = 0; k < 3; k++) {
				Ap[i][j] += matr[i][k]*L[k][j];
			}
		}
	}

}

double _max_(double a, double b) 
{
	if (a>b) return a;
	return b;
}

double minmod(double a, double b) 
{

	if (a*b < 0.0) return 0;
	if (fabs(a) < fabs(b)) return a;
	return b;
}

int main(int argc, char** argv) 
{
	memAlloc();
	Solver* S = SolverFactory::create(SOLVER_TYPE_ZEIDEL);
	S->init(N);
	for (int i = 0; i < N; i++) {
		double x = XMIN+i*h;
		//// Sod
		//if (x < 0.0) {
		//	r[i] = 1.0;
		//	p[i] = 1.0;
		//	u[i] = 0.0;
		//} else {
		//	r[i] = 0.125;
		//	p[i] = 0.1;
		//	u[i] = 0.0;
		//}
		// Lax
		if (x < 0.0) {
			r[i] = 0.445;
			p[i] = 3.528;
			u[i] = 0.698;
		} else {
			r[i] = 0.5;
			p[i] = 0.571;
			u[i] = 0.0;
		}
		//primToCons(r[i], p[i], u[i], ro[i], ru[i], re[i]);
		double dt = CFL*h/(fabs(u[i])+sqrt(GAM*p[i]/r[i]));
		if (dt < tau) tau = dt;
	}
	printf("TAU = %25.16e\n\n", tau);

	double t = 0.0;
	int step = 0;
	while (t < TMAX) {
		t += tau;
		step++;
		S->zero();
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
		
		int maxIter = MAX_ITER;
		//S->printMatr();
		S->solve(EPS, maxIter);
		
		for (int iCell = 0, ind = 0; iCell < N; iCell++, ind += 9) {
			ro[iCell][0] += S->x[ind + 0];
			ro[iCell][1] += S->x[ind + 1];
			ro[iCell][2] += S->x[ind + 2];
			ru[iCell][0] += S->x[ind + 3];
			ru[iCell][1] += S->x[ind + 4];
			ru[iCell][2] += S->x[ind + 5];
			re[iCell][0] += S->x[ind + 6];
			re[iCell][1] += S->x[ind + 7];
			re[iCell][2] += S->x[ind + 8];
			//consToPrim(r[iCell], p[iCell], u[iCell], ro[iCell], ru[iCell], re[iCell]);
		}
		
		printf("%10d | ITER: %4d | RO: ... | RU: ... | RE ... \n", step, maxIter);
		if (step % SAVE_STEP == 0) {
			char str[50];
			sprintf(str, "res_%010d.csv", step);
			FILE * fp = fopen(str, "w");
			for (int i = 0; i < N; i++) {
				fprintf(fp, "%25.15e, %25.15e, %25.15e, %25.15e\n", XMIN+h*i, r[i], p[i], u[i]);
			}
			fclose(fp);
		}
	}

	delete S;
	memFree();
	return 0;
}






























