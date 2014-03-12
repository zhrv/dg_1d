#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "Solver.h"

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

double *ro, *ru, *re;
double *r, *u, *e, *p; 
double *dr, *du, *de;

double **matr, *vect, **matr1;
double **A, **Ap, **Am, **L, **R;

double h	=	(XMAX-XMIN)/N;
double tau	=	1.0e-4;

void roe_orig(double& RI, double& EI, double& PI, double& UI, double& VI, double& WI,
         double RB, double PB, double UB, double VB, double WB,
         double RE, double PE, double UE, double VE, double WE, double GAM) {

	double AGAM  =  (GAM-1.0);

	// Ñץולא ROE
	
	double fG = GAM;

	double fSB = sqrt( RB );
	double fSE = sqrt( RE );
	double fS_ = 1.0 / ( fSB + fSE );

	RI = fSB * fSE;
		
	UI = ( fSB * UB + fSE * UE ) * fS_;
	VI = ( fSB * VB + fSE * VE ) * fS_;
	WI = ( fSB * WB + fSE * WE ) * fS_;

	double EB = PB/(RB*AGAM);
	double EE = PE/(RE*AGAM);
	EI = ( fSB * EB + fSE * EE ) * fS_;
	//TI = ( fSB * TB + fSE * TE ) * fS_;
		
		
	double HB = EB + (UB*UB+VB*VB+WB*WB)*0.5 + PB / RB;
	double HE = EE + (UE*UE+VE*VE+WE*WE)*0.5 + PE / RE;
		
	double HI = ( fSB * HB + fSE * HE ) * fS_;

	PI = ( HI - (UI*UI+VI*VI+WI*WI)*0.5 ) * RI * ( fG - 1.0 ) / fG;

} 

void rim_orig(double& RI, double& EI, double& PI, double& UI, double& VI, double& WI,
         double RB, double PB, double UB, double VB, double WB,
         double RE, double PE, double UE, double VE, double WE, double GAM) {

	double AGAM  =  (GAM-1.0);
	double BGAM  =  (2.0*sqrt(GAM/AGAM));
	double CGAM  =  (1.0/GAM);
	double DGAM  =  (2.0/AGAM);
	double EGAM  =  (AGAM/(GAM+1.0));
	double GGAM  =  (sqrt(GAM*AGAM));
	double HGAM  =  (AGAM/2.0);
	double FGAM  =  (3.0*GAM-1.0);
	double OGAM  =  (AGAM/(2.0*GAM));
	double QGAM  =  (GAM+1.0);
	double PGAM  =  (QGAM/(2.0*GAM));
	double RGAM  =  (4.0*GAM);
	double SGAM  =  (GAM*AGAM);
	double TGAM  =  (QGAM/2.0);
	double UGAM  =  (sqrt(AGAM/GAM));

	double RF,RS,EF,ES,SBL,SFL,SSL,SEL,D;
	double    eps=1.0e-5;
	double    CB= sqrt(GAM*PB/RB);
	double    CE= sqrt(GAM*PE/RE);
	double    EB= CB*CB/SGAM;
	double    EE= CE*CE/SGAM;
	double    RCB=RB*CB;
	double    RCE=RE*CE;
	double    DU=UB-UE;
      if (DU < -2.0*(CB+CE)/AGAM) {
         //printf(" ATTENTION!!!  VACUUM \n");
         RF=0.0;
         RS=0.0;
         EF=0.0;
         ES=0.0;
         SBL=UB-CB;
         SFL=UB+2.0*CB/AGAM;
         SSL=UE-2.0*CE/AGAM;
         SEL=UE+CE;
         goto lbl9;
      }
      double P=(PB*RCE+PE*RCB+DU*RCB*RCE)/(RCB+RCE);
    lbl5:
           if (P<eps) P=eps;

           double PPB=P/PB;
           if(PB>P) goto lbl1;
           double PKB=PGAM*PPB+OGAM;
           double ZNB=RCB*sqrt(PKB);
           double F1=(P-PB)/ZNB ;
           double FS1=(QGAM*PPB+FGAM)/(RGAM*ZNB*PKB);
           goto lbl2;
    lbl1:      
           double ZFB=CB*exp(log(PPB)*OGAM) ;
           F1=DGAM*(ZFB-CB);
           FS1=ZFB/(GAM*P);
    lbl2:      
           double PPE=P/PE;
           if(PE>P) goto lbl3;
           double PKE=PGAM*PPE+OGAM;
           double ZNE=RCE*sqrt(PKE);
           double F2=(P-PE)/ZNE;
           double FS2=(QGAM*PPE+FGAM)/(RGAM*ZNE*PKE);
           goto lbl4;
    lbl3:      
           double ZFE=CE*exp(log(PPE)*OGAM);
           F2=DGAM*(ZFE-CE);
           FS2=ZFE/(GAM*P);
    lbl4:      
           double DP=(DU-F1-F2)/(FS1+FS2);
           P=P+DP;
       if(fabs(DU-F1-F2)>eps) goto lbl5;
 

      PPB=P/PB;
      PPE=P/PE;

//       ZFB=CB*PPB**OGAM;
//       ZFE=CE*PPE**OGAM;
      ZFB=CB*exp(log(PPB)*OGAM);
      ZFE=CE*exp(log(PPE)*OGAM);
      if (PB>P) goto lbl6 ;
      D=UB-sqrt((TGAM*P+HGAM*PB)/RB);
      double UBD=UB-D;
      double RUBD=RB*UBD;
      RF=RUBD*RUBD/(PB-P+RUBD*UBD);
      double UF=D+RUBD/RF;
      EF=P/(AGAM*RF);
      SBL=D;
      SFL=D;
      goto lbl7;
    lbl6: 
      EF=ZFB*ZFB/SGAM;
      UF=UB+DGAM*(CB-ZFB);
      RF=P/(AGAM*EF);
      SBL=UB-CB;
      SFL=UF-ZFB;
    lbl7: 
      if (PE>P) goto lbl8;
      D=UE+sqrt((TGAM*P+HGAM*PE)/RE);
      double UED=UE-D;
      double RUED=RE*UED;
      RS=RUED*RUED/(PE-P+RUED*UED);
      double US=D+RUED/RS;
      ES=P/(AGAM*RS);
      SEL=D;
      SSL=D;
      goto lbl9;
    lbl8: 
      ES=ZFE*ZFE/SGAM;
      US=UE-DGAM*(CE-ZFE);
      RS=P/(AGAM*ES);
      SSL=US+ZFE ;
      SEL=UE+CE ;
    lbl9: 
// 
// C     compute the interpolation value
      if (SEL<=0.0) {
         RI= RE;
         EI= EE;
         UI= UE;
         VI= VE;
         WI= WE;
         goto lbl157;
      }

      if (SBL>=0.0) {
         RI= RB;
         EI= EB;
         UI= UB;
         VI= VB;
         WI= WB;
         goto lbl157;
      }

      if ((SSL>=0.0)&&(SFL<=0.0)) {
         if (US>=0.0) {
            RI= RF;
            EI= EF;
            UI= UF;
            VI= VB;
            WI= WB;
         } else {
            RI= RS;
            EI= ES;
            UI= US;
            VI= VE;
            WI= WE;
         }
         goto lbl157;
      }

      if (SFL>0.0) {
         UI= (UB+DGAM*GGAM*sqrt(EB))/(1+DGAM);
         VI= VB;
		 WI= WB;
		 EI= (UI*UI)/SGAM;
         RI= RB*exp(log(EI/EB)*(1/AGAM));
       } else {
         UI= (UE-DGAM*GGAM*sqrt(EE))/(1+DGAM);
         VI= VE;
		 WI= WE;
		 EI= (UI*UI)/SGAM;
         RI= RE*exp(log(EI/EE)*(1/AGAM)) ;
      }

  lbl157:
      PI= AGAM*EI*RI;
 
} 

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
		primToCons(r[i], p[i], u[i], ro[i], ru[i], re[i]);
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
		
		for (int iCell = 0, ind = 0; iCell < N; iCell++, ind += 3) {
			ro[iCell] += S->x[ind+0];
			ru[iCell] += S->x[ind+1];
			re[iCell] += S->x[ind+2];
			consToPrim(r[iCell], p[iCell], u[iCell], ro[iCell], ru[iCell], re[iCell]);
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































void calcMatrA(double c2, double u, double** A) {
    double cz = sqrt(c2);
    double kk = GAM-1.0;

	A[0][0] = 0.0;
	A[0][1] = 1.0;
	A[0][2] = 0.0;
	
	A[1][0] = (GAM-3.0)*u*u/2.0;
	A[1][1] = (3.0-GAM)*u;
	A[1][2] = kk;
	
	A[2][0] = kk*(u*u*u)-u*(c2/kk+GAM*u*u/2.0);
	A[2][1] = c2/kk+GAM*u*u/2.0-1.5*kk*u*u;
	A[2][2] = GAM*u;
}

void calcMatrAPM__(double c2, double u, double** Am, double** Ap) {
    double cz = sqrt(c2);
	double ll = fabs(u)+cz;
	
	calcMatrA(c2, u, Am);
	calcMatrA(c2, u, Ap);
	
	Am[0][0] -= ll;
	Am[1][1] -= ll;
	Am[2][2] -= ll;
	
	Ap[0][0] += ll;
	Ap[1][1] += ll;
	Ap[2][2] += ll;
	
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			Am[i][j] *= 0.5;
			Ap[i][j] *= 0.5;
		}
	}
}






int main__(int argc, char** argv) 
{
	memAlloc();
	Solver* S = SolverFactory::create(SOLVER_TYPE_ZEIDEL);
	S->init(N);
	for (int i = 0; i < N; i++) {
		double x = XMIN+i*h;
		if (x < 0.0) {
			r[i] = 1.0;
			p[i] = 1.0;
			u[i] = 0.0;
		} else {
			r[i] = 0.125;
			p[i] = 0.1;
			u[i] = 0.0;
		}
		primToCons(r[i], p[i], u[i], ro[i], ru[i], re[i]);
		double dt = CFL*h/(fabs(u[i])+sqrt(GAM*p[i]/r[i]));
		if (dt < tau) tau = dt;
	}
	printf("TAU = %25.16e\n\n");

	double t = 0.0;
	int step = 0;
	while (t < TMAX) {
		t += tau;
		step++;
		S->zero();
		for (int iCell = 1; iCell < N; iCell++) {
			double ri,ei,pi,ui,vi,wi;
			double rb = r[iCell-1];
			double pb = p[iCell-1];
			double ub = u[iCell-1];
			double vb = 0;
			double wb = 0;
			double re = r[iCell];
			double pe = p[iCell];
			double ue = u[iCell];
			double ve = 0;
			double we = 0;
			rim_orig(ri,ei,pi,ui,vi,wi, 
					 rb, pb, ub, vb, wb, 
					 re, pe, ue, ve, we, GAM);
			double c2 = GAM*pi/ri;
			
		
			calcMatrAPM(c2, ui,  Am,  Ap);
		
			// i-1/2
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					matr[i][j] = -Ap[i][j]/h;
				}
			}
			S->addMatrElement(iCell, iCell-1, matr);
		
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					matr[i][j] = -Am[i][j]/h;
				}
			}
			S->addMatrElement(iCell, iCell, matr);
		
			// i+1/2
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					matr[i][j] = Am[i][j]/h;
				}
			}
			S->addMatrElement(iCell-1, iCell, matr);
		
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					matr[i][j] = Ap[i][j]/h;
				}
			}
			S->addMatrElement(iCell-1, iCell-1, matr);

			vect[0] = ri*ui/h;
			vect[1] = (ri*ui*ui+pi)/h;
			vect[2] = (ri*(ei+ui*ui*0.5)+pi)*ui/h;
			S->addRightElement(iCell, vect);
		
			vect[0] *= -1.0;
			vect[1] *= -1.0;
			vect[2] *= -1.0;
			S->addRightElement(iCell-1, vect);
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
			rim_orig(ri,ei,pi,ui,vi,wi, 
					 rb, pb, ub, vb, wb, 
					 re, pe, ue, ve, we, GAM);
			double c2 = GAM*pi/ri;
			
		
			calcMatrAPM(c2, ui, Am,  Ap);
		
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					matr[i][j] = -Am[i][j]/h;
				}
			}
			S->addMatrElement(iCell, iCell, matr);
		
			vect[0] = ri*ui/h;
			vect[1] = (ri*ui*ui+pi)/h;
			vect[2] = (ri*(ei+ui*ui*0.5)+pi)*ui/h;
			S->addRightElement(iCell, vect);
		
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
			rim_orig(ri,ei,pi,ui,vi,wi, 
					 rb, pb, ub, vb, wb, 
					 re, pe, ue, ve, we, GAM);
			double c2 = GAM*pi/ri;
			 
		
			calcMatrAPM(c2, ui, Am, Ap);
		
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					matr[i][j] = Ap[i][j]/h;
				}
			}
			S->addMatrElement(iCell-1, iCell-1, matr);

			vect[0] = -ri*ui/h;
			vect[1] = -(ri*ui*ui+pi)/h;
			vect[2] = -(ri*(ei+ui*ui*0.5)+pi)*ui/h;
			S->addRightElement(iCell-1, vect);
		}

		for (int iCell = 0; iCell < N; iCell++) {
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					matr[i][j] = (i==j ? 1.0/tau : 0.0);
				}
			}
			S->addMatrElement(iCell, iCell, matr);
		}
		int maxIter = MAX_ITER;
		//S->printMatr();
		S->solve(EPS, maxIter);
		for (int iCell = 0, ind = 0; iCell < N; iCell++, ind += 3) {
			ro[iCell] += S->x[ind+0];
			ru[iCell] += S->x[ind+1];
			re[iCell] += S->x[ind+2];
			consToPrim(r[iCell], p[iCell], u[iCell], ro[iCell], ru[iCell], re[iCell]);
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

