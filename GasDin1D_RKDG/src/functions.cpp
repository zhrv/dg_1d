#include <cstdio>
#include <cmath>

double **matr, *vect, **matr1, **matr9, ***matrM;
double **A, **Ap, **Am, **L, **R, ***mA;

double _max_(double a, double b) {
	if (a > b) return a;
	return b;
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

void addMatr3ToMatr9(double **m9, double **m3, int i, int j, int m3size) {

	int ii = i * m3size;
	int jj = j * m3size;
	for (int i1 = 0; i1 < m3size; i1++) {
		for (int j1 = 0; j1 < m3size; j1++) {
			m9[ii + i1][jj + j1] = m9[ii + i1][jj + j1] + m3[i1][j1];
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

/*
// функция расчёта конвективных потоков методом КИРа
REAL Galerkin::FLUX_KIR(Param &rP, Param &rE, const REAL fFieldP, const REAL fFieldE, const VECTOR &vN, const INDEX idField, const POINT &pt)
{
	REAL fFlux;

	static Param rB;
	static VECTOR_5 vQFP, vQFE;
	static VECTOR_5  vUP, vUE;

	rP.UN = rP.U * vN.x + rP.V * vN.y + rP.W * vN.z;
	rE.UN = rE.U * vN.x + rE.V * vN.y + rE.W * vN.z;

	{{ // Схема ROE
		rB.GAMA = rP.GAMA;

		REAL fG = rP.GAMA;

		REAL fSP = sqrt(rP.R);
		REAL fSE = sqrt(rE.R);
		REAL fS_ = 1.0 / (fSP + fSE);

		rB.R = fSP * fSE;

		rB.U = (fSP * rP.U + fSE * rE.U) * fS_;
		rB.V = (fSP * rP.V + fSE * rE.V) * fS_;
		rB.W = (fSP * rP.W + fSE * rE.W) * fS_;

		rB.E = (fSP * rP.E + fSE * rE.E) * fS_;
		rB.T = (fSP * rP.T + fSE * rE.T) * fS_;

		rP.H = rP.E + rP.EK() + rP.P / rP.R;
		rE.H = rE.E + rE.EK() + rE.P / rE.R;

		rB.H = (fSP * rP.H + fSE * rE.H) * fS_;

		rB.P = (rB.H - rB.EK()) * rB.R * (fG - 1.0) / fG;

		rB.CZ = sqrt(fG * rB.P / rB.R);

		rB.UN = rB.U * vN.x + rB.V * vN.y + rB.W * vN.z;
	}}

	static MATRIX_5 mA;
	CalcMatrixA(rB, vN, mA);

	vUP(0) = 1.0;
	vUP(1) = rP.U;
	vUP(2) = rP.V;
	vUP(3) = rP.W;
	vUP(4) = rP.E + rP.EK();
	vUP *= rP.R;

	vUE(0) = 1.0;
	vUE(1) = rE.U;
	vUE(2) = rE.V;
	vUE(3) = rE.W;
	vUE(4) = rE.E + rE.EK();
	vUE *= rE.R;

	vQFP = vUP;
	vQFP *= rP.UN;
	vQFP(1) += rP.P * vN.x;
	vQFP(2) += rP.P * vN.y;
	vQFP(3) += rP.P * vN.z;
	vQFP(4) += rP.P * rP.UN;

	vQFE = vUE;
	vQFE *= rE.UN;
	vQFE(1) += rE.P * vN.x;
	vQFE(2) += rE.P * vN.y;
	vQFE(3) += rE.P * vN.z;
	vQFE(4) += rE.P * rE.UN;

	static VECTOR_5  vQF;
	// vQF = 0.5 * ( vQFE + vQFP + mA * ( vUP - vUE ) );
	vQF = vUP;
	vQF -= vUE;
	vQF = mA * vQF;

	vQF += vQFE;
	vQF += vQFP;
	vQF *= 0.5;

	fFlux = vQF(idField - 1);
	return fFlux;

} // FLUX_KIR
*/


void flux_cir( double &fr, double &fu, double &fe,
	double rb, double pb, double ub, double re, double pe, double ue, double GAM ) 
{
	double ri, ei, pi, ui, vi, wi, vb = 0, ve = 0, wb = 0, we = 0;
	double AGAM = GAM - 1.0;
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

	ri = re-rb;
	ui = re*ue-rb*ub;
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
	double ri, ei, pi, ui, vi, wi, VB=0, VE=0, WB=0, WE=0;
	double AGAM = GAM - 1.0;
	rim_orig(ri, ei, pi, ui, vi, wi,
		RB, PB, UB, VB, WB,
		RE, PE, UE, VE, WE, GAM);
	fr = ri*ui;
	fu = ri*ui*ui + pi;
	fe = (pi / AGAM + ri*ui*ui*0.5 + pi)*ui;

}



