#pragma only

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


extern void addMatr3ToMatr9(double **m9, double **m3, int i, int j, int m3size);
extern void URS(int iCase, double& r, double& p, double& e, double gam);
extern double _max_(double a, double b);

void calcMatrAPM(double c2, double u, double GAM, double** Am, double** Ap);
void calcMatrA(double c2, double u, double GAM, double** A);
void calcMatrL(double c2, double u, double GAM, double** R);
void calcMatrR(double c2, double u, double GAM, double** R);



extern double **matr, *vect, **matr1, **matr9, ***matrM;
extern double **A, **Ap, **Am, **L, **R, ***mA;

