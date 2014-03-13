#pragma only

extern void roe_orig(double& RI, double& EI, double& PI, double& UI, double& VI, double& WI,
	double RB, double PB, double UB, double VB, double WB,
	double RE, double PE, double UE, double VE, double WE, double GAM);

extern void rim_orig(double& RI, double& EI, double& PI, double& UI, double& VI, double& WI,
	double RB, double PB, double UB, double VB, double WB,
	double RE, double PE, double UE, double VE, double WE, double GAM);

extern void addMatr3ToMatr9(double **m9, double **m3, int i, int j);
extern void URS(int iCase, double& r, double& p, double& e, double gam);