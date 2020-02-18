

/*******************************************************/
// this function codes the ordinary differential equations for the coupled system
void BvdP(double t, double y[], double ydot[]);
double * Pertrubation(double* y, double dkick);
double *  Make_step(double *y);
double * init (int nosc_, double epsilon_, double frrms_);
double Calc_mfx(double* y);
double Calc_mfy(double* y);
