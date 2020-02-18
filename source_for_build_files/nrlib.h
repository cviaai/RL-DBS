/**************************************************************
************  supplementary routines from nrutil.c  ***********
**************************************************************/


void nrerror(char error_text[]);
/* Numerical Recipes standard error handler */

float *vector(long nl, long nh);
/* allocate a float vector with subscript range v[nl..nh] */

double *dvector(long nl, long nh);
/* allocate a double vector with subscript range v[nl..nh] */

int *ivector(long nl, long nh);
/* allocate an integer vector with subscript range v[nl..nh] */

unsigned long *lvector(long nl, long nh);
/* allocate an unsigned long vector with subscript range v[nl..nh] */

double **dmatrix(long nrl, long nrh, long ncl, long nch);
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */

void free_vector(float *v, long nl, long nh);
/* free a float vector allocated with vector() */

void free_dvector(double *v, long nl, long nh);
/* free a double vector allocated with dvector() */

void free_ivector(int *v, long nl, long nh);
/* free an int vector allocated with ivector() */

void free_lvector(unsigned long *v, long nl, long nh);
/* free an unsigned long vector allocated with lvector() */

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
/* free a double matrix allocated by dmatrix() */


/**************************************************************
************  routines for AR model and AR-spectrum  **********
**************************************************************/

void memcof(double data[], int n, int m, double *xms, double d[]);

double evlmem(double fdt, double d[], int m, double xms);


/************************************************************
************   random number generators   *******************
************************************************************/
double ran1(long *idum);
double gasdev(long *idum);
void  gasdev2(long *idum, double *gd1, double *gd2);

/*************************************************************
******************   Fourier Transform   *********************
**************************************************************/

void four1(double data[], unsigned long nn, int isign);

void realft(double data[], unsigned long n, int isign);

/*************************************************************
****************   Savitzky-Golay filter   *******************
*************************************************************/


void savgol(double c[], int np, int nl, int nr, int ld, int m);

/*************************************************************
*********************   Interpolation   **********************
*************************************************************/

void spline(double x[], double y[], int n, double yp1, 
	    double ypn, double y2[]);

void splint(double xa[], double ya[], double y2a[], int n, 
	    double x, double *y);


/*************************************************************
********************   ODE Integration   *********************
*******************
******************************************/
//==== Simple RK Integrator, version modified by Arkady ========
void rk(double y[],int n,double &x,double h,
     void (*derivs)(double, double[], double[]), int nt);

      //  Extrapolation integrator (Bulirsch-Stoer method)
void bsstep(double y[], double dydx[], int nv, double *xx, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []));
     // Runge-Kutta integrator
void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []));


/*********************** Linear regression ************************/
void fit(double x[], double y[], int ndata, double sig[], 
	 int mwt, double *a, double *b, double *siga, 
	 double *sigb, double *chi2, double *q);

void simple_fit(double x_0, double xstep, 
		double y[], int ndata,
		double *a, double *b 
		//	double *siga, double *sigb, double *chi2
		);

/**************************** Matrices operations ***************************/ 

void balanc(double **a, int n);
void elmhes(double **a, int n);
void hqr(double **a, int n, double wr[], double wi[]);

/********************** Special functions   **************************/

double bessi0(double x);

/********************** Sorting   **************************/
void sort(unsigned long n, double arr[]);

