#ifndef _NRLIB_
#define _NRLIB_

#include <stdio.h>
#include <stddef.h> 
#include <stdlib.h> 
#include <iostream> 
#include <math.h>
#include "my_defs.h"
using namespace std;

//**************************************************************
//***********  supplementary routines from nrutil.c  ***********
//**************************************************************

#define NR_END 1
#define FREE_ARG char*

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	char msg[]="allocation failure in vector()";
	if (!v) nrerror(msg);
	return v-nl+NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	char msg[]="allocation failure in dvector()";
	if (!v) nrerror(msg);
	return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	char msg[]="allocation failure in ivector()";
	if (!v) nrerror(msg);
        return v-nl+NR_END;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
        unsigned long *v;

        v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
        if (!v) nrerror("allocation failure in lvector()");
        return v-nl+NR_END;
}


double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        double **m;

        /* allocate pointers to rows */
        m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
        char msg[]="allocation failure 1 in matrix()";
	if (!m) nrerror(msg);
        m += NR_END;
        m -= nrl;

        /* allocate rows and set pointers to them */
        m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	char msg1[]="allocation failure 2 in matrix()";
        if (!m[nrl]) nrerror(msg1);
        m[nrl] += NR_END;
        m[nrl] -= ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}


void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}


void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}

//**************************************************************
//***********  routines for AR model and AR-spectrum  **********
//**************************************************************

void memcof(double data[], int n, int m, double *xms, double d[])
{
	int k,j,i;
	double p=0.0,*wk1,*wk2,*wkm;

	wk1=dvector(1,n);
	wk2=dvector(1,n);
	wkm=dvector(1,m);
	for (j=1;j<=n;j++) p += SQR(data[j]);
	*xms=p/n;
	wk1[1]=data[1];
	wk2[n-1]=data[n];
	for (j=2;j<=n-1;j++) {
		wk1[j]=data[j];
		wk2[j-1]=data[j];
	}
	for (k=1;k<=m;k++) {
		double num=0.0,denom=0.0;
		for (j=1;j<=(n-k);j++) {
			num += wk1[j]*wk2[j];
			denom += SQR(wk1[j])+SQR(wk2[j]);
		}
		d[k]=2.0*num/denom;
		*xms *= (1.0-SQR(d[k]));
		for (i=1;i<=(k-1);i++)
			d[i]=wkm[i]-d[k]*wkm[k-i];
		if (k == m) {
			free_dvector(wkm,1,m);
			free_dvector(wk2,1,n);
			free_dvector(wk1,1,n);
			return;
		}
		for (i=1;i<=k;i++) wkm[i]=d[i];
		for (j=1;j<=(n-k-1);j++) {
			wk1[j] -= wkm[k]*wk2[j];
			wk2[j]=wk2[j+1]-wkm[k]*wk1[j+1];
		}
	}
	char msg[]="never get here in memcof.";
	nrerror(msg);
}
/* (C) Copr. 1986-92 Numerical Recipes Software )!0,". */

double evlmem(double fdt, double d[], int m, double xms)
{
	int i;
	double sumr=1.0,sumi=0.0;
	double wr=1.0,wi=0.0,wpr,wpi,wtemp,theta;

	theta=6.28318530717959*fdt;
	wpr=cos(theta);
	wpi=sin(theta);
	for (i=1;i<=m;i++) {
		wr=(wtemp=wr)*wpr-wi*wpi;
		wi=wi*wpr+wtemp*wpi;
		sumr -= d[i]*wr;
		sumi -= d[i]*wi;
	}
	return xms/(sumr*sumr+sumi*sumi);
}
/* (C) Copr. 1986-92 Numerical Recipes Software )!0,". */

//*************************************************************
//************   random number generators   *******************
//*************************************************************
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
/* (C) Copr. 1986-92 Numerical Recipes Software )!0,". */

double gasdev(long *idum)
{
	double ran1(long *idum);
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;

	if (*idum < 0) iset=0;
	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}


/* (C) Copr. 1986-92 Numerical Recipes Software )!0,". */

/******* program gives out 2 Gaussian deviates **********/
/*******             of unit rms               **********/
void  gasdev2(long *idum, double *gd1, double *gd2)
{
	double ran1(long *idum);
	double fac,rsq,v1,v2;

	do {
	  v1=2.0*ran1(idum)-1.0;
	  v2=2.0*ran1(idum)-1.0;
	  rsq=v1*v1+v2*v2;
	} while (rsq >= 1.0 || rsq == 0.0);
	fac=sqrt(-2.0*log(rsq)/rsq);
	*gd1=v1*fac;
	*gd2=v2*fac;
}

//*************************************************************
//******************   Fourier Transform   ********************
//*************************************************************

void four1(double data[], unsigned long nn, int isign)
{
        unsigned long n,mmax,m,j,istep,i;
        double wtemp,wr,wpr,wpi,wi,theta;
        double tempr,tempi;

        n=nn << 1;
        j=1;
        for (i=1;i<n;i+=2) {
                if (j > i) {
                        SWAP(data[j],data[i]);
                        SWAP(data[j+1],data[i+1]);
                }
                m=n >> 1;
                while (m >= 2 && j > m) {
                        j -= m;
                        m >>= 1;
                }
                j += m;
        }
        mmax=2;
        while (n > mmax) {
                istep=mmax << 1;
                theta=isign*(6.28318530717959/mmax);
                wtemp=sin(0.5*theta);
                wpr = -2.0*wtemp*wtemp;
                wpi=sin(theta);
                wr=1.0;
                wi=0.0;
                for (m=1;m<mmax;m+=2) {
                        for (i=m;i<=n;i+=istep) {
                                j=i+mmax;
                                tempr=wr*data[j]-wi*data[j+1];
                                tempi=wr*data[j+1]+wi*data[j];
                                data[j]=data[i]-tempr;
                                data[j+1]=data[i+1]-tempi;
                                data[i] += tempr;
                                data[i+1] += tempi;
                        }
                        wr=(wtemp=wr)*wpr-wi*wpi+wr;
                        wi=wi*wpr+wtemp*wpi+wi;
                }
                mmax=istep;
        }
}

void realft(double data[], unsigned long n, int isign)
{
	unsigned long i,i1,i2,i3,i4,np3;
	double c1=0.5,c2,h1r,h1i,h2r,h2i;
	double wr,wi,wpr,wpi,wtemp,theta;

	theta=3.141592653589793/(double) (n>>1);
	if (isign == 1) {
		c2 = -0.5;
		four1(data,n>>1,1);
	} else {
		c2=0.5;
		theta = -theta;
	}
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	np3=n+3;
	for (i=2;i<=(n>>2);i++) {
		i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
		h1r=c1*(data[i1]+data[i3]);
		h1i=c1*(data[i2]-data[i4]);
		h2r = -c2*(data[i2]+data[i4]);
		h2i=c2*(data[i1]-data[i3]);
		data[i1]=h1r+wr*h2r-wi*h2i;
		data[i2]=h1i+wr*h2i+wi*h2r;
		data[i3]=h1r-wr*h2r+wi*h2i;
		data[i4] = -h1i+wr*h2i+wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1) {
		data[1] = (h1r=data[1])+data[2];
		data[2] = h1r-data[2];
	} else {
		data[1]=c1*((h1r=data[1])+data[2]);
		data[2]=c1*(h1r-data[2]);
		four1(data,n>>1,-1);
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software )!0,". */

//*************************************************************
//****************   Savitzky-Golay filter   ******************
//*************************************************************
void lubksb(double **a, int n, int *indx, double b[])
{
	int i,ii=0,ip,j;
	double sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software )!0,". */


#define TINY 1.0e-20;

void ludcmp(double **a, int n, int *indx, double *d)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;

	vv=dvector(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		char msg[]="Singular matrix in routine ludcmp";
		if (big == 0.0) nrerror(msg);
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_dvector(vv,1,n);
}
#undef TINY
/* (C) Copr. 1986-92 Numerical Recipes Software )!0,". */

static float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

void savgol(double c[], int np, int nl, int nr, int ld, int m)
{
	int imj,ipj,j,k,kk,mm,*indx;
	double d,fac,sum,**a,*b;

	char msg[]="bad args in savgol";
	if (np < nl+nr+1 || nl < 0 || nr < 0 || ld > m || nl+nr < m)
	nrerror(msg);
	indx=ivector(1,m+1);
	a=dmatrix(1,m+1,1,m+1);
	b=dvector(1,m+1);
	for (ipj=0;ipj<=(m << 1);ipj++) {
		sum=(ipj ? 0.0 : 1.0);
		for (k=1;k<=nr;k++) sum += pow((double)k,(double)ipj);
		for (k=1;k<=nl;k++) sum += pow((double)-k,(double)ipj);
		mm=FMIN(ipj,2*m-ipj);
		for (imj = -mm;imj<=mm;imj+=2) a[1+(ipj+imj)/2][1+(ipj-imj)/2]=sum;
	}
	ludcmp(a,m+1,indx,&d);
	for (j=1;j<=m+1;j++) b[j]=0.0;
	b[ld+1]=1.0;
	lubksb(a,m+1,indx,b);
	for (kk=1;kk<=np;kk++) c[kk]=0.0;
	for (k = -nl;k<=nr;k++) {
		sum=b[1];
		fac=1.0;
		for (mm=1;mm<=m;mm++) sum += b[mm+1]*(fac *= k);
		kk=((np-k) % np)+1;
		c[kk]=sum;
	}
	free_dvector(b,1,m+1);
	free_dmatrix(a,1,m+1,1,m+1);
	free_ivector(indx,1,m+1);
}

/* (C) Copr. 1986-92 Numerical Recipes Software )!0,". */

/*************************************************************
*********************   Interpolation   **********************
*************************************************************/

void spline(double x[], double y[], int n, double yp1, 
	    double ypn, double y2[])
{
	int i,k;
	double p,qn,sig,un,*u;

	u=dvector(1,n-1);
	if (yp1 > 0.99e30)
		y2[1]=u[1]=0.0;
	else {
		y2[1] = -0.5;
		u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
	}
	for (i=2;i<=n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
	}
	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	for (k=n-1;k>=1;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
	free_dvector(u,1,n-1);
}

/* (C) Copr. 1986-92 Numerical Recipes Software )!0,". */

void splint(double xa[], double ya[], double y2a[], int n, 
	    double x, double *y)
{
	void nrerror(char error_text[]);
	int klo,khi,k;
	double h,b,a;

	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	char msg[]="Bad xa input to routine splint";
	if (h == 0.0) nrerror(msg);
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}
/* (C) Copr. 1986-92 Numerical Recipes Software )!0,". */


/*************************************************************
********************   ODE Integration   *********************
*************************************************************/
//==== Simple RK Integrator, version modified by Arkady ========

void rk(double y[],int n,double &x,double h,
     void (*derivs)(double, double[], double[]), int nt)
{
  int it,i;
  double hh,h6,xh,*yt,*dy0,*dyt,*dym;
  yt=dvector(1,n);    dy0=dvector(1,n);
  dyt=dvector(1,n);   dym=dvector(1,n);

  hh=h*0.5;    h6=h/6.;
  for (it=1;it<=nt;it++)
    {
      xh=x+hh;
      (*derivs)(x,y,dy0);
      for(i=1;i<=n;i++) yt[i]=y[i]+hh*dy0[i];
      (*derivs)(xh,yt,dyt);
      for(i=1;i<=n;i++) yt[i]=y[i]+hh*dyt[i];
      (*derivs)(xh,yt,dym);
      for(i=1;i<=n;i++)
        {
          yt[i]=y[i]+h*dym[i];
          dym[i]=dyt[i]+dym[i];
        }
      x+=h;
      (*derivs)(x,yt,dyt);
      for(i=1;i<=n;i++) y[i]+=h6*(dy0[i]+2.*dym[i]+dyt[i]);
    }
  free_dvector(yt,1,n);     free_dvector(dyt,1,n);
  free_dvector(dy0,1,n);    free_dvector(dym,1,n);
}
//==================================================================

/********************** extrapolation step integrator ***********************/

#define KMAXX 8
#define IMAXX (KMAXX+1)
#define SAFE1 0.25
#define SAFE2 0.7
#define REDMAX 1.0e-5
#define REDMIN 0.7
#define TINY 1.0e-30
#define SCALMX 0.1


void mmid(double y[], double dydx[], int nvar, double xs, double htot, int nstep,
	double yout[], void (*derivs)(double, double[], double[]))
{
	int n,i;
	double x,swap,h2,h,*ym,*yn;

	ym=dvector(1,nvar);
	yn=dvector(1,nvar);
	h=htot/nstep;
	for (i=1;i<=nvar;i++) {
		ym[i]=y[i];
		yn[i]=y[i]+h*dydx[i];
	}
	x=xs+h;
	(*derivs)(x,yn,yout);
	h2=2.0*h;
	for (n=2;n<=nstep;n++) {
		for (i=1;i<=nvar;i++) {
			swap=ym[i]+h2*yout[i];
			ym[i]=yn[i];
			yn[i]=swap;
		}
		x += h;
		(*derivs)(x,yn,yout);
	}
	for (i=1;i<=nvar;i++)
		yout[i]=0.5*(ym[i]+yn[i]+h*yout[i]);
	free_dvector(yn,1,nvar);
	free_dvector(ym,1,nvar);
}

double **d_bsstep,*x_bsstep; // renamed by us on 26.02.14

void pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv)
{
	int k1,j;
	double q,f2,f1,delta,*c;

	c=dvector(1,nv);
	x_bsstep[iest]=xest;
	for (j=1;j<=nv;j++) dy[j]=yz[j]=yest[j];
	if (iest == 1) {
	  for (j=1;j<=nv;j++)
	    d_bsstep[j][1]=yest[j];
	} else {
		for (j=1;j<=nv;j++) c[j]=yest[j];
		for (k1=1;k1<iest;k1++) {
			delta=1.0/(x_bsstep[iest-k1]-xest);
			f1=xest*delta;
			f2=x_bsstep[iest-k1]*delta;
			for (j=1;j<=nv;j++) {
				q=d_bsstep[j][k1];
				d_bsstep[j][k1]=dy[j];
				delta=c[j]-q;
				dy[j]=f1*delta;
				c[j]=f2*delta;
				yz[j] += dy[j];
			}
		}
		for (j=1;j<=nv;j++) d_bsstep[j][iest]=dy[j];
	}
	free_dvector(c,1,nv);
}

void bsstep(double y[], double dydx[], int nv, double *xx, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []))
{
	void mmid(double y[], double dydx[], int nvar, double xs, double htot,
		int nstep, double yout[], void (*derivs)(double, double[], double[]));
	void pzextr(int iest, double xest, double yest[], double yz[], double dy[],
		int nv);
	int i,iq,k,kk,km;
	static int first=1,kmax,kopt;
	static double epsold = -1.0,xnew;
	double eps1,errmax,fact,h,red,scale,work,wrkmin,xest;
	double *err,*yerr,*ysav,*yseq;
	static double a[IMAXX+1];
	static double alf[KMAXX+1][KMAXX+1];
	static int nseq[IMAXX+1]={0,2,4,6,8,10,12,14,16,18};
	int reduct,exitflag=0;

	d_bsstep=dmatrix(1,nv,1,KMAXX);
	err=dvector(1,KMAXX);
	x_bsstep=dvector(1,KMAXX);
	yerr=dvector(1,nv);
	ysav=dvector(1,nv);
	yseq=dvector(1,nv);
	if (eps != epsold) {
		*hnext = xnew = -1.0e29;
		eps1=SAFE1*eps;
		a[1]=nseq[1]+1;
		for (k=1;k<=KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
		for (iq=2;iq<=KMAXX;iq++) {
			for (k=1;k<iq;k++)
				alf[k][iq]=pow(eps1,(a[k+1]-a[iq+1])/
					((a[iq+1]-a[1]+1.0)*(2*k+1)));
		}
		epsold=eps;
		for (kopt=2;kopt<KMAXX;kopt++)
			if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
		kmax=kopt;
	}
	h=htry;
	for (i=1;i<=nv;i++) ysav[i]=y[i];
	if (*xx != xnew || h != (*hnext)) {
		first=1;
		kopt=kmax;
	}
	reduct=0;
	for (;;) {
		for (k=1;k<=kmax;k++) {
			xnew=(*xx)+h;
			char msg[]="step size underflow in bsstep";
			if (xnew == (*xx)) nrerror(msg);
			mmid(ysav,dydx,nv,*xx,h,nseq[k],yseq,derivs);
			xest=SQR(h/nseq[k]);
			pzextr(k,xest,yseq,y,yerr,nv);
			if (k != 1) {
				errmax=TINY;
				for (i=1;i<=nv;i++) errmax=MAX(errmax,fabs(yerr[i]/yscal[i]));
				errmax /= eps;
				km=k-1;
				err[km]=pow(errmax/SAFE1,1.0/(2*km+1));
			}
			if (k != 1 && (k >= kopt-1 || first)) {
				if (errmax < 1.0) {
					exitflag=1;
					break;
				}
				if (k == kmax || k == kopt+1) {
					red=SAFE2/err[km];
					break;
				}
				else if (k == kopt && alf[kopt-1][kopt] < err[km]) {
						red=1.0/err[km];
						break;
					}
				else if (kopt == kmax && alf[km][kmax-1] < err[km]) {
						red=alf[km][kmax-1]*SAFE2/err[km];
						break;
					}
				else if (alf[km][kopt] < err[km]) {
					red=alf[km][kopt-1]/err[km];
					break;
				}
			}
		}
		if (exitflag) break;
		red=MIN(red,REDMIN);
		red=MAX(red,REDMAX);
		h *= red;
		reduct=1;
	}
	*xx=xnew;
	*hdid=h;
	first=0;
	wrkmin=1.0e35;
	for (kk=1;kk<=km;kk++) {
		fact=MAX(err[kk],SCALMX);
		work=fact*a[kk+1];
		if (work < wrkmin) {
			scale=fact;
			wrkmin=work;
			kopt=kk+1;
		}
	}
	*hnext=h/scale;
	if (kopt >= k && kopt != kmax && !reduct) {
		fact=MAX(scale/alf[kopt-1][kopt],SCALMX);
		if (a[kopt+1]*fact <= wrkmin) {
			*hnext=h/fact;
			kopt++;
		}
	}
	free_dvector(yseq,1,nv);
	free_dvector(ysav,1,nv);
	free_dvector(yerr,1,nv);
	free_dvector(x_bsstep,1,KMAXX);
	free_dvector(err,1,KMAXX);
	free_dmatrix(d_bsstep,1,KMAXX,1,KMAXX);
}
#undef KMAXX
#undef IMAXX
#undef SAFE1
#undef SAFE2
#undef REDMAX
#undef REDMIN
#undef TINY
#undef SCALMX
/* (C) Copr. 1986-92 Numerical Recipes Software )!0,". */


/************************* Runge-Kutta integrator ***************************/

#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4


void rkck(double y[], double dydx[], int n, double x, double h, double yout[],
	double yerr[], void (*derivs)(double, double [], double []))
{
	int i;
	static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
		b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
		b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
		b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
		b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
		c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
		dc5 = -277.0/14336.0;
	double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
		dc4=c4-13525.0/55296.0,dc6=c6-0.25;
	double *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;

	ak2=dvector(1,n);
	ak3=dvector(1,n);
	ak4=dvector(1,n);
	ak5=dvector(1,n);
	ak6=dvector(1,n);
	ytemp=dvector(1,n);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+b21*h*dydx[i];
	(*derivs)(x+a2*h,ytemp,ak2);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
	(*derivs)(x+a3*h,ytemp,ak3);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
	(*derivs)(x+a4*h,ytemp,ak4);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
	(*derivs)(x+a5*h,ytemp,ak5);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
	(*derivs)(x+a6*h,ytemp,ak6);
	for (i=1;i<=n;i++)
		yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
	for (i=1;i<=n;i++)
		yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
	free_dvector(ytemp,1,n);
	free_dvector(ak6,1,n);
	free_dvector(ak5,1,n);
	free_dvector(ak4,1,n);
	free_dvector(ak3,1,n);
	free_dvector(ak2,1,n);
}


void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []))
{
	void rkck(double y[], double dydx[], int n, double x, double h,
		double yout[], double yerr[], void (*derivs)(double, double [], double []));
	int i;
	double errmax,h,xnew,*yerr,*ytemp;

	yerr=dvector(1,n);
	ytemp=dvector(1,n);
	h=htry;
	for (;;) {
		rkck(y,dydx,n,*x,h,ytemp,yerr,derivs);
		errmax=0.0;
		for (i=1;i<=n;i++) errmax=MAX(errmax,fabs(yerr[i]/yscal[i]));
		errmax /= eps;
		if (errmax > 1.0) {
			h=SAFETY*h*pow(errmax,PSHRNK);
			if (h < 0.1*h) h *= 0.1;
			xnew=(*x)+h;
			char msg[]="stepsize underflow in rkqs";
			if (xnew == *x) nrerror(msg);
			continue;
		} else {
			if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
			else *hnext=5.0*h;
			*x += (*hdid=h);
			for (i=1;i<=n;i++) y[i]=ytemp[i];
			break;
		}
	}
	free_dvector(ytemp,1,n);
	free_dvector(yerr,1,n);
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON
/* (C) Copr. 1986-92 Numerical Recipes Software )!0,". */

/*********************** Linear regression ************************/

double gammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}
//--------------------------------------------------------------

#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

void gcf(double *gammcf, double a, double x, double *gln)
{
	double gammln(double xx);
	void nrerror(char error_text[]);
	int i;
	double an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	char msg[]="a too large, ITMAX too small in gcf";
	if (i > ITMAX) nrerror(msg);
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}
#undef ITMAX
#undef EPS
#undef FPMIN


#define ITMAX 100
#define EPS 3.0e-7

void gser(double *gamser, double a, double x, double *gln)
{
	double gammln(double xx);
	void nrerror(char error_text[]);
	int n;
	double sum,del,ap;

	char msg1[]="x less than 0 in routine gser";
	char msg2[]="a too large, ITMAX too small in routine gser";
	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) nrerror(msg1);
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		nrerror(msg2);
		return;
	}
}
#undef ITMAX
#undef EPS
//------------------------------------------------------

double gammq(double a, double x)
{
	void gcf(double *gammcf, double a, double x, double *gln);
	void gser(double *gamser, double a, double x, double *gln);
	void nrerror(char error_text[]);
	double gamser,gammcf,gln;

	char msg[]="Invalid arguments in routine gammq";
	if (x < 0.0 || a <= 0.0) 
	  nrerror(msg);
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}
//------------------------------------------------------

void fit(double x[], double y[], int ndata, double sig[], 
	 int mwt, double *a, double *b, double *siga, 
	 double *sigb, double *chi2, double *q)
{
	double gammq(double a, double x);
	int i;
	double wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;

	*b=0.0;
	if (mwt) {
		ss=0.0;
		for (i=1;i<=ndata;i++) {
			wt=1.0/SQR(sig[i]);
			ss += wt;
			sx += x[i]*wt;
			sy += y[i]*wt;
		}
	} else {
		for (i=1;i<=ndata;i++) {
			sx += x[i];
			sy += y[i];
		}
		ss=ndata;
	}
	sxoss=sx/ss;
	if (mwt) {
		for (i=1;i<=ndata;i++) {
			t=(x[i]-sxoss)/sig[i];
			st2 += t*t;
			*b += t*y[i]/sig[i];
		}
	} else {
		for (i=1;i<=ndata;i++) {
			t=x[i]-sxoss;
			st2 += t*t;
			*b += t*y[i];
		}
	}
	*b /= st2;
	*a=(sy-sx*(*b))/ss;
	*siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
	*sigb=sqrt(1.0/st2);
	*chi2=0.0;
	if (mwt == 0) {
		for (i=1;i<=ndata;i++)
			*chi2 += SQR(y[i]-(*a)-(*b)*x[i]);
		*q=1.0;
		sigdat=sqrt((*chi2)/(ndata-2));
		*siga *= sigdat;
		*sigb *= sigdat;
	} else {
		for (i=1;i<=ndata;i++)
			*chi2 += SQR((y[i]-(*a)-(*b)*x[i])/sig[i]);
		*q=gammq(0.5*(ndata-2),0.5*(*chi2));
	}
}

// simple variant of fit routine for equidistant data. 

void simple_fit(double x_0, double xstep, 
		double y[], int ndata,
		double *a, double *b 
		//	double *siga, double *sigb, double *chi2
		)
{
	double gammq(double a, double x);
	int i;
	double wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;

	*b=0.0;
	for (i=1;i<=ndata;i++) 
	  {
	    sx += x_0+xstep*(i-1);
	    sy += y[i];
	  }
	ss=ndata;
	sxoss=sx/ss;

	for (i=1;i<=ndata;i++) 
	  {
	    t=x_0+xstep*(i-1)-sxoss;
	    st2 += t*t;
	    *b += t*y[i];
	  }

	*b /= st2;
	*a=(sy-sx*(*b))/ss;
	/*
	*siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
	*sigb=sqrt(1.0/st2);
	*chi2=0.0;

	for (i=1;i<=ndata;i++)
	  *chi2 += SQR(y[i]-(*a)-(*b)*(x_0+xstep*(i-1)-sxoss));
	sigdat=sqrt((*chi2)/(ndata-2));
	*siga *= sigdat;
	*sigb *= sigdat;
	*/
}

/**************************** Matrices operations ***************************/ 


#define RADIX 2.0

void balanc(double **a, int n)
{
	int last,j,i;
	double s,r,g,f,c,sqrdx;

	sqrdx=RADIX*RADIX;
	last=0;
	while (last == 0) {
		last=1;
		for (i=1;i<=n;i++) {
			r=c=0.0;
			for (j=1;j<=n;j++)
				if (j != i) {
					c += fabs(a[j][i]);
					r += fabs(a[i][j]);
				}
			if (c && r) {
				g=r/RADIX;
				f=1.0;
				s=c+r;
				while (c<g) {
					f *= RADIX;
					c *= sqrdx;
				}
				g=r*RADIX;
				while (c>g) {
					f /= RADIX;
					c /= sqrdx;
				}
				if ((c+r)/f < 0.95*s) {
					last=0;
					g=1.0/f;
					for (j=1;j<=n;j++) a[i][j] *= g;
					for (j=1;j<=n;j++) a[j][i] *= f;
				}
			}
		}
	}
}
#undef RADIX

#define SWAP(g,h) {y=(g);(g)=(h);(h)=y;}

void elmhes(double **a, int n)
{
	int m,j,i;
	double y,x;

	for (m=2;m<n;m++) {
		x=0.0;
		i=m;
		for (j=m;j<=n;j++) {
			if (fabs(a[j][m-1]) > fabs(x)) {
				x=a[j][m-1];
				i=j;
			}
		}
		if (i != m) {
			for (j=m-1;j<=n;j++) SWAP(a[i][j],a[m][j])
			for (j=1;j<=n;j++) SWAP(a[j][i],a[j][m])
		}
		if (x) {
			for (i=m+1;i<=n;i++) {
				if ((y=a[i][m-1]) != 0.0) {
					y /= x;
					a[i][m-1]=y;
					for (j=m;j<=n;j++)
						a[i][j] -= y*a[m][j];
					for (j=1;j<=n;j++)
						a[j][m] += y*a[j][i];
				}
			}
		}
	}
}
#undef SWAP

#define NRANSI

void hqr(double **a, int n, double wr[], double wi[])
{
	int nn,m,l,k,j,its,i,mmin;
	double z,y,x,w,v,u,t,s,r,q,p,anorm;

	char msg[]="Too many iterations in hqr";
	anorm=0.0;
	for (i=1;i<=n;i++)
		for (j=IMAX(i-1,1);j<=n;j++)
			anorm += fabs(a[i][j]);
	nn=n;
	t=0.0;
	while (nn >= 1) {
		its=0;
		do {
			for (l=nn;l>=2;l--) {
				s=fabs(a[l-1][l-1])+fabs(a[l][l]);
				if (s == 0.0) s=anorm;
				if ((double)(fabs(a[l][l-1]) + s) == s) break;
			}
			x=a[nn][nn];
			if (l == nn) {
				wr[nn]=x+t;
				wi[nn--]=0.0;
			} else {
				y=a[nn-1][nn-1];
				w=a[nn][nn-1]*a[nn-1][nn];
				if (l == (nn-1)) {
					p=0.5*(y-x);
					q=p*p+w;
					z=sqrt(fabs(q));
					x += t;
					if (q >= 0.0) {
						z=p+SIGN(z,p);
						wr[nn-1]=wr[nn]=x+z;
						if (z) wr[nn]=x-w/z;
						wi[nn-1]=wi[nn]=0.0;
					} else {
						wr[nn-1]=wr[nn]=x+p;
						wi[nn-1]= -(wi[nn]=z);
					}
					nn -= 2;
				} else {
					if (its == 30) nrerror(msg);
					if (its == 10 || its == 20) {
						t += x;
						for (i=1;i<=nn;i++) a[i][i] -= x;
						s=fabs(a[nn][nn-1])+fabs(a[nn-1][nn-2]);
						y=x=0.75*s;
						w = -0.4375*s*s;
					}
					++its;
					for (m=(nn-2);m>=l;m--) {
						z=a[m][m];
						r=x-z;
						s=y-z;
						p=(r*s-w)/a[m+1][m]+a[m][m+1];
						q=a[m+1][m+1]-z-r-s;
						r=a[m+2][m+1];
						s=fabs(p)+fabs(q)+fabs(r);
						p /= s;
						q /= s;
						r /= s;
						if (m == l) break;
						u=fabs(a[m][m-1])*(fabs(q)+fabs(r));
						v=fabs(p)*(fabs(a[m-1][m-1])+fabs(z)+fabs(a[m+1][m+1]));
						if ((double)(u+v) == v) break;
					}
					for (i=m+2;i<=nn;i++) {
						a[i][i-2]=0.0;
						if (i != (m+2)) a[i][i-3]=0.0;
					}
					for (k=m;k<=nn-1;k++) {
						if (k != m) {
							p=a[k][k-1];
							q=a[k+1][k-1];
							r=0.0;
							if (k != (nn-1)) r=a[k+2][k-1];
							if ((x=fabs(p)+fabs(q)+fabs(r)) != 0.0) {
								p /= x;
								q /= x;
								r /= x;
							}
						}
						if ((s=SIGN(sqrt(p*p+q*q+r*r),p)) != 0.0) {
							if (k == m) {
								if (l != m)
								a[k][k-1] = -a[k][k-1];
							} else
								a[k][k-1] = -s*x;
							p += s;
							x=p/s;
							y=q/s;
							z=r/s;
							q /= p;
							r /= p;
							for (j=k;j<=nn;j++) {
								p=a[k][j]+q*a[k+1][j];
								if (k != (nn-1)) {
									p += r*a[k+2][j];
									a[k+2][j] -= p*z;
								}
								a[k+1][j] -= p*y;
								a[k][j] -= p*x;
							}
							mmin = nn<k+3 ? nn : k+3;
							for (i=l;i<=mmin;i++) {
								p=x*a[i][k]+y*a[i][k+1];
								if (k != (nn-1)) {
									p += z*a[i][k+2];
									a[i][k+2] -= p*r;
								}
								a[i][k+1] -= p*q;
								a[i][k] -= p;
							}
						}
					}
				}
			}
		} while (l < nn-1);
	}
}
#undef NRANSI

/********************** Special functions   **************************/

double bessi0(double x)
{
	double ax,ans;
	double y;

	if ((ax=fabs(x)) < 3.75) {
		y=x/3.75;
		y*=y;
		ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
			+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
	} else {
		y=3.75/ax;
		ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
			+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
			+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
			+y*0.392377e-2))))))));
	}
	return ans;
}
/********************** Sorting   **************************/
#define NRANSI
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50

void sort(unsigned long n, double arr[])
{
	unsigned long i,ir=n,j,k,l=1,*istack;
	int jstack=0;
	double a,temp;

	istack=lvector(1,NSTACK);
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j];
				for (i=j-1;i>=l;i--) {
					if (arr[i] <= a) break;
					arr[i+1]=arr[i];
				}
				arr[i+1]=a;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(arr[k],arr[l+1])
			if (arr[l] > arr[ir]) {
				SWAP(arr[l],arr[ir])
			}
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir])
			}
			if (arr[l] > arr[l+1]) {
				SWAP(arr[l],arr[l+1])
			}
			i=l+1;
			j=ir;
			a=arr[l+1];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j]);
			}
			arr[l+1]=arr[j];
			arr[j]=a;
			jstack += 2;
			if (jstack > NSTACK) nrerror("NSTACK too small in sort.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free_lvector(istack,1,NSTACK);
}
#undef M
#undef NSTACK
#undef SWAP
#undef NRANSI



#endif /* _NRLIB_ */
