#include <iostream> 
#include <fstream> 
#include <math.h>
#include "my_defs.h"
#include "nrlib.h"

using namespace std;
/*************  important parameters  **************/
int nosc = 1000;           // number of oscillators
int neqn = 2*nosc;     // number of equations
double  epsilon=0.03;         // coupling parameter

//                         here it is just taken from the blue sky, but it 
//                         exactly this quantity that shall be determined by learning
//                         (along with the instants when the kick shall be applied) 
double frrms=0.1;    // width of the distribution of natural frequencies 
const int ntrans=25000, npt = 5000; // length of the transient and of the main integration
const double tstep=0.2;    // integration step
double *I;
double mfx,mfy, t;

void BvdP(double t, double y[], double ydot[])
{  
  double mf;
  int i;
  
  for (i=1, mf=0. ; i <= nosc; i++) mf+=y[2*i-1]; mf/=nosc;  // mean field is just an average of the first variable
  for (i=1; i<=nosc; i++) 
    {
      ydot[2*i-1]=y[2*i-1]-SQR(y[2*i-1])*y[2*i-1]/3.-y[2*i]+I[i] + epsilon*mf;
      ydot[2*i]=0.1*(y[2*i-1]-0.8*y[2*i]+0.7);
    }   
}

  // now comes the perturbation
   // reset for all x-variables
double *  Pertrubation(double* y, double dkick)
{  
  int k;
  //Make simple pertu15000bation for all x-values
  for(k=1; k<=nosc; k++) 
    y[2*k-1] += dkick; 
  return y;

}

//Calculate mean_field x
 double Calc_mfx(double* y)
{  
  int k;
  for (k=1, mfx=0; k<=nosc; k++)  mfx+=y[2*k-1];  
  return mfx/nosc;

}

//Calculate mean_field y
 double Calc_mfy(double* y)
{  
  int k;
  for (k=1, mfy=0; k<=nosc; k++)  mfy+=y[2*k];  
  return mfy/nosc;

}
//Make_step via rk
double *  Make_step(double *y)
{
  
  int i, k;
  rk(y,neqn,t,tstep,BvdP,1); 
  return y;
}
//Init params
double * init (int nosc_, double epsilon_, double frrms_ ) 
{
  nosc = nosc_;           // number of oscillators
  neqn = 2*nosc;     // number of equations
  epsilon = epsilon_;    // coupling parameter
  frrms = frrms_; // width of the distribution of natural frequencies 
  double *y;
  int i, k;
  long iseed1=-51;  // initializati10000on of the random number generator

  y=dvector(1,neqn);  I=dvector(1,nosc);  // reserve space for arrays
  // here I use routines from Numerical Receipes to have arrays with the index
  // running from 1 to neqn or from 1 to nosc 
  for (i=1; i<=nosc; i++) I[i]=0.6 + frrms*gasdev(&iseed1);  // Neuron frequencies are determined
  //by the external current; here it is taken as 0.6 + Gaussian distribution
  for (i=1; i<=neqn; i++) y[i]=2*ran1(&iseed1)-1.;          // Initial conditions for doff equations: just random 
	
  t=0.; // rk(y,neqn,t,tstep,BvdP,ntrans);     // transient: this is time before our system comes to its limit cycle,
  // i.e. all oscillators come to their limit cyces and synchronize
  // Function rk is the Runge-Kutta integrator
 

  return y;
}
  
