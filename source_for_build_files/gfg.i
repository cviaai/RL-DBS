/* file : gfg.i */
  
/* name of module to use*/
%module oscillator_cpp 

%{ 
    /* Every thing in this file is being copied in  
     wrapper file. We include the C header file necessary 
     to compile the interface */
    #include "gfg.h" 
    #include "nrlib.h"
    #include "nrlib.cxx"
    /* variable declaration*/
    double myvar; 
%} 
%include "gfg.h"

%include "carrays.i"
%include "cpointer.i"
%array_class(double,doubleArray);
%pointer_functions(double ,doubleP);


/* explicitly list functions and variables to be interfaced */
//double myvar; 
//long long int fact(long long int n1); 
//int my_mod(int m, int n); 

/* or if we want to interface all functions then we can simply 
   include header file like this -  
   %include "gfg.h" 
*/
