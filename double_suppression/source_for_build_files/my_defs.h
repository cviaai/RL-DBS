#define ON 1
#define OFF 0
#define TRUE 1

const double pi = 3.14159265358979323846;
const double pi2 = pi+pi;

//*****************************************************************
//**************** inline functions definitions *******************

      
template <class X> inline X MIN(X a, X b) { return (a<b ? a : b);}
   
template <class X> inline X MAX(X a, X b) { return (a>b ? a : b);}

template <class X> inline X ABS(X a) { return (a>0 ? a : -a);}

template <class X> inline X SQR(X a) { return (a*a);}

template <class X> inline void SWAP(X &a,X &b) 
  { X temp; temp=a; a=b; b=temp;}

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
