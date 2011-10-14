#ifndef BESSEL_HEADER
#define BESSEL_HEADER

#include<complex>

using namespace std;

//#ifndef bessH
//#define bessH
//#include <complex.h>



//====================================================================
// Collection of Bessel functions and the like, implemented
// in C++ by C. Bond (http://www.crbond.com/), using the
// algorithms in  Zhang and Jin's book "Computation of Special Functions".
// John Wiley and Sons, 1996.
// This is based on the version dated 06/04 which appears to have
// corrected errors in complex Bessel functions.
//====================================================================
namespace CRBond_Bessel
{


 // These used to be defined via #defines! Very bad!
 extern double eps; // 1e-15
 extern double el; // 0.5772156649015329

extern int msta1(double x,int mp);
extern int msta2(double x,int n,int mp);
extern int bessjy01a(double x,double &j0,double &j1,double &y0,double &y1,
    double &j0p,double &j1p,double &y0p,double &y1p);
extern int bessjy01b(double x,double &j0,double &j1,double &y0,double &y1,
    double &j0p,double &j1p,double &y0p,double &y1p);
extern int bessjyna(int n,double x,int &nm,double *jn,double *yn,
    double *jnp,double *ynp);
extern int bessjynb(int n,double x,int &nm,double *jn,double *yn,
    double *jnp,double *ynp);
extern int bessjyv(double v,double x,double &vm,double *jv,double *yv,
    double *jvp,double *yvp);
extern int bessik01a(double x,double &i0,double &i1,double &k0,double &k1,
    double &i0p,double &i1p,double &k0p,double &k1p);
extern int bessik01b(double x,double &i0,double &i1,double &k0,double &k1,
    double &i0p,double &i1p,double &k0p,double &k1p);
extern int bessikna(int n,double x,int &nm,double *in,double *kn,
    double *inp,double *knp);
extern int bessiknb(int n,double x,int &nm,double *in,double *kn,
    double *inp,double *knp);
extern int bessikv(double v,double x,double &vm,double *iv,double *kv,
    double *ivp,double *kvp);
extern int cbessjy01(complex<double> z,complex<double> &cj0,complex<double> &cj1,
    complex<double> &cy0,complex<double> &cy1,complex<double> &cj0p,
    complex<double> &cj1p,complex<double> &cy0p,complex<double> &cy1p);
extern int cbessjyna(int n,complex<double> z,int &nm,complex<double> *cj,
    complex<double> *cy,complex<double> *cjp,complex<double> *cyp);
extern int cbessjynb(int n,complex<double> z,int &nm,complex<double> *cj,
    complex<double> *cy,complex<double> *cjp,complex<double> *cyp);
extern int cbessik01(complex<double>z,complex<double>&ci0,complex<double>&ci1,
    complex<double>&ck0,complex<double>&ck1,complex<double>&ci0p,
    complex<double>&ci1p,complex<double>&ck0p,complex<double>&ck1p);
extern int cbessikna(int n,complex<double> z,int &nm,complex<double> *ci,
    complex<double> *ck,complex<double> *cip,complex<double> *ckp);
extern int cbessiknb(int n,complex<double> z,int &nm,complex<double> *ci,
    complex<double> *ck,complex<double> *cip,complex<double> *ckp);
extern int cbessjyva(double v,complex<double> z,double &vm,complex<double>*cjv,
    complex<double>*cyv,complex<double>*cjvp,complex<double>*cyvp);
extern int cbessikv(double v,complex<double>z,double &vm,complex<double> *civ,
    complex<double> *ckv,complex<double> *civp,complex<double> *ckvp);


}
#endif


