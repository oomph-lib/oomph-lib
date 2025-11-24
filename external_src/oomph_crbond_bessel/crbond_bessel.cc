
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
// Collection of Bessel functions and the like
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


#include "crbond_bessel.h"
#include <cstdlib>

//====================================================================
/// Collection of Bessel functions and the like, implemented
/// in C++ by C. Bond (http://www.crbond.com/), using the
/// algorithms in  Zhang and Jin's book "Computation of Special Functions".
/// John Wiley and Sons, 1996.
/// This is based on the version dated 06/04 which appears to have
/// corrected errors in complex Bessel functions.
//====================================================================
namespace CRBond_Bessel
{

 // Small number
 double eps=1.0e-15;

 // Euler's constant
 double el=0.5772156649015329;

//  bessik.cpp -- computation of modified Bessel functions In, Kn
//      and their derivatives. Algorithms and coefficient values from
//      "Computation of Special Functions", Zhang and Jin, John
//      Wiley and Sons, 1996.
//
//  (C) 2003, C. Bond. All rights reserved.
//
//MH#include <math.h>
//MH#include "bessel.h"

double gamma(double x);

int bessik01a(double x,double &i0,double &i1,double &k0,double &k1,
    double &i0p,double &i1p,double &k0p,double &k1p)
{
    double r,x2,ca,cb,ct,ww,w0,xr,xr2;
    int k,kz;
    static double a[] = {
        0.125,
        7.03125e-2,
        7.32421875e-2,
        1.1215209960938e-1,
        2.2710800170898e-1,
        5.7250142097473e-1,
        1.7277275025845,
        6.0740420012735,
        2.4380529699556e1,
        1.1001714026925e2,
        5.5133589612202e2,
        3.0380905109224e3};
    static double b[] = {
        -0.375,
        -1.171875e-1,
        -1.025390625e-1,
        -1.4419555664063e-1,
        -2.7757644653320e-1,
        -6.7659258842468e-1,
        -1.9935317337513,
        -6.8839142681099,
        -2.7248827311269e1,
        -1.2159789187654e2,
        -6.0384407670507e2,
        -3.3022722944809e3};
    static double a1[] = {
        0.125,
        0.2109375,
        1.0986328125,
        1.1775970458984e1,
        2.1461706161499e2,
        5.9511522710323e3,
        2.3347645606175e5,
        1.2312234987631e7};

    if (x < 0.0) return 1;
    if (x == 0.0) {
        i0 = 1.0;
        i1 = 0.0;
        k0 = 1e308;
        k1 = 1e308;
        i0p = 0.0;
        i1p = 0.5;
        k0p = -1e308;
        k1p = -1e308;
        return 0;
    }
    x2 = x*x;
    if (x <= 18.0) {
        i0 = 1.0;
        r = 1.0;
        for (k=1;k<=50;k++) {
            r *= 0.25*x2/(k*k);
            i0 += r;
            if (fabs(r/i0) < eps) break;
        }
        i1 = 1.0;
        r = 1.0;
        for (k=1;k<=50;k++) {
            r *= 0.25*x2/(k*(k+1));
            i1 += r;
            if (fabs(r/i1) < eps) break;
        }
        i1 *= 0.5*x;
    }
    else {
        if (x >= 50.0) kz = 7;
        else if (x >= 35.0) kz = 9;
        else kz = 12;
        ca = exp(x)/sqrt(2.0*M_PI*x);
        i0 = 1.0;
        xr = 1.0/x;
        for (k=0;k<kz;k++) {
            i0 += a[k]*pow(xr,k+1);
        }
        i0 *= ca;
        i1 = 1.0;
        for (k=0;k<kz;k++) {
            i1 += b[k]*pow(xr,k+1);
        }
        i1 *= ca;
    }
    if (x <= 9.0) {
        ct = -(log(0.5*x)+el);
        k0 = 0.0;
        w0 = 0.0;
        r = 1.0;
        ww = 0.0;
        for (k=1;k<=50;k++) {
            w0 += 1.0/k;
            r *= 0.25*x2/(k*k);
            k0 += r*(w0+ct);
            if (fabs((k0-ww)/k0) < eps) break;
            ww = k0;
        }
        k0 += ct;
    }
    else {
        cb = 0.5/x;
        xr2 = 1.0/x2;
        k0 = 1.0;
        for (k=0;k<8;k++) {
            k0 += a1[k]*pow(xr2,k+1);
        }
        k0 *= cb/i0;
    }
    k1 = (1.0/x - i1*k0)/i0;
    i0p = i1;
    i1p = i0-i1/x;
    k0p = -k1;
    k1p = -k0-k1/x;
    return 0;
}

int bessik01b(double x,double &i0,double &i1,double &k0,double &k1,
    double &i0p,double &i1p,double &k0p,double &k1p)
{
    double t,t2,dtmp,dtmp1;

    if (x < 0.0) return 1;
    if (x == 0.0) {
        i0 = 1.0;
        i1 = 0.0;
        k0 = 1e308;
        k1 = 1e308;
        i0p = 0.0;
        i1p = 0.5;
        k0p = -1e308;
        k1p = -1e308;
        return 0;
    }
    if (x < 3.75) {
        t = x/3.75;
        t2 = t*t;
        i0 = (((((0.0045813*t2+0.0360768)*t2+0.2659732)*t2+
             1.2067492)*t2+3.0899424)*t2+3.5156229)*t2+1.0;
        i1 = x*(((((0.00032411*t2+0.00301532)*t2+0.02658733*t2+
             0.15084934)*t2+0.51498869)*t2+0.87890594)*t2+0.5);
    }
    else {
        t = 3.75/x;
        dtmp1 = exp(x)/sqrt(x);
        dtmp = (((((((0.00392377*t-0.01647633)*t+0.026355537)*t-0.02057706)*t+
          0.00916281)*t-0.00157565)*t+0.00225319)*t+0.01328592)*t+0.39894228;
        i0 = dtmp*dtmp1;
        dtmp = (((((((-0.00420059*t+0.01787654)*t-0.02895312)*t+0.02282967)*t-
          0.01031555)*t+0.00163801)*t-0.00362018)*t-0.03988024)*t+0.39894228;
        i1 = dtmp*dtmp1;
    }
    if (x < 2.0) {
        t = 0.5*x;
        t2 = t*t;       // already calculated above
        dtmp = (((((0.0000074*t2+0.0001075)*t2+0.00262698)*t2+0.0348859)*t2+
          0.23069756)*t2+0.4227842)*t2-0.57721566;
        k0 = dtmp - i0*log(t);
        dtmp = (((((-0.00004686*t2-0.00110404)*t2-0.01919402)*t2-
          0.18156897)*t2-0.67278578)*t2+0.15443144)*t2+1.0;
        k1 = dtmp/x + i1*log(t);
    }
    else {
        t = 2.0/x;
        dtmp1 = exp(-x)/sqrt(x);
        dtmp = (((((0.00053208*t-0.0025154)*t+0.00587872)*t-0.01062446)*t+
          0.02189568)*t-0.07832358)*t+1.25331414;
        k0 = dtmp*dtmp1;
        dtmp = (((((-0.00068245*t+0.00325614)*t-0.00780353)*t+0.01504268)*t-
          0.0365562)*t+0.23498619)*t+1.25331414;
        k1 = dtmp*dtmp1;
    }
    i0p = i1;
    i1p = i0 - i1/x;
    k0p = -k1;
    k1p = -k0 - k1/x;
    return 0;
}
int bessikna(int n,double x,int &nm,double *in,double *kn,
    double *inp,double *knp)
{
    double bi0,bi1,bk0,bk1,g,g0,g1,f,f0,f1,h,h0,h1,s0;
    int k,m,ecode;

    if ((x < 0.0) || (n < 0)) return 1;
    if (x < eps) {
        for (k=0;k<=n;k++) {
            in[k] = 0.0;
            kn[k] = 1e308;
            inp[k] = 0.0;
            knp[k] = -1e308;
        }
        in[0] = 1.0;
        inp[1] = 0.5;
        return 0;
    }
    nm = n;
    ecode = bessik01a(x,in[0],in[1],kn[0],kn[1],inp[0],inp[1],knp[0],knp[1]);
    if (n < 2) return 0;
    bi0 = in[0];
    bi1 = in[1];
    bk0 = kn[0];
    bk1 = kn[1];
    if ((x > 40.0) && (n < (int)(0.25*x))) {
        h0 = bi0;
        h1 = bi1;
        for (k=2;k<=n;k++) {
            h = -2.0*(k-1.0)*h1/x+h0;
            in[k] = h;
            h0 = h1;
            h1 = h;
        }
    }
    else {
        m = msta1(x,200);
        if (m < n) nm = m;
        else m = msta2(x,n,15);
        f0 = 0.0;
        f1 = 1.0e-100;
        for (k=m;k>=0;k--) {
            f = 2.0*(k+1.0)*f1/x+f0;
            if (x <= nm) in[k] = f;
            f0 = f1;
            f1 = f;
        }
        s0 = bi0/f;
        for (k=0;k<=m;k++) {
            in[k] *= s0;
        }
    }
    g0 = bk0;
    g1 = bk1;
    for (k=2;k<=nm;k++) {
        g = 2.0*(k-1.0)*g1/x+g0;
        kn[k] = g;
        g0 = g1;
        g1 = g;
    }
    for (k=2;k<=nm;k++) {
        inp[k] = in[k-1]-k*in[k]/x;
        knp[k] = -kn[k-1]-k*kn[k]/x;
    }
    return 0;
}
int bessiknb(int n,double x,int &nm,double *in,double *kn,
    double *inp,double *knp)
{
    double s0,bs,f,f0,f1,sk0,a0,bkl,vt,r,g,g0,g1;
    int k,kz,m,l;

    if ((x < 0.0) || (n < 0)) return 1;
    if (x < eps) {
        for (k=0;k<=n;k++) {
            in[k] = 0.0;
            kn[k] = 1e308;
            inp[k] = 0.0;
            knp[k] = -1e308;
        }
        in[0] = 1.0;
        inp[1] = 0.5;
        return 0;
    }
    nm = n;
    if (n == 0) nm = 1;
    m = msta1(x,200);
    if (m < nm) nm = m;
    else m = msta2(x,nm,15);
    bs = 0.0;
    sk0 = 0.0;
    f0 = 0.0;
    f1 = 1.0e-100;
    for (k=m;k>=0;k--) {
        f = 2.0*(k+1.0)*f1/x+f0;
        if (k <= nm) in[k] = f;
        if ((k != 0) && (k == 2*(int)(k/2))) {
            sk0 += 4.0*f/k;
        }
        bs += 2.0*f;
        f0 = f1;
        f1 = f;
    }
    s0 = exp(x)/(bs-f);
    for (k=0;k<=nm;k++) {
        in[k] *= s0;
    }
    if (x <= 8.0) {
        kn[0] = -(log(0.5*x)+el)*in[0]+s0*sk0;
        kn[1] = (1.0/x-in[1]*kn[0])/in[0];
    }
    else {
        a0 = sqrt(M_PI_2/x)*exp(-x);
        if (x >= 200.0) kz = 6;
        else if (x >= 80.0) kz = 8;
        else if (x >= 25.0) kz = 10;
        else kz = 16;
        for (l=0;l<2;l++) {
            bkl = 1.0;
            vt = 4.0*l;
            r = 1.0;
            for (k=1;k<=kz;k++) {
                r *= 0.125*(vt-pow(2.0*k-1.0,2))/(k*x);
                bkl += r;
            }
            kn[l] = a0*bkl;
        }
    }
    g0 = kn[0];
    g1 = kn[1];
    for (k=2;k<=nm;k++) {
        g = 2.0*(k-1.0)*g1/x+g0;
        kn[k] = g;
        g0 = g1;
        g1 = g;
    }
    inp[0] = in[1];
    knp[0] = -kn[1];
    for (k=1;k<=nm;k++) {
        inp[k] = in[k-1]-k*in[k]/x;
        knp[k] = -kn[k-1]-k*kn[k]/x;
    }
    return 0;
}

//  The following program computes the modified Bessel functions
//  Iv(x) and Kv(x) for arbitrary positive order. For negative
//  order use:
//
//          I-v(x) = Iv(x) + 2/pi sin(v pi) Kv(x)
//          K-v(x) = Kv(x)
//
int bessikv(double v,double x,double &vm,double *iv,double *kv,
    double *ivp,double *kvp)
{
    double x2,v0,piv,vt,a1,v0p,gap,r,bi0,ca,sum;
    double f,f1,f2,ct,cs,wa,gan,ww,w0,v0n;
    double r1,r2,bk0,bk1,bk2,a2,cb;
    int n,k,kz,m;

    if ((v < 0.0) || (x < 0.0)) return 1;
    x2 = x*x;
    n = (int)v;
    v0 = v-n;
    if (n == 0) n = 1;
    if (x == 0.0) {
        for (k=0;k<=n;k++) {
            iv[k] = 0.0;
            kv[k] = -1e308;
            ivp[k] = 0.0;
            kvp[k] = 1e308;
        }
        if (v0 == 0.0) {
            iv[0] = 1.0;
            ivp[1] = 0.5;
        }
        vm = v;
        return 0;
    }
    piv = M_PI*v0;
    vt = 4.0*v0*v0;
    if (v0 == 0.0) {
        a1 = 1.0;
    }
    else {
        v0p = 1.0+v0;
        gap = gamma(v0p);
        a1 = pow(0.5*x,v0)/gap;
    }
    if (x >= 50.0) kz = 8;
    else if (x >= 35.0) kz = 10;
    else kz = 14;
    if (x <= 18.0) {
        bi0 = 1.0;
        r = 1.0;
        for (k=1;k<=30;k++) {
            r *= 0.25*x2/(k*(k+v0));
            bi0 += r;
            if (fabs(r/bi0) < eps) break;
        }
        bi0 *= a1;
    }
    else {
        ca = exp(x)/sqrt(2.0*M_PI*x);
        sum = 1.0;
        r = 1.0;
        for (k=1;k<=kz;k++) {
            r *= -0.125*(vt-pow(2.0*k-1.0,2.0))/(k*x);
            sum += r;
        }
        bi0 = ca*sum;
    }
    m = msta1(x,200);
    if (m < n) n = m;
    else m = msta2(x,n,15);
    f2 = 0.0;
    f1 = 1.0e-100;
    for (k=m;k>=0;k--) {
        f = 2.0*(v0+k+1.0)*f1/x+f2;
        if (k <= n) iv[k] = f;
        f2 = f1;
        f1 = f;
    }
    cs = bi0/f;
    for (k=0;k<=n;k++) {
        iv[k] *= cs;
    }
    ivp[0] = v0*iv[0]/x+iv[1];
    for (k=1;k<=n;k++) {
        ivp[k] = -(k+v0)*iv[k]/x+iv[k-1];
    }
    ww = 0.0;
    if (x <= 9.0) {
        if (v0 == 0.0) {
            ct = -log(0.5*x)-el;
            cs = 0.0;
            w0 = 0.0;
            r = 1.0;
            for (k=1;k<=50;k++) {
                w0 += 1.0/k;
                r *= 0.25*x2/(k*k);
                cs += r*(w0+ct);
                wa = fabs(cs);
                if (fabs((wa-ww)/wa) < eps) break;
                ww = wa;
            }
            bk0 = ct+cs;
        }
        else {
            v0n = 1.0-v0;
            gan = gamma(v0n);
            a2 = 1.0/(gan*pow(0.5*x,v0));
            a1 = pow(0.5*x,v0)/gap;
            sum = a2-a1;
            r1 = 1.0;
            r2 = 1.0;
            for (k=1;k<=120;k++) {
                r1 *= 0.25*x2/(k*(k-v0));
                r2 *= 0.25*x2/(k*(k+v0));
                sum += a2*r1-a1*r2;
                wa = fabs(sum);
                if (fabs((wa-ww)/wa) < eps) break;
                ww = wa;
            }
            bk0 = M_PI_2*sum/sin(piv);
        }
    }
    else {
        cb = exp(-x)*sqrt(M_PI_2/x);
        sum = 1.0;
        r = 1.0;
        for (k=1;k<=kz;k++) {
            r *= 0.125*(vt-pow(2.0*k-1.0,2.0))/(k*x);
            sum += r;
        }
        bk0 = cb*sum;
    }
    bk1 = (1.0/x-iv[1]*bk0)/iv[0];
    kv[0] = bk0;
    kv[1] = bk1;
    for (k=2;k<=n;k++) {
        bk2 = 2.0*(v0+k-1.0)*bk1/x+bk0;
        kv[k] = bk2;
        bk0 = bk1;
        bk1 = bk2;
    }
    kvp[0] = v0*kv[0]/x-kv[1];
    for (k=1;k<=n;k++) {
        kvp[k] = -(k+v0)*kv[k]/x-kv[k-1];
    }
    vm = n+v0;
    return 0;
}
//  bessjy.cpp -- computation of Bessel functions Jn, Yn and their
//      derivatives. Algorithms and coefficient values from
//      "Computation of Special Functions", Zhang and Jin, John
//      Wiley and Sons, 1996.
//
//  (C) 2003, C. Bond. All rights reserved.
//
// Note that 'math.h' provides (or should provide) values for:
//      pi      M_PI
//      2/pi    M_2_PI
//      pi/4    M_PI_4
//      pi/2    M_PI_2
//
//MH#include <math.h>
//MH#include "bessel.h"

double gamma(double x);
//
//  INPUT:
//      double x    -- argument of Bessel function
//
//  OUTPUT (via address pointers):
//      double j0   -- Bessel function of 1st kind, 0th order
//      double j1   -- Bessel function of 1st kind, 1st order
//      double y0   -- Bessel function of 2nd kind, 0th order
//      double y1   -- Bessel function of 2nd kind, 1st order
//      double j0p  -- derivative of Bessel function of 1st kind, 0th order
//      double j1p  -- derivative of Bessel function of 1st kind, 1st order
//      double y0p  -- derivative of Bessel function of 2nd kind, 0th order
//      double y1p  -- derivative of Bessel function of 2nd kind, 1st order
//
//  RETURN:
//      int error code: 0 = OK, 1 = error
//
//  This algorithm computes the above functions using series expansions.
//
int bessjy01a(double x,double &j0,double &j1,double &y0,double &y1,
    double &j0p,double &j1p,double &y0p,double &y1p)
{
    double x2,r,ec,w0,w1,r0,r1,cs0,cs1;
    double cu,p0,q0,p1,q1,t1,t2;
    int k,kz;
    static double a[] = {
        -7.03125e-2,
         0.112152099609375,
        -0.5725014209747314,
         6.074042001273483,
        -1.100171402692467e2,
         3.038090510922384e3,
        -1.188384262567832e5,
         6.252951493434797e6,
        -4.259392165047669e8,
         3.646840080706556e10,
        -3.833534661393944e12,
         4.854014686852901e14,
        -7.286857349377656e16,
         1.279721941975975e19};
    static double b[] = {
         7.32421875e-2,
        -0.2271080017089844,
         1.727727502584457,
        -2.438052969955606e1,
         5.513358961220206e2,
        -1.825775547429318e4,
         8.328593040162893e5,
        -5.006958953198893e7,
         3.836255180230433e9,
        -3.649010818849833e11,
         4.218971570284096e13,
        -5.827244631566907e15,
         9.476288099260110e17,
        -1.792162323051699e20};
    static double a1[] = {
         0.1171875,
        -0.1441955566406250,
         0.6765925884246826,
        -6.883914268109947,
         1.215978918765359e2,
        -3.302272294480852e3,
         1.276412726461746e5,
        -6.656367718817688e6,
         4.502786003050393e8,
        -3.833857520742790e10,
         4.011838599133198e12,
        -5.060568503314727e14,
         7.572616461117958e16,
        -1.326257285320556e19};
    static double b1[] = {
        -0.1025390625,
         0.2775764465332031,
        -1.993531733751297,
         2.724882731126854e1,
        -6.038440767050702e2,
         1.971837591223663e4,
        -8.902978767070678e5,
         5.310411010968522e7,
        -4.043620325107754e9,
         3.827011346598605e11,
        -4.406481417852278e13,
         6.065091351222699e15,
        -9.833883876590679e17,
         1.855045211579828e20};

    if (x < 0.0) return 1;
    if (x == 0.0) {
        j0 = 1.0;
        j1 = 0.0;
        y0 = -1e308;
        y1 = -1e308;
        j0p = 0.0;
        j1p = 0.5;
        y0p = 1e308;
        y1p = 1e308;
        return 0;
    }
    x2 = x*x;
    if (x <= 12.0) {
        j0 = 1.0;
        r = 1.0;
        for (k=1;k<=30;k++) {
            r *= -0.25*x2/(k*k);
            j0 += r;
            if (fabs(r) < fabs(j0)*1e-15) break;
        }
        j1 = 1.0;
        r = 1.0;
        for (k=1;k<=30;k++) {
            r *= -0.25*x2/(k*(k+1));
            j1 += r;
            if (fabs(r) < fabs(j1)*1e-15) break;
        }
        j1 *= 0.5*x;
        ec = log(0.5*x)+el;
        cs0 = 0.0;
        w0 = 0.0;
        r0 = 1.0;
        for (k=1;k<=30;k++) {
            w0 += 1.0/k;
            r0 *= -0.25*x2/(k*k);
            r = r0 * w0;
            cs0 += r;
            if (fabs(r) < fabs(cs0)*1e-15) break;
        }
        y0 = M_2_PI*(ec*j0-cs0);
        cs1 = 1.0;
        w1 = 0.0;
        r1 = 1.0;
        for (k=1;k<=30;k++) {
            w1 += 1.0/k;
            r1 *= -0.25*x2/(k*(k+1));
            r = r1*(2.0*w1+1.0/(k+1));
            cs1 += r;
            if (fabs(r) < fabs(cs1)*1e-15) break;
        }
        y1 = M_2_PI * (ec*j1-1.0/x-0.25*x*cs1);
    }
    else {
        if (x >= 50.0) kz = 8;          // Can be changed to 10
        else if (x >= 35.0) kz = 10;    //  "       "        12
        else kz = 12;                   //  "       "        14
        t1 = x-M_PI_4;
        p0 = 1.0;
        q0 = -0.125/x;
        for (k=0;k<kz;k++) {
            p0 += a[k]*pow(x,-2*k-2);
            q0 += b[k]*pow(x,-2*k-3);
        }
        cu = sqrt(M_2_PI/x);
        j0 = cu*(p0*cos(t1)-q0*sin(t1));
        y0 = cu*(p0*sin(t1)+q0*cos(t1));
        t2 = x-0.75*M_PI;
        p1 = 1.0;
        q1 = 0.375/x;
        for (k=0;k<kz;k++) {
            p1 += a1[k]*pow(x,-2*k-2);
            q1 += b1[k]*pow(x,-2*k-3);
        }
        j1 = cu*(p1*cos(t2)-q1*sin(t2));
        y1 = cu*(p1*sin(t2)+q1*cos(t2));
    }
    j0p = -j1;
    j1p = j0-j1/x;
    y0p = -y1;
    y1p = y0-y1/x;
    return 0;
}
//
//  INPUT:
//      double x    -- argument of Bessel function
//
//  OUTPUT:
//      double j0   -- Bessel function of 1st kind, 0th order
//      double j1   -- Bessel function of 1st kind, 1st order
//      double y0   -- Bessel function of 2nd kind, 0th order
//      double y1   -- Bessel function of 2nd kind, 1st order
//      double j0p  -- derivative of Bessel function of 1st kind, 0th order
//      double j1p  -- derivative of Bessel function of 1st kind, 1st order
//      double y0p  -- derivative of Bessel function of 2nd kind, 0th order
//      double y1p  -- derivative of Bessel function of 2nd kind, 1st order
//
//  RETURN:
//      int error code: 0 = OK, 1 = error
//
//  This algorithm computes the functions using polynomial approximations.
//
int bessjy01b(double x,double &j0,double &j1,double &y0,double &y1,
    double &j0p,double &j1p,double &y0p,double &y1p)
{
    double t,t2,dtmp,a0,p0,q0,p1,q1,ta0,ta1;
    if (x < 0.0) return 1;
    if (x == 0.0) {
        j0 = 1.0;
        j1 = 0.0;
        y0 = -1e308;
        y1 = -1e308;
        j0p = 0.0;
        j1p = 0.5;
        y0p = 1e308;
        y1p = 1e308;
        return 0;
    }
    if(x <= 4.0) {
        t = x/4.0;
        t2 = t*t;
        j0 = ((((((-0.5014415e-3*t2+0.76771853e-2)*t2-0.0709253492)*t2+
            0.4443584263)*t2-1.7777560599)*t2+3.9999973021)*t2
            -3.9999998721)*t2+1.0;
        j1 = t*(((((((-0.1289769e-3*t2+0.22069155e-2)*t2-0.0236616773)*t2+
            0.1777582922)*t2-0.8888839649)*t2+2.6666660544)*t2-
            3.999999971)*t2+1.9999999998);
        dtmp = (((((((-0.567433e-4*t2+0.859977e-3)*t2-0.94855882e-2)*t2+
            0.0772975809)*t2-0.4261737419)*t2+1.4216421221)*t2-
            2.3498519931)*t2+1.0766115157)*t2+0.3674669052;
        y0 = M_2_PI*log(0.5*x)*j0+dtmp;
        dtmp = (((((((0.6535773e-3*t2-0.0108175626)*t2+0.107657607)*t2-
            0.7268945577)*t2+3.1261399273)*t2-7.3980241381)*t2+
            6.8529236342)*t2+0.3932562018)*t2-0.6366197726;
        y1 = M_2_PI*log(0.5*x)*j1+dtmp/x;
    }
    else {
        t = 4.0/x;
        t2 = t*t;
        a0 = sqrt(M_2_PI/x);
        p0 = ((((-0.9285e-5*t2+0.43506e-4)*t2-0.122226e-3)*t2+
             0.434725e-3)*t2-0.4394275e-2)*t2+0.999999997;
        q0 = t*(((((0.8099e-5*t2-0.35614e-4)*t2+0.85844e-4)*t2-
            0.218024e-3)*t2+0.1144106e-2)*t2-0.031249995);
        ta0 = x-M_PI_4;
        j0 = a0*(p0*cos(ta0)-q0*sin(ta0));
        y0 = a0*(p0*sin(ta0)+q0*cos(ta0));
        p1 = ((((0.10632e-4*t2-0.50363e-4)*t2+0.145575e-3)*t2
            -0.559487e-3)*t2+0.7323931e-2)*t2+1.000000004;
        q1 = t*(((((-0.9173e-5*t2+0.40658e-4)*t2-0.99941e-4)*t2
            +0.266891e-3)*t2-0.1601836e-2)*t2+0.093749994);
        ta1 = x-0.75*M_PI;
        j1 = a0*(p1*cos(ta1)-q1*sin(ta1));
        y1 = a0*(p1*sin(ta1)+q1*cos(ta1));
    }
    j0p = -j1;
    j1p = j0-j1/x;
    y0p = -y1;
    y1p = y0-y1/x;
    return 0;
}
int msta1(double x,int mp)
{
    double a0,f0,f1,f;
    int i,n0,n1,nn;

    a0 = fabs(x);
    n0 = (int)(1.1*a0)+1;
    f0 = 0.5*log10(6.28*n0)-n0*log10(1.36*a0/n0)-mp;
    n1 = n0+5;
    f1 = 0.5*log10(6.28*n1)-n1*log10(1.36*a0/n1)-mp;
    for (i=0;i<20;i++) {
        nn = n1-(n1-n0)/(1.0-f0/f1);
        f = 0.5*log10(6.28*nn)-nn*log10(1.36*a0/nn)-mp;
        int aux=int(std::abs(float(nn-n1)));
        //if (int(std::abs(nn-n1)) < 1) break;
        if (aux < 1) break;
        n0 = n1;
        f0 = f1;
        n1 = nn;
        f1 = f;
    }
    return nn;
}
int msta2(double x,int n,int mp)
{
    double a0,ejn,hmp,f0,f1,f,obj;
    int i,n0,n1,nn;

    a0 = fabs(x);
    hmp = 0.5*mp;
    ejn = 0.5*log10(6.28*n)-n*log10(1.36*a0/n);
    if (ejn <= hmp) {
        obj = mp;
        n0 = (int)(1.1*a0);
        if (n0 < 1) n0 = 1;
    }
    else {
        obj = hmp+ejn;
        n0 = n;
    }
    f0 = 0.5*log10(6.28*n0)-n0*log10(1.36*a0/n0)-obj;
    n1 = n0+5;
    f1 = 0.5*log10(6.28*n1)-n1*log10(1.36*a0/n1)-obj;
    for (i=0;i<20;i++) {
        nn = n1-(n1-n0)/(1.0-f0/f1);
        f = 0.5*log10(6.28*nn)-nn*log10(1.36*a0/nn)-obj;
        int aux=int(std::abs(float(nn-n1)));
        if (aux < 1) break;
        n0 = n1;
        f0 = f1;
        n1 = nn;
        f1 = f;
    }
    return nn+10;
}
//
//  INPUT:
//  double x    -- argument of Bessel function of 1st and 2nd kind.
//  int n       -- order
//
//  OUPUT:
//
//  int nm      -- highest order actually computed (nm <= n)
//  double jn[] -- Bessel function of 1st kind, orders from 0 to nm
//  double yn[] -- Bessel function of 2nd kind, orders from 0 to nm
//  double j'n[]-- derivative of Bessel function of 1st kind,
//                      orders from 0 to nm
//  double y'n[]-- derivative of Bessel function of 2nd kind,
//                      orders from 0 to nm
//
//  Computes Bessel functions of all order up to 'n' using recurrence
//  relations. If 'nm' < 'n' only 'nm' orders are returned.
//
int bessjyna(int n,double x,int &nm,double *jn,double *yn,
    double *jnp,double *ynp)
{
    double bj0,bj1,f,f0,f1,f2,cs;
    int i,k,m,ecode;

    nm = n;
    if ((x < 0.0) || (n < 0)) return 1;
    if (x < 1e-15) {
        for (i=0;i<=n;i++) {
            jn[i] = 0.0;
            yn[i] = -1e308;
            jnp[i] = 0.0;
            ynp[i] = 1e308;
        }
        jn[0] = 1.0;
        jnp[1] = 0.5;
        return 0;
    }
    ecode = bessjy01a(x,jn[0],jn[1],yn[0],yn[1],jnp[0],jnp[1],ynp[0],ynp[1]);
    if (n < 2) return 0;
    bj0 = jn[0];
    bj1 = jn[1];
    if (n < (int)0.9*x) {
        for (k=2;k<=n;k++) {
            jn[k] = 2.0*(k-1.0)*bj1/x-bj0;
            bj0 = bj1;
            bj1 = jn[k];
        }
    }
    else {
        m = msta1(x,200);
        if (m < n) nm = m;
        else m = msta2(x,n,15);
        f2 = 0.0;
        f1 = 1.0e-100;
        for (k=m;k>=0;k--) {
            f = 2.0*(k+1.0)/x*f1-f2;
            if (k <= nm) jn[k] = f;
            f2 = f1;
            f1 = f;
        }
        if (fabs(bj0) > fabs(bj1)) cs = bj0/f;
        else cs = bj1/f2;
        for (k=0;k<=nm;k++) {
            jn[k] *= cs;
        }
    }
    for (k=2;k<=nm;k++) {
        jnp[k] = jn[k-1]-k*jn[k]/x;
    }
    f0 = yn[0];
    f1 = yn[1];
    for (k=2;k<=nm;k++) {
        f = 2.0*(k-1.0)*f1/x-f0;
        yn[k] = f;
        f0 = f1;
        f1 = f;
    }
    for (k=2;k<=nm;k++) {
        ynp[k] = yn[k-1]-k*yn[k]/x;
    }
    return 0;
}
//
//  Same input and output conventions as above. Different recurrence
//  relations used for 'x' < 300.
//
int bessjynb(int n,double x,int &nm,double *jn,double *yn,
    double *jnp,double *ynp)
{
    double t1,t2,f,f1,f2,bj0,bj1,bjk,by0,by1,cu,s0,su,sv;
    double ec,bs,byk,p0,p1,q0,q1;
    static double a[] = {
        -0.7031250000000000e-1,
         0.1121520996093750,
        -0.5725014209747314,
         6.074042001273483};
    static double b[] = {
         0.7324218750000000e-1,
        -0.2271080017089844,
         1.727727502584457,
        -2.438052969955606e1};
    static double a1[] = {
         0.1171875,
        -0.1441955566406250,
         0.6765925884246826,
        -6.883914268109947};
    static double b1[] = {
       -0.1025390625,
        0.2775764465332031,
       -1.993531733751297,
        2.724882731126854e1};

    int i,k,m;
    nm = n;
    if ((x < 0.0) || (n < 0)) return 1;
    if (x < 1e-15) {
        for (i=0;i<=n;i++) {
            jn[i] = 0.0;
            yn[i] = -1e308;
            jnp[i] = 0.0;
            ynp[i] = 1e308;
        }
        jn[0] = 1.0;
        jnp[1] = 0.5;
        return 0;
    }
    if (x <= 300.0 || n > (int)(0.9*x)) {
        if (n == 0) nm = 1;
        m = msta1(x,200);
        if (m < nm) nm = m;
        else m = msta2(x,nm,15);
        bs = 0.0;
        su = 0.0;
        sv = 0.0;
        f2 = 0.0;
        f1 = 1.0e-100;
        for (k = m;k>=0;k--) {
            f = 2.0*(k+1.0)/x*f1 - f2;
            if (k <= nm) jn[k] = f;
            if ((k == 2*(int)(k/2)) && (k != 0)) {
                bs += 2.0*f;
//                su += pow(-1,k>>1)*f/(double)k;
                su += (-1)*((k & 2)-1)*f/(double)k;
            }
            else if (k > 1) {
//                sv += pow(-1,k>>1)*k*f/(k*k-1.0);
                sv += (-1)*((k & 2)-1)*(double)k*f/(k*k-1.0);
            }
            f2 = f1;
            f1 = f;
        }
        s0 = bs+f;
        for (k=0;k<=nm;k++) {
            jn[k] /= s0;
        }
        ec = log(0.5*x) +0.5772156649015329;
        by0 = M_2_PI*(ec*jn[0]-4.0*su/s0);
        yn[0] = by0;
        by1 = M_2_PI*((ec-1.0)*jn[1]-jn[0]/x-4.0*sv/s0);
        yn[1] = by1;
    }
    else {
        t1 = x-M_PI_4;
        p0 = 1.0;
        q0 = -0.125/x;
        for (k=0;k<4;k++) {
            p0 += a[k]*pow(x,-2*k-2);
            q0 += b[k]*pow(x,-2*k-3);
        }
        cu = sqrt(M_2_PI/x);
        bj0 = cu*(p0*cos(t1)-q0*sin(t1));
        by0 = cu*(p0*sin(t1)+q0*cos(t1));
        jn[0] = bj0;
        yn[0] = by0;
        t2 = x-0.75*M_PI;
        p1 = 1.0;
        q1 = 0.375/x;
        for (k=0;k<4;k++) {
            p1 += a1[k]*pow(x,-2*k-2);
            q1 += b1[k]*pow(x,-2*k-3);
        }
        bj1 = cu*(p1*cos(t2)-q1*sin(t2));
        by1 = cu*(p1*sin(t2)+q1*cos(t2));
        jn[1] = bj1;
        yn[1] = by1;
        for (k=2;k<=nm;k++) {
            bjk = 2.0*(k-1.0)*bj1/x-bj0;
            jn[k] = bjk;
            bj0 = bj1;
            bj1 = bjk;
        }
    }
    jnp[0] = -jn[1];
    for (k=1;k<=nm;k++) {
        jnp[k] = jn[k-1]-k*jn[k]/x;
    }
    for (k=2;k<=nm;k++) {
        byk = 2.0*(k-1.0)*by1/x-by0;
        yn[k] = byk;
        by0 = by1;
        by1 = byk;
    }
    ynp[0] = -yn[1];
    for (k=1;k<=nm;k++) {
        ynp[k] = yn[k-1]-k*yn[k]/x;
    }
    return 0;

}

//  The following routine computes Bessel Jv(x) and Yv(x) for
//  arbitrary positive order (v). For negative order, use:
//
//      J-v(x) = Jv(x)cos(v pi) - Yv(x)sin(v pi)
//      Y-v(x) = Jv(x)sin(v pi) + Yv(x)cos(v pi)
//
int bessjyv(double v,double x,double &vm,double *jv,double *yv,
    double *djv,double *dyv)
{
    double v0,vl,vg,vv,a,a0,r,x2,bjv0,bjv1,bjvl,f,f0,f1,f2;
    double r0,r1,ck,cs,cs0,cs1,sk,qx,px,byv0,byv1,rp,xk,rq;
    double b,ec,w0,w1,bju0,bju1,pv0,pv1,byvk;
    int j,k,l,m,n,kz;

    x2 = x*x;
    n = (int)v;
    v0 = v-n;
    if ((x < 0.0) || (v < 0.0)) return 1;
    if (x < 1e-15) {
        for (k=0;k<=n;k++) {
            jv[k] = 0.0;
            yv[k] = -1e308;
            djv[k] = 0.0;
            dyv[k] = 1e308;
            if (v0 == 0.0) {
                jv[0] = 1.0;
                djv[1] = 0.5;
            }
            else djv[0] = 1e308;
        }
        vm = v;
        return 0;
    }
    if (x <= 12.0) {
        for (l=0;l<2;l++) {
            vl = v0 + l;
            bjvl = 1.0;
            r = 1.0;
            for (k=1;k<=40;k++) {
                r *= -0.25*x2/(k*(k+vl));
                bjvl += r;
                if (fabs(r) < fabs(bjvl)*1e-15) break;
            }
            vg = 1.0 + vl;
            a = pow(0.5*x,vl)/gamma(vg);
            if (l == 0) bjv0 = bjvl*a;
            else bjv1 = bjvl*a;
        }
    }
    else {
        if (x >= 50.0) kz = 8;
        else if (x >= 35.0) kz = 10;
        else kz = 11;
        for (j=0;j<2;j++) {
            vv = 4.0*(j+v0)*(j+v0);
            px = 1.0;
            rp = 1.0;
            for (k=1;k<=kz;k++) {
                rp *= (-0.78125e-2)*(vv-pow(4.0*k-3.0,2.0))*
                    (vv-pow(4.0*k-1.0,2.0))/(k*(2.0*k-1.0)*x2);
                px += rp;
            }
            qx = 1.0;
            rq = 1.0;
            for (k=1;k<=kz;k++) {
                rq *= (-0.78125e-2)*(vv-pow(4.0*k-1.0,2.0))*
                    (vv-pow(4.0*k+1.0,2.0))/(k*(2.0*k+1.0)*x2);
                qx += rq;
            }
            qx *= 0.125*(vv-1.0)/x;
            xk = x-(0.5*(j+v0)+0.25)*M_PI;
            a0 = sqrt(M_2_PI/x);
            ck = cos(xk);
            sk = sin(xk);

            if (j == 0) {
                bjv0 = a0*(px*ck-qx*sk);
                byv0 = a0*(px*sk+qx*ck);
            }
            else {
                bjv1 = a0*(px*ck-qx*sk);
                byv1 = a0*(px*sk+qx*ck);
            }
        }
    }
    jv[0] = bjv0;
    jv[1] = bjv1;
    djv[0] = v0*jv[0]/x-jv[1];
    djv[1] = -(1.0+v0)*jv[1]/x+jv[0];
    if ((n >= 2) && (n <= (int)(0.9*x))) {
        f0 = bjv0;
        f1 = bjv1;
        for (k=2;k<=n;k++) {
            f = 2.0*(k+v0-1.0)*f1/x-f0;
            jv[k] = f;
            f0 = f1;
            f1 = f;
        }
    }
    else if (n >= 2) {
        m = msta1(x,200);
        if (m < n) n = m;
        else m = msta2(x,n,15);
        f2 = 0.0;
        f1 = 1.0e-100;
        for (k=m;k>=0;k--) {
            f = 2.0*(v0+k+1.0)/x*f1-f2;
            if (k <= n) jv[k] = f;
            f2 = f1;
            f1 = f;
        }
        if (fabs(bjv0) > fabs(bjv1)) cs = bjv0/f;
        else cs = bjv1/f2;
        for (k=0;k<=n;k++) {
            jv[k] *= cs;
        }
    }
    for (k=2;k<=n;k++) {
        djv[k] = -(k+v0)*jv[k]/x+jv[k-1];
    }
    if (x <= 12.0) {
        if (v0 != 0.0) {
            for (l=0;l<2;l++) {
                vl = v0 +l;
                bjvl = 1.0;
                r = 1.0;
                for (k=1;k<=40;k++) {
                    r *= -0.25*x2/(k*(k-vl));
                    bjvl += r;
                    if (fabs(r) < fabs(bjvl)*1e-15) break;
                }
                vg = 1.0-vl;
                b = pow(2.0/x,vl)/gamma(vg);
                if (l == 0) bju0 = bjvl*b;
                else bju1 = bjvl*b;
            }
            pv0 = M_PI*v0;
            pv1 = M_PI*(1.0+v0);
            byv0 = (bjv0*cos(pv0)-bju0)/sin(pv0);
            byv1 = (bjv1*cos(pv1)-bju1)/sin(pv1);
        }
        else {
            ec = log(0.5*x)+el;
            cs0 = 0.0;
            w0 = 0.0;
            r0 = 1.0;
            for (k=1;k<=30;k++) {
                w0 += 1.0/k;
                r0 *= -0.25*x2/(k*k);
                cs0 += r0*w0;
            }
            byv0 = M_2_PI*(ec*bjv0-cs0);
            cs1 = 1.0;
            w1 = 0.0;
            r1 = 1.0;
            for (k=1;k<=30;k++) {
                w1 += 1.0/k;
                r1 *= -0.25*x2/(k*(k+1));
                cs1 += r1*(2.0*w1+1.0/(k+1.0));
            }
            byv1 = M_2_PI*(ec*bjv1-1.0/x-0.25*x*cs1);
        }
    }
    yv[0] = byv0;
    yv[1] = byv1;
    for (k=2;k<=n;k++) {
        byvk = 2.0*(v0+k-1.0)*byv1/x-byv0;
        yv[k] = byvk;
        byv0 = byv1;
        byv1 = byvk;
    }
    dyv[0] = v0*yv[0]/x-yv[1];
    for (k=1;k<=n;k++) {
        dyv[k] = -(k+v0)*yv[k]/x+yv[k-1];
    }
    vm = n + v0;
    return 0;
}
//  cbessik.cpp -- complex modified Bessel functions.
//  Algorithms and coefficient values from "Computation of Special
//  Functions", Zhang and Jin, John Wiley and Sons, 1996.
//
//  (C) 2003, C. Bond. All rights reserved.
//
//MH#include <complex.h>
//MH#include "bessel.h"

static complex<double> cii(0.0,1.0);
static complex<double> czero(0.0,0.0);
static complex<double> cone(1.0,0.0);

double gamma(double x);

int cbessik01(complex<double>z,complex<double>&ci0,complex<double>&ci1,
    complex<double>&ck0,complex<double>&ck1,complex<double>&ci0p,
    complex<double>&ci1p,complex<double>&ck0p,complex<double>&ck1p)
{
    complex<double> z1,z2,zr,zr2,cr,ca,cb,cs,ct,cw;
    double a0,w0;
    int k,kz;
    static double a[] = {
        0.125,
        7.03125e-2,
        7.32421875e-2,
        1.1215209960938e-1,
        2.2710800170898e-1,
        5.7250142097473e-1,
        1.7277275025845,
        6.0740420012735,
        2.4380529699556e1,
        1.1001714026925e2,
        5.5133589612202e2,
        3.0380905109224e3};
    static double b[] = {
        -0.375,
        -1.171875e-1,
        -1.025390625e-1,
        -1.4419555664063e-1,
        -2.7757644653320e-1,
        -6.7659258842468e-1,
        -1.9935317337513,
        -6.8839142681099,
        -2.7248827311269e1,
        -1.2159789187654e2,
        -6.0384407670507e2,
        -3.3022722944809e3};
    static double a1[] = {
        0.125,
        0.2109375,
        1.0986328125,
        1.1775970458984e1,
        2.1461706161499e2,
        5.9511522710323e3,
        2.3347645606175e5,
        1.2312234987631e7,
        8.401390346421e08,
        7.2031420482627e10};

    a0 = abs(z);
    z2 = z*z;
    z1 = z;
    if (a0 == 0.0) {
        ci0 = cone;
        ci1 = czero;
        ck0 = complex<double> (1e308,0);
        ck1 = complex<double> (1e308,0);
        ci0p = czero;
        ci1p = complex<double>(0.5,0.0);
        ck0p = complex<double>(-1e308,0);
        ck1p = complex<double>(-1e308,0);
        return 0;
    }
    if (real(z) < 0.0) z1 = -z;
    if (a0 <= 18.0) {
        ci0 = cone;
        cr = cone;
        for (k=1;k<=50;k++) {
            cr *= 0.25*z2/(double)(k*k);
            ci0 += cr;
            if (abs(cr/ci0) < eps) break;
        }
        ci1 = cone;
        cr = cone;
        for (k=1;k<=50;k++) {
            cr *= 0.25*z2/(double)(k*(k+1.0));
            ci1 += cr;
            if (abs(cr/ci1) < eps) break;
        }
        ci1 *= 0.5*z1;
    }
    else {
        if (a0 >= 50.0) kz = 7;
        else if (a0 >= 35.0) kz = 9;
        else kz = 12;
        ca = exp(z1)/sqrt(2.0*M_PI*z1);
        ci0 = cone;
        zr = 1.0/z1;
        for (k=0;k<kz;k++) {
            ci0 += a[k]*pow(zr,k+1.0);
        }
        ci0 *= ca;
        ci1 = cone;
        for (k=0;k<kz;k++) {
            ci1 += b[k]*pow(zr,k+1.0);
        }
        ci1 *= ca;
    }
    if (a0 <= 9.0) {
        cs = czero;
        ct = -log(0.5*z1)-el;
        w0 = 0.0;
        cr = cone;
        for (k=1;k<=50;k++) {
            w0 += 1.0/k;
            cr *= 0.25*z2/(double)(k*k);
            cs += cr*(w0+ct);
            if (abs((cs-cw)/cs) < eps) break;
            cw = cs;
        }
        ck0 = ct+cs;
    }
    else {
        cb = 0.5/z1;
        zr2 = 1.0/z2;
        ck0 = cone;
        for (k=0;k<10;k++) {
            ck0 += a1[k]*pow(zr2,k+1.0);
        }
        ck0 *= cb/ci0;
    }
    ck1 = (1.0/z1 - ci1*ck0)/ci0;
    if (real(z) < 0.0) {
        if (imag(z) < 0.0) {
            ck0 += cii*M_PI*ci0;
            ck1 = -ck1+cii*M_PI*ci1;
        }
        else if (imag(z) > 0.0) {
            ck0 -= cii*M_PI*ci0;
            ck1 = -ck1-cii*M_PI*ci1;
        }
        ci1 = -ci1;
    }
    ci0p = ci1;
    ci1p = ci0-1.0*ci1/z;
    ck0p = -ck1;
    ck1p = -ck0-1.0*ck1/z;
    return 0;
}
int cbessikna(int n,complex<double> z,int &nm,complex<double> *ci,
    complex<double> *ck,complex<double> *cip,complex<double> *ckp)
{
    complex<double> ci0,ci1,ck0,ck1,ckk,cf,cf1,cf2,cs;
    double a0;
    int k,m,ecode;
    a0 = abs(z);
    nm = n;
    if (a0 < 1.0e-100) {
        for (k=0;k<=n;k++) {
            ci[k] = czero;
            ck[k] = complex<double>(-1e308,0);
            cip[k] = czero;
            ckp[k] = complex<double>(1e308,0);
        }
        ci[0] = complex<double>(1e308,0);
        cip[1] = complex<double>(0.5,0.0);
        return 0;
    }
    ecode = cbessik01(z,ci[0],ci[1],ck[0],ck[1],cip[0],cip[1],ckp[0],ckp[1]);
    if (n < 2) return 0;
    ci0 = ci[0];
    ci1 = ci[1];
    ck0 = ck[0];
    ck1 = ck[1];
    m = msta1(a0,200);
    if (m < n) nm = m;
    else m = msta2(a0,n,15);
    cf2 = czero;
    cf1 = complex<double>(1.0e-100,0.0);
    for (k=m;k>=0;k--) {
        cf = 2.0*(k+1.0)*cf1/z+cf2;
        if (k <= nm) ci[k] = cf;
        cf2 = cf1;
        cf1 = cf;
    }
    cs = ci0/cf;
    for (k=0;k<=nm;k++) {
        ci[k] *= cs;
    }
    for (k=2;k<=nm;k++) {
        if (abs(ci[k-1]) > abs(ci[k-2])) {
            ckk = (1.0/z-ci[k]*ck[k-1])/ci[k-1];
        }
        else {
            ckk = (ci[k]*ck[k-2]+2.0*(k-1.0)/(z*z))/ci[k-2];
        }
        ck[k] = ckk;
    }
    for (k=2;k<=nm;k++) {
        cip[k] = ci[k-1]-(double)k*ci[k]/z;
        ckp[k] = -ck[k-1]-(double)k*ck[k]/z;
    }
    return 0;
}
int cbessiknb(int n,complex<double> z,int &nm,complex<double> *ci,
    complex<double> *ck,complex<double> *cip,complex<double> *ckp)
{
    complex<double> z1,cbs,csk0,cf,cf0,cf1,ca0,cbkl;
    complex<double> cg,cg0,cg1,cs0,cs,cr;
    double a0,vt,fac;
    int k,kz,l,m;

    a0 = abs(z);
    nm = n;
    if (a0 < 1.0e-100) {
        for (k=0;k<=n;k++) {
            ci[k] = czero;
            ck[k] = complex<double>(1e308,0);
            cip[k] = czero;
            ckp[k] = complex<double>(-1e308,0);
        }
        ci[0] = complex<double>(1.0,0.0);
        cip[1] = complex<double>(0.5,0.0);
        return 0;
    }
    z1 = z;
    if (real(z) < 0.0) z1 = -z;
    if (n == 0) nm = 1;
    m = msta1(a0,200);
    if (m < nm) nm = m;
    else m = msta2(a0,nm,15);
    cbs = czero;
    csk0 = czero;
    cf0 = czero;
    cf1 = complex<double>(1.0e-100,0.0);
    for (k=m;k>=0;k--) {
        cf = 2.0*(k+1.0)*cf1/z1+cf0;
        if (k <=nm) ci[k] = cf;
        if ((k != 0) && (k == 2*(k>>1)))  csk0 += 4.0*cf/(double)k;
        cbs += 2.0*cf;
        cf0 = cf1;
        cf1 = cf;
    }
    cs0 = exp(z1)/(cbs-cf);
    for (k=0;k<=nm;k++) {
        ci[k] *= cs0;
    }
    if (a0 <= 9.0) {
        ck[0] = -(log(0.5*z1)+el)*ci[0]+cs0*csk0;
        ck[1] = (1.0/z1-ci[1]*ck[0])/ci[0];
    }
    else {
        ca0 = sqrt(M_PI_2/z1)*exp(-z1);
        if (a0 >= 200.0) kz = 6;
        else if (a0 >= 80.0) kz = 8;
        else if (a0 >= 25.0) kz = 10;
        else kz = 16;
        for (l=0;l<2;l++) {
            cbkl = cone;
            vt = 4.0*l;
            cr = cone;
            for (k=1;k<=kz;k++) {
                cr *= 0.125*(vt-pow(2.0*k-1.0,2.0))/((double)k*z1);
                cbkl += cr;
            }
            ck[l] = ca0*cbkl;
        }
    }
    cg0 = ck[0];
    cg1 = ck[1];
    for (k=2;k<=nm;k++) {
        cg = 2.0*(k-1.0)*cg1/z1+cg0;
        ck[k] = cg;
        cg0 = cg1;
        cg1 = cg;
    }
    if (real(z) < 0.0) {
        fac = 1.0;
        for (k=0;k<=nm;k++) {
            if (imag(z) < 0.0) {
                ck[k] = fac*ck[k]+cii*M_PI*ci[k];
            }
            else {
                ck[k] = fac*ck[k]-cii*M_PI*ci[k];
            }
            ci[k] *= fac;
            fac = -fac;
        }
    }
    cip[0] = ci[1];
    ckp[0] = -ck[1];
    for (k=1;k<=nm;k++) {
        cip[k] = ci[k-1]-(double)k*ci[k]/z;
        ckp[k] = -ck[k-1]-(double)k*ck[k]/z;
    }
    return 0;
}
int cbessikv(double v,complex<double>z,double &vm,complex<double> *civ,
    complex<double> *ckv,complex<double> *civp,complex<double> *ckvp)
{
    complex<double> z1,z2,ca1,ca,cs,cr,ci0,cbi0,cf,cf1,cf2;
    complex<double> ct,cp,cbk0,ca2,cr1,cr2,csu,cws,cb;
    complex<double> cg0,cg1,cgk,cbk1,cvk;
    double a0,v0,v0p,v0n,vt,w0,piv,gap,gan;
    int m,n,k,kz;

    a0 = abs(z);
    z1 = z;
    z2 = z*z;
    n = (int)v;
    v0 = v-n;
    piv = M_PI*v0;
    vt = 4.0*v0*v0;
    if (n == 0) n = 1;
    if (a0 < 1e-100) {
        for (k=0;k<=n;k++) {
            civ[k] = czero;
            ckv[k] = complex<double>(-1e308,0);
            civp[k] = czero;
            ckvp[k] = complex<double>(1e308,0);
        }
        if (v0 == 0.0) {
            civ[0] = cone;
            civp[1] = complex<double> (0.5,0.0);
        }
        vm = v;
        return 0;
    }
    if (a0 >= 50.0) kz = 8;
    else if (a0 >= 35.0) kz = 10;
    else kz = 14;
    if (real(z) <= 0.0) z1 = -z;
    if (a0 < 18.0) {
        if (v0 == 0.0) {
            ca1 = cone;
        }
        else {
            v0p = 1.0+v0;
            gap = gamma(v0p);
            ca1 = pow(0.5*z1,v0)/gap;
        }
        ci0 = cone;
        cr = cone;
        for (k=1;k<=50;k++) {
            cr *= 0.25*z2/(k*(k+v0));
            ci0 += cr;
            if (abs(cr/ci0) < eps) break;
        }
        cbi0 = ci0*ca1;
    }
    else {
        ca = exp(z1)/sqrt(2.0*M_PI*z1);
        cs = cone;
        cr = cone;
        for (k=1;k<=kz;k++) {
            cr *= -0.125*(vt-pow(2.0*k-1.0,2.0))/((double)k*z1);
            cs += cr;
        }
        cbi0 = ca*cs;
    }
    m = msta1(a0,200);
    if (m < n) n = m;
    else m = msta2(a0,n,15);
    cf2 = czero;
    cf1 = complex<double>(1.0e-100,0.0);
    for (k=m;k>=0;k--) {
        cf = 2.0*(v0+k+1.0)*cf1/z1+cf2;
        if (k <= n) civ[k] = cf;
        cf2 = cf1;
        cf1 = cf;
    }
    cs = cbi0/cf;
    for (k=0;k<=n;k++) {
        civ[k] *= cs;
    }
    if (a0 <= 9.0) {
        if (v0 == 0.0) {
            ct = -log(0.5*z1)-el;
            cs = czero;
            w0 = 0.0;
            cr = cone;
            for (k=1;k<=50;k++) {
                w0 += 1.0/k;
                cr *= 0.25*z2/(double)(k*k);
                cp = cr*(w0+ct);
                cs += cp;
                if ((k >= 10) && (abs(cp/cs) < eps)) break;
            }
            cbk0 = ct+cs;
        }
        else {
            v0n = 1.0-v0;
            gan = gamma(v0n);
            ca2 = 1.0/(gan*pow(0.5*z1,v0));
            ca1 = pow(0.5*z1,v0)/gap;
            csu = ca2-ca1;
            cr1 = cone;
            cr2 = cone;
            cws = czero;
            for (k=1;k<=50;k++) {
                cr1 *= 0.25*z2/(k*(k-v0));
                cr2 *= 0.25*z2/(k*(k+v0));
                csu += ca2*cr1-ca1*cr2;
                if ((k >= 10) && (abs((cws-csu)/csu) < eps)) break;
                cws = csu;
            }
            cbk0 = csu*M_PI_2/sin(piv);
        }
    }
    else {
        cb = exp(-z1)*sqrt(M_PI_2/z1);
        cs = cone;
        cr = cone;
        for (k=1;k<=kz;k++) {
            cr *= 0.125*(vt-pow(2.0*k-1.0,2.0))/((double)k*z1);
            cs += cr;
        }
        cbk0 = cb*cs;
    }
    cbk1 = (1.0/z1-civ[1]*cbk0)/civ[0];
    ckv[0] = cbk0;
    ckv[1] = cbk1;
    cg0 = cbk0;
    cg1 = cbk1;
    for (k=2;k<=n;k++) {
        cgk = 2.0*(v0+k-1.0)*cg1/z1+cg0;
        ckv[k] = cgk;
        cg0 = cg1;
        cg1 = cgk;
    }
    if (real(z) < 0.0) {
        for (k=0;k<=n;k++) {
            cvk = exp((k+v0)*M_PI*cii);
            if (imag(z) < 0.0) {
                ckv[k] = cvk*ckv[k]+M_PI*cii*civ[k];
                civ[k] /= cvk;
            }
            else if (imag(z) > 0.0) {
                ckv[k] = ckv[k]/cvk-M_PI*cii*civ[k];
                civ[k] *= cvk;
            }
        }
    }
    civp[0] = v0*civ[0]/z+civ[1];
    ckvp[0] = v0*ckv[0]/z-ckv[1];
    for (k=1;k<=n;k++) {
        civp[k] = -(k+v0)*civ[k]/z+civ[k-1];
        ckvp[k] = -(k+v0)*ckv[k]/z-ckv[k-1];
    }
    vm = n+v0;
    return 0;
}
//  cbessjy.cpp -- complex Bessel functions.
//  Algorithms and coefficient values from "Computation of Special
//  Functions", Zhang and Jin, John Wiley and Sons, 1996.
//
//  (C) 2003, C. Bond. All rights reserved.
//
//MH#include <complex.h>
//MH#include "bessel.h"
double gamma(double);

//MH static complex<double> cii(0.0,1.0);
//MH static complex<double> cone(1.0,0.0);
//MH static complex<double> czero(0.0,0.0);

int cbessjy01(complex<double> z,complex<double> &cj0,complex<double> &cj1,
    complex<double> &cy0,complex<double> &cy1,complex<double> &cj0p,
    complex<double> &cj1p,complex<double> &cy0p,complex<double> &cy1p)
{
    complex<double> z1,z2,cr,cp,cs,cp0,cq0,cp1,cq1,ct1,ct2,cu;
    double a0,w0,w1;
    int k,kz;

    static double a[] = {
        -7.03125e-2,
         0.112152099609375,
        -0.5725014209747314,
         6.074042001273483,
        -1.100171402692467e2,
         3.038090510922384e3,
        -1.188384262567832e5,
         6.252951493434797e6,
        -4.259392165047669e8,
         3.646840080706556e10,
        -3.833534661393944e12,
         4.854014686852901e14,
        -7.286857349377656e16,
         1.279721941975975e19};
    static double b[] = {
         7.32421875e-2,
        -0.2271080017089844,
         1.727727502584457,
        -2.438052969955606e1,
         5.513358961220206e2,
        -1.825775547429318e4,
         8.328593040162893e5,
        -5.006958953198893e7,
         3.836255180230433e9,
        -3.649010818849833e11,
         4.218971570284096e13,
        -5.827244631566907e15,
         9.476288099260110e17,
        -1.792162323051699e20};
    static double a1[] = {
         0.1171875,
        -0.1441955566406250,
         0.6765925884246826,
        -6.883914268109947,
         1.215978918765359e2,
        -3.302272294480852e3,
         1.276412726461746e5,
        -6.656367718817688e6,
         4.502786003050393e8,
        -3.833857520742790e10,
         4.011838599133198e12,
        -5.060568503314727e14,
         7.572616461117958e16,
        -1.326257285320556e19};
    static double b1[] = {
        -0.1025390625,
         0.2775764465332031,
        -1.993531733751297,
         2.724882731126854e1,
        -6.038440767050702e2,
         1.971837591223663e4,
        -8.902978767070678e5,
         5.310411010968522e7,
        -4.043620325107754e9,
         3.827011346598605e11,
        -4.406481417852278e13,
         6.065091351222699e15,
        -9.833883876590679e17,
         1.855045211579828e20};

    a0 = abs(z);
    z2 = z*z;
    z1 = z;
    if (a0 == 0.0) {
        cj0 = cone;
        cj1 = czero;
        cy0 = complex<double>(-1e308,0);
        cy1 = complex<double>(-1e308,0);
        cj0p = czero;
        cj1p = complex<double>(0.5,0.0);
        cy0p = complex<double>(1e308,0);
        cy1p = complex<double>(1e308,0);
        return 0;
    }
    if (real(z) < 0.0) z1 = -z;
    if (a0 <= 12.0) {
        cj0 = cone;
        cr = cone;
        for (k=1;k<=40;k++) {
            cr *= -0.25*z2/(double)(k*k);
            cj0 += cr;
            if (abs(cr) < abs(cj0)*eps) break;
        }
        cj1 = cone;
        cr = cone;
        for (k=1;k<=40;k++) {
            cr *= -0.25*z2/(k*(k+1.0));
            cj1 += cr;
            if (abs(cr) < abs(cj1)*eps) break;
        }
        cj1 *= 0.5*z1;
        w0 = 0.0;
        cr = cone;
        cs = czero;
        for (k=1;k<=40;k++) {
            w0 += 1.0/k;
            cr *= -0.25*z2/(double)(k*k);
            cp = cr*w0;
            cs += cp;
            if (abs(cp) < abs(cs)*eps) break;
        }
        cy0 = M_2_PI*((log(0.5*z1)+el)*cj0-cs);
        w1 = 0.0;
        cr = cone;
        cs = cone;
        for (k=1;k<=40;k++) {
            w1 += 1.0/k;
            cr *= -0.25*z2/(k*(k+1.0));
            cp = cr*(2.0*w1+1.0/(k+1.0));
            cs += cp;
            if (abs(cp) < abs(cs)*eps) break;
        }
        cy1 = M_2_PI*((log(0.5*z1)+el)*cj1-1.0/z1-0.25*z1*cs);
    }
    else {
        if (a0 >= 50.0) kz = 8;         // can be changed to 10
        else if (a0 >= 35.0) kz = 10;   //   "      "     "  12
        else kz = 12;                   //   "      "     "  14
        ct1 = z1 - M_PI_4;
        cp0 = cone;
        for (k=0;k<kz;k++) {
            cp0 += a[k]*pow(z1,-2.0*k-2.0);
        }
        cq0 = -0.125/z1;
        for (k=0;k<kz;k++) {
            cq0 += b[k]*pow(z1,-2.0*k-3.0);
        }
        cu = sqrt(M_2_PI/z1);
        cj0 = cu*(cp0*cos(ct1)-cq0*sin(ct1));
        cy0 = cu*(cp0*sin(ct1)+cq0*cos(ct1));
        ct2 = z1 - 0.75*M_PI;
        cp1 = cone;
        for (k=0;k<kz;k++) {
            cp1 += a1[k]*pow(z1,-2.0*k-2.0);
        }
        cq1 = 0.375/z1;
        for (k=0;k<kz;k++) {
            cq1 += b1[k]*pow(z1,-2.0*k-3.0);
        }
        cj1 = cu*(cp1*cos(ct2)-cq1*sin(ct2));
        cy1 = cu*(cp1*sin(ct2)+cq1*cos(ct2));
    }
    if (real(z) < 0.0) {
        if (imag(z) < 0.0) {
            cy0 -= 2.0*cii*cj0;
            cy1 = -(cy1-2.0*cii*cj1);
        }
        else if (imag(z) > 0.0) {
            cy0 += 2.0*cii*cj0;
            cy1 = -(cy1+2.0*cii*cj1);
        }
        cj1 = -cj1;
    }
    cj0p = -cj1;
    cj1p = cj0-cj1/z;
    cy0p = -cy1;
    cy1p = cy0-cy1/z;
    return 0;
}

int cbessjyna(int n,complex<double> z,int &nm,complex<double> *cj,
    complex<double> *cy,complex<double> *cjp,complex<double> *cyp)
{
    complex<double> cbj0,cbj1,cby0,cby1,cj0,cjk,cj1,cf,cf1,cf2;
    complex<double> cs,cg0,cg1,cyk,cyl1,cyl2,cylk,cp11,cp12,cp21,cp22;
    complex<double> ch0,ch1,ch2;
    double a0,yak,ya1,ya0,wa;
    int m,k,lb,lb0;

    if (n < 0) return 1;
    a0 = abs(z);
    nm = n;
    if (a0 < 1.0e-100) {
        for (k=0;k<=n;k++) {
            cj[k] = czero;
            cy[k] = complex<double> (-1e308,0);
            cjp[k] = czero;
            cyp[k] = complex<double>(1e308,0);
        }
        cj[0] = cone;
        cjp[1] = complex<double>(0.5,0.0);
        return 0;
    }
    cbessjy01(z,cj[0],cj[1],cy[0],cy[1],cjp[0],cjp[1],cyp[0],cyp[1]);
    cbj0 = cj[0];
    cbj1 = cj[1];
    cby0 = cy[0];
    cby1 = cy[1];
    if (n <= 1) return 0;
    if (n < (int)0.25*a0) {
        cj0 = cbj0;
        cj1 = cbj1;
        for (k=2;k<=n;k++) {
            cjk = 2.0*(k-1.0)*cj1/z-cj0;
            cj[k] = cjk;
            cj0 = cj1;
            cj1 = cjk;
        }
    }
    else {
        m = msta1(a0,200);
        if (m < n) nm = m;
        else m = msta2(a0,n,15);
        cf2 = czero;
        cf1 = complex<double> (1.0e-100,0.0);
        for (k=m;k>=0;k--) {
            cf = 2.0*(k+1.0)*cf1/z-cf2;
            if (k <=nm) cj[k] = cf;
            cf2 = cf1;
            cf1 = cf;
        }
        if (abs(cbj0) > abs(cbj1)) cs = cbj0/cf;
        else cs = cbj1/cf2;
        for (k=0;k<=nm;k++) {
            cj[k] *= cs;
        }
    }
    for (k=2;k<=nm;k++) {
        cjp[k] = cj[k-1]-(double)k*cj[k]/z;
    }
    ya0 = abs(cby0);
    lb = 0;
    cg0 = cby0;
    cg1 = cby1;
    for (k=2;k<=nm;k++) {
        cyk = 2.0*(k-1.0)*cg1/z-cg0;
        yak = abs(cyk);
        ya1 = abs(cg0);
        if ((yak < ya0) && (yak < ya1)) lb = k;
        cy[k] = cyk;
        cg0 = cg1;
        cg1 = cyk;
    }
    lb0 = 0;
    if ((lb > 4) && (imag(z) != 0.0)) {
        while (lb != lb0) {
            ch2 = cone;
            ch1 = czero;
            lb0 = lb;
            for (k=lb;k>=1;k--) {
                ch0 = 2.0*k*ch1/z-ch2;
                ch2 = ch1;
                ch1 = ch0;
            }
            cp12 = ch0;
            cp22 = ch2;
            ch2 = czero;
            ch1 = cone;
            for (k=lb;k>=1;k--) {
                ch0 = 2.0*k*ch1/z-ch2;
                ch2 = ch1;
                ch1 = ch0;
            }
            cp11 = ch0;
            cp21 = ch2;
            if (lb == nm)
                cj[lb+1] = 2.0*lb*cj[lb]/z-cj[lb-1];
            if (abs(cj[0]) > abs(cj[1])) {
                cy[lb+1] = (cj[lb+1]*cby0-2.0*cp11/(M_PI*z))/cj[0];
                cy[lb] = (cj[lb]*cby0+2.0*cp12/(M_PI*z))/cj[0];
            }
            else {
                cy[lb+1] = (cj[lb+1]*cby1-2.0*cp21/(M_PI*z))/cj[1];
                cy[lb] = (cj[lb]*cby1+2.0*cp22/(M_PI*z))/cj[1];
            }
            cyl2 = cy[lb+1];
            cyl1 = cy[lb];
            for (k=lb-1;k>=0;k--) {
                cylk = 2.0*(k+1.0)*cyl1/z-cyl2;
                cy[k] = cylk;
                cyl2 = cyl1;
                cyl1 = cylk;
            }
            cyl1 = cy[lb];
            cyl2 = cy[lb+1];
            for (k=lb+1;k<n;k++) {
                cylk = 2.0*k*cyl2/z-cyl1;
                cy[k+1] = cylk;
                cyl1 = cyl2;
                cyl2 = cylk;
            }
            for (k=2;k<=nm;k++) {
                wa = abs(cy[k]);
                if (wa < abs(cy[k-1])) lb = k;
            }
        }
    }
    for (k=2;k<=nm;k++) {
        cyp[k] = cy[k-1]-(double)k*cy[k]/z;
    }
    return 0;
}

int cbessjynb(int n,complex<double> z,int &nm,complex<double> *cj,
    complex<double> *cy,complex<double> *cjp,complex<double> *cyp)
{
    complex<double> cf,cf0,cf1,cf2,cbs,csu,csv,cs0,ce;
    complex<double> ct1,cp0,cq0,cp1,cq1,cu,cbj0,cby0,cbj1,cby1;
    complex<double> cyy,cbjk,ct2;
    double a0,y0;
    int k,m;
    static double a[] = {
        -0.7031250000000000e-1,
         0.1121520996093750,
        -0.5725014209747314,
         6.074042001273483};
    static double b[] = {
         0.7324218750000000e-1,
        -0.2271080017089844,
         1.727727502584457,
        -2.438052969955606e1};
    static double a1[] = {
         0.1171875,
        -0.1441955566406250,
         0.6765925884246826,
        -6.883914268109947};
    static double b1[] = {
       -0.1025390625,
        0.2775764465332031,
       -1.993531733751297,
        2.724882731126854e1};

    y0 = abs(imag(z));
    a0 = abs(z);
    nm = n;
    if (a0 < 1.0e-100) {
        for (k=0;k<=n;k++) {
            cj[k] = czero;
            cy[k] = complex<double> (-1e308,0);
            cjp[k] = czero;
            cyp[k] = complex<double>(1e308,0);
        }
        cj[0] = cone;
        cjp[1] = complex<double>(0.5,0.0);
        return 0;
    }
    if ((a0 <= 300.0) || (n > (int)(0.25*a0))) {
        if (n == 0) nm = 1;
        m = msta1(a0,200);
        if (m < nm) nm = m;
        else m = msta2(a0,nm,15);
        cbs = czero;
        csu = czero;
        csv = czero;
        cf2 = czero;
        cf1 = complex<double> (1.0e-100,0.0);
        for (k=m;k>=0;k--) {
            cf = 2.0*(k+1.0)*cf1/z-cf2;
            if (k <= nm) cj[k] = cf;
            if (((k & 1) == 0) && (k != 0)) {
                if (y0 <= 1.0) {
                    cbs += 2.0*cf;
                }
                else {
                    cbs += (-1)*((k & 2)-1)*2.0*cf;
                }
                csu += (double)((-1)*((k & 2)-1))*cf/(double)k;
            }
            else if (k > 1) {
                csv += (double)((-1)*((k & 2)-1)*k)*cf/(double)(k*k-1.0);
            }
            cf2 = cf1;
            cf1 = cf;
        }
        if (y0 <= 1.0) cs0 = cbs+cf;
        else cs0 = (cbs+cf)/cos(z);
        for (k=0;k<=nm;k++) {
            cj[k] /= cs0;
        }
        ce = log(0.5*z)+el;
        cy[0] = M_2_PI*(ce*cj[0]-4.0*csu/cs0);
        cy[1] = M_2_PI*(-cj[0]/z+(ce-1.0)*cj[1]-4.0*csv/cs0);
    }
    else {
        ct1 = z-M_PI_4;
        cp0 = cone;
        for (k=0;k<4;k++) {
            cp0 += a[k]*pow(z,-2.0*k-2.0);
        }
        cq0 = -0.125/z;
        for (k=0;k<4;k++) {
            cq0 += b[k] *pow(z,-2.0*k-3.0);
        }
        cu = sqrt(M_2_PI/z);
        cbj0 = cu*(cp0*cos(ct1)-cq0*sin(ct1));
        cby0 = cu*(cp0*sin(ct1)+cq0*cos(ct1));
        cj[0] = cbj0;
        cy[0] = cby0;
        ct2 = z-0.75*M_PI;
        cp1 = cone;
        for (k=0;k<4;k++) {
            cp1 += a1[k]*pow(z,-2.0*k-2.0);
        }
        cq1 = 0.375/z;
        for (k=0;k<4;k++) {
            cq1 += b1[k]*pow(z,-2.0*k-3.0);
        }
        cbj1 = cu*(cp1*cos(ct2)-cq1*sin(ct2));
        cby1 = cu*(cp1*sin(ct2)+cq1*sin(ct2));
        cj[1] = cbj1;
        cy[1] = cby1;
        for (k=2;k<=n;k++) {
            cbjk = 2.0*(k-1.0)*cbj1/z-cbj0;
            cj[k] = cbjk;
            cbj0 = cbj1;
            cbj1 = cbjk;
        }
    }
    cjp[0] = -cj[1];
    for (k=1;k<=nm;k++) {
        cjp[k] = cj[k-1]-(double)k*cj[k]/z;
    }
    if (abs(cj[0]) > 1.0)
        cy[1] = (cj[1]*cy[0]-2.0/(M_PI*z))/cj[0];
    for (k=2;k<=nm;k++) {
        if (abs(cj[k-1]) >= abs(cj[k-2]))
            cyy = (cj[k]*cy[k-1]-2.0/(M_PI*z))/cj[k-1];
        else
            cyy = (cj[k]*cy[k-2]-4.0*(k-1.0)/(M_PI*z*z))/cj[k-2];
        cy[k] = cyy;
    }
    cyp[0] = -cy[1];
    for (k=1;k<=nm;k++) {
        cyp[k] = cy[k-1]-(double)k*cy[k]/z;
    }

    return 0;
}

int cbessjyva(double v,complex<double> z,double &vm,complex<double>*cjv,
    complex<double>*cyv,complex<double>*cjvp,complex<double>*cyvp)
{
    complex<double> z1,z2,zk,cjvl,cr,ca,cjv0,cjv1,cpz,crp;
    complex<double> cqz,crq,ca0,cck,csk,cyv0,cyv1,cju0,cju1,cb;
    complex<double> cs,cs0,cr0,cs1,cr1,cec,cf,cf0,cf1,cf2;
    complex<double> cfac0,cfac1,cg0,cg1,cyk,cp11,cp12,cp21,cp22;
    complex<double> ch0,ch1,ch2,cyl1,cyl2,cylk;

    double a0,v0,pv0,pv1,vl,ga,gb,vg,vv,w0,w1,ya0,yak,ya1,wa;
    int j,n,k,kz,l,lb,lb0,m;

    a0 = abs(z);
    z1 = z;
    z2 = z*z;
    n = (int)v;


    v0 = v-n;

    pv0 = M_PI*v0;
    pv1 = M_PI*(1.0+v0);
    if (a0 < 1.0e-100) {
        for (k=0;k<=n;k++) {
            cjv[k] = czero;
            cyv[k] = complex<double> (-1e308,0);
            cjvp[k] = czero;
            cyvp[k] = complex<double> (1e308,0);

        }
        if (v0 == 0.0) {
            cjv[0] = cone;
            cjvp[1] = complex<double> (0.5,0.0);
        }
        else {
            cjv[0] = complex<double> (1e308,0);
        }
        vm = v;
        return 0;
    }
    if (real(z1) < 0.0) z1 = -z;
    if (a0 <= 12.0) {
        for (l=0;l<2;l++) {
            vl = v0+l;
            cjvl = cone;
            cr = cone;
            for (k=1;k<=40;k++) {
                cr *= -0.25*z2/(k*(k+vl));
                cjvl += cr;
                if (abs(cr) < abs(cjvl)*eps) break;
            }
           vg = 1.0 + vl;
           ga = gamma(vg);
           ca = pow(0.5*z1,vl)/ga;
           if (l == 0) cjv0 = cjvl*ca;
           else cjv1 = cjvl*ca;
        }
    }
    else {
        if (a0 >= 50.0) kz = 8;
        else if (a0 >= 35.0) kz = 10;
        else kz = 11;
        for (j=0;j<2;j++) {
            vv = 4.0*(j+v0)*(j+v0);
            cpz = cone;
            crp = cone;
            for (k=1;k<=kz;k++) {
                crp = -0.78125e-2*crp*(vv-pow(4.0*k-3.0,2.0))*
                    (vv-pow(4.0*k-1.0,2.0))/(k*(2.0*k-1.0)*z2);
                cpz += crp;
            }
            cqz = cone;
            crq = cone;
            for (k=1;k<=kz;k++) {
                crq = -0.78125e-2*crq*(vv-pow(4.0*k-1.0,2.0))*
                    (vv-pow(4.0*k+1.0,2.0))/(k*(2.0*k+1.0)*z2);
                cqz += crq;
            }
            cqz *= 0.125*(vv-1.0)/z1;
            zk = z1-(0.5*(j+v0)+0.25)*M_PI;
            ca0 = sqrt(M_2_PI/z1);
            cck = cos(zk);
            csk = sin(zk);
            if (j == 0) {
                cjv0 = ca0*(cpz*cck-cqz*csk);
                cyv0 = ca0*(cpz*csk+cqz+cck);
            }
            else {
                cjv1 = ca0*(cpz*cck-cqz*csk);
                cyv1 = ca0*(cpz*csk+cqz*cck);
            }
        }
    }
    if (a0 <= 12.0) {
        if (v0 != 0.0) {
            for (l=0;l<2;l++) {
                vl = v0+l;
                cjvl = cone;
                cr = cone;
                for (k=1;k<=40;k++) {
                    cr *= -0.25*z2/(k*(k-vl));
                    cjvl += cr;
                    if (abs(cr) < abs(cjvl)*eps) break;
                }
                vg = 1.0-vl;
                gb = gamma(vg);
                cb = pow(2.0/z1,vl)/gb;
                if (l == 0) cju0 = cjvl*cb;
                else cju1 = cjvl*cb;
            }
            cyv0 = (cjv0*cos(pv0)-cju0)/sin(pv0);
            cyv1 = (cjv1*cos(pv1)-cju1)/sin(pv1);
        }
        else {
            cec = log(0.5*z1)+el;
            cs0 = czero;
            w0 = 0.0;
            cr0 = cone;
            for (k=1;k<=30;k++) {
                w0 += 1.0/k;
                cr0 *= -0.25*z2/(double)(k*k);
                cs0 += cr0*w0;
            }
            cyv0 = M_2_PI*(cec*cjv0-cs0);
            cs1 = cone;
            w1 = 0.0;
            cr1 = cone;
            for (k=1;k<=30;k++) {
                w1 += 1.0/k;
                cr1 *= -0.25*z2/(k*(k+1.0));
                cs1 += cr1*(2.0*w1+1.0/(k+1.0));
            }
            cyv1 = M_2_PI*(cec*cjv1-1.0/z1-0.25*z1*cs1);
        }
    }
    if (real(z) < 0.0) {
        cfac0 = exp(pv0*cii);
        cfac1 = exp(pv1*cii);
        if (imag(z) < 0.0) {
            cyv0 = cfac0*cyv0-2.0*cii*cos(pv0)*cjv0;
            cyv1 = cfac1*cyv1-2.0*cii*cos(pv1)*cjv1;
            cjv0 /= cfac0;
            cjv1 /= cfac1;
        }
        else if (imag(z) > 0.0) {
            cyv0 = cyv0/cfac0+2.0*cii*cos(pv0)*cjv0;
            cyv1 = cyv1/cfac1+2.0*cii*cos(pv1)*cjv1;
            cjv0 *= cfac0;
            cjv1 *= cfac1;
        }
    }
    cjv[0] = cjv0;
    cjv[1] = cjv1;
    if ((n >= 2) && (n <= (int)(0.25*a0))) {
        cf0 = cjv0;
        cf1 = cjv1;
        for (k=2;k<= n;k++) {
            cf = 2.0*(k+v0-1.0)*cf1/z-cf0;
            cjv[k] = cf;
            cf0 = cf1;
            cf1 = cf;
        }
    }
    else if (n >= 2) {
        m = msta1(a0,200);
        if (m < n) n = m;
        else  m = msta2(a0,n,15);
        cf2 = czero;
        cf1 = complex<double>(1.0e-100,0.0);
        for (k=m;k>=0;k--) {
            cf = 2.0*(v0+k+1.0)*cf1/z-cf2;
            if (k <= n) cjv[k] = cf;
            cf2 = cf1;
            cf1 = cf;
        }
        if (abs(cjv0) > abs(cjv1)) cs = cjv0/cf;
        else cs = cjv1/cf2;
        for (k=0;k<=n;k++) {
            cjv[k] *= cs;
        }
    }
    cjvp[0] = v0*cjv[0]/z-cjv[1];
    for (k=1;k<=n;k++) {
        cjvp[k] = -(k+v0)*cjv[k]/z+cjv[k-1];
    }
    cyv[0] = cyv0;
    cyv[1] = cyv1;
    ya0 = abs(cyv0);
    lb = 0;
    cg0 = cyv0;
    cg1 = cyv1;
    for (k=2;k<=n;k++) {
        cyk = 2.0*(v0+k-1.0)*cg1/z-cg0;
        yak = abs(cyk);
        ya1 = abs(cg0);
        if ((yak < ya0) && (yak< ya1)) lb = k;
        cyv[k] = cyk;
        cg0 = cg1;
        cg1 = cyk;
    }
    lb0 = 0;
    if ((lb > 4) && (imag(z) != 0.0)) {
        while(lb != lb0) {
            ch2 = cone;
            ch1 = czero;
            lb0 = lb;
            for (k=lb;k>=1;k--) {
                ch0 = 2.0*(k+v0)*ch1/z-ch2;
                ch2 = ch1;
                ch1 = ch0;
            }
            cp12 = ch0;
            cp22 = ch2;
            ch2 = czero;
            ch1 = cone;
            for (k=lb;k>=1;k--) {
                ch0 = 2.0*(k+v0)*ch1/z-ch2;
                ch2 = ch1;
                ch1 = ch0;
            }
            cp11 = ch0;
            cp21 = ch2;
            if (lb == n)
                cjv[lb+1] = 2.0*(lb+v0)*cjv[lb]/z-cjv[lb-1];
            if (abs(cjv[0]) > abs(cjv[1])) {
                cyv[lb+1] = (cjv[lb+1]*cyv0-2.0*cp11/(M_PI*z))/cjv[0];
                cyv[lb] = (cjv[lb]*cyv0+2.0*cp12/(M_PI*z))/cjv[0];
            }
            else {
                cyv[lb+1] = (cjv[lb+1]*cyv1-2.0*cp21/(M_PI*z))/cjv[1];
                cyv[lb] = (cjv[lb]*cyv1+2.0*cp22/(M_PI*z))/cjv[1];
            }
            cyl2 = cyv[lb+1];
            cyl1 = cyv[lb];
            for (k=lb-1;k>=0;k--) {
                cylk = 2.0*(k+v0+1.0)*cyl1/z-cyl2;
                cyv[k] = cylk;
                cyl2 = cyl1;
                cyl1 = cylk;
            }
            cyl1 = cyv[lb];
            cyl2 = cyv[lb+1];
            for (k=lb+1;k<n;k++) {
                cylk = 2.0*(k+v0)*cyl2/z-cyl1;
                cyv[k+1] = cylk;
                cyl1 = cyl2;
                cyl2 = cylk;
            }
            for (k=2;k<=n;k++) {
                wa = abs(cyv[k]);
                if (wa < abs(cyv[k-1])) lb = k;
            }
        }
    }
    cyvp[0] = v0*cyv[0]/z-cyv[1];
    for (k=1;k<=n;k++) {
        cyvp[k] = cyv[k-1]-(k+v0)*cyv[k]/z;
    }
    vm = n+v0;
    return 0;
}


//  gamma.cpp -- computation of gamma function.
//      Algorithms and coefficient values from "Computation of Special
//      Functions", Zhang and Jin, John Wiley and Sons, 1996.
//
//  (C) 2003, C. Bond. All rights reserved.
//
// Returns gamma function of argument 'x'.
//
// NOTE: Returns 1e308 if argument is a negative integer or 0,
//      or if argument exceeds 171.
//
//MH#include <math.h>
double gamma(double x)
{
    int i,k,m;
    double ga,gr,r,z;

    static double g[] = {
        1.0,
        0.5772156649015329,
       -0.6558780715202538,
       -0.420026350340952e-1,
        0.1665386113822915,
       -0.421977345555443e-1,
       -0.9621971527877e-2,
        0.7218943246663e-2,
       -0.11651675918591e-2,
       -0.2152416741149e-3,
        0.1280502823882e-3,
       -0.201348547807e-4,
       -0.12504934821e-5,
        0.1133027232e-5,
       -0.2056338417e-6,
        0.6116095e-8,
        0.50020075e-8,
       -0.11812746e-8,
        0.1043427e-9,
        0.77823e-11,
       -0.36968e-11,
        0.51e-12,
       -0.206e-13,
       -0.54e-14,
        0.14e-14};

    if (x > 171.0) return 1e308;    // This value is an overflow flag.
    if (x == (int)x) {
        if (x > 0.0) {
            ga = 1.0;               // use factorial
            for (i=2;i<x;i++) {
               ga *= i;
            }
         }
         else
            ga = 1e308;
     }
     else {
        if (fabs(x) > 1.0) {
            z = fabs(x);
            m = (int)z;
            r = 1.0;
            for (k=1;k<=m;k++) {
                r *= (z-k);
            }
            z -= m;
        }
        else
            z = x;
        gr = g[24];
        for (k=23;k>=0;k--) {
            gr = gr*z+g[k];
        }
        ga = 1.0/(gr*z);
        if (fabs(x) > 1.0) {
            ga *= r;
            if (x < 0.0) {
                ga = -M_PI/(x*ga*sin(M_PI*x));
            }
        }
    }
    return ga;
}

}
