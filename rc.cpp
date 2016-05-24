#include <iostream>
#include <fstream>
#include <cmath>

#include <gsl/gsl_integration.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>


#include "rc.h"

using namespace std;




// ---------------------------------------------------------------------------

float func_n(float x) {

  extern double A, APH, BTA, EP;


  float f = band_func(x,A,APH,BTA,EP);
  return f;
}

// ---------------------------------------------------------------------------

float func_f(float x) {

  extern double A, APH, BTA, EP;

  float f = band_flux(x,A,APH,BTA,EP);
  return f;
}


// ---------------------------------------------------------------------------


double band_func(double x,double A,double a,double b,double Ep){

  
  double F;

  double Ec = (a-b)*Ep/(2+a);

  if(x < Ec)
    F = A*pow((x/100.),a)*exp(-x*(2+a)/Ep);
  else
    F = A*pow(((a-b)*Ep/100./(2+a)),(a-b))
        * exp(b-a)*pow((x/100.),b);

  return F;
}


// -----------------------------------------------------------------------


double band_flux(double x,double A,double a,double b,double Ep){

  
  double F;

  double Ec = (a-b)*Ep/(2+a);

  if(x < Ec)
    F = x*A*pow((x/100.),a)*exp(-x*(2+a)/Ep);
  else
    F = x*A*pow(((a-b)*Ep/100./(2+a)),(a-b))
        * exp(b-a)*pow((x/100.),b);

  return F;
}


// -----------------------------------------------------------------------


double band_dydA(double x,double A,double a,double b,double Ep){

  double F;
  A = 1;

  double Ec = (a-b)*Ep/(2+a);

  double p = pow((1./100.),a);

  if(x < Ec)
    F = p * exp((2+a)*x/Ep) * pow(x,a);
  else
    F = p * exp(b-a) * pow(((a-b)*Ep/(2+a)),a-b) * pow(x,b);

  return F;
}


// -----------------------------------------------------------------------


double band_dyda(double x,double A,double a,double b,double Ep){

  double F;

  double Ec = (a-b)*Ep/(2+a);

  double p  = pow((1./100.),a);
  double e  = exp((2+a)*x/Ep);

  double ex = exp(b-a);
  double pw = pow(((a-b)*Ep/(2+a)),(a-b));

  if(x < Ec)
    F = p * A * e * pow(x,(1+a)) / Ep - p * A * e * pow(x,a) * log(100.) +
        p * A * e * pow(x,a) * log(x);
  else
    F = - p * A * ex * pw * pow(x,b) -
          p * A * ex * pw * pow(x,b) * log(100.) +
          p * A * ex * pw * pow(x,b) *
      ( ( (2+a) * (Ep/(2+a) - (a-b)*Ep/pow(2+a,2)) ) / Ep + log((a-b)*Ep/(2+a)));

  return F;
}


// -----------------------------------------------------------------------


double band_dydb(double x,double A,double a,double b,double Ep){

  double F;

  double Ec = (a-b)*Ep/(2+a);

  double p  = pow((1./100.),a);
  double pw = pow(((a-b)*Ep/(2+a)),(a-b));
  double ex = exp(b-a);
  double xb = pow(x,b);

  if(x < Ec)
    F = 0;
  else
    F = p * A * ex + pw * xb +
        p * A * ex * pw * xb * (-1 - log((a-b)*Ep)/(2+a)) +
        p * A * ex * pw * xb * log(x);

  return F;
}


// -----------------------------------------------------------------------


double band_dydE(double x,double A,double a,double b,double Ep){
  
  double F;

  double Ec = (a-b)*Ep/(2+a);

  double p  = pow((1./100.),a);
  double ex = exp(b-a);
  double e  = exp((2+a)*x/Ep);
  double xb = pow(x,b);
  double fr = (a-b)*Ep/(2+a);


  if(x < Ec)
    F = - p * (2+a) * A * e * pow(x,1+a)/Ep/Ep;
  else
    F = p * A * pow(a-b,2) * ex * pow(fr,a-b-1) * xb / (2+a);

  return F;
}


// -----------------------------------------------------------------------
