#ifndef RC_H
#define RC_H

#define E1  50.
#define E2  300.

#define MAX 500000

// ---------------------------------------------------------------


struct constant{
  double A;
  double APH;
  double BTA;
  double EP;
};

// ---------------------------------------------------------------

float func_n(float); 
float func_f(float); 

double band_func(double,double,double,double,double);
double band_dydA(double,double,double,double,double);
double band_dyda(double,double,double,double,double);
double band_dydb(double,double,double,double,double);
double band_dydE(double,double,double,double,double);
double band_flux(double,double,double,double,double);

#endif
