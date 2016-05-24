#ifndef PAR_H
#define PAR_H

#include "rc-Int.h"

class Par: Int{

 public:
  std::string trig;
  double t_start;
  double t_end;
  double amp;
  double amp_err;
  double epeak;
  double ep_err;
  double ebreak;
  double eb_err;
  double alpha;
  double al_err;
  double beta;
  double b_err;
  double pflux;
  double pfx_err;
  double pflnc;
  double pfc_err;
  double eflux;
  double efx_err;
  double eflnc;
  double efc_err;
  double redchi2;
  int    dof;

  Int itg;

  Par(){}
  ~Par(){}

  Par& operator=(const Par&);
  bool operator==(const Par&)const;

  // metoda wypisujaca triggery
  std::string& trigger(void);

  int integration(void);

  Par rand(void)const;

  // operatory wejscia wyjscia
  friend std::istream& operator>>(std::istream&,Par&);
  friend std::ostream& operator<<(std::ostream&,const Par&); 

  friend double rand_par(double,double);
};

double inver_gauss(double,double);
double gauss(double,double);
#endif
