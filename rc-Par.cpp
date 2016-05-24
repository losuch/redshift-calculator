#include <iostream>
#include <fstream>
#include <cmath>

#include <cstdlib>
#include <ctime>

#include <gsl/gsl_integration.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include <string>
#include <vector>

#include "rc.h"
#include "rc-Par.h"
#include "rc-Burst.h"
#include "rc-Int.h"
#include "trapzd.c"

using namespace std;

double A;
double APH;
double BTA;
double EP;
double Ec;


// ************************* PAR ***********************************

Par& Par::operator=(const Par &p){
  
  this->trig    = p.trig;
  this->t_start = p.t_start;
  this->t_end   = p.t_end;
  this->amp     = p.amp;
  this->amp_err = p.amp_err;
  this->epeak   = p.epeak;
  this->ep_err  = p.ep_err;
  this->ebreak  = p.ebreak;
  this->eb_err  = p.eb_err;
  this->alpha   = p.alpha;
  this->al_err  = p.al_err;
  this->beta    = p.beta;
  this->b_err   = p.b_err;
  this->pflux   = p.pflux;
  this->pfx_err = p.pfx_err;
  this->pflnc   = p.pflnc;
  this->pfc_err = p.pfc_err;
  this->eflux   = p.eflux;
  this->efx_err = p.efx_err;
  this->eflnc   = p.eflnc;
  this->efc_err = p.efc_err;
  this->redchi2 = p.redchi2;
  this->dof     = p.dof;

  this->itg     = p.itg;

  return *this;
}


//------------------------------------------------------------------


bool Par::operator==(const Par& p)const{

  bool war;

  if(this->trig == p.trig)
    war = true;
  else
    war = false;

  return war;
}


//------------------------------------------------------------------


string& Par::trigger(void){

  return this->trig;
}


//------------------------------------------------------------------


int Par::integration(void){
  
  double I_res; //, I_res_err;
  double F_res; //, F_res_err;

  /*
  struct constant con = {this->amp , this->alpha , this->beta , this->epeak};

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

  gsl_function F1,F2;

  F1.function = &func_n;
  F1.params   = &con;

  F2.function = &func_f;
  F2.params   = &con;

  gsl_integration_qags (&F1, E1, E2, 0, 1e-5,1000 ,w, &I_res, &I_res_err);
  gsl_integration_qags (&F2, E1, E2, 0, 1e-5,1000 ,w, &F_res, &F_res_err);
  */

  A   = this -> amp;
  APH = this -> alpha;
  BTA = this -> beta;
  EP  = this -> epeak;
  Ec  = (APH-BTA)*EP/(2+APH);
  

  I_res = 2*trapzd(func_n,E1,E2,9);
  F_res = 2*trapzd(func_f,E1,E2,9);

  this -> itg.I = I_res;
  this -> itg.F = F_res;

  

  
  
  return EXIT_SUCCESS;
}


//------------------------------------------------------------------

Par Par::rand(void)const{

  Par p;
  int i = 0;

  p = *this;
  do{
   

    p.amp   = rand_par(this->amp,this->amp_err);
    p.alpha = rand_par(this->alpha,this->al_err);
    p.beta  = rand_par(this->beta,this->b_err);
    p.epeak = rand_par(this->epeak,this->ep_err);

    if(i>1000){
      cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
      exit(EXIT_FAILURE);
    }

    i++;
  }while(p.amp <= 0 || /* p.alpha < -2 || p.beta >= -2 ||*/ p.epeak <= 0);
 
  return p;
}

//------------------------------------------------------------------


istream& operator>>(istream& is,Par &p){

  is >> p.trig 
     >> p.t_start
     >> p.t_end
     >> p.amp
     >> p.amp_err
     >> p.epeak
     >> p.ep_err
     >> p.ebreak
     >> p.eb_err
     >> p.alpha
     >> p.al_err
     >> p.beta
     >> p.b_err
     >> p.pflux
     >> p.pfx_err
     >> p.pflnc
     >> p.pfc_err
     >> p.eflux
     >> p.efx_err
     >> p.eflnc
     >> p.efc_err
     >> p.redchi2
     >> p.dof;

  return is;
}


//------------------------------------------------------------------


ostream& operator<<(ostream &os,const Par &p){

  os << p.trig;

  return os;
}


//------------------------------------------------------------------


double rand_par(double x,double er){

  double par;

  //srand((unsigned int)time(NULL));         //- generator liczb pseudo losowych
  double l = (double)(rand() % 1000000 + 1); //- losuje liczbe z przedzialu <1,1000>
  
  double s = l/1000000.;

  par = x + (2*s-1)*er;

  return par;
}

