#include <iostream>
#include <fstream>
#include <cmath>

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
#include "rc-Est.h"

using namespace std;



dane& dane::operator=(const dane &d){

  this->x = d.x;
  this->y = d.y;
  this->s = d.s;
  return *this;
}


string& Est::trigger(void){
  
  return this->trig;
}

//------------------------------------------------------------------


double Est::bf_fit_data(int nr,string &stat){

  const gsl_multifit_fdfsolver_type *T;

  gsl_multifit_fdfsolver *s;
  
  int status;
  size_t i, iter = 0;

  const size_t n = nr;
  const size_t p = 4;

  gsl_matrix *covar = gsl_matrix_alloc (p, p);

  double t[nr], y[nr], sigma[nr];

  struct data d = { nr, t, y, sigma};

  gsl_multifit_function_fdf f;

  //double x_init[4] = {fpar->amp , fpar->alpha , fpar->beta , fpar->epeak };
  double x_init[4] = {this->amp_s , this->alpha_s , this->beta_s , this->epeak_s };

  gsl_vector_view x = gsl_vector_view_array (x_init, p);

  const gsl_rng_type *type;

  gsl_rng *r;

  gsl_rng_env_setup();

  type = gsl_rng_default;
  r    = gsl_rng_alloc(type);

  
  f.f      = &band_f;
  f.df     = &band_df;
  f.fdf    = &band_fdf;
  f.n      = n;
  f.p      = p;
  f.params = &d;

  
  for(i = 0; i < n ; i++){

    t[i]     = this->v[i].x;
    y[i]     = this->v[i].y;
    sigma[i] = this->v[i].s;
  }

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc(T, n, p);

  gsl_multifit_fdfsolver_set(s, &f, &x.vector);

  do{
    iter++;
    status = gsl_multifit_fdfsolver_iterate(s);
    
    if(status)
      break;
      
    status = gsl_multifit_test_delta (s->dx, s->x,
					1e-5, 1e-5);
  }while(status == GSL_CONTINUE && iter < 500);


  gsl_multifit_covar(s->J, 0.0, covar);

  double chi = gsl_blas_dnrm2(s->f);
  double dof = n - p;
  double c   = GSL_MAX_DBL(1, chi / sqrt(dof)); 


  this -> amp_f     = FIT(0);
  this -> amp_f_err = c*ERR(0);
  this -> alpha_f   = FIT(1);
  this -> al_f_err  = c*ERR(1);
  this -> beta_f    = FIT(2);
  this -> b_f_err   = c*ERR(2);
  this -> epeak_f   = FIT(3);
  this -> ep_f_err  = c*ERR(3);

 
  gsl_matrix_free(covar);
  stat = gsl_strerror(status);

  gsl_multifit_fdfsolver_free(s);
  //gsl_multifit_fdfsolver_type_free(T);
  //free(T);
  //free(type);
  free(r);

  return (pow(chi, 2.0));
}


//------------------------------------------------------------------

bool Est::operator==(const string &str)const{
  
  bool war = false;

  if(this->trig == str)
    war = true;

  return war;
}


//------------------------------------------------------------------


Est& Est::operator=(const Est &e){
  
  this->trig    = e.trig;
  this->amp_s   = e.amp_s;
  this->alpha_s = e.alpha_s;
  this->beta_s  = e.beta_s;
  this->epeak_s = e.epeak_s;
  
  
  return *this;
}

//------------------------------------------------------------------


istream& operator>>(istream &is,Est &p){
  
  is >> p.trig >> p.amp_s >> p.alpha_s >> p.beta_s >> p.epeak_s;

  return is;
}

ostream& operator<<(ostream &os,const Est &p){
  
  os << p.amp_s << "   " << p.alpha_s << "   " << p.beta_s << "   " << p.epeak_s;

  return os;
}


// -----------------------------------------------------------------------

int band_f(const gsl_vector *x, void *data, gsl_vector *f){


  size_t  n     = ((struct data *)data) -> n;
  double *t     = ((struct data *)data) -> t;
  double *y     = ((struct data *)data) -> y;
  double *sigma = ((struct data *)data) -> sigma;
     
  double A      = gsl_vector_get (x, 0);
  double alpha  = gsl_vector_get (x, 1);
  double beta   = gsl_vector_get (x, 2);
  double Ep     = gsl_vector_get (x, 3);

  double Yi;
  
  size_t i;
  
  for (i = 0; i < n; i++){
    /* BAND FUNCTION */

    Yi = band_func(t[i],A,alpha,beta,Ep);

    gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
  }

  //free(t);
  //free(y);
  //free(sigma);

  
  return GSL_SUCCESS;
}


// -----------------------------------------------------------------------


int band_df (const gsl_vector *x, void *data, gsl_matrix *J){


  size_t  n     = ((struct data *)data) -> n;
  double *t     = ((struct data *)data) -> t;
  double *sigma = ((struct data *)data) -> sigma;
  
  double A      = gsl_vector_get (x, 0);
  double alpha  = gsl_vector_get (x, 1);
  double beta   = gsl_vector_get (x, 2);
  double Ep     = gsl_vector_get (x, 3);
  
  size_t i;
  
  for (i = 0; i < n; i++){

      /* Jacobian matrix J(i,j) = dfi / dxj,              */
      /* where fi = (Yi - yi)/sigma[i],                   */
      /*       Yi = band_function(x; A, alpha, beta, Ep)  */
      /* and the xj are the parameters (A,alpha,beta,Ep)  */


      double s = sigma[i];
      
      double dydA = band_dydA(t[i],A,alpha,beta,Ep);
      double dyda = band_dyda(t[i],A,alpha,beta,Ep);
      double dydb = band_dydb(t[i],A,alpha,beta,Ep);
      double dydE = band_dydE(t[i],A,alpha,beta,Ep);

      gsl_matrix_set (J, i, 0, dydA/s); 
      gsl_matrix_set (J, i, 1, dyda/s);
      gsl_matrix_set (J, i, 2, dydb/s);
      gsl_matrix_set (J, i, 3, dydE/s);
    }

  //free(t);
  //free(sigma);

  return GSL_SUCCESS;
}


// -----------------------------------------------------------------------


int band_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J){


  band_f(x, data, f);
  band_df(x, data, J);
  
  return GSL_SUCCESS;
}

