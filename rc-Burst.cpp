#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <gsl/gsl_integration.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include <vector>
#include <string>

#include "rc.h"
#include "rc-Par.h"
#include "rc-Burst.h"
#include "rc-Est.h"

using namespace std;

// ************************** BURST ********************************


int Burst::read_burst(const Par &p, const vector<Par> &v){

  this->fpar = p;

  for(int i = 0 ; i < (int)v.size() ; i++){
     
    if(p == v[i])
      this->bpar.push_back(v[i]);
  }

  return (this->bpar.size());
}


//------------------------------------------------------------------

double Burst::integration(void){

  fpar.integration();

  for(int i = 0 ; i < (int)(this->bpar.size()) ; i++)
    this->bpar[i].integration();
  
  return this->fpar.itg.F;
}


// -----------------------------------------------------------------

string& Burst::trigger(void){

  return this->fpar.trig;
}

// -----------------------------------------------------------------

double Burst::compute_redshift(void){
  
  int nf, k;
  vector<double> nph, tph;
  double t_start, t_end;
 
  nf = this->mapping(nph,tph);
 


  this->N15 = this->n15_estimate(nph,tph,k);

  

  t_start = tph[k];
  t_end   = tph[k+1500];
  
  this->est.t_start = t_start;
  this->est.t_end   = t_end;
 
  this->generate_spectrum(this->bpar,this->est.v,t_start,t_end);

  string Status;
  double Chi2;
  Chi2 = this->est.bf_fit_data((int)this->est.v.size(),Status);

  this->est.status = Status;

  this->Ep = this->est.epeak_f;

  this->X = (this->N15)/(this->Ep);

  if(Chi2 <= 26.296)
    this->status = "success";

  this->chi2 = Chi2;

  //nph.empty();
  //tph.empty();

  return this->X;
}


//------------------------------------------------------------------


int Burst::mapping(vector<double> &nph,vector<double> &tph){

  int i = 0;
  int j = 0;

  double np;

  for(double t = this->fpar.t_start; t < this->fpar.t_end; t += 0.01, i++){

    if(t > this->bpar[j].t_end)
      j++;

    np = this->bpar[j].itg.I;
    
    nph.push_back(np);
    tph.push_back(t);
  }

  return (--i);
}


//------------------------------------------------------------------

  
double Burst::n15_estimate(vector<double> &nph,vector<double> &tph,int &h){
  
  double N15 = 0;
  double n15 = 0;
  int i, j;
  int k = 0;
  
  for( i = 1500 ; i < (int)nph.size() ; i++){
    
    n15 = 0;
    for(j = (i-1500) ; j < i ; j++){
      n15 += 0.01*nph[j];
    }
    if(n15>N15){
      N15 = n15;
      k   = i - 1500;
    }
  }
  
  h = k;
  
  return N15;
}


//------------------------------------------------------------------

  
void Burst::generate_spectrum(vector<Par> &par,vector<dane> &v,double t_start,double t_end){

  int i, j;

  double E;
  double dE = (E2 - E1)/ 20.;

  double sum, sum_err;

  double time;

  double dNdE;
  double dNdA, dNda, dNdb, dNdEp;

  double ts = t_start;
  double te = t_end;

  double s1,s2,s3,s4,sc;

  double tp,tk;

  
  for(E = E1, j = 0 ; E < E2 ; E += dE, j++){
    
    sum     = 0.0;
    sum_err = 0.0;

    for( i = 0 ; i < (int)par.size(); i++){

      tp = par[i].t_start;
      tk = par[i].t_end;

      if(tp <= ts && tk >= ts && tk <= te) //1

	time = tk - ts;

      else if(tp >= ts && tk <= te) //2

	time = tk - tp;

      else if(tp >= ts && tp <= te && tk >= te) //3

	time = te - tp;

      else if(tk <= ts || tp >= te) //4

	time = 0;

      else if(tp <= ts && tk >= te) //5

	time = te - ts;

      else{
	time = 0;
	fprintf(stderr,"generate_spectrum: Unexpected exception.\n");
	exit(EXIT_FAILURE);
      }

      sum += band_func(E,par[i].amp,par[i].alpha,par[i].beta,par[i].epeak)*time;

      
      dNdA  = band_dydA(E,par[i].amp,par[i].alpha,par[i].beta,par[i].epeak);
      dNda  = band_dyda(E,par[i].amp,par[i].alpha,par[i].beta,par[i].epeak);
      dNdb  = band_dydb(E,par[i].amp,par[i].alpha,par[i].beta,par[i].epeak);
      dNdEp = band_dydE(E,par[i].amp,par[i].alpha,par[i].beta,par[i].epeak);

      s1 = dNdA  * par[i].amp_err;
      s2 = dNda  * par[i].al_err;
      s3 = dNdb  * par[i].b_err;
      s4 = dNdEp * par[i].ep_err;
      
      sc = sqrt( s1*s1 + s2*s2 + s3*s3 + s4*s4);
      
      sum_err += sc*time;
    }

    dNdE = sum / 15.;
    
    dane w;
    
    w.x = E;
    w.y = dNdE;
    //w.s = sum_err/15.;
    w.s = 0.01;

    v.push_back(w);

  }
}

//------------------------------------------------------------------

Burst& Burst::operator=(const Burst &b){

  this->fpar = b.fpar;
  this->est  = b.est;
  
  for(int i = 0 ; i < (int)b.bpar.size() ; i++){
    Par p = b.bpar[i];
    this->bpar.push_back(p);
  }

  return *this;
}

//------------------------------------------------------------------


bool Burst::operator==(const string &str)const{

  bool war;

  if(this->fpar.trig == str)
    war = true;
  else
    war = false;
  return war;
}


//------------------------------------------------------------------


void Burst::print_burst(FILE *out){

  fprintf(out,"%s\t%.6f\t%.6f\t%.6f\t%s\n",this->fpar.trig.c_str(),this->X,this->X_err,this->dX,this->status.c_str());
  //fprintf(out,"%s\n",this->fpar.trig.c_str());
}


//------------------------------------------------------------------

Burst Burst::rand(void)const{

  Burst r;

  r.fpar = this->fpar.rand();

  for(int i = 0 ; i < (int)this->bpar.size() ; i++){
    r.bpar.push_back(this->bpar[i].rand());
  }

  r.est = this->est;
  
  return r;
}


//------------------------------------------------------------------


ostream& operator<<(ostream &os,const Burst &b){
  
  os << b.fpar.trig << "\t" << b.X << "\t" << b.X_err << "\t" << b.dX << "\t" << b.status;
  
  return os;
}


//------------------------------------------------------------------


void compute_redshift_error(Burst &b){

  vector<double> v2;
  srand((unsigned int)time(NULL));
  for(int i = 0 ; i < 1000 ; i++){

    if(i == 0){
      fprintf(stderr,"%s\t\t0.00",b.trigger().c_str());
    }
    else if(i<999){
      if(!(i % 50)){
	fprintf(stderr,"\b\b\b\b%.2f",(double)i/1000.);
      }
    }
    else{
      fprintf(stderr,"\b\b\b\b1.00\tOK\n");
    }
    
    Burst mc;
    double rds;

    mc = b.rand(); // losowanie
    
    
    
    mc.integration();
    rds = mc.compute_redshift();


    v2.push_back(rds);

  }

  sort(v2.begin(),v2.end());

  double Xerr = v2[842] - v2[158];

  b.X_err = Xerr;

  b.dX = sqrt((b.X - v2[500])*(b.X - v2[500]));
  
  
}

//------------------------------------------------------------------

void Burst::calibrate(void){

  this->npz     = x2z(this->X);
  this->npz_err = x2z_err_syst(this->X,this->X_err);
}


//------------------------------------------------------------------

double x2z(double x){

  //double f = A*log10(x)*log10(x) + B*log10(x) + C;

  double Log = log(10);
  double f   = _C + _B*log(x)/Log + _A*log(x)*log(x)/Log/Log;
  double g   = pow(10,f);

  return g;
}

// --------------------------------------------------------------------

double x2z_err(double x, double x_err){

  double Log = log(10.);
  double f   = x2z(x) * Log * ( _B/x/Log + 2*_A*log(x)/x/Log/Log ) * x_err;

  f = fabs(f);
  return f;
}

// --------------------------------------------------------------------

double da(double x){
  
  return fabs((x2z(x)*log(x)*log(x)/log(10.))*_DA);
}


double db(double x){

  return fabs(x2z(x)*log(x)*_DB);
}


double dc(double x){

  return fabs(x2z(x)*log(10.)*_DC);
}

// --------------------------------------------------------------------

double x2z_err_syst(double x, double x_err){

  double syst2 = da(x)*da(x) + db(x)*db(x) + dc(x)*dc(x); 
  
  double f     = sqrt( x2z_err(x,x_err)*x2z_err(x,x_err) + syst2 );

  return f;
}

// --------------------------------------------------------------------
