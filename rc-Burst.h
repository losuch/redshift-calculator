#ifndef BURST_H
#define BURST_H

#include "rc-Par.h"
#include "rc-Est.h"

// parametry estymacji funkcji: z(X) = 10^( A*log10(X)^2 + B*log10(X) + C )

#define _A   -0.055462
#define _B   -0.373666
#define _C    0.266277

#define _DA   0.06893
#define _DB   0.09794
#define _DC   0.03562


class Burst: Par, Est{

 private:

  Par fpar;
  std::vector<Par> bpar;

  int    mapping(std::vector<double>&,std::vector<double>&);
  double n15_estimate(std::vector<double>&,std::vector<double>&,int&);
  void   generate_spectrum(std::vector<Par>&,std::vector<dane>&,double,double);

 public:

  Est est;

  double chi2;

  double N15;
  double Ep;

  double X;
  double X_err;
  double dX;

  double npz, npz_err;

  std::string status;

  Burst(): chi2(0), N15(0), Ep(0), X(0), X_err(0), dX(0), npz(0), npz_err(0), status("failure"){}// constructor
  ~Burst(){
    //bpar.empty();
  }      // destructor

  // metoda wczytujace dane do klasy z pliku
  int read_burst(const Par&,const std::vector<Par>&);

  std::string& trigger(void);
  
  // metoda calkujaca funkcje banda 
  // (wyznaczenie strumienia fotonow i energji)
  double integration(void);
  
  // metoda obliczajaca przesuniecie ku czerwieni
  double compute_redshift(void);

  Burst rand(void)const;

  void print_burst(FILE*);

  bool operator==(const std::string&)const;
  Burst& operator=(const Burst&);
  

  friend std::ostream& operator<<(std::ostream&,const Burst&);

  friend void compute_redshift_error(Burst&);



  void calibrate(void);
};

void compute_redshift_error(Burst&);

double x2z(double);
double x2z_err(double,double);
double da(double);
double db(double);
double dc(double);
double x2z_err_syst(double,double);

#endif
