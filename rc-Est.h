#ifndef EST_H
#define EST_H

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))


struct data{

  size_t  n;
  double *t;
  double *y;
  double *sigma;
};

class dane{

 public:
  double x;
  double y;
  double s;

  dane& operator=(const dane&);  
};

class Est{

 public:
  
  std::string trig;
  
  double amp_s, alpha_s, beta_s, epeak_s;

  double amp_f, amp_f_err;
  double alpha_f, al_f_err;
  double beta_f, b_f_err;
  double epeak_f, ep_f_err;

  double t_start, t_end;

  std::vector<dane> v;

  std::string status;

  Est():amp_f(0), amp_f_err(0), alpha_f(0), al_f_err(0), beta_f(0), b_f_err(0),epeak_f(0), ep_f_err(0){}
    ~Est(){}
    
  std::string& trigger(void);

  double bf_fit_data(int,std::string&);

  bool operator==(const std::string&)const;

  Est& operator=(const Est&);

  friend std::istream& operator>>(std::istream&,Est&); 
  friend std::ostream& operator<<(std::ostream&,const Est&);
};

int band_f   (const gsl_vector*,void*,gsl_vector*);
int band_df  (const gsl_vector*, void*, gsl_matrix*);
int band_fdf (const gsl_vector*, void*, gsl_vector*, gsl_matrix*);


#endif
