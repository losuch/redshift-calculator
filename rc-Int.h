#ifndef INT_H
#define INT_H

class Int{

 public:
  double I, I_err;
  double F, F_err;

  Int& operator=(const Int &i){
    this->I     = i.I;
    this->I_err = i.I_err;
    this->F     = i.F;
    this->F_err = i.F_err;
    return *this;
  }

  Int(){}
  ~Int(){}
};

#endif
