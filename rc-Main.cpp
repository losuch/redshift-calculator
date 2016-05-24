#include <iostream>
#include <fstream>
#include <cmath>

#include <vector>
#include <string>
#include <algorithm>

#include <gsl/gsl_integration.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>


#include "rc.h"
#include "rc-Par.h"
#include "rc-Burst.h"
#include "rc-Est.h"


using namespace std;

const string FPAR = "data/FPAR.dat";
const string BPAR = "data/BPAR.dat";
const string LIST = "data/LIST.dat";
const string STAR = "data/STAR.dat";


int main(int argc, char **argv){


  
  if(argc<2){
    cerr << "Try `rc --help' for more information." << endl;
    exit(EXIT_FAILURE);
  }


  // HELP
  if(!strcmp(argv[1],"--help")){
    cout << "\n\n\t\t *** Redshift Calculator *** \n\t\t\t\t\t(c)2007 LO\n\n\n" 
	 << "\trc [ --help | -list | -trig 'trigger_number' ]\n\nOPTION:\n"
	 << "   --help\n\tshow this help.\n\n"
	 << "   -list\n\tlist the trigger numbers.\n\n"
	 << "   -trig 'trigger_number'\n\tcalculate redshift for 'trigger_number'\n"
	 << endl;
    exit(EXIT_SUCCESS);
  }


  
  vector<string>::iterator st;

  vector<string> list;

  ifstream in1(LIST.c_str());

  while(true){

    string str;

    in1 >> str;

    if(!in1) break;
   
    list.push_back(str);
  }

  in1.close();

  if(!strcmp(argv[1],"-list")){
    for(int i = 0 ; i < (int)list.size() ; i++){
      cout << list[i] << endl;
    }
    exit(EXIT_SUCCESS);
  }

  if(argc<3 || strcmp(argv[1],"-trig")){

    cerr << "Podales zla liczbe lub niewlasciwe argumenty" << endl;
    exit(EXIT_FAILURE);
  }


  
  vector<Par> pfpar,pbpar;

  ifstream in2,in3;

  in2.open(FPAR.c_str());

  while(true){
    
    Par p;

    in2 >> p;

    if(!in2) break;

    if(argc == 3){
      if(!strcmp("-trig",argv[1])){
	if(strcmp(argv[2],p.trigger().c_str())) continue;
      }
    }
    
    st = find(list.begin(),list.end(),p.trigger());
    if(st == list.end()) continue;
    
    if((p.t_end - p.t_start)<15.) continue;

    pfpar.push_back(p);
  }

  in2.close();
  

  in3.open(BPAR.c_str());
  
  while(true){

    Par p;

    in3 >> p;

    if(!in3) break;

    pbpar.push_back(p);

  
  }


  in3.close();

  
  vector<Burst> redsh;
  vector<Burst>::iterator rt;

  for(int i = 0 ; i < (int) pfpar.size() ; i++){
    
    Burst b;

    if(b.read_burst(pfpar[i],pbpar)){
      redsh.push_back(b);
    }
  }

  if((int)redsh.size() == 0){
    cerr << "Nie znaleziono blysku o podanym numerze trigger\nzobacz -list" << endl;
    exit(EXIT_FAILURE);
  }

  // wczytywanie warunkow poczatkowych do fitowania -----------
  ifstream in;
  in.open(STAR.c_str());

  vector<Est> vest;
  vector<Est>::iterator et;

  while(true){
    
    Est p;

    in >> p;

    if(!in) break;

    vest.push_back(p);
  }
  in.close();


  
  for(rt = redsh.begin() ; rt != redsh.end() ; rt++){
    
    et = find(vest.begin(),vest.end(),rt->trigger());

    if(et == vest.end()){
      cout << "No trigger: " << (rt->trigger()) << endl;
      exit(EXIT_FAILURE);
    }
    rt->est = *et;
  }
  
  for(rt = redsh.begin() ; rt != redsh.end() ; rt++)
    rt -> integration();


  
  for(int i = 0 ; i < (int)redsh.size() ; i++){
    
    redsh[i].compute_redshift();
    compute_redshift_error(redsh[i]);
    redsh[i].calibrate();
    printf("%s\t%.6f\t%.6f\t%s\n",redsh[i].trigger().c_str(),redsh[i].npz,redsh[i].npz_err,redsh[i].status.c_str());
  }



  return EXIT_SUCCESS;
}
