#include "MBUtils.hh"

std::string itos(int val) {
  char temp[32];
  sprintf(temp,"%i",val);
  return std::string(temp);
};

std::string ftos(double val) {
  char temp[32];
  sprintf(temp,"%.1f",val);
  return std::string(temp);
};

double power(double num, int pow) {
  double result = 1.;
  if (pow==0) return result;
  else if (pow>0) {
    for (int i=0; i<pow; i++) result = result*num;
    return result;
  }
  else if (pow<0) {
    for (int i=0; i>pow; i--) result = result/num;   
    return result;
  }
  else throw "Error in power somewhere";
};
