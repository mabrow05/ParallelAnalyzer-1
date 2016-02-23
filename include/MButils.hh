#include <cstdlib>



std::string itos(int val) {
  char temp[32];
  sprintf(temp,"%i",val);
  return std::string(temp);
};

std::string ftos(double val) {
  char temp[32];
  sprintf(temp,"%f",val);
  return std::string(temp);
};
