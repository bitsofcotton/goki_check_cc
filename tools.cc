#include <cstdio>
#include <iostream>
#include <Eigen/Core>
#include "ppm2eigen.hh"
#include "enlarge.hh"

using namespace std;

void usage() {
  cout << "Usage: tools (enlarge|enlargeds|bump) <input filename>.p[gp]m <output filename>.p[gp]m";
  return;
}

int main(int argc, const char* argv[]) {
  if(argc < 4) {
    usage();
    return 0;
  }
  int mode = - 1;
  if(strcmp(argv[1], "enlarge") == 0)
    mode = 0;
  else if(strcmp(argv[1], "enlargeds") == 0)
    mode = 1;
  else if(strcmp(argv[1], "bump") == 0)
    mode = 2;
  if(mode < 0) {
    usage();
    return - 1;
  }
  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> data[3];
  if(!loadp2or3<float>(data, argv[2]))
    return - 2;
  switch(mode) {
  case 0:
    {
      enlarger2ex<float, complex<float> > enlarger;
      for(int i = 0; i < 3; i ++) {
        data[i] = enlarger.enlarge2(data[i], enlarger2ex<float, complex<float> >::ENLARGE_BOTH);
      }
    }
    break;
  case 1:
    {
      enlarger2exds<float, complex<float> > enlarger;
      for(int i = 0; i < 3; i ++)
        data[i] = enlarger.enlarge2ds(data[i], enlarger2exds<float, complex<float> >::ENLARGE_BOTH);
    }
    break;
  case 2:
    cerr << "not implemented now." << endl;
    break;
  default:
    break;
  }
  if(!savep2or3<float>(argv[3], data, ! true))
    return - 3;
  return 0;
}

