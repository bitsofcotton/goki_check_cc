#if !defined(_OBJ2VECTOR_)

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using std::getline;
using std::istringstream;
using std::vector;

template <typename T> bool loadobj(vector<Eigen::Matrix<T, 3, 1> >& data, const char* filename) {
  ifstream input;
  input.open(filename);
  if(input.is_open()) {
    std::string work;
    while(getline(input, work)) {
      int i = 0;
      for(; i < work.size() && work[i] == ' '; i ++);
      if(i < work.size() && work[i] == 'v' &&
         i + 1 < work.size() && work[i + 1] != 'f' && work[i + 1] != 'n') {
        std::stringstream sub(work.substr(i + 2, work.size() - (i + 2)));
        Eigen::Matrix<T, 3, 1> buf;
        sub >> buf[0];
        sub >> buf[1];
        sub >> buf[2];
        data.push_back(buf);
      }
    }
  } else {
    cerr << "Unable to open file: " << filename << endl;
    return false;
  }
  return true;
 
}

#define _OBJ2VECTOR_
#endif

