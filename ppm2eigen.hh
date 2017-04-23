#if !defined(_PPM2EIGEN_)

#include <Eigen/Core>
#include <iostream>
#include <fstream>

using namespace std;

template <typename T> bool loadstub(ifstream& input, const int& nmax, const int& ncolor, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>* datas) {
  int i = 0, j = 0, k = 0;
  char buf;
  int  work = 0;
  bool mode = false;
  while(input.get(buf) && j < datas[0].rows()) {
    if('0' <= buf && buf <= '9') {
      work *= 10;
      work += buf - '0';
      mode  = true;
      continue;
    } else if(mode) {
      mode = false;
      datas[k](j, i) = T(work) / nmax;
      work = 0;
      if(++ k >= ncolor) {
        i ++;
        k = 0;
      }
      if(i >= datas[0].cols()) {
        j ++;
        i = 0;
      }
    }
  }
  return true;
}

template <typename T> bool loadp2or3(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> data[3], const char* filename) {
  using namespace std;
  string line;
  string line2;
  string line3;
  ifstream input;
  input.open(filename);
  if(input.is_open()) {
    try {
      getline(input, line);
      while(line[0] == '#' && getline(input, line)) ;
      getline(input, line2);
      while(line2[0] == '#' && getline(input, line2)) ;
      getline(input, line3);
      while(line3[0] == '#' && getline(input, line3)) ;
      istringstream iline2(line2);
      int w, h;
      iline2 >> w;
      iline2 >> h;
      if(w <= 0 || h <= 0) {
        cerr << "unknown size." << endl;
        input.close();
        return false;
      } 
      istringstream iline3(line3);
      int nmax;
      iline3 >> nmax;
      if(line[0] == 'P') {
        if(line[1] == '2') {
          data[0] = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(h, w);
          loadstub<T>(input, nmax, 1, data);
        } else if(line[1] == '3') {
          for(int i = 0; i < 3; i ++)
            data[i] = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(h, w);
          loadstub<T>(input, nmax, 3, data);
        } else {
          cerr << "unknown file type." << endl;
          input.close();
          return false;
        }
      } else {
        cerr << "unknown file type." << endl;
        input.close();
        return false;
      }
    } catch (const ifstream::failure& e) {
      cerr << "Exception while reading." << endl;
    } catch (const istringstream::failure& e) {
      cerr << "Exception while reading." << endl;
    }
    input.close();
  } else {
    cerr << "Unable to open file: " << filename << endl;
    return false;
  }
  return true;
}

template <typename T> bool savep2or3(const char* filename, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> data[3], const bool& gray) {
  using namespace std;
  ofstream output;
  output.open(filename);
  if(output.is_open()) {
    try {
      if(gray)
        output << "P2" << endl;
      else
        output << "P3" << endl;
      output << data[0].cols() << " " << data[0].rows() << endl;
      output << 255 << endl;
      for(int i = 0; i < data[0].rows(); i ++)
        for(int j = 0; j < data[0].cols(); j ++)
          if(gray)
            output << int(data[0](i, j) * 255) << endl;
          else
            for(int k = 0; k < 3; k ++)
              output << int(data[k](i, j) * 255) << endl;
    } catch (const ofstream::failure& e) {
      cerr << "An error has occured while writing file." << endl;
    }
    output.close();
  } else {
    cerr << "Unable to open file: " << filename << endl;
    return false;
  }
  return true;
}

template <typename T> void normalize(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> data[3], const T& upper) {
  T MM(0), mm(0);
  for(int k = 0; k < 3; k ++)
    for(int i = 0; i < data[k].rows(); i ++)
      for(int j = 0; j < data[k].cols(); j ++) {
        MM = std::max(MM, data[k](i, j));
        mm = std::min(mm, data[k](i, j));
      }
  if(MM == mm)
    return;
  for(int k = 0; k < 3; k ++) {
    for(int i = 0; i < data[k].rows(); i ++)
      for(int j = 0; j < data[k].cols(); j ++)
        data[k](i, j) -= mm;
    data[k] *= upper / (MM - mm);
  }
}

#define _PPM2EIGEN_
#endif

