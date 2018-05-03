/* BSD 3-Clause License:
 * Copyright (c) 2018, bitsofcotton.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 *    Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *    Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation or other materials provided with the distribution.
 *    Neither the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#if !defined(_PPM2EIGEN_)

#include <Eigen/Core>
#include <iostream>
#include <fstream>

using std::max;
using std::min;
using std::cerr;
using std::endl;
using std::string;
using std::ifstream;
using std::getline;
using std::istringstream;
using std::ofstream;

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> rgb2l(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> rgb[3]) {
  // CIE 1931 XYZ from wikipedia.org
  auto X(.49000/.17697 * rgb[0] + .31000/.17697  * rgb[1] + .30000/.17697  * rgb[2]);
  auto Y(.17697/.17697 * rgb[0] + .81240/.17697  * rgb[1] + .010630/.17697 * rgb[2]);
  auto Z(                         .010000/.17697 * rgb[1] + .99000/.17697  * rgb[2]);
  return Y;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result(rgb[0].rows(), rgb[0].cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int j = 0; j < rgb[0].rows(); j ++)
    for(int k = 0; k < rgb[0].cols(); k ++)
      result(j, k) = sqrt(X(j, k) * X(j, k) + Y(j, k) * Y(j, k) + Z(j, k) * Z(j, k));
  return result;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tilt45(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& in, const bool& invert, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& orig = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>()) {
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> res;
  if(invert) {
    assert(orig.rows() + orig.cols() == in.rows() &&
           orig.rows() + orig.cols() == in.cols());
    res = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(orig.rows(), orig.cols());
    for(int i = 0; i < orig.rows(); i ++)
      for(int j = 0; j < orig.cols(); j ++)
        res(i, j) = in(i + j, orig.rows() - i + j);
  } else {
    res = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(in.rows() + in.cols(), in.cols() + in.rows());
    for(int i = - in.rows(); i < 2 * in.rows(); i ++)
      for(int j = - in.cols(); j < 2 * in.cols(); j ++) {
              int y(i + j);
        const int x(in.rows() - i + j);
        if(0 <= y && y < res.rows() && 0 <= x && x < res.cols())
          res(y, x) = in(abs(abs((i + in.rows()) % (2 * in.rows())) - in.rows()) % in.rows(),
                         abs(abs((j + in.cols()) % (2 * in.cols())) - in.cols()) % in.cols());
        y ++;
        if(0 <= y && y < res.rows() && 0 <= x && x < res.cols())
          res(y, x) = in(abs(abs((i + in.rows()) % (2 * in.rows())) - in.rows()) % in.rows(),
                         abs(abs((j + in.cols()) % (2 * in.cols())) - in.cols()) % in.cols());
      }
  }
  return res;
}

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
        if(++ i >= datas[0].cols()) {
          j ++;
          i = 0;
        }
        k = 0;
      }
    }
  }
  return true;
}

template <typename T> bool loadp2or3(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> data[3], const char* filename) {
  string line;
  string line2;
  string line3;
  ifstream input;
  input.open(filename);
  if(input.is_open()) {
    try {
      getline(input, line);
      while(line[0] == '#' && getline(input, line) && !input.eof() && !input.bad()) ;
      getline(input, line2);
      while(line2[0] == '#' && getline(input, line2) && !input.eof() && !input.bad()) ;
      getline(input, line3);
      while(line3[0] == '#' && getline(input, line3) && !input.eof() && !input.bad()) ;
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
  ofstream output;
  output.open(filename);
  if(output.is_open()) {
    try {
      if(gray)
        output << "P2" << "\n";
      else
        output << "P3" << "\n";
      output << data[0].cols() << " " << data[0].rows() << "\n";
      output << 255 << "\n";
      for(int i = 0; i < data[0].rows(); i ++)
        for(int j = 0; j < data[0].cols(); j ++)
          if(gray)
            output << int(data[0](i, j) * 255) << "\n";
          else
            for(int k = 0; k < 3; k ++)
              output << int(data[k](i, j) * 255) << "\n";
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
  T MM(data[0](0, 0)), mm(data[0](0, 0));
  for(int k = 0; k < 3; k ++)
    for(int i = 0; i < data[k].rows(); i ++)
      for(int j = 0; j < data[k].cols(); j ++) {
        MM = max(MM, data[k](i, j));
        mm = min(mm, data[k](i, j));
      }
  if(MM == mm)
    return;
  for(int k = 0; k < 3; k ++) {
    for(int i = 0; i < data[k].rows(); i ++)
      for(int j = 0; j < data[k].cols(); j ++)
        data[k](i, j) -= mm;
    data[k] *= upper / (MM - mm);
  }
  return;
}

#define _PPM2EIGEN_
#endif

