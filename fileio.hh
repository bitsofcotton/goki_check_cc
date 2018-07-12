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

#if !defined(_FILEIO_GOKI_)

#include <Eigen/Core>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using std::max;
using std::min;
using std::cerr;
using std::endl;
using std::string;
using std::ifstream;
using std::getline;
using std::istringstream;
using std::stringstream;
using std::ofstream;
using std::sqrt;
using std::vector;
using std::pair;
using std::make_pair;

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

template <typename T> bool saveobj(const vector<Eigen::Matrix<T, 3, 1> >& data, const vector<Eigen::Matrix<int, 3, 1> >& polys, const char* filename, const bool& arout = false, const bool& addstand = false, const vector<vector<int> >& edges = vector<vector<int> >(), const T& zs = T(2)) {
  ofstream output;
  output.open(filename, std::ios::out);
  if(output.is_open()) {
    double Mh(1), Mw(1), lz(0);
    int lfs(0);
    for(int fslash(0) ; filename[fslash]; fslash ++)
      if(filename[fslash] == '/')
        lfs = fslash;
    if(lfs) lfs ++;
    output << "mtllib " << &filename[lfs] << ".mtl" << endl;
    for(int i = 0; i < data.size(); i ++) {
      Mh = max(data[i][1], Mh);
      Mw = max(data[i][0], Mw);
      lz = min(- data[i][2], lz);
    }
    if(arout) {
      assert(!addstand);
      for(int i = 0; i < data.size(); i ++)
        output << "v " << data[i][1] / Mh - T(.5) << " " << T(.5) - data[i][0] / Mw << " " << - (data[i][2] + lz / T(2)) / sqrt(Mw * Mh) << endl;
    } else if(addstand) {
      for(int i = 0; i < data.size(); i ++)
        output << "v " << data[i][1] << " " << Mw - data[i][0] << " " << - data[i][2] + zs - lz << endl;
      for(int i = 0; i < data.size(); i ++)
        output << "v " << data[i][1] << " " << Mw - data[i][0] << " " << 0 << endl;
    } else
      for(int i = 0; i < data.size(); i ++)
        output << "v " << data[i][1] << " " << - data[i][0] << " " << - data[i][2] << endl;
    for(int i = 0; i < data.size(); i ++)
      output << "vt " << data[i][1] / Mh << " " << 1. - data[i][0] / Mw << endl;
    output << "usemtl material0" << endl;
    // xchg with clockwise/counter clockwise.
    for(int i = 0; i < polys.size(); i ++) {
      output << "f " << polys[i][0] + 1 << "/" << polys[i][0] + 1 << "/" << polys[i][0] + 1;
      output << " "  << polys[i][1] + 1 << "/" << polys[i][1] + 1 << "/" << polys[i][1] + 1;
      output << " "  << polys[i][2] + 1 << "/" << polys[i][2] + 1 << "/" << polys[i][2] + 1 << endl;
      if(addstand) {
        output << "f " << data.size() + polys[i][0] + 1;
        output << " "  << data.size() + polys[i][2] + 1;
        output << " "  << data.size() + polys[i][1] + 1 << endl;
      }
    }
    if(addstand && edges.size())
      cerr << edges.size() << "parts found." << endl;
      for(int ii = 0; ii < edges.size(); ii ++) if(edges[ii].size()) {
        const vector<int>& outer(edges[ii]);
        for(int i = 0; i < outer.size() - 1; i ++) {
          const int& i0(i);
          const int  i1(i + 1);
          output << "f " << data.size() + outer[i0] + 1;
          output << " "  << outer[i1] + 1;
          output << " "  << data.size() + outer[i1] + 1 << endl;
          output << "f " << data.size() + outer[i0] + 1;
          output << " "  << outer[i0] + 1;
          output << " "  << outer[i1] + 1 << endl;
        }
      }
    output.close();
  } else {
    cerr << "Unable to open file: " << filename << endl;
    return false;
  }
  return true;
}

template <typename T> bool loadobj(vector<Eigen::Matrix<T, 3, 1> >& data, vector<Eigen::Matrix<int, 3, 1> >& polys, const char* filename) {
  ifstream input;
  input.open(filename);
  if(input.is_open()) {
    string work;
    while(getline(input, work) && !input.eof() && !input.bad()) {
      int i = 0;
      for(; i < work.size() && work[i] == ' '; i ++);
      if(i + 1 < work.size() && work[i] == 'v' && work[i + 1] == ' ') {
        stringstream sub(work.substr(i + 2, work.size() - (i + 2)));
        Eigen::Matrix<T, 3, 1> buf;
        sub >> buf[1];
        sub >> buf[0];
        sub >> buf[2];
        buf[0] = - buf[0];
        buf[2] = - buf[2];
        data.push_back(buf);
      } else if(i + 1 < work.size() && work[i] == 'f' && work[i + 1] == ' ') {
        stringstream sub(work.substr(i + 2, work.size() - (i + 2)));
        Eigen::Matrix<int, 3, 1> wbuf;
        int  widx(0);
        bool flag(false);
        while(!sub.eof() && !sub.bad()) {
          sub >> wbuf[widx];
          if(wbuf[widx] >= 0)
            wbuf[widx] --;
          widx ++;
          if(widx > 2)
            flag = true;
          if(flag)
            polys.push_back(wbuf);
          widx %= 3;
          if(sub.eof() || sub.bad())
            break;
          sub.ignore(20, ' ');
        }
      }
    }
    for(int i = 0; i < polys.size(); i ++)
      for(int j = 0; j < polys[i].size(); j ++) {
        while(polys[i][j] < 0) polys[i][j] += data.size();
        polys[i][j] = polys[i][j] % data.size();
      }
    input.close();
  } else {
    cerr << "Unable to open file: " << filename << endl;
    return false;
  }
  return true;
}

#define _FILEIO_GOKI_
#endif

