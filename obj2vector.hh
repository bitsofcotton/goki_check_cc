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

#if !defined(_OBJ2VECTOR_)

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "tilt.hh"

using std::cerr;
using std::flush;
using std::endl;
using std::string;
using std::ifstream;
using std::ofstream;
using std::getline;
using std::stringstream;
using std::vector;
using std::binary_search;
using std::abs;
using std::pair;
using std::make_pair;
using std::unique;

template <typename T> void maskVectors(vector<Eigen::Matrix<T, 3, 1> >& points, vector<Eigen::Matrix<int, 3, 1> >& polys, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& data) {
  std::vector<int> elim, elimp, after;
  for(int i = 0, ii = 0; i < points.size(); i ++) {
    const int y(std::max(std::min(int(points[i][0]), int(data.rows() - 1)), 0));
    const int x(std::max(std::min(int(points[i][1]), int(data.cols() - 1)), 0));
    if(data(y, x) > .5) {
      elim.push_back(i);
      after.push_back(- 1);
    } else
      after.push_back(ii ++);
  }
  for(int i = 0; i < polys.size(); i ++)
    if(std::binary_search(elim.begin(), elim.end(), polys[i][0]) ||
       std::binary_search(elim.begin(), elim.end(), polys[i][1]) ||
       std::binary_search(elim.begin(), elim.end(), polys[i][2]))
      elimp.push_back(i);
  for(int i = 0, j = 0; i < elim.size(); i ++) {
    points.erase(points.begin() + (elim[i] - j), points.begin() + (elim[i] - j + 1));
    j ++;
  }
  for(int i = 0; i < polys.size(); i ++)
    for(int j = 0; j < polys[i].size(); j ++)
      polys[i][j] = after[polys[i][j]];
  for(int i = 0, j = 0; i < elimp.size(); i ++) {
    polys.erase(polys.begin() + (elimp[i] - j), polys.begin() + (elimp[i] - j + 1));
    j ++;
  }
  return;
}

template <typename T> void floodfill(Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& checked, vector<pair<int, int> >& store, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mask, const int& y, const int& x) {
  assert(mask.rows() == checked.rows() && mask.cols() == checked.cols());
  vector<pair<int, int> > tries;
  tries.push_back(make_pair(+ 1,   0));
  tries.push_back(make_pair(  0, + 1));
  tries.push_back(make_pair(- 1,   0));
  tries.push_back(make_pair(  0, - 1));
  vector<pair<int, int> > stack;
  stack.push_back(make_pair(y, x));
  while(stack.size()) {
    const auto pop(stack[stack.size() - 1]);
    stack.pop_back();
    const int& yy(pop.first);
    const int& xx(pop.second);
    if(! (0 <= yy && yy < checked.rows() && 0 <= xx && xx < checked.cols()) )
      store.push_back(pop);
    else if(!checked(yy, xx)) {
      checked(yy, xx) = true;
      if(T(.5) < mask(yy, xx))
        store.push_back(pop);
      else
        for(int i = 0; i < tries.size(); i ++)
          stack.push_back(make_pair(yy + tries[i].first, xx + tries[i].second));
    }
  }
  return;
}

// N.B. mask 2 vector but this takes almost bruteforce.
// XXX there must exists good libraries, but I cannot find with License and Eigen STL pair.
// XXX but this is terrible slow.
template <typename T> vector<vector<int> > getEdges(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mask, const vector<Eigen::Matrix<T, 3, 1> >& points, const int& vbox) {
  cerr << "getEdges" << flush;
  Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> checked(mask.rows(), mask.cols());
  for(int i = 0; i < checked.rows(); i ++)
    for(int j = 0; j < checked.cols(); j ++)
      checked(i, j) = false;
  vector<pair<int, int> > store;
  for(int i = 0; i < checked.rows(); i ++)
    for(int j = 0; j < checked.cols(); j ++)
      if(mask(i, j) < T(.5))
        floodfill(checked, store, mask, i, j);
  sort(store.begin(), store.end());
  store.erase(unique(store.begin(), store.end()), store.end());
  cerr << " with " << store.size() << " edge points " << flush;
  
  vector<vector<int> > result;
  // stored indices.
  vector<int>          se;
  for( ; se.size() < store.size(); ) {
    // tree index.
    vector<int> e;
    int         i(0);
    for( ; i < store.size(); i ++)
      if(!binary_search(se.begin(), se.end(), i))
        break;
    if(store.size() <= i)
      break;
    // get edge point tree.
    for( ; i < store.size(); ) {
      // get nearest points.
      vector<int> si;
      for(int j = 0; j < store.size(); j ++)
        if(!binary_search(se.begin(), se.end(), j) &&
           abs(store[i].first  - store[j].first)  <= 1 &&
           abs(store[i].second - store[j].second) <= 1)
          si.push_back(j);
      if(!si.size())
        break;
      // normal vector direction.
      sort(si.begin(), si.end());
      int j(0);
      for( ; j < si.size(); j ++) {
        const auto& sti(store[i]);
        const auto& stj(store[si[j]]);
        const T y(stj.first  - sti.first);
        const T x(stj.second - stj.second);
        const T x2(x - y + sti.second);
        const T y2(x + y + sti.first);
        if(0 <= x2 && x2 < mask.cols() && 0 <= y2 && y2 < mask.rows() &&
           mask(y2, x2) < T(.5))
          break;
      }
      if(si.size() <= j)
        j = 0;
      // store.
      i = si[j];
      e.push_back(i);
      se.push_back(i);
      sort(se.begin(), se.end());
    }
    if(1 < e.size()) {
      result.push_back(vector<int>());
      // apply to points index.
      vector<int> pj;
      for(int i = 0; i < e.size(); i ++) {
        const auto& s(store[e[i]]);
        vector<pair<T, int> > distances;
        for(int j = 0; j < points.size(); j ++)
          distances.push_back(make_pair(
                                pow(T(s.first  - points[j][0]), T(2)) +
                                  pow(T(s.second - points[j][1]), T(2)),
                                j));
        sort(distances.begin(), distances.end());
        if(distances.size() && !binary_search(pj.begin(), pj.end(), distances[0].second)) {
          result[result.size() - 1].push_back(distances[0].second);
          pj.push_back(distances[0].second);
        }
      }
      result[result.size() - 1].push_back(result[result.size() - 1][0]);
    } else
      se.push_back(i);
    cerr << "." << flush;
  }
  return result;
}

template <typename T> bool saveobj(const vector<Eigen::Matrix<T, 3, 1> >& data, const vector<Eigen::Matrix<int, 3, 1> >& polys, const char* filename, const bool& addstand = false, const vector<vector<int> >& edges = vector<vector<int> >(), const T& zs = T(2)) {
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
    if(addstand) {
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

#define _OBJ2VECTOR_
#endif

