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

template <typename T> int clockwise(const Eigen::Matrix<T, 3, 1>& p0, const Eigen::Matrix<T, 3, 1>& p1, const Eigen::Matrix<T, 3, 1>& p2, const T& epsilon = T(1e-4)) {
  Eigen::Matrix<T, 3, 3> dc0;
  dc0(0, 0) = T(1);
  dc0(0, 1) = p0[0];
  dc0(0, 2) = p0[1];
  dc0(1, 0) = T(1);
  dc0(1, 1) = p1[0];
  dc0(1, 2) = p1[1];
  dc0(2, 0) = T(1);
  dc0(2, 1) = p2[0];
  dc0(2, 2) = p2[1];
  const T det(dc0.determinant());
  if(abs(det) <= pow(epsilon, T(3)))
    return 0;
  return det < T(0) ? - 1 : 1;
}

template <typename T> bool isCrossing(const Eigen::Matrix<T, 3, 1>& p0, const Eigen::Matrix<T, 3, 1>& p1, const Eigen::Matrix<T, 3, 1>& q0, const Eigen::Matrix<T, 3, 1>& q1, const T& err = T(1e-4)) {
  // t * p0 + (1 - t) * p1 == s * q0 + (1 - s) * q1
  // <=> p1 + (p0 - p1) t == q1 + (q0 - q1) s
  // <=> [(p0 - p1), (q1 - q0)][t, s] == q1 - p1.
  // <=> Ax==b.
  Eigen::Matrix<T, 2, 2> A;
  Eigen::Matrix<T, 2, 1> b;
  A(0, 0) = p0[0] - p1[0];
  A(1, 0) = p0[1] - p1[1];
  A(0, 1) = q1[0] - q0[0];
  A(1, 1) = q1[1] - q0[1];
  b[0]    = q1[0] - p1[0];
  b[1]    = q1[1] - p1[1];
  if(abs(A.determinant()) <= pow(err, T(2)))
    return false;
  auto x(A.inverse() * b);
  return (err <= x[0] && x[0] <= T(1) - err) &&
         (err <= x[1] && x[1] <= T(1) - err);
}

template <typename T> bool isSameLine2(const Eigen::Matrix<T, 3, 1>& a, const Eigen::Matrix<T, 3, 1>& b, const Eigen::Matrix<T, 3, 1>& c, const T& rerr = T(1e-6)) {
  Eigen::Matrix<T, 3, 1> bcn(b - c);
  bcn[2] = T(0);
  Eigen::Matrix<T, 3, 1> err(a);
  err[2] = T(0);
  err   -= err.dot(bcn) * bcn / bcn.dot(bcn);
  return err.dot(err) <= rerr;
}

// N.B. delaunay trianglation algorithm best works but this is brute force.
// XXX: terrible inefficient temporary stub.
template <typename T> vector<Eigen::Matrix<int, 3, 1> > loadBumpSimpleMesh(const vector<Eigen::Matrix<T, 3, 1> >& dst, const vector<int>& dstpoints, const T epsilon = T(1e-5)) {
  vector<Eigen::Matrix<int, 3, 1> > res;
  tilter<T> tilt;
  cerr << "Delaunay(" << dstpoints.size() << ")";
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < dstpoints.size(); i ++) {
    cerr << "." << flush;
    for(int j = i + 1; j < dstpoints.size(); j ++)
      for(int k = j + 1; k < dstpoints.size(); k ++) {
        if(isSameLine2<T>(dst[dstpoints[i]], dst[dstpoints[j]], dst[dstpoints[k]]))
          continue;
        Eigen::Matrix<T, 4, 4> dc;
        auto g((dst[dstpoints[i]] +
                dst[dstpoints[j]] +
                dst[dstpoints[k]]) / T(3));
        int idxs[3];
        idxs[0] = i;
        idxs[1] = j;
        idxs[2] = k;
        for(int l = 0; l < 3; l ++) {
          dc(l, 0) = T(1);
          dc(l, 1) = dst[dstpoints[idxs[l]]][0] - g[0];
          dc(l, 2) = dst[dstpoints[idxs[l]]][1] - g[1];
          dc(l, 3) = dc(l, 1) * dc(l, 1) + dc(l, 2) * dc(l, 2);
        }
        dc(3, 0) = T(1);
        int cw(clockwise<T>(dst[dstpoints[i]] - g,
                            dst[dstpoints[j]] - g,
                            dst[dstpoints[k]] - g, epsilon));
        bool flag(true);
        for(int l = 0; l < dstpoints.size(); l ++) {
          dc(3, 1) = dst[dstpoints[l]][0] - g[0];
          dc(3, 2) = dst[dstpoints[l]][1] - g[1];
          dc(3, 3) = dc(3, 1) * dc(3, 1) + dc(3, 2) * dc(3, 2);
          // if one of them not meets delaunay or inside triangle, break;
          if(cw * dc.determinant() < - epsilon ||
             (tilt.sameSide2(dst[dstpoints[i]],
                             dst[dstpoints[j]],
                             dst[dstpoints[k]],
                             dst[dstpoints[l]], false) &&
              tilt.sameSide2(dst[dstpoints[j]],
                             dst[dstpoints[k]],
                             dst[dstpoints[i]],
                             dst[dstpoints[l]], false) &&
              tilt.sameSide2(dst[dstpoints[k]],
                             dst[dstpoints[i]],
                             dst[dstpoints[j]],
                             dst[dstpoints[l]], false)) ) {
            flag = false;
            break;
          }
        }
        if(flag) {
          Eigen::Matrix<int, 3, 1> idx;
          idx[0] = i;
          if(cw < 0) {
            idx[1] = k;
            idx[2] = j;
          } else {
            idx[1] = j;
            idx[2] = k;
          }
#if defined(_OPENMP)
#pragma omp critical
#endif
          {
            for(int i = 0; i < res.size(); i ++)
              for(int j = 0; j < res[i].size(); j ++)
                for(int k = 0; k < idx.size(); k ++)
                  if(isCrossing(dst[dstpoints[res[i][(j + 0) % 3]]],
                                dst[dstpoints[res[i][(j + 1) % 3]]],
                                dst[dstpoints[idx[(k + 0) % 3]]],
                                dst[dstpoints[idx[(k + 1) % 3]]]))
                    goto fixnext;
            res.push_back(idx);
           fixnext:
            ;
          }
        }
      }
  }
  return res;
}

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
      else
        checked(i, j) = true;
  if(!store.size()) {
    for(int i = 0; i < checked.rows(); i ++) {
      store.push_back(make_pair(i, 0));
      store.push_back(make_pair(i, checked.cols() - 1));
    }
    for(int i = 0; i < checked.cols(); i ++) {
      store.push_back(make_pair(0,                  i));
      store.push_back(make_pair(checked.rows() - 1, i));
    }
  }
  sort(store.begin(), store.end());
  store.erase(unique(store.begin(), store.end()), store.end());
  cerr << " with " << store.size() << " edge points " << flush;
  
  vector<vector<int> > result;
  // stored indices.
  vector<int>          se;
  for( ; se.size() < store.size(); ) {
    // tree index.
    vector<int> e;
    int i(0);
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
           abs(store[i].second - store[j].second) <= 1 &&
           !(store[i].first  == store[j].first &&
             store[i].second == store[j].second))
          si.push_back(j);
      if(!si.size())
        break;
      // normal vector direction.
      sort(si.begin(), si.end());
      int j(0);
      for( ; j < si.size(); j ++) {
        const auto& sti(store[i]);
        const auto& stj(store[si[j]]);
        if(0 < abs(sti.first - stj.first)) {
          if(0 < abs(sti.second - stj.second)) {
            // 45 degree.
            if((sti.first  - stj.first) *
               (sti.second - stj.second) > 0) {
              if(0 <= sti.first  && sti.first  < mask.rows() &&
                 0 <= sti.second && sti.second < mask.cols() &&
                 mask(sti.first, stj.second) < T(.5))
                break;
            } else if(0 <= sti.first  && sti.first  < mask.rows() &&
                      0 <= sti.second && sti.second < mask.cols() &&
                      mask(stj.first, sti.second) < T(.5))
              break;
          } else if(0 <= stj.first && stj.first < mask.rows() &&
                    0 <= sti.second + (stj.first - sti.first) &&
                         sti.second + (stj.first - sti.first) < mask.cols() &&
                    mask(stj.first,
                         sti.second + (stj.first - sti.first)) < T(.5))
            break;
          // it is guaranteed 0 < abs(sti.second - stj.second) from addition.
        } else if(0 <= stj.second && stj.second < mask.cols() &&
                  0 <= sti.first - (stj.second - sti.second) &&
                       sti.first - (stj.second - sti.second) < mask.rows() &&
                  mask(sti.first - (stj.second - sti.second),
                       stj.second) < T(.5))
          break;
      }
      if(si.size() <= j)
        j = 0;
      // store.
      e.push_back(si[j]);
      se.push_back(si[j]);
      i = si[j];
      sort(se.begin(), se.end());
    }
    if(1 < e.size()) {
      result.push_back(vector<int>());
      // apply to points index.
      for(int i = 0; i < e.size(); i ++) {
        const auto& s(store[e[i]]);
        vector<pair<T, int> > distances;
        for(int j = 0; j < points.size(); j ++)
          distances.push_back(make_pair(
                                sqrt(pow(T(s.first  - points[j][0]), T(2)) +
                                     pow(T(s.second - points[j][1]), T(2))),
                                j));
        sort(distances.begin(), distances.end());
        result[result.size() - 1].push_back(distances[0].second);
      }
      const auto& head(points[result[result.size() - 1][0]]);
      const auto& tail(points[result[result.size() - 1][result[result.size() - 1].size() - 1]]);
      const auto  delta(head - tail);
      if(sqrt(delta[0] * delta[0] + delta[1] * delta[1]) <= T(2))
        result.push_back(result[0]);
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
          output << " "  << data.size() + outer[i1] + 1;
          output << " "  << outer[i1] + 1 << endl;
          output << "f " << data.size() + outer[i0] + 1;
          output << " "  << outer[i1] + 1;
          output << " "  << outer[i0] + 1 << endl;
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
      for(int j = 0; j < polys[i].size(); j ++)
        polys[i][j] = abs(int((polys[i][j] + data.size()) % (2 * data.size())) - int(data.size()));
    input.close();
  } else {
    cerr << "Unable to open file: " << filename << endl;
    return false;
  }
  return true;
}

#define _OBJ2VECTOR_
#endif

