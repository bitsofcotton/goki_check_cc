#if !defined(_OBJ2VECTOR_)

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "tilt.hh"

using std::cerr;
using std::cout;
using std::endl;
using std::isfinite;
using std::string;
using std::ifstream;
using std::ofstream;
using std::getline;
using std::istringstream;
using std::stringstream;
using std::vector;
using std::floor;
using std::min;
using std::max;
using std::abs;

template <typename T> int clockwise(const Eigen::Matrix<T, 3, 1>& p0, const Eigen::Matrix<T, 3, 1>& p1, const Eigen::Matrix<T, 3, 1>& p2, const T& epsilon = T(0)) {
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
  if(abs(det) <= epsilon)
    return 0;
  return det < T(0) ? - 1 : 1;
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
    cerr << ".";
    fflush(stderr);
    for(int j = i + 1; j < dstpoints.size(); j ++)
      for(int k = j + 1; k < dstpoints.size(); k ++) {
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
            res.push_back(idx);
          }
        }
      }
  }
  return res;
}

template <typename T> bool saveobj(const vector<Eigen::Matrix<T, 3, 1> >& data, const vector<Eigen::Matrix<int, 3, 1> >& polys, const char* filename, const T& zr = T(20)) {
  ofstream output;
  output.open(filename, std::ios::out);
  if(output.is_open()) {
    for(int i = 0; i < data.size(); i ++)
      output << "v " << data[i][1] << " " << - data[i][0] << " " << data[i][2] * zr << endl;
    // xchg with clockwise/counter clockwise.
    for(int i = 0; i < polys.size(); i ++)
      output << "f " << polys[i][0] + 1 << " " << polys[i][1] + 1 << " " << polys[i][2] + 1 << endl;
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
        sub >> buf[0];
        sub >> buf[1];
        sub >> buf[2];
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
          if(sub.eof() || sub.bad())
            break;
          if(widx > 2)
            flag = true;
          if(flag)
            polys.push_back(wbuf);
          widx %= 3;
          sub.ignore(20, ' ');
        }
      }
    }
    for(int i = 0; i < polys.size(); i ++)
      for(int j = 0; j < polys[i].size(); j ++)
        polys[i][j] = abs(int(polys[i][j] % (2 * data.size())) - int(data.size()));
    input.close();
  } else {
    cerr << "Unable to open file: " << filename << endl;
    return false;
  }
  return true;
}

#define _OBJ2VECTOR_
#endif

