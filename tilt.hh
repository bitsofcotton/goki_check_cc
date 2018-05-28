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

#if !defined(_TILT_)

#include <Eigen/Core>
#include <Eigen/LU>
#include <vector>
#include <cmath>

using namespace Eigen;
using std::min;
using std::max;
using std::ceil;
using std::sqrt;
using std::cerr;
using std::endl;
using std::flush;
using std::vector;
using std::abs;
using std::isfinite;

template <typename T> class triangles_t {
public:
  Eigen::Matrix<T, 3, 3> p;
  Eigen::Matrix<T, 3, 1> n;
  T                      c;
  T                      z;
  triangles_t<T>& rotate(const Eigen::Matrix<T, 3, 3>& R, const Eigen::Matrix<T, 3, 1>& origin) {
    for(int i = 0; i < 3; i ++)
      p.col(i) = R * (p.col(i) - origin) + origin;
    return *this;
  }
  triangles_t<T>& solveN() {
    const auto pq(p.col(1) - p.col(0));
    const auto pr(p.col(2) - p.col(0));
    n[0] =   (pq[1] * pr[2] - pq[2] * pr[1]);
    n[1] = - (pq[0] * pr[2] - pq[2] * pr[0]);
    n[2] =   (pq[0] * pr[1] - pq[1] * pr[0]);
    if(n.dot(n) > 0)
      n /= sqrt(n.dot(n));
    z = n.dot(p.col(0));
    return *this;
  }
};

template <typename T> class match_t;
template <typename T> class tilter {
public:
  typedef Matrix<T, Dynamic, Dynamic> Mat;
  typedef Matrix<T, 3, 3>             Mat3x3;
  typedef Matrix<T, 3, 1>             Vec3;
  typedef Matrix<T, 2, 1>             Vec2;
  typedef triangles_t<T>              Triangles;
  
  tilter();
  ~tilter();
  void initialize(const T& z_ratio);
  
  Mat  tilt(const Mat& in, const Mat& bump, const int& idx, const int& samples, const T& psi);
  Mat  tilt(const Mat& in, const Mat& bump, const match_t<T>& m);
  Mat  tilt(const Mat& in, const vector<Triangles>& triangles0, const match_t<T>& m);
  Mat  tilt(const Mat& in, const vector<Triangles>& triangles);
  Triangles makeTriangle(const int& u, const int& v, const Mat& in, const Mat& bump, const int& flg);
  bool sameSide2(const Vec2& p0, const Vec2& p1, const Vec2& p, const Vec2& q, const bool& extend = true, const T& err = T(1e-5));
  bool sameSide2(const Vec3& p0, const Vec3& p1, const Vec3& p, const Vec3& q, const bool& extend = true, const T& err = T(1e-5));
  
private:
  T    sgn(const T& x);
  bool onTriangle(T& z, const Triangles& tri, const Vec2& geom);
  T    Pi;
  T    z_ratio;
  T    thresh;
};

template <typename T> tilter<T>::tilter() {
  initialize(1.);
  return;
}

template <typename T> tilter<T>::~tilter() {
  ;
}

template <typename T> void tilter<T>::initialize(const T& z_ratio) {
  this->z_ratio = z_ratio;
  Pi            = atan2(T(1), T(1)) * T(4);
  thresh        = 1e-8;
  return;
}

template <typename T> T tilter<T>::sgn(const T& x) {
  if(x < T(0))
    return - T(1);
  if(x > T(0))
    return   T(1);
  return T(0);
}

template <typename T> triangles_t<T> tilter<T>::makeTriangle(const int& u, const int& v, const Mat& in, const Mat& bump, const int& flg) {
  Triangles work;
  if(flg) {
    work.p(0, 0) = u;
    work.p(1, 0) = v;
    work.p(0, 1) = u + 1;
    work.p(1, 1) = v;
    work.p(0, 2) = u + 1;
    work.p(1, 2) = v + 1;
  } else {
    work.p(0, 0) = u;
    work.p(1, 0) = v;
    work.p(0, 1) = u;
    work.p(1, 1) = v + 1;
    work.p(0, 2) = u + 1;
    work.p(1, 2) = v + 1;
  }
  work.c = T(0);
  for(int i = 0; i < 3;  i ++) {
    work.p(2, i) = bump(int(work.p(0, i)), int(work.p(1, i))) * z_ratio;
    work.c      += in(int(work.p(0, i)), int(work.p(1, i)));
  }
  work.c /= T(3);
  return work.solveN();
}

template <typename T> bool tilter<T>::sameSide2(const Vec2& p0, const Vec2& p1, const Vec2& p2, const Vec2& q, const bool& extend, const T& err) {
  const Vec2 dlt(p1 - p0);
  Vec2 dp(p2 - p0);
  dp -= dlt.dot(dp) * dlt / dlt.dot(dlt);
  // N.B. dp.dot(p1 - p0) >> 0.
  return dp.dot(q - p0) >= (extend ? - T(1) : T(1)) * (abs(dp[0]) + abs(dp[1])) * err;
}

template <typename T> bool tilter<T>::sameSide2(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& q, const bool& extend, const T& err) {
  return sameSide2(Vec2(p0[0], p0[1]), Vec2(p1[0], p1[1]),
                   Vec2(p2[0], p2[1]), Vec2(q[0],  q[1]),
                   extend, err);
}

// <[x, y, t], triangle.n> == triangle.z
template <typename T> bool tilter<T>::onTriangle(T& z, const Triangles& tri, const Vec2& geom) {
  Vec3 v0;
  Vec3 camera;
  v0[0] = 0;
  v0[1] = 0;
  v0[2] = 1;
  camera[0] = geom[0];
  camera[1] = geom[1];
  camera[2] = 0;
  // <v0 t + camera, v4> = tri.
  const T t((tri.z - tri.n.dot(camera)) / (tri.n.dot(v0)));
  z = camera[2] + v0[2] * t;
  Vec3 geom3;
  geom3[0] = geom[0];
  geom3[1] = geom[1];
  geom3[2] = T(0);
  return (sameSide2(tri.p.col(0), tri.p.col(1), tri.p.col(2), geom3, true, T(.125)) &&
          sameSide2(tri.p.col(1), tri.p.col(2), tri.p.col(0), geom3, true, T(.125)) &&
          sameSide2(tri.p.col(2), tri.p.col(0), tri.p.col(1), geom3, true, T(.125)));
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tilter<T>::tilt(const Mat& in, const Mat& bump, const int& idx, const int& samples, const T& psi) {
  const T theta(2. * Pi * idx / samples);
  const T lpsi(Pi * psi);
  Mat3x3 R0;
  Mat3x3 R1;
  R0(0, 0) =   cos(theta);
  R0(0, 1) = - sin(theta);
  R0(0, 2) = 0.;
  R0(1, 0) =   sin(theta);
  R0(1, 1) =   cos(theta);
  R0(1, 2) = 0.;
  R0(2, 0) = 0.;
  R0(2, 1) = 0.;
  R0(2, 2) = 1.;
  R1(0, 0) = 1.;
  R1(0, 1) = 0.;
  R1(0, 2) = 0.;
  R1(1, 0) = 0.;
  R1(1, 1) =   cos(lpsi);
  R1(1, 2) = - sin(lpsi);
  R1(2, 0) = 0.;
  R1(2, 1) =   sin(lpsi);
  R1(2, 2) =   cos(lpsi);
  Vec3 pcenter;
  pcenter[0] = T(in.rows() - 1.) / 2;
  pcenter[1] = T(in.cols() - 1.) / 2;
  pcenter[2] = T(.5);
  match_t<T> m;
  m.rot    = R0.transpose() * R1 * R0;
  m.offset = pcenter - m.rot * pcenter;
  m.ratio  = T(1);
  return tilt(in, bump, m);
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tilter<T>::tilt(const Mat& in, const Mat& bump, const match_t<T>& m) {
  assert(in.rows() == bump.rows() && in.cols() == bump.cols());
  vector<Triangles> triangles;
  triangles.reserve((in.rows() - 1) * (in.cols() - 1) * 2);
  for(int i = 0; i < in.rows() - 1; i ++)
    for(int j = 0; j < in.cols() - 1; j ++) {
      triangles.push_back(makeTriangle(i, j, in, bump, false));
      triangles.push_back(makeTriangle(i, j, in, bump, true));
    }
  return tilt(in, triangles, m);
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tilter<T>::tilt(const Mat& in, const vector<Triangles>& triangles0, const match_t<T>& m) {
  vector<Triangles> triangles(triangles0);
  for(int j = 0; j < triangles.size(); j ++) {
    for(int k = 0; k < 3; k ++)
      triangles[j].p.col(k) = m.transform(triangles[j].p.col(k));
    triangles[j].solveN();
  }
  return tilt(in, triangles);
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tilter<T>::tilt(const Mat& in, const vector<Triangles>& triangles) {
  Mat result(in.rows(), in.cols());
  for(int i = 0; i < in.rows(); i ++)
    for(int j = 0; j < in.cols(); j ++)
      result(i, j) = 0.;
  cerr << "t" << flush;
  Mat zb(in.rows(), in.cols());
  for(int j = 0; j < zb.rows(); j ++)
    for(int k = 0; k < zb.cols(); k ++)
      zb(j, k) = - T(1e8);
  // able to boost with divide and conquer.
  for(int j = 0; j < triangles.size(); j ++) {
    const Triangles& tri(triangles[j]);
    Vec2 gs[3];
    for(int k = 0; k < 3; k ++) {
      gs[k][0] = tri.p(0, k);
      gs[k][1] = tri.p(1, k);
    }
    int ll = int( min(min(gs[0][0], gs[1][0]), gs[2][0]));
    int rr = ceil(max(max(gs[0][0], gs[1][0]), gs[2][0])) + 1;
    int bb = int( min(min(gs[0][1], gs[1][1]), gs[2][1]));
    int tt = ceil(max(max(gs[0][1], gs[1][1]), gs[2][1])) + 1;
    for(int y = max(0, ll); y < min(rr, int(in.rows())); y ++)
      for(int x = max(0, bb); x < min(tt, int(in.cols())); x ++) {
        T z;
        Vec2 midgeom;
        midgeom[0] = y;
        midgeom[1] = x;
        if(onTriangle(z, tri, midgeom) && isfinite(z) && zb(y, x) < z) {
          result(y, x) = tri.c;
          zb(y, x)     = z;
        }
      }
  }
  return result;
}

#define _TILT_
#endif

