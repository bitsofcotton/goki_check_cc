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
using std::fflush;
using std::vector;

template <typename T> class tilter {
public:
  typedef Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Matrix<T, 3, 3>                           Mat3x3;
  typedef Matrix<T, Eigen::Dynamic, 1>              Vec;
  typedef Matrix<T, 3, 1>                           Vec3;
  typedef Matrix<T, 2, 1>                           Vec2;
  typedef Matrix<T, 3, 5>                           Triangles;
  
  tilter();
  ~tilter();
  void initialize(const T& z_ratio);
  
  Mat tilt(const Mat& in, const Mat& bump, const int& idx, const int& samples, const T& psi);
  Mat tilt(const Mat& in, const Mat& bump, const Mat3x3& rot, const Mat3x3& rotrev, const Vec3& moveto, const T& rto, const Vec3& origin0);
  Vec3 solveN(const Vec3& p, const Vec3& q, const Vec3& r);
  Eigen::Matrix<T, 3, 5> makeTriangle(const int& u, const int& v, const Mat& in, const Mat& bump, const int& flg);
  Mat tiltsub(const Mat& in, const vector<Triangles>& triangles, const Mat3x3& rot, const Mat3x3& rotrev, const Vec3& origin, const Vec3& moveto, const T& rto);
  bool sameSide2(const Vec2& p0, const Vec2& p1, const Vec2& p, const Vec2& q);
private:
  T    sgn(const T& x);
  Vec3 rotate0(const Vec3& work, const Mat3x3& rot, const Vec3& origin);
  Eigen::Matrix<T, 3, 5> rotate(const Eigen::Matrix<T, 3, 5>& triangle, const Mat3x3& rot, const Vec3& origin, const T& rr);
  bool onTriangle(T& z, const Triangles& tri, const Vec2& geom);
  bool scale(Mat3x3& A, Vec3& b, const Mat3x3& vorig, const Mat3x3& vto);
  T   Pi;
  T   z_ratio;
  T   thresh;
};

template <typename T> tilter<T>::tilter() {
  initialize(.75);
  return;
}

template <typename T> tilter<T>::~tilter() {
  ;
}

template <typename T> void tilter<T>::initialize(const T& z_ratio) {
  this->z_ratio = z_ratio;
  Pi            = atan2(T(1.), T(1.)) * T(4.);
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

// normal vector to plane(p, q, r).
template <typename T> Eigen::Matrix<T, 3, 1> tilter<T>::solveN(const Vec3& p, const Vec3& q, const Vec3& r) {
  const Vec3 pq(q - p);
  const Vec3 pr(r - p);
  const Vec3 prorth(pr - pr.dot(pq) * pq / pq.dot(pq));
  Vec3 nn;
  for(int i = 0; i < nn.size(); i ++) {
    for(int j = 0; j < nn.size(); j ++)
      nn[j] = (i == j ? T(1) : T(0));
    nn -= pq.dot(nn) * pq / pq.dot(pq) + prorth.dot(nn) * prorth / prorth.dot(prorth);
    if(nn.dot(nn) > 0)
      break;
  }
  if(nn.dot(nn) == 0) {
    cerr << "internal error." << endl;
    return nn;
  }
  return nn / sqrt(nn.dot(nn));
}

template <typename T> Eigen::Matrix<T, 3, 5> tilter<T>::makeTriangle(const int& u, const int& v, const Mat& in, const Mat& bump, const int& flg) {
  Eigen::Matrix<T, 3, 5> work;
  if(flg) {
    work(0, 0) = u;
    work(1, 0) = v;
    work(0, 1) = u + 1;
    work(1, 1) = v;
    work(0, 2) = u + 1;
    work(1, 2) = v + 1;
  } else {
    work(0, 0) = u;
    work(1, 0) = v;
    work(0, 1) = u;
    work(1, 1) = v + 1;
    work(0, 2) = u + 1;
    work(1, 2) = v + 1;
  }
  for(int i = 0; i < 3; i ++)
    work(i, 3) = 0.;
  for(int i = 0; i < 3;  i ++) {
    work(2, i)  = bump(int(work(0, i)), int(work(1, i))) * sqrt(T(in.rows() * in.cols())) * z_ratio;
    work(0, 3) += in(int(work(0, i)), int(work(1, i)));
  }
  work(0, 3) /= 3.;
  work.col(4) = solveN(work.col(0), work.col(1), work.col(2));
  work(1, 3)  = work.col(4).dot(work.col(0));
  return work;
}

template <typename T> Eigen::Matrix<T, 3, 1> tilter<T>::rotate0(const Vec3& work, const Mat3x3& rot, const Vec3& origin) {
  return rot * (work - origin) + origin;
}

template <typename T> bool tilter<T>::sameSide2(const Vec2& p0, const Vec2& p1, const Vec2& p2, const Vec2& q) {
  const Vec2 dlt(p1 - p0);
  Vec2 dp(p2 - p0);
  dp -= dlt.dot(dp) * dlt / dlt.dot(dlt);
  // XXX refixed. dp.dot(p2 - p0) >> 0.
  // XXX refixed. another place uses.
  return dp.dot(q - p0) + max(abs(dp[0]) + abs(dp[1]), T(0)) / T(2) / max(max(sqrt((p0 - p1).dot(p0 - p1)), sqrt((p1 - p2).dot(p1 - p2))), sqrt((p2 - p0).dot(p2 - p0))) >= T(0);
}

template <typename T> Eigen::Matrix<T, 3, 5> tilter<T>::rotate(const Eigen::Matrix<T, 3, 5>& triangle, const Mat3x3& rot, const Vec3& origin, const T& rr) {
  Triangles res;
  Vec3      zero;
  for(int i = 0; i < 3; i ++) {
    res.col(i)  = rotate0(triangle.col(i), rot, origin);
    res.col(i) *= rr;
    zero[i]     = 0.;
  }
  res.col(3) = triangle.col(3);
  res.col(4) = rotate0(triangle.col(4), rot, zero);
  res(1, 3)  = res.col(4).dot(res.col(0));
  return res;
}

// <[x, y, t], triangle.n> == triangle.k
template <typename T> bool tilter<T>::onTriangle(T& z, const Eigen::Matrix<T, 3, 5>& tri, const Vec2& geom) {
  Vec3 v0;
  Vec3 camera;
  v0[0] = 0;
  v0[1] = 0;
  v0[2] = 1;
  camera[0] = geom[0];
  camera[1] = geom[1];
  camera[2] = 0;
  // <v0 t + camera, v4> = tri.
  const T t((tri(1, 3) - tri.col(4).dot(camera)) / (tri.col(4).dot(v0)));
  z = camera[2] + v0[2] * t;
  Eigen::Matrix<T, 2, 3> tritri;
  for(int i = 0; i < 3; i ++) {
    tritri(0, i) = tri(0, i);
    tritri(1, i) = tri(1, i);
  }
  // XXX: don't now why, but logic is inverted in calculation.
  return !(sameSide2(tritri.col(0), tritri.col(1), tritri.col(2), geom) &&
           sameSide2(tritri.col(1), tritri.col(2), tritri.col(0), geom) &&
           sameSide2(tritri.col(2), tritri.col(0), tritri.col(1), geom));
}

template <typename T> bool tilter<T>::scale(Mat3x3& A, Vec3& b, const Mat3x3& vorig, const Mat3x3& vto) {
  b = vto.col(0) - vorig.col(0);
  Mat3x3 vdorig, vdto;
  vdorig.col(0) = vorig.col(1) - vorig.col(0);
  vdorig.col(1) = vorig.col(2) - vorig.col(0);
  vdorig(0, 2)  = 0.;
  vdorig(1, 2)  = 0.;
  vdorig(2, 2)  = 1.;
  vdto.col(0)   = vto.col(1)   - vto.col(0);
  vdto.col(1)   = vto.col(2)   - vto.col(0);
  vdto(0, 2)    = 0.;
  vdto(1, 2)    = 0.;
  vdto(2, 2)    = 1.;
  A = vdto * vdorig.inverse();
  return true;
} 

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tilter<T>::tilt(const Mat& in, const Mat& bump, const int& idx, const int& samples, const T& psi) {
  const T theta(2. * Pi * idx / samples);
  const T lpsi(Pi / 2. - psi * Pi / 2.);
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
  Vec3 zero3;
  zero3[0]   = T(0);
  zero3[1]   = T(0);
  zero3[2]   = T(0);
  Vec3 pcenter;
  pcenter[0] = in.rows() / 2.;
  pcenter[1] = in.cols() / 2.;
  pcenter[2] = 0.;
  return tilt(in, bump, R1 * R0, R0.transpose(), zero3, T(1), pcenter);
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tilter<T>::tilt(const Mat& in, const Mat& bump, const Mat3x3& rot, const Mat3x3& rotrev, const Vec3& moveto, const T& rto, const Vec3& origin0) {
  if(in.rows() != bump.rows() || in.cols() != bump.cols()) {
    cerr << "tilt: size mismatch..." << endl;
    return Mat();
  }
  cerr << " making triangles";
  fflush(stderr);
  vector<Triangles> triangles;
  for(int i = 0; i < in.rows() - 1; i ++)
    for(int j = 0; j < in.cols() - 1; j ++) {
      triangles.push_back(makeTriangle(i, j, in, bump, false));
      triangles.push_back(makeTriangle(i, j, in, bump, true));
    }
  return tiltsub(in, triangles, rot, rotrev, origin0, moveto, rto);
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tilter<T>::tiltsub(const Mat& in, const vector<Triangles>& triangles, const Mat3x3& rot, const Mat3x3& rotrev, const Vec3& origin, const Vec3& moveto, const T& rto) {
  Mat result(in.rows(), in.cols());
  for(int i = 0; i < in.rows(); i ++)
    for(int j = 0; j < in.cols(); j ++)
      result(i, j) = 0.;
  vector<Triangles> rotriangles;
  for(int j = 0; j < triangles.size(); j ++) {
    rotriangles.push_back(rotate(rotate(triangles[j], rot, origin, rto), rotrev, origin, T(1)));
    rotriangles[j].col(0) += moveto;
    rotriangles[j].col(1) += moveto;
    rotriangles[j].col(2) += moveto;
  }
  cerr << " draw";
  fflush(stderr);
  Mat zb(in.rows(), in.cols());
  for(int j = 0; j < zb.rows(); j ++)
    for(int k = 0; k < zb.cols(); k ++)
      zb(j, k) = - T(200000);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for
#endif
  for(int j = 0; j < rotriangles.size(); j ++) {
    Triangles& tri = rotriangles[j];
    Vec2 gs[3];
    for(int k = 0; k < 3; k ++) {
      gs[k][0] = tri(0, k);
      gs[k][1] = tri(1, k);
    }
    int ll = int(min(min(gs[0][0], gs[1][0]), gs[2][0]));
    int rr = ceil(max(max(gs[0][0], gs[1][0]), gs[2][0])) + 1;
    int bb = int(min(min(gs[0][1], gs[1][1]), gs[2][1]));
    int tt = ceil(max(max(gs[0][1], gs[1][1]), gs[2][1])) + 1;
    for(int y = max(0, ll); y < min(rr, int(in.rows())); y ++)
      for(int x = max(0, bb); x < min(tt, int(in.cols())); x ++) {
        T z;
        Vec2 midgeom;
        midgeom[0] = y;
        midgeom[1] = x;
        if(!onTriangle(z, tri, midgeom)
           && zb(y, x) < z) {
            result(y, x)  = tri(0, 3);
            zb(y, x)      = z;
        }
      }
  }
  return result;
}

#define _TILT_
#endif

