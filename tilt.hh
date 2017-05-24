#if !defined(_TILT_)

#include <Eigen/Core>
#include <Eigen/LU>
#include <vector>
#include <cmath>

using namespace Eigen;

template <typename T> class tilter {
public:
  typedef Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Matrix<T, 3, 3>                           Mat3x3;
  typedef Matrix<T, 2, 2>                           Mat2x2;
  typedef Matrix<T, Eigen::Dynamic, 1>              Vec;
  typedef Matrix<T, 3, 1>                           Vec3;
  typedef Matrix<T, 2, 1>                           Vec2;
  typedef Matrix<T, 3, 5>                           Triangles;
  
  tilter();
  ~tilter();
  void initialize(const T& z_ratio, const T& psi, const int& samples);
  
  Mat tilt(const Mat& in, const Mat& bump, const int& idx);
private:
  T    sgn(const T& x);
  Vec3 solveN(const Vec3& p, const Vec3& q, const int& cnt);
  Eigen::Matrix<T, 3, 5> makeTriangle(const int& u, const int& v, const Mat& in, const Mat& bump, const int& flg);
  Vec3 rotate0(const Vec3& work, const Vec3& origin, const T& theta, const T& psi0);
  Eigen::Matrix<T, 3, 5> rotate(const Eigen::Matrix<T, 3, 5>& triangle, const Eigen::Matrix<T, 3, 1>& origin, const T& theta, const T& psi);
  Vec2 rotate2d(const Vec2& work, const Vec2& origin, const T& theta);
  bool sameSide2(const Vec2& p0, const Vec2& p1, const Vec2& p, const Vec2& q);
  bool onTriangle(T& z, const Triangles& tri, const Vec2& geom);
  bool scale(Mat3x3& A, Vec3& b, const Mat3x3& vorig, const Mat3x3& vto);
  T   Pi;
  T   z_ratio;
  T   c_ratio;
  T   psi;
  T   thresh;
  int samples;
};

template <typename T> tilter<T>::tilter() {
  initialize(.08, .99, 4);
  return;
}

template <typename T> tilter<T>::~tilter() {
  ;
}

template <typename T> void tilter<T>::initialize(const T& z_ratio, const T& psi, const int& samples) {
  this->z_ratio = z_ratio;
  this->samples = samples;
  Pi            = atan2(T(1.), T(1.)) * T(4.);
  this->psi     = psi * Pi / 2.;
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

// <p, n> = 0, <q, n> = 0.
template <typename T> Eigen::Matrix<T, 3, 1> tilter<T>::solveN(const Eigen::Matrix<T, 3, 1>& p, const Eigen::Matrix<T, 3, 1>& q, const int& cnt) {
  const T det(p[1] * q[2] - p[2] * q[1]);
  Vec3 nn;
  if(std::abs(det) <= thresh) {
    if(cnt > 2)
      for(int i = 0; i < nn.size(); i ++)
        nn[i] = 0;
    else {
      Vec3 p2;
      Vec3 q2;
      p2[0] = p[1];
      p2[1] = p[2];
      p2[2] = p[0];
      q2[0] = q[1];
      q2[1] = q[2];
      q2[2] = q[0];
      Vec3 res(solveN(p2, q2, cnt + 1));
      nn[0] = res[2];
      nn[1] = res[0];
      nn[2] = res[1];
    }
  } else {
    Vec2 b;
    b[0]  = - p[2];
    b[1]  = - q[2];
    nn[0] = (  q[2] * b[0] - p[2] * b[1]) / det;
    nn[1] = (- p[1] * b[0] + q[1] * b[1]) / det;
    nn[2] = 1.;
  }
  if(nn.dot(nn) == 0)
    return nn;
  return nn / std::sqrt(nn.dot(nn));
}

template <typename T> Eigen::Matrix<T, 3, 5> tilter<T>::makeTriangle(const int& u, const int& v, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& in, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& bump, const int& flg) {
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
    work(2, i)  = bump(int(work(0, i)), int(work(1, i))) * std::sqrt(T(in.rows() * in.cols())) * z_ratio;
    work(0, 3) += in(int(work(0, i)), int(work(1, i)));
  }
  work(0, 3) /= 3.;
  work.col(4)  = solveN(work.col(0) - work.col(1), work.col(2) - work.col(1), 0);
  for(int i = 0; i < 3; i ++)
    work(1, 3) += work.col(4).dot(work.col(i));
  work(1, 3)   /= 3.;
  return work;
}

template <typename T> Eigen::Matrix<T, 3, 1> tilter<T>::rotate0(const Eigen::Matrix<T, 3, 1>& work, const Eigen::Matrix<T, 3, 1>& origin, const T& theta, const T& psi0) {
  const T psi(Pi / 2. - psi0);
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
  R1(1, 1) =   cos(psi);
  R1(1, 2) = - sin(psi);
  R1(2, 0) = 0.;
  R1(2, 1) =   sin(psi);
  R1(2, 2) =   cos(psi);
  // glitch : rotate in 2D, then tilt.
  return R1 * (R0 * (work - origin)) + origin;
}

template <typename T> Eigen::Matrix<T, 2, 1> tilter<T>::rotate2d(const Eigen::Matrix<T, 2, 1>& work, const Eigen::Matrix<T, 2, 1>& origin, const T& theta) {
  Mat2x2 R;
  R(0, 0) =   cos(theta);
  R(0, 1) = - sin(theta);
  R(1, 0) =   sin(theta);
  R(1, 1) =   cos(theta);
  return R * (work - origin) + origin;
}

template <typename T> bool tilter<T>::sameSide2(const Eigen::Matrix<T, 2, 1>& p0, const Eigen::Matrix<T, 2, 1>& p1, const Eigen::Matrix<T, 2, 1>& p2, const Eigen::Matrix<T, 2, 1>& q) {
  const Vec2 dlt(p1 - p0);
  Vec2 dp(p2 - p0);
  dp -= dlt.dot(dp) * dlt / dlt.dot(dlt);
  return dp.dot(q) * dp.dot(p0) >= 0;
}

template <typename T> Eigen::Matrix<T, 3, 5> tilter<T>::rotate(const Eigen::Matrix<T, 3, 5>& triangle, const Eigen::Matrix<T, 3, 1>& origin, const T& theta, const T& psi) {
  Triangles res;
  Vec3      zero;
  for(int i = 0; i < 3; i ++) {
    res.col(i) = rotate0(triangle.col(i), origin, theta, psi);
    zero[i]    = 0.;
  }
  res.col(3) = triangle.col(3);
  res.col(4) = rotate0(triangle.col(4), zero, theta, psi);
  res(1, 3)  = 0.;
  for(int i = 0; i < 3; i ++)
    res(1, 3) += res.col(4).dot(res.col(i));
  res(1, 3)   /= 3.;
  return res;
}

// <[x, y, t], triangle.n> == triangle.k
template <typename T> bool tilter<T>::onTriangle(T& z, const Eigen::Matrix<T, 3, 5>& tri, const Eigen::Matrix<T, 2, 1>& geom) {
  Vec3 v0;
  Vec3 camera;
  v0[0] = 0;
  v0[1] = 0;
  v0[2] = 1;
  camera[0] = geom[0];
  camera[1] = geom[1];
  camera[2] = 0;
  const T t((tri(1, 3) - tri.col(4).dot(camera)) / (tri.col(4).dot(v0)));
  Eigen::Matrix<T, 2, 3> tritri;
  for(int i = 0; i < 3; i ++) {
    tritri(0, i) = tri(0, i);
    tritri(1, i) = tri(1, i);
  }
  bool f = true;
  f = f && sameSide2(tritri.col(0), tritri.col(1), tritri.col(2), geom);
  f = f && sameSide2(tritri.col(1), tritri.col(2), tritri.col(0), geom);
  f = f && sameSide2(tritri.col(2), tritri.col(0), tritri.col(1), geom);
  z = camera[2] + v0[2] * t;
  return f;
}

template <typename T> bool tilter<T>::scale(Eigen::Matrix<T, 3, 3>& A, Eigen::Matrix<T, 3, 1>& b, const Eigen::Matrix<T, 3, 3>& vorig, const Eigen::Matrix<T, 3, 3>& vto) {
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

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tilter<T>::tilt(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& in, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& bump, const int& idx) {
  Mat result(in.rows(), in.cols());
  for(int i = 0; i < in.rows(); i ++)
    for(int j = 0; j < in.cols(); j ++)
      result(i, j) = 0.;
  if(in.rows() != bump.rows() || in.cols() != bump.cols()) {
    std::cerr << "tilt: size mismatch..." << std::endl;
    return result;
  }
  std::cerr << "making triangles" << std::endl;
  std::vector<Triangles> triangles;
  for(int i = 0; i < in.rows() - 1; i ++)
    for(int j = 0; j < in.cols() - 1; j ++) {
      triangles.push_back(makeTriangle(i, j, in, bump, false));
      triangles.push_back(makeTriangle(i, j, in, bump, true));
    }
  Mat zb(in.rows(), in.cols());
  for(int i = 0; i < zb.rows(); i ++)
    for(int j = 0; j < zb.cols(); j ++)
      zb(i, j) = 0.;
  const T radius(std::max(T(in.rows()), T(in.cols())) * 2.);
  Vec2 zero2;
  Vec3 pcenter;
  Vec3 zero3;
  zero2[0]   = 0.;
  zero2[1]   = 0.;
  pcenter[0] = in.rows() / 2.;
  pcenter[1] = in.cols() / 2.;
  pcenter[2] = radius;
  zero3[0]   = 0.;
  zero3[1]   = 0.;
  zero3[2]   = 0.;
  const int i = idx % samples;
  {
    std::cerr << "rotate: " << i << "/" << samples << std::endl;
    const T theta(2. * Pi * i / samples);
    std::vector<Triangles> rotriangles;
    for(int j = 0; j < triangles.size(); j ++)
      rotriangles.push_back(rotate(rotate(triangles[j], pcenter, theta, psi), pcenter, - theta, Pi / 2.));
    std::cerr << "draw: " << i << "/" << samples << std::endl;
    for(int j = 0; j < zb.rows(); j ++)
      for(int k = 0; k < zb.cols(); k ++)
        zb(j, k) = - 20000;
    Mat zb0(zb);
    Mat result0(result);
    for(int j = 0; j < rotriangles.size(); j ++) {
      Triangles& tri = rotriangles[j];
      Vec2 gs[3];
      for(int k = 0; k < 3; k ++) {
        gs[k][0] = tri(0, k);
        gs[k][1] = tri(1, k);
      }
      int ll = int(std::min(std::min(gs[0][0], gs[1][0]), gs[2][0]));
      int rr = std::ceil(std::max(std::max(gs[0][0], gs[1][0]), gs[2][0]));
      int bb = int(std::min(std::min(gs[0][1], gs[1][1]), gs[2][1]));
      int tt = std::ceil(std::max(std::max(gs[0][1], gs[1][1]), gs[2][1]));
      for(int y = std::max(0, ll); y < std::min(rr, int(in.rows())); y ++)
        for(int x = std::max(0, bb); x < std::min(tt, int(in.cols())); x ++) {
          T z;
          Vec2 midgeom;
          midgeom[0] = x + .5;
          midgeom[1] = y + .5;
          if(onTriangle(z, tri, midgeom) && zb(y, x) < z) {
            result(y, x)  = tri(0, 3);
            zb(y, x)      = z;
          }
          if(zb0(y, x) < tri(1, 3)) {
            result0(y, x) = tri(0, 3);
            zb0(y, x)     = tri(1, 3);
          }
        }
    }
    for(int y = 0; y < result.rows(); y ++)
      for(int x = 0; x < result.cols(); x ++)
        if(zb(y, x) < 0)
          result(y, x) = result0(y, x);
  }
  return result;
}

#define _TILT_
#endif

