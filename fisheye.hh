#if !defined(_2D3D_PSEUDO_)

#include <cstdio>
#include <cmath>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "edgedetect.hh"

using std::complex;
using std::abs;
using std::vector;
using std::sort;
using std::max;
using std::min;
using std::isfinite;

template <typename T> class PseudoBump {
public:
  int    z_max;
  int    stp;
  int    rstp;
  T      roff;
  T      cdist;
  T      rdist;
  T      cthresh;
  typedef complex<T> U;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> MatU;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1>              Vec;
  typedef Eigen::Matrix<U, Eigen::Dynamic, 1>              VecU;
  
  PseudoBump();
  void initialize(const int& z_max, const int& stp, const T& cthresh);
  ~PseudoBump();
  
  Mat getPseudoBumpSub(const Mat& work, const int& rstp);
  Mat getPseudoBump(const Mat& input, const bool& y_only);
  Mat rgb2l(const Mat rgb[3]);
  Vec complementLine(const Vec& line, const T& rratio = T(.8));
private:
  T zz(const T& t);
  T sgn(const T& x);
  Vec getLineAxis(Vec p, Vec c, const int& w, const int& h);
  Eigen::Matrix<Mat, Eigen::Dynamic, Eigen::Dynamic> prepareLineAxis(const Vec& p0, const Vec& p1, const int& z0, const int& rstp);
  T   getImgPt(const Mat& img, const T& y, const T& x);
  Vec indiv(const Vec& p0, const Vec& p1, const T& pt);
  
  int ww;
  int hh;
  Mat Dop;
  T   Pi;
};

template <typename T> PseudoBump<T>::PseudoBump() {
  initialize(20, 20, .5);
}

template <typename T> PseudoBump<T>::~PseudoBump() {
  ;
}

template <typename T> void PseudoBump<T>::initialize(const int& z_max, const int& stp, const T& cthresh) {
  this->z_max   = z_max;
  this->stp     = stp;
  this->rstp    = stp * 2;
  this->roff    = 1. / 6.;
  this->cdist   = - 1.;
  this->rdist   =    .5;
  this->cthresh = cthresh;
  Pi            = 4. * atan2(T(1.), T(1.));
  return;
};

template <typename T> T PseudoBump<T>::sgn(const T& x) {
  if(x < T(0))
    return - T(1);
  if(x > T(0))
    return T(1);
  return T(0);
}

template <typename T> T PseudoBump<T>::zz(const T& t) {
  // XXX select me:
  // return (t + 1) / z_max;
  return (t + 1 + z_max) / z_max / T(3);
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::rgb2l(const Mat rgb[3]) {
  Mat Y(T( .299  ) * rgb[0] + T( .587  ) * rgb[1] + T(.114   ) * rgb[2]);
  Mat U(T(-.14713) * rgb[0] + T(-.28886) * rgb[1] + T(.436   ) * rgb[2]);
  Mat V(T( .615  ) * rgb[0] + T(-.51499) * rgb[1] + T(-.10001) * rgb[2]);
  Mat result(rgb[0].rows(), rgb[0].cols());
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for(int j = 0; j < rgb[0].rows(); j ++)
    for(int k = 0; k < rgb[0].cols(); k ++)
      result(j, k) = sqrt(Y(j, k) * Y(j, k) + U(j, k) * U(j, k) + V(j, k) * V(j, k));
  return result;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::getPseudoBumpSub(const Mat& work, const int& rstp) {
  Mat result(work.rows(), work.cols());
  ww = work.cols();
  hh = work.rows();
  Vec p0(3), p1(3);
  p0[0] = 0;
  p0[1] = work.cols() / 2;
  p0[2] = 0;
  p1[0] = work.rows();
  p1[1] = work.cols() / 2;
  p1[2] = 0;
  cerr << " bump";
  cerr.flush();
  Eigen::Matrix<Mat, Eigen::Dynamic, Eigen::Dynamic> lrf(prepareLineAxis(p0, p1, z_max, rstp));
  for(int i = 0; i < lrf.rows(); i ++)
    for(int j = 0; j < lrf(i, 0).cols(); j ++) {
      lrf(i, 0)(1, j) = 0;
      lrf(i, 1)(1, j) = 0;
    }
  const int& delta(result.rows());
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for(int i = 0; i < result.cols(); i ++) {
    Vec p0(3), p1(3);
    p0[0] = 0;
    p0[1] = i;
    p0[2] = 0;
    p1[0] = work.rows();
    p1[1] = i;
    p1[2] = 0;
    cerr << ".";
    cerr.flush();
    for(int j = 0; j < result.rows(); j ++)
      result(j, i) = - T(1.);
    Vec zval(result.rows());
    for(int j = 0; j < zval.size(); j ++)
      zval[j] = T(0);
    for(int s = 0; s < delta; s ++) {
      Vec pt(indiv(p0, p1, s / T(delta)));
      for(int zz = 0; zz < lrf.rows(); zz ++) {
        Vec c(lrf(zz, 0).cols());
        for(int u = 0; u < c.size(); u ++) {
          c[u]  = getImgPt(work, lrf(zz, 0)(0, u) + pt[0],
                                 lrf(zz, 0)(1, u) + pt[1]);
          c[u] -= getImgPt(work, lrf(zz, 1)(0, u) + pt[0],
                                 lrf(zz, 1)(1, u) + pt[1]);
        }
        Vec lc(c.size() / 2 + 1);
        Vec rc(c.size() / 2 + 1);
        for(int u = 0; u < lc.size(); u ++)
          lc[u] = c[u];
        for(int u = 0; u < rc.size(); u ++)
          rc[u] = c[u - rc.size() + c.size()];
        const Vec msl(Dop * lc);
        const Vec msr(Dop * rc);
        const T   lr(msl[msl.size() - 1] * msl[msl.size() - 1] / (msl.dot(msl) - msl[msl.size() - 1] * msl[msl.size() - 1]));
        const T   rr(msr[0] * msr[0] / (msr.dot(msr) - msr[0] * msr[0]));
        const T   n2(abs(msl[msl.size() - 1] - msr[0]));
        if(isfinite(n2) && zval[s] < n2 &&
           cthresh / msl.size() <= lr &&
           cthresh / msr.size() <= rr) {
          result(s, i) = zz / T(z_max);
          zval[s]      = n2;
        }
      }
    }
    result.col(i) = complementLine(result.col(i));
  }
  return result;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::getPseudoBump(const Mat& input, const bool& y_only) {
  Mat result(input.rows(), input.cols());
  if(y_only)
    result = getPseudoBumpSub(input, rstp);
  else
    result = getPseudoBumpSub(input.transpose(), rstp).transpose() +
             getPseudoBumpSub(input, rstp);
  // bump map is in the logic exchanged convex part and concave part.
  return - result;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> PseudoBump<T>::getLineAxis(Vec p, Vec c, const int& w, const int& h) {
  // suppose straight ray from infinitely distant 1 x 1 size.
  // camera geometry with c, lookup on p[2] z-distance
  // virtual images.
  p[0]  = (p[0] - h / 2.) / (h / 2.);
  p[1]  = (p[1] - w / 2.) / (w / 2.);
  // XXX: virtually : z = p[2], p[2] = 0.
  // <c + (p - c) * t, [0, 0, 1]> = z
  T   t((p[2] - c[2]) / (- c[2]));
  Vec work((p - c) * t + c);
  work[0] = (work[0] * h / 2. + h / 2.);
  work[1] = (work[1] * w / 2. + w / 2.);
  return work;
}

template <typename T> Eigen::Matrix<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::prepareLineAxis(const Vec& p0, const Vec& p1, const int& z0, const int& rstp) {
  Vec dir(p1 - p0);
  dir /= sqrt(dir.dot(dir));
  Vec camera0(3);
  camera0[0] = dir[0] * roff;
  camera0[1] = dir[1] * roff;
  camera0[2] = cdist;
  Vec camera1(- camera0);
  camera1[2] = cdist;
  const Vec center(indiv(p0, p1, T(1) / T(2)));
  
  Eigen::Matrix<Mat, Eigen::Dynamic, Eigen::Dynamic> result(z0, 2);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for
#endif
  for(int zi = 0; zi < z0; zi ++) {
    result(zi, 0) = Mat(3, stp);
    result(zi, 1) = Mat(3, stp);
    for(int s = 0; s < stp; s ++) {
      Vec cpoint(3);
      cpoint[0] = (s / (stp - 1.) - 1. / 2.) * rstp + center[0];
      cpoint[1] = (s / (stp - 1.) - 1. / 2.) * rstp + center[1];
      cpoint[2] = zz(zi + 1) * rdist;
      result(zi, 0).col(s) = getLineAxis(cpoint, camera0, T(ww), T(hh)) - center;
      result(zi, 1).col(s) = getLineAxis(cpoint, camera1, T(ww), T(hh)) - center;
    }
  }
  MatU Dopb(stp / 2 + 1, stp / 2 + 1);
  MatU Iopb(Dopb.rows(), Dopb.cols());
  U I(sqrt(complex<T>(- 1)));
#if defined(_OPENMP)
#pragma omp for
#endif
  for(int i = 0; i < Dopb.rows(); i ++)
    for(int j = 0; j < Dopb.cols(); j ++) {
      Dopb(i, j) = exp(complex<T>(- 2) * Pi * I * complex<T>(i * j) / T(Dopb.rows()));
      Iopb(i, j) = exp(complex<T>(  2) * Pi * I * complex<T>(i * j) / T(Iopb.rows()));
    }
#if defined(_OPENMP)
#pragma omp for
#endif
  for(int i = 0; i < Dopb.rows(); i ++)
    Dopb.row(i) *= complex<T>(- 2.) * Pi * I * T(i) / T(Dopb.rows());
  Dop = (Iopb * Dopb).real();
  return result;
}

template <typename T> T PseudoBump<T>::getImgPt(const Mat& img, const T& y, const T& x) {
  const int& w(img.cols());
  const int& h(img.rows());
  const int xx(abs((int(x + .5) + 3 * w) % (2 * w) - w) % w);
  const int yy(abs((int(y + .5) + 3 * h) % (2 * h) - h) % h);
  return img(yy, xx);
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> PseudoBump<T>::indiv(const Vec& p0, const Vec& p1, const T& pt) {
  return (p1 - p0) * pt + p0;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> PseudoBump<T>::complementLine(const Vec& line, const T& rratio) {
  vector<int> ptsi;
  vector<T>   pts;
  bool flag = true;
  for(int i = 0; i < line.size(); i ++)
    if(T(0) <= line[i]) {
      if(!i)
        flag = false;
      else if(flag) {
        ptsi.push_back(- i);
        pts.push_back(line[i]);
        flag = false;
      }
      ptsi.push_back(i);
      pts.push_back(line[i]);
    }
  Vec result(line);
  if(ptsi.size() <= 0) {
    for(int i = 0; i < line.size(); i ++)
      result[i] = T(.5);
    return result;
  }
  ptsi.push_back(ptsi[ptsi.size() - 1] + 2. * (line.size() - ptsi[ptsi.size() - 1]));
  pts.push_back(pts[pts.size() - 1]);
  int rng[3];
  rng[0] = rng[1] = rng[2] = 0;
  for(int i = 0; i < line.size(); i ++) {
    if(result[i] >= T(0))
      continue;
    for(; rng[0] < ptsi.size() - 1; rng[0] ++)
      if(i <= ptsi[rng[0] + 1])
        break;
    for(; rng[1] < ptsi.size(); rng[1] ++)
      if(i <= ptsi[rng[1]] && ptsi[rng[1]] <= i)
        break;
    for(; rng[2] < ptsi.size(); rng[2] ++)
      if(i < ptsi[rng[2]])
        break;
    if(rng[2] < rng[1])
      rng[1] = (rng[2] + rng[0]) / 2;
    if(rng[1] == rng[0] && 0 < rng[0])
      rng[0] --;
    if(rng[1] == rng[2] && rng[2] < ptsi.size() - 1)
      rng[2] ++;
    const T ratio(abs((ptsi[rng[2]] - ptsi[rng[1]] + 1.) /
                      (ptsi[rng[1]] - ptsi[rng[0]] + 1.)));
    if(ratio < rratio || T(1) / ratio < rratio) {
      if(abs(ptsi[rng[2]] - ptsi[rng[1]]) >
         abs(ptsi[rng[1]] - ptsi[rng[0]])) {
        for(; rng[0] > 0; rng[0] --)
          if(ptsi[rng[1]] - ptsi[rng[0]] > (ptsi[rng[2]] - ptsi[rng[1]]) * rratio)
            break;
      } else {
        for(; rng[2] < ptsi.size() - 1; rng[2] ++)
          if(ptsi[rng[2]] - ptsi[rng[1]] > (ptsi[rng[1]] - ptsi[rng[0]]) * rratio)
            break;
      }
    }
    result[i] = 0.;
    for(int ii = 0; ii < 3; ii ++) {
      T work(1);
      for(int jj = 0; jj < 3; jj ++)
        if(ptsi[rng[ii]] != ptsi[rng[jj]])
          work *= (T(i) - ptsi[rng[jj]]) / T(ptsi[rng[ii]] - ptsi[rng[jj]]);
      result[i] += work * pts[rng[ii]];
    }
    result[i] = max(min(result[i], T(1)), T(0));
  }
  return result;
}

#define _2D3D_PSEUDO_
#endif

