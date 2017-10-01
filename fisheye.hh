#if !defined(_2D3D_PSEUDO_)

#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "enlarge.hh"
#include "edgedetect.hh"
#include "obj2vector.hh"
#include "tilt.hh"

using std::complex;
using std::abs;
using std::vector;
using std::sort;
using std::max;
using std::min;
using std::isfinite;

template <typename T> int cmpbump(const Eigen::Matrix<T, 3, 1>& x0, const Eigen::Matrix<T, 3, 1>& x1) {
  return x0[2] < x1[2];
}

template <typename T> class PseudoBump {
public:
  int    z_max;
  int    stp;
  int    rstp;
  T      roff;
  T      cdist;
  T      rdist;
  T      cutoff;
  T      cutz;
  T      cthresh;
  int    crowd;
  int    vmax;
  T      sthresh;
  typedef complex<T> U;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> MatU;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1>              Vec;
  typedef Eigen::Matrix<U, Eigen::Dynamic, 1>              VecU;
  typedef Eigen::Matrix<T, 2, 1> Vec2;
  typedef Eigen::Matrix<T, 3, 1> Vec3;
  
  PseudoBump();
  ~PseudoBump();
  void initialize(const int& z_max, const int& stp, const int& crowd, const int& vmax, const T& rdist);
  
  Mat getPseudoBumpSub(const Mat& work, const int& rstp, const bool& elim);
  Mat getPseudoBump(const Mat& input, const bool& y_only = false);
  Mat getPseudoBumpVec(const Mat& input, vector<Vec3>& points, vector<Eigen::Matrix<int, 3, 1> >& delaunays, const bool& y_only = false);
private:
  T sgn(const T& x);
  Vec getLineAxis(Vec p, Vec c, const int& w, const int& h);
  Eigen::Matrix<Mat, Eigen::Dynamic, Eigen::Dynamic> prepareLineAxis(const Vec& p0, const Vec& p1, const int& z0, const int& rstp);
  T   getImgPt(const Mat& img, const T& y, const T& x);
  Vec indiv(const Vec& p0, const Vec& p1, const T& pt);
  Mat complement(const Mat& in, const int& crowd, const int& vmax, vector<Vec3>& geoms, vector<Eigen::Matrix<int, 3, 1> >& delaunay);
  Mat integrate(const Mat& input);
  
  int ww;
  int hh;
  Mat Dop;
  T   Pi;
};

template <typename T> PseudoBump<T>::PseudoBump() {
  initialize(20, 20, 20, 800, T(.75));
}

template <typename T> PseudoBump<T>::~PseudoBump() {
  ;
}

template <typename T> void PseudoBump<T>::initialize(const int& z_max, const int& stp, const int& crowd, const int& vmax, const T& rdist) {
  this->z_max   = z_max;
  this->stp     = stp;
  this->rstp    = stp * 2;
  this->crowd   = crowd;
  this->vmax    = vmax;
  this->roff    = 1. / 6.;
  this->cdist   = - 1.;
  this->rdist   = rdist * (- this->cdist);
  this->cutoff  = T(.5);
  this->cutz    = T(.125);
  this->cthresh = T(1);
  this->sthresh = T(1e-8) / T(256);
  Pi = 4. * atan2(T(1.), T(1.));
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
  return;
};

template <typename T> T PseudoBump<T>::sgn(const T& x) {
  if(x < T(0))
    return - T(1);
  if(x > T(0))
    return T(1);
  return T(0);
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::getPseudoBumpSub(const Mat& work0, const int& rstp, const bool& elim) {
  // N.B. we assume bump map operation as a linear, but it isn't.
  // N.B. integrate for local to global operation.
  const Mat work(integrate(work0));
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
      result(j, i) = T(.5);
    Vec zval(result.rows());
    for(int j = 0; j < zval.size(); j ++)
      zval[j] = T(0);
    for(int s = 0; s < result.rows(); s ++) {
      Vec pt(indiv(p0, p1, s / T(result.rows())));
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
        const T   lr(msl[msl.size() - 1] * msl[msl.size() - 1] / (msl.dot(msl) -
 msl[msl.size() - 1] * msl[msl.size() - 1]));
        const T   rr(msr[0] * msr[0] / (msr.dot(msr) - msr[0] * msr[0]));
        const T   n2(abs(msl[msl.size() - 1] - msr[0]));
        if(isfinite(n2) && zval[s] < n2 &&
           (cthresh / msl.size() <= lr ||
            cthresh / msr.size() <= rr) ) {
          result(s, i) = zz / T(lrf.rows());
          zval[s]      = n2;
        }
      }
    }
    // XXX: unfortunately, too near or too far is not stable.
    for(int s = 0; s < result.rows(); s ++)
      if(result(s, i) <= cutz || T(1) - cutz <= result(s, i))
        result(s, i) = T(.5);
  }
  if(elim) {
    edgedetect<T> detect;
    auto edge(detect.detect(detect.detect(result, detect.DETECT_Y), detect.COLLECT_Y));
    T avg(0);
    for(int i = 0; i < result.cols(); i ++)
      avg = max(avg, edge.col(i).dot(edge.col(i)));
    avg = sqrt(avg / result.rows());
    for(int i = 0; i < result.rows(); i ++)
      for(int j = 0; j < result.cols(); j ++)
        if(edge(i, j) <= cutoff * avg)
          result(i, j) = - T(4);
  }
  return result;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::getPseudoBump(const Mat& input, const bool& y_only) {
  Mat result(input.rows(), input.cols());
  if(y_only)
    result = getPseudoBumpSub(input, rstp, false);
  else
    result = getPseudoBumpSub(input.transpose(), rstp, false).transpose() +
             getPseudoBumpSub(input, rstp, false);
  // bump map is in the logic exchanged convex part and concave part.
  return - result;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::getPseudoBumpVec(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& input, vector<Eigen::Matrix<T, 3, 1> >& points, vector<Eigen::Matrix<int, 3, 1> >& delaunays, const bool& y_only) {
  Mat result(input.rows(), input.cols());
  if(y_only)
    result = getPseudoBumpSub(input, rstp, true);
  else
    result = getPseudoBumpSub(input.transpose(), rstp, true).transpose() +
             getPseudoBumpSub(input, rstp, true);
  return - complement(result, crowd, vmax, points, delaunays);
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> PseudoBump<T>::getLineAxis(Vec p, Vec c, const int& w, const int& h) {
  // suppose straight ray from infinitely distant 1 x 1 size.
  // camera geometry with c, lookup on p[2] z-distance
  // virtual images.
  p[0]  = (p[0] - h / 2.) / (h / 2.);
  p[1]  = (p[1] - w / 2.) / (w / 2.);
  // virtually : z = p[2], p[2] = 0.
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
      // XXX checkme:
      cpoint[2] = (zi + 1) / T(z0) * rdist;
      result(zi, 0).col(s) = getLineAxis(cpoint, camera0, T(ww), T(hh)) - center;
      result(zi, 1).col(s) = getLineAxis(cpoint, camera1, T(ww), T(hh)) - center;
    }
  }
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

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::complement(const Mat& in, const int& crowd, const int& vmax, vector<Vec3>& geoms, vector<Eigen::Matrix<int, 3, 1> >& delaunay) {
  Mat result(in);
  T   avg(0);
  int cavg(0);
  geoms = vector<Vec3>();
  for(int i = 0; i < in.rows() / crowd; i ++)
    for(int j = 0; j < in.cols() / crowd; j ++) {
      Vec work(3);
      work[0] = 0;
      work[1] = 0;
      work[2] = 0;
      int count(0);
      for(int ii = i * crowd; ii < min((i + 1) * crowd, int(in.rows())); ii ++)
        for(int jj = j * crowd; jj < min((j + 1) * crowd, int(in.cols())); jj ++) {
          if(0 <= in(ii, jj)) {
            work[0] += ii;
            work[1] += jj;
            work[2] += in(ii, jj);
            count ++;
            avg += in(ii, jj);
            cavg ++;
          }
        }
      if(count)
        geoms.push_back(work / count);
    }
  sort(geoms.begin(), geoms.end(), cmpbump<T>);
  int start(max(0, int(geoms.size() - vmax) / 2));
  int end(min(int(geoms.size()), int(geoms.size() + vmax) / 2));
  geoms.erase(geoms.begin(), geoms.begin() + start);
  geoms.resize(end - start);
  Vec work(3);
  work[0] = 0;
  work[1] = 0;
  if(cavg)
    work[2] = avg / cavg;
  else
    work[2] = 0;
  geoms.push_back(work);
  work[0] = in.rows();
  geoms.push_back(work);
  work[1] = in.cols();
  geoms.push_back(work);
  work[0] = 0;
  geoms.push_back(work);
  vector<int> idxs;
  for(int i = 0; i < geoms.size(); i ++)
    idxs.push_back(i);
  delaunay = loadBumpSimpleMesh<T>(geoms, idxs);
  tilter<T> tilt;
  result *= T(0);
  for(int i = 0; i < delaunay.size(); i ++) {
    const Vec3& p0(geoms[delaunay[i][0]]);
    const Vec3& p1(geoms[delaunay[i][1]]);
    const Vec3& p2(geoms[delaunay[i][2]]);
    for(int j = 0; j < in.rows(); j ++)
      for(int k = 0; k < in.cols(); k ++) {
        Vec3 q;
        q[0] = j;
        q[1] = k;
        q[2] = 0;
        if(tilt.sameSide2(p0, p1, p2, q) &&
           tilt.sameSide2(p1, p2, p0, q) &&
           tilt.sameSide2(p2, p0, p1, q))
          result(j, k) = (geoms[delaunay[i][0]][2] + 
                          geoms[delaunay[i][1]][2] +
                          geoms[delaunay[i][2]][2]) / T(3);
      }
  }
  return result;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::integrate(const Mat& input) {
  MatU Iop(input.rows(), input.rows()), A(input.rows(), input.rows());
  U    I(sqrt(complex<T>(- 1)));
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for
#endif
  for(int i = 0; i < Iop.rows(); i ++)
    for(int j = 0; j < Iop.cols(); j ++) {
      Iop(i, j) = exp(complex<T>(- 2) * Pi * I * complex<T>(i * j) / T(Iop.rows()));
      A(i, j)   = exp(complex<T>(  2) * Pi * I * complex<T>(i * j) / T(Iop.rows()));
    }
  Iop.row(0) *= T(0);
  for(int i = 1; i < Iop.rows(); i ++)
    Iop.row(i) /= complex<T>(- 2.) * Pi * I * T(i) / T(Iop.rows());
  Iop = A * Iop;
  
  Mat centerized(input.rows(), input.cols());
  Vec avg(input.cols());
#if defined(_OPENMP)
#pragma omp for
#endif
  for(int i = 0; i < input.cols(); i ++) {
    avg[i] = T(0);
    for(int j = 0; j < input.rows(); j ++)
      avg[i] += input(j, i);
    avg[i] /= input.rows();
    for(int j = 0; j < input.rows(); j ++)
      centerized(j, i) = input(j, i) - avg[i];
  }
  centerized = Iop.real().template cast<T>() * centerized;
#if defined(_OPENMP)
#pragma omp for
#endif
  for(int i = 0; i < input.cols(); i ++)
    for(int j = 0; j < input.rows(); j ++)
      centerized(j, i) += avg[i] * j / input.rows();
  return centerized;
}

#define _2D3D_PSEUDO_
#endif

