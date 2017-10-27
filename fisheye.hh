#if !defined(_2D3D_PSEUDO_)

#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "tilt.hh"
#include "obj2vector.hh"

using std::complex;
using std::abs;
using std::vector;
using std::sort;
using std::max;
using std::min;
using std::cos;
using std::sqrt;
using std::log;
using std::isfinite;
using std::cerr;
using std::flush;

template <typename T> int cmpbump(const Eigen::Matrix<T, 3, 1>& x0, const Eigen::Matrix<T, 3, 1>& x1) {
  return x0[2] < x1[2];
}

template <typename T> class PseudoBump {
public:
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1>              Vec;
  typedef Eigen::Matrix<T, 3, 1> Vec3;
  
  PseudoBump();
  ~PseudoBump();
  void initialize(const int& z_max, const int& stp, const int& rstp, const int& vmax);
  Mat  getPseudoBumpVec(const Mat& in, vector<Vec3>& geoms, vector<Eigen::Matrix<int, 3, 1> >& delaunay, Mat& bumps);
  Mat  getPseudoBumpSub(const Mat& work);
  
  int vmax;

private:
  Vec  getLineAxis(Vec p, Vec c, const int& w, const int& h);
  Eigen::Matrix<Mat, Eigen::Dynamic, Eigen::Dynamic> prepareLineAxis(const Vec& p0, const Vec& p1, const int& ww, const int& hh);
  T    getImgPt(const Mat& img, const T& y, const T& x);
  Vec  indiv(const Vec& p0, const Vec& p1, const T& pt);
  void complement(const Mat& in, const int& vmax, const int& guard, vector<Vec3>& geoms, vector<Eigen::Matrix<int, 3, 1> >& delaunay);
  Vec  complementLine(const Vec& line, const T& rratio = T(.5));
  void autoLevel(Mat& data, int npad = - 1);
  
  int z_max;
  int stp;
  T   rstp;
  
  int guard;
  
  T   roff;
  T   cdist;
  T   rdist;
  
  T   Pi;
};

template <typename T> PseudoBump<T>::PseudoBump() {
  initialize(30, 16, 4, 800);
}

template <typename T> PseudoBump<T>::~PseudoBump() {
  ;
}

template <typename T> void PseudoBump<T>::initialize(const int& z_max, const int& stp, const int& rstp, const int& vmax) {
  this->Pi    = T(4) * atan2(T(1), T(1));
  this->z_max = z_max;
  this->stp   = stp;
  this->rstp  = stp * rstp;
  this->vmax  = vmax;
  this->guard = rstp / 2;
  this->roff  = T(1) / T(6);
  // N.B. ray is from infinite far, so same side of these.
  this->cdist = T(2);
  this->rdist = T(1);
  return;
};

// bump it with abs(d/dt local color) / abs(d/dt color) >> 0.
template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::getPseudoBumpSub(const Mat& work) {
  Mat result(work.rows(), work.cols());
  Mat workv(work.rows() * stp, work.cols());
  Mat workh(work.rows(), work.cols() * stp);
  for(int i = 0; i < workv.rows(); i ++)
    if(i % stp == stp / 2)
      workv.row(i) = work.row(i / stp);
    else
      for(int j = 0; j < workv.cols(); j ++)
        workv(i, j) = - T(1);
  for(int i = 0; i < workh.cols(); i ++)
    if(i % stp == stp / 2)
      workh.col(i) = work.col(i / stp);
    else
      for(int j = 0; j < workh.rows(); j ++)
        workh(j, i) = - T(1);
  for(int i = 0; i < workv.cols(); i ++)
    workv.col(i) = complementLine(workv.col(i));
  for(int i = 0; i < workh.rows(); i ++)
    workh.row(i) = complementLine(workh.row(i));
  Vec p0(3), p1(3), p0x(3), p1x(3);
  p0[0]  = 0;
  p0[1]  = workv.cols() / 2;
  p0[2]  = 0;
  p1[0]  = workv.rows();
  p1[1]  = workv.cols() / 2;
  p1[2]  = 0;
  p0x[0] = workh.rows() / 2;
  p0x[1] = 0;
  p0x[2] = 0;
  p1x[0] = workh.rows() / 2;
  p1x[1] = workh.cols();
  p1x[2] = 0;
  cerr << " bump" << flush;
  Eigen::Matrix<Mat, Eigen::Dynamic, Eigen::Dynamic> lrf( prepareLineAxis(p0,  p1, workv.cols(), workv.rows()));
  Eigen::Matrix<Mat, Eigen::Dynamic, Eigen::Dynamic> lrfx(prepareLineAxis(p0x, p1x, workh.cols(), workh.rows()));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < result.cols(); i ++) {
    Vec p0(3), p1(3);
    p0[0] = 0;
    p0[1] = i;
    p0[2] = 0;
    p1[0] = work.rows();
    p1[1] = i;
    p1[2] = 0;
    cerr << "." << flush;
    for(int j = 0; j < result.rows(); j ++)
      result(j, i) = - T(1);
    Vec zval(result.rows());
    for(int j = 0; j < zval.size(); j ++)
      zval[j] = T(0);
    for(int s = 0; s < result.rows(); s ++) {
      Vec pt(indiv(p0, p1, s / T(result.rows())));
      for(int zz = 0; zz < lrf.rows(); zz ++) {
        Vec c(lrf(zz, 0).cols()), cx(lrfx(zz, 0).cols());
        for(int u = 0; u < c.size(); u ++) {
          c[u]   = getImgPt(workv, lrf( zz, 0)(0, u) + pt[0] * stp + stp / 2,
                                   lrf( zz, 0)(1, u) + pt[1]);
          c[u]  -= getImgPt(workv, lrf( zz, 1)(0, u) + pt[0] * stp + stp / 2,
                                   lrf( zz, 1)(1, u) + pt[1]);
          cx[u]  = getImgPt(workh, lrfx(zz, 0)(0, u) + pt[0],
                                   lrfx(zz, 0)(1, u) + pt[1] * stp + stp / 2);
          cx[u] -= getImgPt(workh, lrfx(zz, 1)(0, u) + pt[0],
                                   lrfx(zz, 1)(1, u) + pt[1] * stp + stp / 2);
        }
        Vec cc(c.size() / 4 * 2 + 1);
        Vec ccx(cc.size());
        for(int u =  c.size() / 2 - c.size() / 4;
                u <= c.size() / 2 + c.size() / 4;
                u ++) {
          const int idx(u - c.size() / 2 + c.size() / 4);
          cc[idx]  = c[u];
          ccx[idx] = cx[u];
        }
        Vec cl(c.size() / 2 + 1);
        Vec cr(cl.size()), cxl(cl.size()), cxr(cl.size());
        for(int u = 0; u < cl.size(); u ++) {
          cl[ u] = c[ u];
          cxl[u] = cx[u];
          cr[ u] = c[ u - cl.size() + c.size()];
          cxr[u] = cx[u - cl.size() + cx.size()];
        }
        Vec ccl(c.size() / 4 + 1);
        Vec ccr(ccl.size()), ccxl(ccl.size()), ccxr(ccl.size());
        for(int u = 0; u < ccl.size(); u ++) {
          ccl[ u] = cc[ u - ccl.size() + cc.size() / 2 + 1];
          ccxl[u] = ccx[u - ccl.size() + cc.size() / 2 + 1];
          ccr[ u] = cc[ u - ccl.size() + cc.size() / 2 + cc.size() / 4 + 2];
          ccxr[u] = ccx[u - ccl.size() + cc.size() / 2 + cc.size() / 4 + 2];
        }
        // N.B. simply take the ratio of local and far difference on two eyes
        //      on left side and right side and center, then ratio it.
        const T n2((ccl.dot( ccl)  / cl.dot( cl)  +
                    ccr.dot( ccr)  / cr.dot( cr)) /
                   (cc.dot( cc)  / c.dot( c)) +
                    cc.dot( cc)  / c.dot( c)  +
                   (ccxl.dot(ccxl) / cxl.dot(cxl) +
                    ccxr.dot(ccxr) / cxr.dot(cxr)) /
                   (ccx.dot(ccx) / cx.dot(cx)) +
                    ccx.dot(ccx) / cx.dot(cx));
        if(isfinite(n2) && zval[s] < n2) {
          // N.B. If increase zz, decrease the distance from camera.
          //      And, zz->0 is treated as distant in tilter.
          result(s, i) = (T(1) + zz) / T(lrf.rows());
          zval[s]      = n2;
        }
      }
    }
  }
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < result.cols(); i ++)
    result.col(i) = complementLine(result.col(i));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < result.rows(); i ++)
    result.row(i) = complementLine(result.row(i));
  autoLevel(result);
  return result;
}

// get bump with multiple scale and vectored result.
template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::getPseudoBumpVec(const Mat& in, vector<Vec3>& geoms, vector<Eigen::Matrix<int, 3, 1> >& delaunay, Mat& bumps) {
  bumps = getPseudoBumpSub(in);
  complement(bumps, vmax, guard, geoms, delaunay);
  Mat result(in.rows(), in.cols());
  for(int i = 0; i < result.rows(); i ++)
    for(int j = 0; j < result.cols(); j ++)
      result(i, j) = T(0);
  tilter<T> tilt;
  vector<Vec3> newgeoms(geoms), newgeomsM(geoms);
  for(int i = 0; i < newgeoms.size(); i ++) {
    newgeoms[i][2]  = - T(1);
    newgeomsM[i][2] =   T(2);
  }
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
          newgeoms[delaunay[i][0]][2] = max(newgeoms[delaunay[i][0]][2], result(j, k));
          newgeoms[delaunay[i][1]][2] = max(newgeoms[delaunay[i][1]][2], result(j, k));
          newgeoms[delaunay[i][2]][2] = max(newgeoms[delaunay[i][2]][2], result(j, k));
          newgeomsM[delaunay[i][0]][2] = min(newgeomsM[delaunay[i][0]][2], result(j, k));
          newgeomsM[delaunay[i][1]][2] = min(newgeomsM[delaunay[i][1]][2], result(j, k));
          newgeomsM[delaunay[i][2]][2] = min(newgeomsM[delaunay[i][2]][2], result(j, k));
      }
  }
  autoLevel(result);
  for(int i = 0; i < geoms.size(); i ++)
    geoms[i][2] = newgeoms[i][2] + newgeomsM[i][2];
  return result;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> PseudoBump<T>::getLineAxis(Vec p, Vec c, const int& w, const int& h) {
  // suppose straight ray from infinitely distant 2 x 2 size.
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

template <typename T> Eigen::Matrix<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::prepareLineAxis(const Vec& p0, const Vec& p1, const int& ww, const int& hh) {
  Vec dir(p1 - p0);
  dir /= sqrt(dir.dot(dir));
  Vec camera0(3);
  camera0[0] = dir[0] * roff;
  camera0[1] = dir[1] * roff;
  camera0[2] = cdist;
  Vec camera1(- camera0);
  camera1[2] = cdist;
  const Vec center(indiv(p0, p1, T(1) / T(2)));
  
  Eigen::Matrix<Mat, Eigen::Dynamic, Eigen::Dynamic> result(z_max, 2);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int zi = 0; zi < z_max; zi ++) {
    result(zi, 0) = Mat(3, stp);
    result(zi, 1) = Mat(3, stp);
    for(int s = 0; s < stp; s ++) {
      Vec cpoint(3);
      cpoint[0] = (s / (stp - 1.) - 1. / 2.) * rstp + center[0];
      cpoint[1] = (s / (stp - 1.) - 1. / 2.) * rstp + center[1];
      cpoint[2] = (zi + 1) / T(z_max) * rdist;
      result(zi, 0).col(s) = getLineAxis(cpoint, camera0, T(ww), T(hh)) - center;
      result(zi, 1).col(s) = getLineAxis(cpoint, camera1, T(ww), T(hh)) - center;
    }
  }
  return result;
}

template <typename T> T PseudoBump<T>::getImgPt(const Mat& img, const T& y, const T& x) {
  const int& w(img.cols());
  const int& h(img.rows());
  const int  xx(abs((int(x + .5) + 3 * w) % (2 * w) - w) % w);
  const int  yy(abs((int(y + .5) + 3 * h) % (2 * h) - h) % h);
  return img(yy, xx);
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> PseudoBump<T>::indiv(const Vec& p0, const Vec& p1, const T& pt) {
  return (p1 - p0) * pt + p0;
}

template <typename T> void PseudoBump<T>::complement(const Mat& in, const int& vmax, const int& guard, vector<Vec3>& geoms, vector<Eigen::Matrix<int, 3, 1> >& delaunay) {
  T avg(0);
  geoms = vector<Vec3>();
  for(int i = 1; i < (in.rows() / guard - 1); i ++)
    for(int j = 1; j < (in.cols() / guard - 1); j ++) {
      Vec work(3);
      work[0] = i * guard;
      work[1] = j * guard;
      work[2] = 0;
      vector<T> wbuf;
      for(int ii = i * guard; ii < min((i + 1) * guard, int(in.rows())); ii ++)
        for(int jj = j * guard; jj < min((j + 1) * guard, int(in.cols())); jj ++)
          wbuf.push_back(in(ii, jj));
      work[2] = (wbuf[0] + wbuf[wbuf.size() - 1]) / T(2);
      geoms.push_back(work);
      avg += work[2];
    }
  avg /= T((in.rows() / guard) * (in.cols() / guard));
  sort(geoms.begin(), geoms.end(), cmpbump<T>);
  // XXX lowpoly is needed.
  if(vmax < geoms.size())
    geoms.erase(geoms.begin() + vmax / 2, geoms.end() - vmax / 2);
  Vec work(3);
  work[0] = 0;
  work[1] = 0;
  work[2] = avg;
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
  return;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> PseudoBump<T>::complementLine(const Vec& line, const T& rratio) {
  vector<int> ptsi;
  vector<T>   pts;
  for(int i = 0; i < line.size(); i ++)
    if(T(0) <= line[i]) {
      ptsi.push_back(i);
      pts.push_back(line[i]);
    }
  Vec result(line);
  if(ptsi.size() < 2)
    return line;
  ptsi.insert(ptsi.begin(), - ptsi[1]);
  pts.insert(pts.begin(), pts[1]);
  ptsi.push_back(line.size() + (line.size() - ptsi[ptsi.size() - 3]));
  pts.push_back(pts[pts.size() - 3]);
  int rng[3];
  rng[0] = rng[1] = rng[2] = 1;
  for(int i = 0; i < line.size(); i ++) {
    if(result[i] >= T(0))
      continue;
    for(; rng[1] < ptsi.size() - 2 && ptsi[rng[1]] < i; rng[1] ++) ;
    rng[0] = rng[1] - 1;
    rng[2] = rng[1] + 1;
    const T ratio((ptsi[rng[2]] - ptsi[rng[1]]) /
                  (ptsi[rng[1]] - ptsi[rng[0]]));
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
    result[i] = max(min(result[i], T(2)), - T(1));
  }
  return result;
}

template <typename T> void PseudoBump<T>::autoLevel(Mat& data, int npad) {
  if(npad <= 0)
    npad = (data.rows() + data.cols()) * 8;
  vector<T> stat;
  for(int j = 0; j < data.rows(); j ++)
    for(int k = 0; k < data.cols(); k ++)
      stat.push_back(data(j, k));
  sort(stat.begin(), stat.end());
  const T mm(stat[npad]);
  T MM(stat[stat.size() - 1 - npad]);
  if(MM == mm)
    MM = mm + 1.;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < data.rows(); i ++)
    for(int j = 0; j < data.cols(); j ++)
      data(i, j) = (max(min(data(i, j), MM), mm) - mm) / (MM - mm);
  return;
}

#define _2D3D_PSEUDO_
#endif

