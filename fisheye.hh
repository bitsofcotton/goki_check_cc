#if !defined(_2D3D_PSEUDO_)

#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "enlarge.hh"
#include "tilt.hh"
#include "obj2vector.hh"
#include "scancontext.hh"

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
  void initialize(const int& z_max, const int& stp, const int& vmax);
  Mat  getPseudoBumpVec(const Mat& in, vector<Vec3>& geoms, vector<Eigen::Matrix<int, 3, 1> >& delaunay, Mat& bumpp);
  
  int vmax;
  Vec  complementLine(const Vec& line, const T& rratio = T(.5));
  T    getImgPt(const Mat& img, const int& y, const int& x);
  void setImgPt(Mat& img, const int& y, const int& x, const T& v);

private:
  Vec  getPseudoBumpSub(const Vec& work, const Eigen::Matrix<Mat, Eigen::Dynamic, 1>& cf);
  Vec  getLineAxis(Vec p, Vec c, const int& w, const int& h);
  Eigen::Matrix<Mat, Eigen::Dynamic, 1> prepareLineAxis(const Vec& p0, const Vec& p1, const int& hh, const int& ww, const int& rstp);
  T    getImgPt(const Vec& img, const T& y);
  Vec  indiv(const Vec& p0, const Vec& p1, const T& pt);
  void autoLevel(Mat& data, int npad = - 1);
  
  int z_max;
  int stp;
  
  T   Pi;
  Mat Dop;
};

template <typename T> PseudoBump<T>::PseudoBump() {
  initialize(20, 15, 800);
}

template <typename T> PseudoBump<T>::~PseudoBump() {
  ;
}

template <typename T> void PseudoBump<T>::initialize(const int& z_max, const int& stp, const int& vmax) {
  this->Pi      = T(4) * atan2(T(1), T(1));
  this->z_max   = z_max;
  this->stp     = stp;
  this->vmax    = vmax;
  Eigen::Matrix<complex<T>, Eigen::Dynamic, Eigen::Dynamic> DFT(stp, stp), IDFT(stp, stp);
  for(int i = 0; i < DFT.rows(); i ++)
    for(int j = 0; j < DFT.cols(); j ++) {
      DFT( i, j) = exp(complex<T>(- 2.) * Pi * sqrt(complex<T>(- 1)) * complex<T>(i * j / T(stp)));
      IDFT(i, j) = exp(complex<T>(  2.) * Pi * sqrt(complex<T>(- 1)) * complex<T>(i * j / T(stp))) / T(stp);
    }
  for(int i = 0; i < DFT.rows(); i ++)
    DFT.row(i) *= complex<T>(2.) * Pi * sqrt(complex<T>(- 1)) * T(i) / T(DFT.rows());
  Dop  = (IDFT * DFT).real();
  return;
};

// bump it with abs(d/dt local color) / ||local color|| >> 0.
template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> PseudoBump<T>::getPseudoBumpSub(const Vec& work, const Eigen::Matrix<Mat, Eigen::Dynamic, 1>& cf) {
  cerr << "." << flush;
  Vec result(work.size());
  Vec workv(work.size() * stp);
  for(int i = 0; i < workv.size(); i ++)
    if(i % stp == stp / 2)
      workv[i] = work[i / stp];
    else
      workv[i] = - T(1);
  workv = complementLine(workv);
  Vec p0(3), p1(3);
  p0[0]  = 0;
  p0[1]  = 0;
  p0[2]  = 0;
  p1[0]  = result.size() * stp;
  p1[1]  = 0;
  p1[2]  = 0;
  for(int j = 0; j < result.size(); j ++)
    result[j] = - T(1);
  Vec zval(result.size());
  for(int j = 0; j < zval.size(); j ++)
    zval[j] = T(0);
  for(int s = 0; s < result.size(); s ++) {
    Vec pt(indiv(p0, p1, s / T(result.size())));
    for(int zz = 0; zz < cf.size(); zz ++) {
      Vec c(cf[zz].cols());
      for(int u = 0; u < c.size(); u ++)
        c[u] = getImgPt(workv, cf[zz](0, u) + pt[0] + stp / 2);
      Vec rc(c);
      for(int u = 0; u < c.size(); u ++)
        rc[u] = c[c.size() - 1 - u];
      const T n2(max(Dop.row(c.size() / 2).dot(c)  / sqrt(c.dot(c)),
                     Dop.row(c.size() / 2).dot(rc) / sqrt(c.dot(c))));
      if(isfinite(n2) && zval[s] < n2) {
        // N.B. If increase zz, decrease the distance from camera.
        //      And, zz->0 is treated as distant in tilter.
        result[s] = (T(1) + zz) / T(cf.size());
        zval[s]   = n2;
      }
    }
  }
  return complementLine(result);
}

// get bump with multiple scale and vectored result.
template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::getPseudoBumpVec(const Mat& in, vector<Vec3>& geoms, vector<Eigen::Matrix<int, 3, 1> >& delaunay, Mat& bumpp) {
  cerr << "bump" << flush;
  vector<Eigen::Matrix<Mat, Eigen::Dynamic, 1> > cf;
  // this is NOT a good idea to get each regions.
  for(int i = 0; !i ||
                 6 < max(in.rows(), in.cols()) / 2 * 6 *
                       stp / pow(2, i) / 18 / stp; i ++) {
    Vec p0(3), p1(3);
    p0[0]  = 0;
    p0[1]  = 0;
    p0[2]  = 0;
    p1[0]  = max(in.rows(), in.cols()) / 2 * 6 * stp / pow(2, i);
    p1[1]  = 0;
    p1[2]  = 0;
    auto cfv(prepareLineAxis(p0, p1, p1[0], 1, p1[0] / 18 / stp));
    cf.push_back(cfv);
  }
  vector<Eigen::Matrix<int, Eigen::Dynamic, 1> > corners;
  for(int i = 0; i < 5; i ++)
    corners.push_back(Eigen::Matrix<int, Eigen::Dynamic, 1>(2));
  corners[0][0] = corners[1][0] = - in.rows() / 2;
  corners[0][1] = corners[1][1] = - in.cols() / 2;
  corners[2][0] = corners[3][0] = in.rows() / 2;
  corners[1][1] = corners[3][1] = in.cols() / 2;
  corners[4][0] = corners[4][1] = 0;
  T avg0(0);
  bumpp = Mat(in * T(0));
  for(int ic = 0; ic < corners.size(); ic ++) {
    Mat bumps(in * T(0));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 1; i <= sqrt(T(2)) * max(bumps.rows(), bumps.cols()); i ++) {
      Vec work(max(bumps.rows(), bumps.cols()) * 6 * 3);
      const int i0(work.size() / 6 / 3);
      const int i1(i0 * 6);
      if(i0 <= 0)
        break;
      for(int j = 0; j < work.size(); j ++)
        work[j] = - T(1);
      for(int j = 0; j < i; j ++) {
        work[(j + i * 0) * i0 / i] = getImgPt(in, corners[ic][0] + in.rows() / 2 - i,     corners[ic][1] + in.cols() / 2 + i / 2 - j);
        work[(j + i * 1) * i0 / i] = getImgPt(in, corners[ic][0] + in.rows() / 2 - i + j, corners[ic][1] + in.cols() / 2 - i / 2 - j / 2);
        work[(j + i * 2) * i0 / i] = getImgPt(in, corners[ic][0] + in.rows() / 2     + j, corners[ic][1] + in.cols() / 2 - i     + j / 2);
        work[(j + i * 3) * i0 / i] = getImgPt(in, corners[ic][0] + in.rows() / 2 + i,     corners[ic][1] + in.cols() / 2 - i / 2 + j);
        work[(j + i * 4) * i0 / i] = getImgPt(in, corners[ic][0] + in.rows() / 2 + i - j, corners[ic][1] + in.cols() / 2 + i / 2 + j / 2);
        work[(j + i * 5) * i0 / i] = getImgPt(in, corners[ic][0] + in.rows() / 2     - j, corners[ic][1] + in.cols() / 2 + i     - j / 2);
      }
      for(int j = i1; j < work.size(); j ++)
        work[j] = work[j % i1];
      work = complementLine(work);
      Vec res(work.size());
      for(int j = 0; j < res.size(); j ++)
        res[j] = T(0);
      for(int j = 0; j < cf.size(); j ++) {
        Vec lwork((work.size() + pow(2, j) - 1) / pow(2, j));
        for(int k = 0; k < lwork.size(); k ++) {
          int cnt(0);
          for(int l = 0; l < pow(2, j) && k * pow(2, j) + l < work.size(); l ++) {
            lwork[k] += work[k * pow(2, j) + l];
            cnt ++;
          }
          if(cnt)
            lwork[k] /= cnt;
          else
            lwork[k]  = T(0);
        }
        lwork = getPseudoBumpSub(lwork, cf[j]);
        Vec buf(work.size());
        for(int k = 0; k < lwork.size(); k ++)
          for(int l = 0; l < pow(2, j) && k * pow(2, j) + l < work.size(); l ++)
            buf[k * pow(2, j) + l] = l == pow(2, j) / 2 ? lwork[k] : - T(1);
        res += complementLine(buf);
      }
      work = res;
      if(avg0 == T(0)) {
        for(int j = 0; j < work.size(); j ++)
          avg0 += work[j];
        avg0 /= work.size();
      } else {
        T avg1(0);
        for(int j = 0; j < work.size(); j ++)
          avg1 += work[j];
        avg1 /= work.size();
        for(int j = 0; j < work.size(); j ++)
          work[j] += avg0 - avg1;
      }
      setImgPt(bumps, corners[ic][0] + in.rows() / 2 - i, corners[ic][1] + in.cols() / 2 + i / 2 + 1, work[(0 +     i * 0) * i0 / i + i1]);
      setImgPt(bumps, corners[ic][0] + in.rows() / 2 + i, corners[ic][1] + in.cols() / 2 - i / 2 - 1, work[(i - 1 + i * 2) * i0 / i + i1]);
      for(int j = 0; j < i; j ++) {
        setImgPt(bumps, corners[ic][0] + in.rows() / 2 - i,     corners[ic][1] + in.cols() / 2 + i / 2 - j    , work[(j + i * 0) * i0 / i + i1]);
        setImgPt(bumps, corners[ic][0] + in.rows() / 2 - i + j, corners[ic][1] + in.cols() / 2 - i / 2 - j / 2, work[(j + i * 1) * i0 / i + i1]);
        // XXX
        setImgPt(bumps, corners[ic][0] + in.rows() / 2 - i + j, corners[ic][1] + in.cols() / 2 - i / 2 - j / 2 + 1, work[(j + i * 1) * i0 / i + i1]);
        setImgPt(bumps, corners[ic][0] + in.rows() / 2     + j, corners[ic][1] + in.cols() / 2 - i     + j / 2, work[(j + i * 2) * i0 / i + i1]);
        setImgPt(bumps, corners[ic][0] + in.rows() / 2 + i,     corners[ic][1] + in.cols() / 2 - i / 2 + j    , work[(j + i * 3) * i0 / i + i1]);
        setImgPt(bumps, corners[ic][0] + in.rows() / 2 + i - j, corners[ic][1] + in.cols() / 2 + i / 2 + j / 2, work[(j + i * 4) * i0 / i + i1]);
        // XXX
        setImgPt(bumps, corners[ic][0] + in.rows() / 2 + i - j, corners[ic][1] + in.cols() / 2 + i / 2 + j / 2 + 1, work[(j + i * 4) * i0 / i + i1]);
        setImgPt(bumps, corners[ic][0] + in.rows() / 2     - j, corners[ic][1] + in.cols() / 2 + i     - j / 2, work[(j + i * 5) * i0 / i + i1]);
      }
    }
    bumpp -= bumps;
  }
  autoLevel(bumpp);
  lowFreq<T> lf;
  geoms = lf.getLowFreq(bumpp, vmax);
  vector<Eigen::Matrix<T, 3, 1> > gwork;
  for(int i = 0; i < geoms.size(); i ++)
    if(geoms[i][2] != T(0))
      gwork.push_back(geoms[i]);
  geoms = gwork;
  vector<int> idxs;
  for(int i = 0; i < geoms.size(); i ++)
    idxs.push_back(i);
  delaunay = loadBumpSimpleMesh<T>(geoms, idxs);
  Mat result(bumpp.rows(), bumpp.cols());
  for(int i = 0; i < result.rows(); i ++)
    for(int j = 0; j < result.cols(); j ++)
      result(i, j) = T(0);
  tilter<T> tilt;
  for(int i = 0; i < delaunay.size(); i ++) {
    const Vec3& p0(geoms[delaunay[i][0]]);
    const Vec3& p1(geoms[delaunay[i][1]]);
    const Vec3& p2(geoms[delaunay[i][2]]);
    for(int j = 0; j < result.rows(); j ++)
      for(int k = 0; k < result.cols(); k ++) {
        Vec3 q;
        q[0] = j;
        q[1] = k;
        q[2] = 0;
        if(tilt.sameSide2(p0, p1, p2, q, true, T(.5)) &&
           tilt.sameSide2(p1, p2, p0, q, true, T(.5)) &&
           tilt.sameSide2(p2, p0, p1, q, true, T(.5)))
          result(j, k) = (geoms[delaunay[i][0]][2] + 
                          geoms[delaunay[i][1]][2] +
                          geoms[delaunay[i][2]][2]) / T(3);
      }
  }
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

template <typename T> Eigen::Matrix<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Dynamic, 1> PseudoBump<T>::prepareLineAxis(const Vec& p0, const Vec& p1, const int& hh, const int& ww, const int& rstp) {
  // N.B. ray is from infinite far, so same side of these.
  const T cdist( 1);
  const T rdist0(0);
  const T rdist1(1);
  
  Vec dir(p1 - p0);
  dir /= sqrt(dir.dot(dir));
  Vec camera(3);
  camera[0] = T(0);
  camera[1] = T(0);
  camera[2] = cdist;
  const Vec center(indiv(p0, p1, T(1) / T(2)));
  
  Eigen::Matrix<Mat, Eigen::Dynamic, 1> result(z_max);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int zi = 0; zi < z_max; zi ++) {
    result[zi] = Mat(3, stp);
    for(int s = 0; s < stp; s ++) {
      Vec cpoint(3);
      cpoint[0] = (s / (stp - 1.) - 1. / 2.) * stp * rstp + center[0];
      cpoint[1] = (s / (stp - 1.) - 1. / 2.) * stp * rstp + center[1];
      cpoint[2] = rdist0 + (zi + 1) / T(z_max) * (rdist1 - rdist0);
      result[zi].col(s) = getLineAxis(cpoint, camera, T(ww), T(hh)) - center;
    }
  }
  return result;
}

template <typename T> T PseudoBump<T>::getImgPt(const Vec& img, const T& y) {
  const int& h(img.size());
  const int  yy(abs((int(y + .5) + 3 * h) % (2 * h) - h) % h);
  return img[yy];
}

template <typename T> T PseudoBump<T>::getImgPt(const Mat& img, const int& y, const int& x) {
  const int& h(img.rows());
  const int& w(img.cols());
  const int  xx(abs((x + 3 * w) % (2 * w) - w) % w);
  const int  yy(abs((y + 3 * h) % (2 * h) - h) % h);
  return img(yy, xx);
}

template <typename T> void PseudoBump<T>::setImgPt(Mat& img, const int& y, const int& x, const T& v) {
  const int& h(img.rows());
  const int& w(img.cols());
  if(0 <= x && x < w && 0 <= y && y < h)
    img(y, x) = v;
  return;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> PseudoBump<T>::indiv(const Vec& p0, const Vec& p1, const T& pt) {
  return (p1 - p0) * pt + p0;
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

