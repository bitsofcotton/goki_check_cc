#if !defined(_2D3D_PSEUDO_)

#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "enlarge.hh"
#include "tilt.hh"
#include "edgedetect.hh"
#include "obj2vector.hh"

using std::complex;
using std::abs;
using std::pow;
using std::vector;
using std::sort;
using std::max;
using std::min;
using std::cos;
using std::sqrt;
using std::log;
using std::isfinite;

template <typename T> int cmpbump(const Eigen::Matrix<T, 3, 1>& x0, const Eigen::Matrix<T, 3, 1>& x1) {
  return x0[2] < x1[2];
}

template <typename T> class PseudoBump {
public:
  typedef complex<T> U;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> MatU;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1>              Vec;
  typedef Eigen::Matrix<U, Eigen::Dynamic, 1>              VecU;
  typedef Eigen::Matrix<T, 2, 1> Vec2;
  typedef Eigen::Matrix<T, 3, 1> Vec3;
  
  PseudoBump();
  ~PseudoBump();
  void initialize(const int& z_max, const int& stp, const int& crowd, const int& vmax, const int& nloop, const int& ndiv, const T& rdist);
  Mat  getPseudoBumpLoop(const Mat& work, const bool& y_only = false, const bool& alevel = true);
  Mat  getPseudoBumpVec(const Mat& input, vector<Vec3>& points, vector<Eigen::Matrix<int, 3, 1> >& delaunays, Mat& bumps, const bool& y_only = false);
  
  int vmax;
  
private:
  void correctSubPseudoBump(Mat& result, const Mat& workl, const Mat& worku, const int& j, const int& k, const T& ry, const T& t0, const T& t1);
  Mat  getPseudoBumpSub(const Mat& work, const int& rstp);
private:
  T    sgn(const T& x);
  Vec  getLineAxis(Vec p, Vec c, const int& w, const int& h);
  Eigen::Matrix<Mat, Eigen::Dynamic, Eigen::Dynamic> prepareLineAxis(const Vec& p0, const Vec& p1, const int& z0, const int& rstp);
  T    getImgPt(const Mat& img, const T& y, const T& x);
  Vec  indiv(const Vec& p0, const Vec& p1, const T& pt);
  void complement(const Mat& in, const int& crowd, const int& vmax, vector<Vec3>& geoms, vector<Eigen::Matrix<int, 3, 1> >& delaunay);
  Vec  complementLine(const Vec& line, const T& rratio = T(.5));
  
  int    z_max;
  int    stp;
  int    rstp;
  int    rrstp;
  T      roff;
  T      cdist;
  T      rdist;
  T      cutoff;
  T      corrnl;
  int    crowd;
  int    stpmax;
  int    nloop;
  int    ndiv;
  
  int ww;
  int hh;
  T   Pi;
};

template <typename T> PseudoBump<T>::PseudoBump() {
  initialize(20, 25, 25, 800, 8, 12, T(2.));
}

template <typename T> PseudoBump<T>::~PseudoBump() {
  ;
}

template <typename T> void PseudoBump<T>::initialize(const int& z_max, const int& stp, const int& crowd, const int& vmax, const int& nloop, const int& ndiv, const T& rdist) {
  this->Pi       = T(4) * atan2(T(1), T(1));
  this->z_max    = z_max;
  this->stp      = stp;
  this->rstp     = stp;
  this->rrstp    = 2;
  assert(1 < this->rrstp);
  this->crowd    = crowd;
  this->stpmax   = int(log(T(crowd) / (rstp / 2) / 2) / log(T(rrstp))) + 1;
  assert(1 <= this->stpmax);
  this->vmax     = vmax;
  this->nloop    = nloop;
  this->ndiv     = ndiv;
  assert(nloop / T(ndiv) < T(1));
  this->roff     = 1. / 6.;
  this->cdist    = - 1.;
  this->rdist    = rdist * (- this->cdist);
  this->corrnl   = T(.8);
  return;
};

template <typename T> T PseudoBump<T>::sgn(const T& x) {
  if(x < T(0))
    return - T(1);
  if(x > T(0))
    return T(1);
  return T(0);
}

// correct sub bump with tilted sub bump.
template <typename T> void PseudoBump<T>::correctSubPseudoBump(Mat& result, const Mat& workl, const Mat& worku, const int& j, const int& k, const T& ry, const T& t0, const T& t1) {
  const int j1((result.rows() / 2 - 1 -  j     ) / ry + result.rows() / 2);
  const int j0((result.rows() / 2 - 1 - (j + 1)) / ry + result.rows() / 2);
  const int jj1(max(int(result.rows() - 1 -  j),      0));
  const int jj0(max(int(result.rows() - 1 - (j + 1)), 0));
  const T x0(max(worku(jj0, k) * t1, workl(jj0, k) * t0));
  const T x1(max(worku(jj1, k) * t1, workl(jj1, k) * t0));
  if(x0 < T(0) || x1 < T(0))
    return;
  for(int l =  max(0, min(j0, int(result.rows() - 1)));
          l <= max(0, min(j1, int(result.rows() - 1)));
          l ++) {
    const T t(j1 == j0 ? T(1) : (l - j0) * (j1 - j0));
    const T buf((x0 * t + x1 * (T(1) - t)) * ry);
    if(result(l, k) < T(0))
      result(l, k) = min(max(buf, T(0)), T(1));
  }
  return;
}

// bump it with abs(d/dt local color) / abs(d/dt color) >> 0.
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
    cerr << ".";
    cerr.flush();
    for(int j = 0; j < result.rows(); j ++)
      result(j, i) = - T(1);
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
        Vec cc(c.size() / 4 * 2 + 1);
        for(int u =  c.size() / 2 - c.size() / 4;
                u <= c.size() / 2 + c.size() / 4;
                u ++)
          cc[u - c.size() / 2 + c.size() / 4] = c[u];
        // N.B. simply take the ratio of local and far difference on two eyes.
        const T n2(cc.dot(cc) / c.dot(c));
        if(isfinite(n2) && (n2 < zval[s] || zz == 0)) {
          // N.B. If increase zz, increase distance from camera.
          //      And, zz->0 is treated as distant in tilter.
          result(s, i) = T(1) - zz / T(lrf.rows());
          zval[s]      = n2;
        }
      }
    }
  }
  return result;
}

// overall bump it, complement with sub sub bump.
template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::getPseudoBumpLoop(const Mat& input, const bool& y_only, const bool& alevel) {
  Mat result;
  if(y_only) {
    result = getPseudoBumpSub(input, rstp);
    for(int i = 0; i < nloop; i ++ ) {
      tilter<T> tilter;
      Mat work(result);
      for(int j = 0; j < work.cols(); j ++)
        work.col(j) = complementLine(work.col(j));
      const Mat worku(getPseudoBumpSub(tilter.tilt(input, work, 1, 4, (i + 1) / T(ndiv) + T(1)), rstp));
      const Mat workl(getPseudoBumpSub(tilter.tilt(input, work, 3, 4, (i + 1) / T(ndiv) + T(1)), rstp));
      const T ry(cos(Pi / T(2) * (i + 1) / T(ndiv)));
      const T dy(result.rows() / 2 - result.rows() / 2 * ry);
      const T t0(abs(T(1) - sqrt(T(2) * (T(1) - ry))));
      const T t1(abs(T(1) + sqrt(T(2) * (T(1) - ry))));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int j = dy; j < int(result.rows()) - int(dy) + 1; j ++)
        for(int k = 0; k < result.cols(); k ++)
          correctSubPseudoBump(result, workl, worku, j, k, ry, t0, t1);
    }
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < result.cols(); i ++)
      result.col(i) = complementLine(result.col(i));
#if defined(_OPENMP)
#prgma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < result.rows(); i ++)
      result.row(i) = complementLine(result.row(i));
  } else {
    const Mat resx(getPseudoBumpLoop(input.transpose(), true, false).transpose());
    const Mat resy(getPseudoBumpLoop(input,             true, false));
    result = (resx + resy) / T(2);
  }
  if(alevel)
    autoLevel<T>(&result, 1);
  return result;
}

// bump with multiple scale and merge, then, get vectored result.
template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::getPseudoBumpVec(const Mat& input, vector<Vec3>& points, vector<Eigen::Matrix<int, 3, 1> >& delaunays, Mat& bumps, const bool& y_only) {
  bumps = getPseudoBumpLoop(input, y_only);
  complement(bumps, crowd, vmax, points, delaunays);
  Mat result(input.rows(), input.cols());
  for(int i = 0; i < result.rows(); i ++)
    for(int j = 0; j < result.cols(); j ++)
      result(i, j) = T(0);
  const vector<Eigen::Matrix<int, 3, 1> >& delaunay(delaunays);
  const vector<Vec3>& geoms(points);
  const Mat& in(input);
        tilter<T> tilt;
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
  autoLevel<T>(&result, 1);
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
#pragma omp parallel for schedule(static, 1)
#endif
  for(int zi = 0; zi < z0; zi ++) {
    result(zi, 0) = Mat(3, stp);
    result(zi, 1) = Mat(3, stp);
    for(int s = 0; s < stp; s ++) {
      Vec cpoint(3);
      cpoint[0] = (s / (stp - 1.) - 1. / 2.) * rstp + center[0];
      cpoint[1] = (s / (stp - 1.) - 1. / 2.) * rstp + center[1];
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
  const int  xx(abs((int(x + .5) + 3 * w) % (2 * w) - w) % w);
  const int  yy(abs((int(y + .5) + 3 * h) % (2 * h) - h) % h);
  return img(yy, xx);
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> PseudoBump<T>::indiv(const Vec& p0, const Vec& p1, const T& pt) {
  return (p1 - p0) * pt + p0;
}

template <typename T> void PseudoBump<T>::complement(const Mat& in, const int& crowd, const int& vmax, vector<Vec3>& geoms, vector<Eigen::Matrix<int, 3, 1> >& delaunay) {
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
  if(vmax < geoms.size()) {
    const int diff(geoms.size() - vmax);
    geoms.erase(geoms.end() - diff / 2, geoms.end());
    geoms.erase(geoms.begin(), geoms.begin() + diff / 2);
  }
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
    result[i] = max(min(result[i], T(1)), T(0));
  }
  return result;
}

#define _2D3D_PSEUDO_
#endif

