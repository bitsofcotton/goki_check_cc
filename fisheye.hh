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
using std::vector;
using std::sort;
using std::max;
using std::min;
using std::cos;
using std::sqrt;
using std::isfinite;

template <typename T> int cmpbump(const Eigen::Matrix<T, 3, 1>& x0, const Eigen::Matrix<T, 3, 1>& x1) {
  return x0[2] < x1[2];
}

template <typename T> class PseudoBump {
public:
  int    z_max;
  int    stp;
  int    rstp;
  int    rrstp;
  T      roff;
  T      cdist;
  T      rdist;
  T      cutoff;
  T      cutz;
  T      cthresh;
  int    crowd;
  int    vmax;
  int    nloop;
  int    ndiv;
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
  
  Mat  getPseudoBumpSub(const Mat& work, const int& rstp);
  Mat  getPseudoBumpLoop(const Mat& work);
  void getPseudoBumpVecSub(const Mat& input, vector<Vec3>& points, vector<Eigen::Matrix<int, 3, 1> >& delaunays, Mat& result, const int& vmax, const bool& y_only = false);
  Mat  getPseudoBumpVec(const Mat& input, vector<Vec3>& points, vector<Eigen::Matrix<int, 3, 1> >& delaunays, Mat& bumps, const bool& y_only = false);
private:
  T sgn(const T& x);
  Vec getLineAxis(Vec p, Vec c, const int& w, const int& h);
  Eigen::Matrix<Mat, Eigen::Dynamic, Eigen::Dynamic> prepareLineAxis(const Vec& p0, const Vec& p1, const int& z0, const int& rstp);
  T    getImgPt(const Mat& img, const T& y, const T& x);
  Vec  indiv(const Vec& p0, const Vec& p1, const T& pt);
  void complement(const Mat& in, const int& crowd, const int& vmax, vector<Vec3>& geoms, vector<Eigen::Matrix<int, 3, 1> >& delaunay);
  Vec  complementLine(const Vec& line, const T& rratio = T(.8));
  Mat  shrink(const Mat& in);
  
  int ww;
  int hh;
  Mat Dop;
  T   Pi;
};

template <typename T> PseudoBump<T>::PseudoBump() {
  initialize(20, 12, 16, 800, 2, 6, T(2.));
}

template <typename T> PseudoBump<T>::~PseudoBump() {
  ;
}

template <typename T> void PseudoBump<T>::initialize(const int& z_max, const int& stp, const int& crowd, const int& vmax, const int& nloop, const int& ndiv, const T& rdist) {
  this->Pi      = 4. * atan2(T(1.), T(1.));
  this->z_max   = z_max;
  this->stp     = stp;
  this->rstp    = stp;
  this->rrstp   = rstp / 4;
  assert(1 < this->rrstp);
  this->crowd   = crowd;
  this->vmax    = vmax;
  this->nloop   = nloop;
  this->ndiv    = ndiv;
  assert(nloop / T(ndiv) < T(1));
  this->roff    = 1. / 6.;
  this->cdist   = - 1.;
  this->rdist   = rdist * (- this->cdist);
  this->cutoff  = T(.5);
  this->cutz    = T(.0125);
  this->cthresh = T(0);
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
    for(int s = 0; s < result.rows(); s ++)
      if(result(s, i) <= cutz || T(1) - cutz <= result(s, i))
        result(s, i) = - T(1);
  }
  return result;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::getPseudoBumpLoop(const Mat& input) {
  Mat result(getPseudoBumpSub(input, rstp));
  for(int i = 0; i < nloop; i ++ ) {
    tilter<T> tilter;
    Mat work(result);
    Mat workc(work);
    for(int j = 0; j < workc.cols(); j ++)
      workc.col(j) = complementLine(workc.col(j));
    const Mat worku(getPseudoBumpSub(tilter.tilt(workc, work, 1, 4, (i + 1) / T(ndiv) + T(1)), rstp));
    const Mat workl(getPseudoBumpSub(tilter.tilt(workc, work, 3, 4, (i + 1) / T(ndiv) + T(1)), rstp));
    const T ry(cos(Pi / T(2) * (i + 1) / T(ndiv)));
    const T dy(result.rows() / 2 - result.rows() / 2 * ry);
    // XXX : inverse of t?
    const T t(T(1) + sqrt(T(2) * (T(1) - ry)));
    for(int j = dy; j < result.rows() - dy + 1; j ++)
      for(int k = 0; k < result.cols(); k ++) {
        const int j1((result.rows() / 2 - 1 -  j     ) / ry + result.rows() / 2);
        const int j0((result.rows() / 2 - 1 - (j + 1)) / ry + result.rows() / 2);
        const int jc1(min(max(j1, 0), int(result.rows() - 1)));
        const int jc0(min(max(j0, 0), int(result.rows() - 1)));
        const int jj1(result.rows() - 1 - j);
        const int jj0(min(result.rows() - 1 - (j + 1), result.rows() - 1));
        T x0(worku(jj0, k) / t);
        if(x0 < T(0) || (T(0) <= workl(jj0, k) * t && workl(jj0, k) * t < x0))
          x0 = t * workl(jj0, k);
        T x1(worku(jj1, k) / t);
        if(x1 < T(0) || (T(0) <= workl(jj1, k) * t && workl(jj1, k) * t < x1))
          x1 = t * workl(jj1, k);
        if(x0 < T(0) || x1 < T(0))
          continue;
        for(int l = jj0; l < jj1; l ++) {
          const T t((l - j0) / T(j1 - j0));
          if(result(l, k) < T(0))
            result(l, k) = max(result(l, k), (x0 * t + x1 * (T(1) - t)) * ry);
        }
      }
  }
  for(int i = 0; i < result.cols(); i ++)
    result.col(i) = complementLine(result.col(i));
  return result;
}

template <typename T> void PseudoBump<T>::getPseudoBumpVecSub(const Mat& input, vector<Vec3>& points, vector<Eigen::Matrix<int, 3, 1> >& delaunays, Mat& result, const int& vmax, const bool& y_only) {
  if(y_only)
    result = getPseudoBumpLoop(input);
  else
    result = getPseudoBumpLoop(input.transpose()).transpose() +
             getPseudoBumpLoop(input);
  Mat sute(result);
  for(int i = 0; i < sute.cols(); i ++)
    sute.col(i) = complementLine(sute.col(i));
  const T med(autoLevel<T>(&sute, 1, 0));
  for(int i = 0; i < result.rows(); i ++)
    for(int j = 0; j < result.cols(); j ++)
      if(T(1) <= sute(i, j) && T(0) <= result(i, j))
        result(i, j) = med;
  sute = Mat(result);
  for(int i = 0; i < sute.cols(); i ++)
    sute.col(i) = complementLine(sute.col(i));
  edgedetect<T> detect;
  auto edge(detect.detect(detect.detect(sute, detect.DETECT_Y), detect.COLLECT_Y));
  T avg(0);
  for(int i = 0; i < result.cols(); i ++) {
    T buf(0);
    int count(0);
    for(int j = 0; j < result.rows(); j ++)
      if(T(0) <= result(j, i)) {
        buf += result(j, i);
        count ++;
      }
    if(count)
      avg = max(avg, buf / count);
  }
  avg = sqrt(avg / result.rows());
  Mat work(result);
  for(int i = 0; i < result.rows(); i ++)
    for(int j = 0; j < result.cols(); j ++)
      if(edge(i, j) <= cutoff * avg)
        work(i, j) = - T(4);
  complement(work, crowd, vmax, points, delaunays);
  return;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::getPseudoBumpVec(const Mat& input, vector<Vec3>& points, vector<Eigen::Matrix<int, 3, 1> >& delaunays, Mat& bumps, const bool& y_only) {
  Mat mats(input);
  getPseudoBumpVecSub(mats, points, delaunays, bumps, vmax, y_only);
  bumps = Mat(- bumps);
  tilter<T> tilt;
  T ratio(1);
  while(rrstp * rrstp * 8 < min(mats.rows(), mats.cols())) {
    mats   = shrink(mats);
    ratio *= rrstp;
    vector<Vec3> workpoints;
    vector<Eigen::Matrix<int, 3, 1> > workdels;
    Mat sute;
    getPseudoBumpVecSub(mats, workpoints, workdels, sute, vmax / ratio, y_only);
    for(int i = 0; i < workdels.size(); i ++) {
      const Vec3& p0(workpoints[workdels[i][0]] * ratio);
      const Vec3& p1(workpoints[workdels[i][1]] * ratio);
      const Vec3& p2(workpoints[workdels[i][2]] * ratio);
      const Vec3  g((p0 + p1 + p2) / T(3));
      for(int j = 0; j < points.size(); j ++) {
        const Vec3& q(points[j]);
        if(tilt.sameSide2(p0, p1, p2, q) &&
           tilt.sameSide2(p1, p2, p0, q) &&
           tilt.sameSide2(p2, p0, p1, q))
          // XXX check me.
          points[j][2] += g[2] / sqrt(ratio);
      }
    }
  }
  vector<Vec3> ppoints(points);
  sort(ppoints.begin(), ppoints.end(), cmpbump<T>);
  int ii(0), jj(ppoints.size() - 1);
  for(int i = 0; i < points.size() - 1; i ++)
    if(ppoints[i][0] < T(2) || ppoints[i][1] < T(2) ||
       input.rows() - 3 < ppoints[i][0] ||
       input.cols() - 3 < ppoints[i][1])
      ii = i + 1;
    else
      break;
  for(int i = ii, i0 = ii; i < points.size() - 1; i ++)
    if(points[i0][2] == points[i][2])
      ii = i + 1;
    else
      break;
  for(int j = ppoints.size() - 1; j > 0; j --)
     if(ppoints[j][0] < T(2) || ppoints[j][1] < T(2) ||
        input.rows() - 3 < ppoints[j][0] ||
        input.cols() - 3 < ppoints[j][1])
       jj = j - 1;
     else
       break;
  for(int j = jj, j0 = jj; j > 0; j --)
    if(ppoints[j0][2] == points[j][2])
      jj = j - 1;
    else
      break;
  const T mm(ppoints[ii][2]);
  const T MM(ppoints[jj][2]);
  for(int i = 0; i < points.size(); i ++)
    points[i][2] = max(min(MM, points[i][2]), mm);
  Mat result(input.rows(), input.cols());
  for(int i = 0; i < result.rows(); i ++)
    for(int j = 0; j < result.cols(); j ++)
      result(i, j) = T(0);
  const vector<Eigen::Matrix<int, 3, 1> >& delaunay(delaunays);
  const vector<Vec3>& geoms(points);
  const Mat& in(input);
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

template <typename T> void PseudoBump<T>::complement(const Mat& in, const int& crowd, const int& vmax, vector<Vec3>& geoms, vector<Eigen::Matrix<int, 3, 1> >& delaunay) {
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
  if(vmax < geoms.size()) {
    const int diff(geoms.size() - vmax);
    const int start((geoms.size() - diff) / 2);
    const int end((geoms.size() + diff) / 2);
    geoms.erase(geoms.begin() + start, geoms.begin() + end);
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

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::shrink(const Mat& in) {
  Mat res(in.rows() / rrstp, in.cols() / rrstp);
  for(int i = 0; i < res.rows(); i ++)
    for(int j = 0; j < res.cols(); j ++) {
      res(i, j) = T(0);
      int count(0);
      for(int ii = i * rrstp; ii < min((i + 1) * rrstp, int(in.rows())); ii ++)
        for(int jj = j * rrstp; jj < min((j + 1) * rrstp, int(in.cols())); jj ++) {
          res(i, j) += in(ii, jj);
          count ++;
        }
      res(i, j) /= count;
    }
  return res;
}

#define _2D3D_PSEUDO_
#endif

