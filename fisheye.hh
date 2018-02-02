#if !defined(_2D3D_PSEUDO_)

#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Dense>

using std::complex;
using std::abs;
using std::vector;
using std::sort;
using std::pair;
using std::make_pair;
using std::max;
using std::min;
using std::cos;
using std::sqrt;
using std::log;
using std::isfinite;
using std::cerr;
using std::flush;

template <typename T> class PseudoBump {
public:
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1>              Vec;
  typedef Eigen::Matrix<T,   3, 1> Vec3;
  typedef Eigen::Matrix<int, 2, 1> Vec2i;
  
  PseudoBump();
  ~PseudoBump();
  void  initialize(const int& z_max, const int& stp, const T& rthresh);
  Mat   getPseudoBumpVec(const Mat& in, vector<Vec3>& geoms, vector<Eigen::Matrix<int, 3, 1> >& delaunay);
  
  int vbox;
  
private:
  Vec   complementLine(const Vec& line, const T& rratio = T(.5), const int& guard = int(1));
  Mat   autoLevel(const Mat& data, int npad = - 1);
  Vec   getPseudoBumpSub(const Vec& work, const Eigen::Matrix<Mat, Eigen::Dynamic, 1>& cf, const int& lower = 0, int upper = - 1);
  Mat   getPseudoBump(const Mat& in);
  Vec   getLineAxis(Vec p, Vec c);
  Eigen::Matrix<Mat, Eigen::Dynamic, 1> prepareLineAxis(const Vec& d, const int& rstp);
  T     getImgPt(const Vec& img, const T& y);
  T     getImgPt(const Mat& img, const int& y, const int& x);
  void  setImgPt(Mat& img, const int& y, const int& x, const T& v);
  Vec   indiv(const Vec& p0, const Vec& p1, const T& pt);
  Vec2i getHexGeom(const Vec2i& center, const int& idx, const int& Midx);
  
  int z_max;
  int stp;
  T   rthresh;
  int lloop;
  T   cdist;
  T   zdist;
  T   rz;
  
  T   Pi;
  Mat Dop;
};

template <typename T> PseudoBump<T>::PseudoBump() {
  initialize(20, 15, 2.);
}

template <typename T> PseudoBump<T>::~PseudoBump() {
  ;
}

template <typename T> void PseudoBump<T>::initialize(const int& z_max, const int& stp, const T& rthresh) {
  this->z_max   = z_max;
  this->stp     = stp;
  this->rthresh = rthresh;
  this->lloop   = 4;
  this->vbox    = 1;
  this->cdist   = T(1);
  this->zdist   = T(1);
  this->rz      = T(1) / T(6);
  this->Pi      = T(4) * atan2(T(1), T(1));
  Eigen::Matrix<complex<T>, Eigen::Dynamic, Eigen::Dynamic> DFT(stp, stp), IDFT(stp, stp);
  for(int i = 0; i < DFT.rows(); i ++)
    for(int j = 0; j < DFT.cols(); j ++) {
      DFT( i, j) = exp(complex<T>(- 2.) * Pi * sqrt(complex<T>(- 1)) * complex<T>(i * j / T(stp)));
      IDFT(i, j) = exp(complex<T>(  2.) * Pi * sqrt(complex<T>(- 1)) * complex<T>(i * j / T(stp))) / T(stp);
    }
  for(int i = 0; i < DFT.rows(); i ++)
    DFT.row(i) *= complex<T>(2.) * Pi * sqrt(complex<T>(- 1)) * T(i) / T(DFT.rows()) * T(8);
  Dop  = (IDFT * DFT).real();
  return;
};

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> PseudoBump<T>::getPseudoBumpSub(const Vec& work, const Eigen::Matrix<Mat, Eigen::Dynamic, 1>& cf, const int& lower, int upper) {
  cerr << "." << flush;
  assert(0 < upper);
  assert(lower < upper);
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
  for(int s = lower; s < upper; s ++) {
    Vec pt(indiv(p0, p1, s / T(result.size())));
    for(int zz = 0; zz < cf.size(); zz ++) {
      // d/dt (local color):
      Vec c(cf[zz].cols());
      for(int u = 0; u < c.size(); u ++)
        c[u] = getImgPt(workv, cf[zz](0, u) + pt[0] + stp / 2);
      Vec rc(c);
      for(int u = 0; u < c.size(); u ++)
        rc[u] = c[c.size() - 1 - u];
      const T n2(abs(Dop.row( c.size() / 2).dot( c)) +
                 abs(Dop.row(rc.size() / 2).dot(rc)) );
      // ||local color|| >> ||near local color||
      // for suppress pseudo detect.
      T ci(0), co(0);
      for(int u = 0; u < c.size() / 3; u ++)
        co += c[u] * c[u] + c[c.size() - u - 1] * c[c.size() - u - 1];
      for(int u = c.size() / 3; u < c.size(); u ++)
        ci += c[u] * c[u];
      const T r2(co / ci);
      if(isfinite(n2) && zval[s] < n2 && isfinite(r2) && r2 < rthresh) {
        // N.B. If increase zz, decrease the distance from camera.
        //      And, zz->0 is treated as distant in tilter.
        result[s] = (T(1) + zz) / T(cf.size());
        zval[s]   = n2;
      }
    }
  }
  return result;
}

// Hextile.
template <typename T> Eigen::Matrix<int, 2, 1> PseudoBump<T>::getHexGeom(const Vec2i& center, const int& idx, const int& Midx) {
  Vec2i result(center);
  const int i(Midx / 6);
  const int j(idx % i);
  switch(idx / i) {
  case 0:
    result[0] += - i;
    result[1] +=   i / 2 - j;
    break;
  case 1:
    result[0] += - i     + j;
    result[1] += - i / 2 - j / 2;
    break;
  case 2:
    result[0] +=           j;
    result[1] += - i     + j / 2;
    break;
  case 3:
    result[0] +=   i;
    result[1] += - i / 2 + j     - 1;
    break;
  case 4:
    result[0] +=   i     - j;
    result[1] +=   i / 2 + j / 2 - 1;
    break;
  case 5:
    result[0] +=         - j;
    result[1] +=   i     - j / 2 - 1;
    break;
  default:
    result[0] = result[1] = - 1;
  }
  return result;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::getPseudoBump(const Mat& in) {
  /*
   * Get hex tile based bumpmap.
   */
  Vec p0(3), p1(3);
  p0[0] = 0;
  p0[1] = 0;
  p0[2] = 0;
  p1[1] = 0;
  p1[2] = 0;
  vector<Eigen::Matrix<Mat, Eigen::Dynamic, 1> > cfs;
  for(int i = 0; i < lloop; i ++) {
    p1[0] = max(in.rows(), in.cols()) / 2 * 6 * 3;
    cfs.push_back(prepareLineAxis(p0 - p1, p1[0] / 18 * stp * (i + 1) / lloop));
  }
  // XXX fixme:
  vector<Vec2i> corners;
  corners.push_back(Vec2i(0, 0));
  
  Mat result(in.rows(), in.cols());
  for(int i = 0; i < result.rows(); i ++)
    for(int j = 0; j < result.cols(); j ++)
      result(i, j) = T(0);
  for(int ic = 0; ic < corners.size(); ic ++) {
    Mat bumps(in * T(0));
    Vec work(int(sqrt(2) * max(bumps.rows(), bumps.cols()) * 6 * 3));
    Vec bwork(work.size());
    for(int i = 0; i < bwork.size(); i ++)
      bwork[i] = T(0);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 1; i < sqrt(T(2)) * max(bumps.rows(), bumps.cols()); i ++) {
      const int i0(work.size() / 6 / 3);
      const int i1(i0 * 6);
      if(i0 <= 0)
        break;
      for(int j = 0; j < work.size(); j ++)
        work[j] = - T(1);
      int ii(- 1), jj(- 1);
      for(int j = 0; j < i * 6; j ++) {
        Vec2i geom(getHexGeom(corners[ic], j, i * 6));
        if(0 <= geom[0] && geom[0] < in.rows() &&
           0 <= geom[1] && geom[1] < in.cols()) {
          if(ii < 0)
            ii = j;
          jj = j;
        }
        work[int(j * i1 / i / 6)] = getImgPt(in, geom[0], geom[1]);
      }
      if(jj < 0 || ii < 0) continue;
      for(int j = i1; j < work.size(); j ++)
        work[j] = work[j % i1];
      work = complementLine(work);
      Vec res(work * T(0));
      for(int k = 0; k < cfs.size(); k ++) {
        const Vec local(complementLine(getPseudoBumpSub(work, cfs[k], i1, i1 * 2)) * (k + 1));
        if(bwork.dot(local) >= T(0)) {
          bwork = local;
          res  += local;
        } else
          res  -= local;
      }
      for(int j = 0; j < i * 6; j ++) {
        Vec2i geom(getHexGeom(corners[ic], j, i * 6));
        T buf(0);
        for(int k = i1 + int(j * i1 / i / 6); k < i1 + int((j + 1) * i1 / i / 6); k ++)
          buf += res[k];
        buf /= int((j + 1) * i1 / i / 6) - int(j * i1 / i / 6);
        setImgPt(bumps, geom[0], geom[1], buf);
        if(j / i == 1)
          setImgPt(bumps, geom[0], geom[1] - 1, buf);
        if(j / i == 4)
          setImgPt(bumps, geom[0], geom[1] + 1, buf);
      }
    }
    result += bumps;
  }
  return result;
}

// get bump with multiple scale and vectored result.
template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::getPseudoBumpVec(const Mat& in, vector<Vec3>& geoms, vector<Eigen::Matrix<int, 3, 1> >& delaunay) {
  cerr << "bump" << flush;
  Mat result(autoLevel(getPseudoBump(autoLevel(in))));
  
  /*
   * Get vector based bumps.
   */
  geoms = vector<Eigen::Matrix<T, 3, 1> >();
  T aavg(0);
  for(int i = 0; i < result.rows(); i ++)
    for(int j = 0; j < result.cols(); j ++)
      aavg += result(i, j);
  aavg /= result.rows() * result.cols();
  for(int i = 0; i < result.rows() / vbox + 1; i ++)
    for(int j = 0; j < result.cols() / vbox + 1; j ++) {
      T   avg(0);
      int cnt(0);
      for(int ii = i * vbox; ii < min((i + 1) * vbox, int(result.rows())); ii ++)
        for(int jj = j * vbox; jj < min((j + 1) * vbox, int(result.cols())); jj ++) {
          avg += result(ii, jj);
          cnt ++;
        }
      if(! cnt)
        geoms.push_back(geoms[geoms.size() - 1]);
      else {
        Eigen::Matrix<T, 3, 1> work;
        work[0]  = i * vbox;
        work[1]  = j * vbox;
        // XXX fixme: sign.
        work[2]  = (avg / cnt - aavg) / T(2) * sqrt(T(result.rows() * result.cols())) * rz;
        geoms.push_back(work);
      }
    }
  Eigen::Matrix<T, 3, 1> avg;
  avg[0] = avg[1] = avg[2] = T(0);
  for(int i = 0; i < geoms.size(); i ++)
    avg += geoms[i];
  avg /= geoms.size();
  for(int i = 0; i < geoms.size(); i ++)
    geoms[i][2] -= avg[2];
  delaunay = vector<Eigen::Matrix<int, 3, 1> >();
  for(int i = 1; i < result.rows() / vbox + 1; i ++)
    for(int j = 0; j < result.cols() / vbox; j ++) {
      Eigen::Matrix<int, 3, 1> work, work2;
      work[0]  = (i - 1) * (result.cols() / vbox + 1) + j;
      work[1]  =  i      * (result.cols() / vbox + 1) + j;
      work[2]  =  i      * (result.cols() / vbox + 1) + j + 1;
      work2[0] = (i - 1) * (result.cols() / vbox + 1) + j;
      work2[2] = (i - 1) * (result.cols() / vbox + 1) + j + 1;
      work2[1] =  i      * (result.cols() / vbox + 1) + j + 1;
      delaunay.push_back(work);
      delaunay.push_back(work2);
    }
  return result;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> PseudoBump<T>::getLineAxis(Vec p, Vec c) {
  // suppose straight ray from infinitely distant 2 x 2 size.
  // camera geometry with c, lookup on p[2] z-distance
  // virtual images.
  // virtually : z = p[2], p[2] = 0.
  // <c + (p - c) * t, [0, 0, 1]> = z
  T   t((p[2] - c[2]) / (- c[2]));
  Vec work((p - c) * t + c);
  return work;
}

template <typename T> Eigen::Matrix<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Dynamic, 1> PseudoBump<T>::prepareLineAxis(const Vec& d, const int& rstp) {
  // N.B. ray is from infinite far, so same side of these.
  const Vec dir(d / sqrt(d.dot(d)));
  Vec camera(3);
  camera[0] = T(0);
  camera[1] = T(0);
  camera[2] = cdist;
  
  Eigen::Matrix<Mat, Eigen::Dynamic, 1> result(z_max);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int zi = 0; zi < z_max; zi ++) {
    result[zi] = Mat(3, stp);
    for(int s = 0; s < stp; s ++) {
      Vec cpoint(3);
      cpoint[0] = (s / (stp - 1.) - 1. / 2.) * stp * dir[0];
      cpoint[1] = (s / (stp - 1.) - 1. / 2.) * stp * dir[1];
      cpoint[2] = (zi + 1) / T(z_max + 1) * zdist;
      result[zi].col(s) = getLineAxis(cpoint, camera) * rstp;
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

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> PseudoBump<T>::complementLine(const Vec& line, const T& rratio, const int& guard) {
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
  int i;
  for(i = 0; i < line.size(); i ++) {
    if(result[i] >= T(0))
      continue;
    for(; rng[1] < ptsi.size() - 2 && ptsi[rng[1]] < i; rng[1] ++) ;
    rng[0] = rng[2] = rng[1];
    while(0 < rng[0] && ptsi[rng[1]] - ptsi[rng[0]] <= guard)
      rng[0] --;
    while(rng[2] < ptsi.size() - 1 && ptsi[rng[2]] - ptsi[rng[1]] <= guard)
      rng[2] ++;
    if(rng[0] == rng[1]) rng[0] --;
    if(rng[2] == rng[1]) rng[2] ++;
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
    for(; i < line.size() && result[i] < T(0); i ++) {
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
  }
  return result;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::autoLevel(const Mat& data, int npad) {
  if(npad <= 0)
    npad = 0;
  vector<T> stat;
  for(int j = 0; j < data.rows(); j ++)
    for(int k = 0; k < data.cols(); k ++)
      stat.push_back(data(j, k));
  sort(stat.begin(), stat.end());
  const T mm(stat[npad]);
  T MM(stat[stat.size() - 1 - npad]);
  if(MM == mm)
    MM = mm + 1.;
  Mat result(data.rows(), data.cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < data.rows(); i ++)
    for(int j = 0; j < data.cols(); j ++)
      result(i, j) = (max(min(data(i, j), MM), mm) - mm) / (MM - mm);
  return result;
}

#define _2D3D_PSEUDO_
#endif

