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
  T      rdist;
  T      brange;
  int    nlevel;
  T      sthresh;
  T      cthresh;
  typedef complex<T> U;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> MatU;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1>              Vec;
  typedef Eigen::Matrix<U, Eigen::Dynamic, 1>              VecU;
  
  PseudoBump();
  void initialize(const int& z_max, const int& stp, const int& nlevel, const T& cthresh);
  ~PseudoBump();
  
  Mat getPseudoBumpSub(const Mat& input, const int& rstp);
  Mat getPseudoBump(const Mat& input, const bool& y_only);
  Mat rgb2l(const Mat rgb[3]);
  Vec complementLine(const Vec& line, const T& rratio = T(.8));
private:
  T zz(const T& t);
  T sgn(const T& x);
  Vec minSquare(const Vec& input);
  Vec getLineAxis(Vec p, Vec c, const int& w, const int& h);
  Eigen::Matrix<Mat, Eigen::Dynamic, Eigen::Dynamic> prepareLineAxis(const Vec& p0, const Vec& p1, const int& z0, const int& rstp);
  T   getImgPt(const Mat& img, const T& y, const T& x);
  Vec indiv(const Vec& p0, const Vec& p1, const T& pt);
  Mat integrate(const Mat& input);
  Mat autoLevel(const Mat& input);
  
  int ww;
  int hh;
  Mat Dop;
  T   Pi;
};

template <typename T> PseudoBump<T>::PseudoBump() {
  initialize(20, 12, 64, .5);
}

template <typename T> PseudoBump<T>::~PseudoBump() {
  ;
}

template <typename T> void PseudoBump<T>::initialize(const int& z_max, const int& stp, const int& nlevel, const T& cthresh) {
  this->z_max   = z_max;
  this->stp     = stp;
  this->rstp    = stp * 2;
  this->roff    = 1. / 6.;
  this->rdist   = 2.;
  this->nlevel  = nlevel;
  this->cthresh = cthresh;
  Pi            = 4. * atan2(T(1.), T(1.));
  sthresh       = T(1e-8) / T(256);
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

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::getPseudoBumpSub(const Mat& input, const int& rstp) {
  Mat work(input);
  Mat result(work.rows(), work.cols());
  const int ppratio = (result.rows() + input.rows() - 1) / input.rows();
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
  Vec ms0(input.cols()), ms1(input.cols());
#if defined(_OPENMP)
#pragma omp for
#endif
  for(int i = 0; i < input.cols(); i ++) {
    Vec ms(minSquare(input.col(i)));
    ms0[i] = ms[0];
    ms1[i] = ms[1];
    for(int j = 0; j < input.rows(); j ++)
      centerized(j, i) = input(j, i) - (ms[0] + ms[1] * j);
  }
  centerized = Iop.real().template cast<T>() * centerized;
#if defined(_OPENMP)
#pragma omp for
#endif
  for(int i = 0; i < input.cols(); i ++)
    for(int j = 0; j < input.rows(); j ++)
      centerized(j, i) += (ms0[i] * j + ms1[i] * j * j / 2.) / input.rows();
  return centerized / input.rows();
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::autoLevel(const Mat& input) {
  Mat res(input.rows(), input.cols());
  vector<T> stat;
  for(int i = 0; i < input.rows(); i ++)
    for(int j = 0; j < input.cols(); j ++)
      stat.push_back(input(i, j));
  sort(stat.begin(), stat.end());
  const T mm(stat[stat.size() / nlevel]);
  T MM(stat[stat.size() * (nlevel - 1) / nlevel]);
  if(MM == mm)
    MM = mm + 1.;
#if defined(_OPENMP)
#pragma omp for
#endif
  for(int i = 0; i < input.rows(); i ++)
    for(int j = 0; j < input.cols(); j ++)
      res(i, j) = (max(min(input(i, j), MM), mm) - mm) / (MM - mm);
  return res;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::getPseudoBump(const Mat& input, const bool& y_only) {
  edgedetect<T> detect;
  Mat result(input);
  result *= T(0);
  if(y_only) {
    // N.B. We may need integrate on local to global operation.
    // N.B. We assume differential/integral of getPseudoBumpSub as a linear,
    //      but in fact, it's not.
    result = detect.detect(getPseudoBumpSub(integrate(input), rstp),
                           edgedetect<T>::DETECT_Y);
  } else {
    const Mat cachex(integrate(input.transpose()));
    const Mat cachey(integrate(input));
    const Mat pbx(getPseudoBumpSub(cachex, rstp).transpose());
    const Mat pby(getPseudoBumpSub(cachey, rstp));
    result = detect.detect(pbx, edgedetect<T>::DETECT_X) +
             detect.detect(pby, edgedetect<T>::DETECT_Y);
  }
  // bump map is in the logic exchanged convex part and concave part.
  return - autoLevel(result);
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> PseudoBump<T>::minSquare(const Vec& input) {
  for(int i = 0; i < input.size(); i ++)
    if(!isfinite(input[i])) {
      Vec res(2);
      res[0] = res[1] = T(0);
      return res;
    }
  Vec avg(2);
  avg[0] = T(0);
  avg[1] = T(0);
  for(int i = 0; i < input.size(); i ++) {
    avg[0] += input[i];
    avg[1] += i;
  }
  avg[0] /= input.size();
  avg[1] /= input.size();
  
  Vec b(input.size());
  for(int i = 0; i < input.size(); i ++)
    b[i] = input[i] - avg[0];
  
  Mat A(input.size(), 2);
  for(int i = 0; i < input.size(); i ++) {
    A(i, 0) = pow(T(input[i] - avg[0]), T(2));
    A(i, 1) = (i - avg[1]) * (input[i] - avg[0]);
  }
  Eigen::JacobiSVD<Mat> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Mat Ut(svd.matrixU().transpose());
  Mat Vt(svd.matrixV());
  Vec w(svd.singularValues());
  for(int i = 0; i < w.size(); i ++)
    if(abs(w[i]) > sthresh)
      w[i] = T(1) / w[i];
  Vec result(Ut * b);
  for(int i = 0; i < w.size(); i ++)
    result[i] *= w[i];
  result = Vt * result;
  result[0] += avg[0];
  result[1] /= input.size();
  return result;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> PseudoBump<T>::getLineAxis(Vec p, Vec c, const int& w, const int& h) {
  // suppose straight ray from infinitely distant 1 x 1 size.
  // camera geometry with c, lookup on p[2] z-distance
  // virtual images.
  p[0]  = (p[0] - h / 2.) / (h / 2.);
  p[1]  = (p[1] - w / 2.) / (w / 2.);
  p[2]  = zz(p[2]);
  c[2]  = zz(c[2]);
  // <c + (p - c) * t, [0, 0, 1]> = z
  // fix z = 0, p_z = zz(z_max).
  T   t((- c[2]) / (p[2] - c[2]));
  Vec work((p - c) * t + c);
  work[0] = (work[0] * h / 2. + h / 2.);
  work[1] = (work[1] * w / 2. + w / 2.);
  return work;
}

template <typename T> Eigen::Matrix<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::prepareLineAxis(const Vec& p0, const Vec& p1, const int& z0, const int& rstp) {
  Vec dir(p1 - p0);
  T ndir(sqrt(dir.dot(dir)));
  dir /= ndir;
  
  Eigen::Matrix<Mat, Eigen::Dynamic, Eigen::Dynamic> result(z0, 2);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for
#endif
  for(int zi = 0; zi < z0; zi ++) {
    result(zi, 0) = Mat(3, stp);
    result(zi, 1) = Mat(3, stp);
    for(int s = 0; s < stp; s ++) {
      Vec cpoint(indiv(p0, p1, 1. / 2.));
      cpoint[0] += (s / (stp - 1.) - 1. / 2.) * rstp * dir[0];
      cpoint[1] += (s / (stp - 1.) - 1. / 2.) * rstp * dir[1];
      cpoint[2]  = z_max;

      const T zzi(zz(zi + 1));
      Vec camera(3);
      camera[0] = dir[0] * roff;
      camera[1] = dir[1] * roff;
      camera[2] = - zzi * rdist * ndir;
      Vec rd(getLineAxis(cpoint, camera, T(ww), T(hh)));
      camera    = - camera;
      camera[2] = - camera[2];
      Vec ld(getLineAxis(cpoint, camera, T(ww), T(hh)));
      rd -= indiv(p0, p1, 1. / 2.);
      ld -= indiv(p0, p1, 1. / 2.);
      result(zi, 0).col(s) = rd;
      result(zi, 1).col(s) = ld;
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
  T xx((int(x + .5) + 2 * w) % (2 * w));
  T yy((int(y + .5) + 2 * h) % (2 * h));
  if(abs(xx) >= w)
    xx = - xx + sgn(xx) * w;
  if(abs(yy) >= h)
    yy = - yy + sgn(yy) * h;
  return T(img((int(abs(yy)) + h) % h, (int(abs(xx)) + w) % w));
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

