#if !defined(_2D3D_PSEUDO_)

#include <cstdio>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "enlarge.hh"

template <typename T> class PseudoBump {
public:
  int    z_max;
  int    stp;
  int    rstp;
  int    enll;
  T      roff;
  T      rdist;
  T      sthresh;
  T      gthresh;
  int    zband;
  int    rband;
  int    pratio;
  T      Pi;
  
  PseudoBump();
  void initialize(const int& z_max, const int& stp, const int& rstp, const int& enll, const T& roff, const T& rdist, const T& sthresh, const T& gthresh, const int& rband, const int& zband);
  ~PseudoBump();
  
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> getPseudoBump(const Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& input, const bool& y_only);
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> rgb2l(const Matrix<T, Eigen::Dynamic, Eigen::Dynamic> rgb[3]);
private:
  T zz(const T& t);
  T sgn(const T& x);
  Eigen::Matrix<T, Eigen::Dynamic, 1> minSquare(const Eigen::Matrix<T, Eigen::Dynamic, 1>& input);
  Eigen::Matrix<T, Eigen::Dynamic, 1> getLineAxis(Eigen::Matrix<T, Eigen::Dynamic, 1> p, Eigen::Matrix<T, Eigen::Dynamic, 1> c, const int& w, const int& h);
  Eigen::Matrix<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Dynamic, Eigen::Dynamic> prepareLineAxis(const Eigen::Matrix<T, Eigen::Dynamic, 1>& p0, const Eigen::Matrix<T, Eigen::Dynamic, 1>& p1, const int& z0);
  T getImgPt(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& img, const T& x, const T& y);
  Eigen::Matrix<T, Eigen::Dynamic, 1> indiv(const Eigen::Matrix<T, Eigen::Dynamic, 1>& p0, const Eigen::Matrix<T, Eigen::Dynamic, 1>& p1, const T& pt);
  enlarger2exds<T, complex<T> > enlarger;
  Eigen::Matrix<T, Eigen::Dynamic, 1> complementLine(const Eigen::Matrix<T, Eigen::Dynamic, 1>& line);
  T complement(const Eigen::Matrix<int, Eigen::Dynamic, 1>& ptsi, const Eigen::Matrix<T, Eigen::Dynamic, 1>& pts, const T& idx);
  
  int ww;
  int hh;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Dop;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Dop2;
};

template <typename T> PseudoBump<T>::PseudoBump() {
  // initialize(32, 16, 16, 3, 100., 1e4, 1e-17, 2 / 256., 16, 4);
  initialize(8, 16, 16, 2, 100., 1e4, 1e-17, 2 / 256., 16, 4);
}

template <typename T> PseudoBump<T>::~PseudoBump() {
  ;
}

template <typename T> void PseudoBump<T>::initialize(const int& z_max, const int& stp, const int& rstp, const int& enll, const T& roff, const T& rdist, const T& sthresh, const T& gthresh, const int& rband, const int& zband) {
  const int enlp(std::max(int(pow(2, enll) - 1), int(1)));
  this->z_max   = z_max;
  this->stp     = stp;
  this->rstp    = rstp * enlp;
  this->enll    = enll;
  this->roff    = roff;
  this->rdist   = rdist;
  this->sthresh = sthresh;
  this->gthresh = gthresh;
  this->rband   = rband;
  this->zband   = zband;
  this->pratio  = enlp;
  this->Pi      = 4. * atan2(T(1.), T(1.));
};

template <typename T> T PseudoBump<T>::sgn(const T& x) {
  if(x < T(0))
    return - T(1);
  if(x > T(0))
    return T(1);
  return T(0);
}

template <typename T> T PseudoBump<T>::zz(const T& t) {
  return (t + T(1)) / z_max;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::rgb2l(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> rgb[3]) {
  return .212639 * rgb[0] + .715169 * rgb[1] + .072192 * rgb[2];
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::getPseudoBump(const Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& input, const bool& y_only) {
  if(y_only) {
    auto work(input);
    for(int i = 0; i < enll; i ++)
      work = enlarger.enlarge2ds(work, enlarger2exds<T, complex<T> >::ENLARGE_BOTH);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result(work.rows(), work.cols());
    ww = work.cols();
    hh = work.rows();
    Eigen::Matrix<T, Eigen::Dynamic, 1> p0(3), p1(3);
    p0[0] = work.cols() / 2;
    p0[1] = 0;
    p0[2] = 0;
    p1[0] = work.cols() / 2;
    p1[1] = work.rows();
    p1[2] = 0;
    cerr << "bump: preparing matrix." << endl;
    auto lrf(prepareLineAxis(p0, p1, z_max));
    for(int i = 0; i < lrf.rows(); i ++)
      for(int j = 0; j < lrf(i, 0).cols(); j ++) {
        lrf(i, 0)(0, j) = 0;
        lrf(i, 1)(0, j) = 0;
      }
    const int& delta(result.rows());
    for(int i = 0; i < result.cols(); i ++) {
      p0[0] = i;
      p1[0] = i;
      cerr << "bump: " << i << "/" << result.cols() << endl;
      for(int j = 0; j < result.rows(); j ++)
        result(j, i) = - T(1);
      Eigen::Matrix<T, Eigen::Dynamic, 1> zval(result.rows());
      for(int j = 0; j < zval.size(); j ++)
        zval[j] = T(0);
      for(int zz = 0; zz < lrf.rows(); zz ++) {
        Eigen::Matrix<T, Eigen::Dynamic, 1> pb(indiv(p0, p1, - 1 / T(delta)));
        for(int s = 0; s < delta; s ++) {
          auto pt(indiv(p0, p1, s / T(delta)));
          auto pf(indiv(p0, p1, (s + 1) / T(delta)));
          if(abs(getImgPt(work, pb[0], pb[1]) - getImgPt(work, pt[0], pt[1])) <= gthresh)
            continue;
          pb = pt;
          Eigen::Matrix<T, Eigen::Dynamic, 1> c(lrf(zz, 0).cols());
          for(int u = 0; u < c.size(); u ++) {
            c[u]  = getImgPt(work, lrf(zz, 0)(0, u) + pt[0], lrf(zz, 0)(1, u) + pt[1]);
            c[u] -= getImgPt(work, lrf(zz, 1)(0, u) + pt[0], lrf(zz, 1)(1, u) + pt[1]);
          }
          Eigen::Matrix<T, Eigen::Dynamic, 1> lc(c.size() / 2 + 1), rc(c.size() / 2 + 1);
          for(int u = 0; u < lc.size(); u ++)
            lc[u] = c[u];
          for(int u = 0; u < rc.size(); u ++)
            rc[u] = c[u - rc.size() + c.size()];
          Eigen::Matrix<T, Eigen::Dynamic, 1> lc2(zband), rc2(zband);
          for(int u = 0; u < lc2.size(); u ++)
            lc2[u] = lc[u - lc2.size() + lc.size()];
          for(int u = 0; u < rc2.size(); u ++)
            rc2[u] = rc[u - rc2.size() + rc.size()];
          auto msl(minSquare(Dop * lc)),    msr(minSquare(Dop * rc));
          auto msl2(minSquare(Dop2 * lc2)), msr2(minSquare(Dop2 * rc2));
          auto ms( abs(msl[0])  + abs(msr[0]) );
          auto ms2(abs(msl2[0]) + abs(msr2[0]));
          T n2(0);
          if(ms2 > ms)
            n2 = abs(ms2 - ms);
          if(zval[s] < n2) {
            result(s, i) = zz / T(z_max);
            zval[s]      = n2;
          }
        }
      }
      result.col(i) = complementLine(result.col(i));
    }
    int ppratio = (result.cols() + input.cols() - 1) / input.cols();
    cerr << input.cols() * ppratio - result.cols() << endl;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> res(input.rows(), result.cols());
    for(int i = 0; i < input.rows(); i ++) {
      res.row(i) = result.row(i * ppratio);
      int j;
      for(j = 1; j < ppratio && i * ppratio + j < result.rows(); j ++)
        res.row(i) += result.row(i * ppratio + j);
      res.row(i) /= j;
    }
    result = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(input.rows(), input.cols());
    for(int i = 0; i < input.cols(); i ++) {
      result.col(i) = res.col(i * ppratio);
      int j;
      for(j = 1; j < ppratio && i * ppratio + j < result.cols(); j ++)
        result.col(i) += res.col(i * ppratio + j);
      result.col(i) /= j;
    }
    return result;
  } else {
    auto pbx(getPseudoBump(input.transpose(), true).transpose());
    auto pby(getPseudoBump(input, true));
    return (pbx + pby) / 2.;
  }
  return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>();
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> PseudoBump<T>::minSquare(const Eigen::Matrix<T, Eigen::Dynamic, 1>& input) {
  for(int i = 0; i < input.size(); i ++)
    if(!isfinite(input[i])) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> res(2);
      res[0] = res[1] = T(0);
      return res;
    }
  Matrix<T, Eigen::Dynamic, 1> avg(2);
  avg[0] = T(0);
  avg[1] = T(0);
  for(int i = 0; i < input.size(); i ++) {
    avg[0] += input[i];
    avg[1] += i;
  }
  avg[0] /= input.size();
  avg[1] /= input.size();
  
  Matrix<T, Eigen::Dynamic, 1> b(input.size());
  for(int i = 0; i < input.size(); i ++)
    b[i] = input[i] - avg[0];
  
  Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A(input.size(), 2);
  for(int i = 0; i < input.size(); i ++) {
    A(i, 0) = pow(input[i] - avg[0], 2.);
    A(i, 1) = (i - avg[1]) * (input[i] - avg[0]);
  }
  Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
  auto U(svd.matrixU().transpose());
  auto Vt(svd.matrixV().transpose());
  auto w(svd.singularValues());
  for(int i = 0; i < w.size(); i ++)
    if(abs(w[i]) > sthresh)
      w[i] = T(1) / w[i];
  Eigen::Matrix<T, Eigen::Dynamic, 1> result(U * b);
  for(int i = 0; i < w.size(); i ++)
    result[i] *= w[i];
  result = Vt * result;
  return result;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> PseudoBump<T>::getLineAxis(Eigen::Matrix<T, Eigen::Dynamic, 1> p, Eigen::Matrix<T, Eigen::Dynamic, 1> c, const int& w, const int& h) {
  // suppose straight ray from infinitely distant 1 x 1 size.
  // camera geometry with c, lookup on p[2] z-distance
  // virtual images.
  p[0]  = (p[0] - w / 2.) / (w / 2.);
  p[1]  = (p[1] - h / 2.) / (h / 2.);
  p[2]  = zz(p[2]);
  c[2]  = zz(c[2]);
  // <c + (p - c) * t, [0, 0, 1]> = z
  // fix z = 0, p_z = zz(z_max).
  T t((- c[2]) / (p[2] - c[2]));
  Eigen::Matrix<T, Eigen::Dynamic, 1> work((p - c) * t + c);
  work[0] = (work[0] * w / 2. + w / 2.);
  while(abs(work[0]) > w) work[0] -= sgn(work[0]) * w;
  work[1] = (work[1] * h / 2. + h / 2.);
  while(abs(work[1]) > h) work[1] -= sgn(work[1]) * h;
  return work;
}

template <typename T> Eigen::Matrix<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::prepareLineAxis(const Eigen::Matrix<T, Eigen::Dynamic, 1>& p0, const Eigen::Matrix<T, Eigen::Dynamic, 1>& p1, const int& z0) {
  Eigen::Matrix<T, Eigen::Dynamic, 1> dir(p1 - p0);
  auto dirn2(sqrt(dir.dot(dir)));
  dir /= dirn2;
  
  Eigen::Matrix<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Dynamic, Eigen::Dynamic> result(z0, 2);
  for(int zi = 0; zi < z0; zi ++) {
    result(zi, 0) = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(3, stp);
    result(zi, 1) = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(3, stp);
    for(int s = 0; s < stp; s ++) {
      auto cpoint(indiv(p0, p1, 1. / 2.));
      cpoint[0] += (s / (stp - 1.) - 1. / 2.) * rstp * dir[0];
      cpoint[1] += (s / (stp - 1.) - 1. / 2.) * rstp * dir[1];
      cpoint[2]  = z_max;

      auto zzi(zi + 1);
      Eigen::Matrix<T, Eigen::Dynamic, 1> camera(3);
      camera[0] = dir[0] * roff;
      camera[1] = dir[1] * roff;
      camera[2] = - zzi / T(z_max) * rdist;
      auto rd(getLineAxis(cpoint, camera, T(ww), T(hh)));
      camera    = - camera;
      camera[2] = - camera[2];
      auto ld(getLineAxis(cpoint, camera, T(ww), T(hh)));
      rd -= indiv(p0, p1, 1. / 2.);
      ld -= indiv(p0, p1, 1. / 2.);
      result(zi, 0).col(s) = rd;
      result(zi, 1).col(s) = ld;
    }
  }
  Eigen::Matrix<complex<T>, Eigen::Dynamic, Eigen::Dynamic> Dopb(stp / 2 + 1, stp / 2 + 1);
  Eigen::Matrix<complex<T>, Eigen::Dynamic, Eigen::Dynamic> Dop2b(zband, zband);
  Eigen::Matrix<complex<T>, Eigen::Dynamic, Eigen::Dynamic> Iopb(Dopb.rows(), Dopb.cols());
  Eigen::Matrix<complex<T>, Eigen::Dynamic, Eigen::Dynamic> Iop2b(Dop2b.rows(), Dop2b.cols());
  complex<T> I(sqrt(complex<T>(- 1)));
  for(int i = 0; i < Dopb.rows(); i ++)
    for(int j = 0; j < Dopb.cols(); j ++) {
      Dopb(i, j) = exp(complex<T>(- 2) * Pi * I * complex<T>(i * j) / T(Dopb.rows()));
      Iopb(i, j) = exp(complex<T>(  2) * Pi * I * complex<T>(i * j) / T(Dopb.rows())) / T(Dopb.rows());
    }
  for(int i = 0; i < Dopb.rows(); i ++)
    Dopb.row(i) *= complex<T>(- 2.) * Pi * I * T(i) / T(Dopb.rows());
  Dop = (Iopb * Dopb).real();
  for(int i = 0; i < Dop2b.rows(); i ++)
    for(int j = 0; j < Dop2b.cols(); j ++) {
      Dop2b(i, j) = exp(complex<T>(- 2) * Pi * I * complex<T>(i * j) / T(Dop2b.rows()));
      Iop2b(i, j) = exp(complex<T>(  2) * Pi * I * complex<T>(i * j) / T(Dop2b.rows())) / T(Dop2b.rows());
    }
  for(int i = 0; i < Dop2b.rows(); i ++)
    Dop2b.row(i) *= complex<T>(- 2.) * Pi * I * T(i) / T(Dop2b.rows());
  Dop2 = (Iop2b * Dop2b).real();
  return result;
}

template <typename T> T PseudoBump<T>::getImgPt(const Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& img, const T& x, const T& y) {
  const int& w(img.cols());
  const int& h(img.rows());
  auto xx((int(x + .5) + 2 * w) % (2 * w));
  auto yy((int(y + .5) + 2 * h) % (2 * h));
  if(abs(xx) >= w)
    xx = - xx + sgn(xx) * w;
  if(abs(yy) >= h)
    yy = - yy + sgn(yy) * h;
  return T(img((int(abs(yy)) + h) % h, (int(abs(xx)) + w) % w));
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> PseudoBump<T>::indiv(const Eigen::Matrix<T, Eigen::Dynamic, 1>& p0, const Eigen::Matrix<T, Eigen::Dynamic, 1>& p1, const T& pt) {
  return (p1 - p0) * pt + p0;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> PseudoBump<T>::complementLine(const Eigen::Matrix<T, Eigen::Dynamic, 1>& line) {
  Eigen::Matrix<int, Eigen::Dynamic, 1> ptsi(line.size() * 2);
  Eigen::Matrix<T,   Eigen::Dynamic, 1> pts(line.size() * 2);
  Eigen::Matrix<T,   Eigen::Dynamic, 1> result(line);
  int idx = 0;
  for(int i = 0; i < line.size(); i ++) {
    if(0 <= line[i]) {
      if(idx < 1 && i != 0) {
        ptsi[idx] = - i;
        pts[idx]  = line[i];
        idx ++;
      }
      ptsi[idx] = i;
      pts[idx]  = line[i];
      idx ++;
    }
  }
  if(idx > 0) {
    ptsi[idx] = line.size() + (line.size() - ptsi[idx - 1]);
    pts[idx]  = pts[idx - 1];
    idx ++;
  }
  for(int i = idx; i < ptsi.size(); i ++) {
    ptsi[i] = - 1;
    pts[i]  = - 1;
  }
  if(idx < 3)
    for(int i = 0; i < result.size(); i ++)
      result[i] = T(1. / 2);
  else
    for(int i = 0; i < result.size(); i ++)
      if(line[i] < T(0))
        result[i] = complement(ptsi, pts, i);
  return result;
}

template <typename T> T PseudoBump<T>::complement(const Eigen::Matrix<int, Eigen::Dynamic, 1>& ptsi, const Eigen::Matrix<T, Eigen::Dynamic, 1>& pts, const T& idx) {
  int size(pts.size() - 1);
  for(int i = 0; i < pts.size(); i ++)
    if(pts[i] < 0) {
      size = i;
      break;
    }
  size ++;
  if(size <= 1) {
    if(size < 1)
      return 0.;
    return pts[0];
  }
  Eigen::Matrix<int, Eigen::Dynamic, 1> rng(3);
  for(int i = 0; i < rng.size(); i ++)
    rng[i] = - 1;
  for(int i = 0; i < size - 1; i ++)
    if(ptsi[i] <= idx - rband && idx - rband <= ptsi[i + 1]) {
      rng[0] = i;
      break;
    }
  for(int i = 0; i < size - 1; i ++)
    if(ptsi[i] <= idx && idx <= ptsi[i + 1]) {
      rng[1] = i;
      break;
    }
  for(int i = 0; i < size - 1; i ++)
    if(ptsi[i] <= idx + rband && idx + rband <= ptsi[i + 1]) {
      rng[2] = i;
      break;
    }
  if((rng[0] == rng[1] && rng[1] == rng[2]) || rng[0] < 0 || rng[1] < 0 || rng[2] < 0) {
    if(rng[0] < 0)
      return 0.;
    return pts[rng[0]];
  }
  // specific complement.
  T res(0);
  for(int i = 0; i < rng.size(); i ++) {
    T work(1);
    for(int j = 0; j < rng.size(); j ++) {
      if(ptsi[rng[i]] == ptsi[rng[j]])
        continue;
      work *= (T(idx) - ptsi[rng[j]]) / T(pts[rng[i]] - pts[rng[j]]);
    }
    res += pts[rng[i]] * work;
  }
  return std::max(T(0), std::min(T(2), res));
}

#define _2D3D_PSEUDO_
#endif

