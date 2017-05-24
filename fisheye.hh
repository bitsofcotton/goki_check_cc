#if !defined(_2D3D_PSEUDO_)

#include <cstdio>
#include <cmath>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "enlarge.hh"
#include "edgedetect.hh"

template <typename T> class PseudoBump {
public:
  int    z_max;
  int    stp;
  int    rstp;
  int    estp;
  int    enll;
  T      roff;
  T      rdist;
  T      gratio;
  int    zband;
  int    nlevel;
  int    nedge;
  int    rband;
  T      sthresh;
  T      Pi;
  
  PseudoBump();
  void initialize(const int& z_max, const int& stp, const int& rstp, const int& estp, const int& enll, const T& roff, const T& rdist, const T& gratio, const int& nlevel, const int& nedge);
  ~PseudoBump();
  
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> getPseudoBumpSub(const Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& input);
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> getPseudoBump(const Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& input, const bool& y_only);
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> rgb2l(const Matrix<T, Eigen::Dynamic, Eigen::Dynamic> rgb[3]);
private:
  T zz(const T& t);
  T sgn(const T& x);
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> getMosaic(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& input, const int& i);
  Eigen::Matrix<T, Eigen::Dynamic, 1> minSquare(const Eigen::Matrix<T, Eigen::Dynamic, 1>& input);
  Eigen::Matrix<T, Eigen::Dynamic, 1> getLineAxis(Eigen::Matrix<T, Eigen::Dynamic, 1> p, Eigen::Matrix<T, Eigen::Dynamic, 1> c, const int& w, const int& h);
  Eigen::Matrix<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Dynamic, Eigen::Dynamic> prepareLineAxis(const Eigen::Matrix<T, Eigen::Dynamic, 1>& p0, const Eigen::Matrix<T, Eigen::Dynamic, 1>& p1, const int& z0);
  T getImgPt(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& img, const T& y, const T& x);
  Eigen::Matrix<T, Eigen::Dynamic, 1> indiv(const Eigen::Matrix<T, Eigen::Dynamic, 1>& p0, const Eigen::Matrix<T, Eigen::Dynamic, 1>& p1, const T& pt);
  enlarger2exds<T, complex<T> > enlarger;
  Eigen::Matrix<T, Eigen::Dynamic, 1> complementLine(const Eigen::Matrix<T, Eigen::Dynamic, 1>& line, const int& rband);
  
  int ww;
  int hh;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Dop;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Dop2;
};

template <typename T> PseudoBump<T>::PseudoBump() {
  initialize(8, 16, 3, 4, 1, .5, 1e3, 0., 8, 32);
}

template <typename T> PseudoBump<T>::~PseudoBump() {
  ;
}

template <typename T> void PseudoBump<T>::initialize(const int& z_max, const int& stp, const int& rstp, const int& estp, const int& enll, const T& roff, const T& rdist, const T& gratio, const int& nlevel, const int& nedge) {
  const int enlp(std::max(int(pow(2, enll) - 1), int(1)));
  this->z_max   = z_max;
  this->stp     = stp;
  this->rstp    = rstp;
  this->estp    = estp;
  this->enll    = enll;
  this->roff    = roff;
  this->rdist   = rdist;
  this->gratio  = gratio;
  this->zband   = int(std::max(1., std::sqrt(stp)));
  this->nlevel  = nlevel;
  this->nedge   = nedge;
  Pi            = 4. * atan2(T(1.), T(1.));
  rband         = this->rstp;
  sthresh       = pow(1. / 256., 6.);
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

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::getPseudoBumpSub(const Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& input) {
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> work(input);
  for(int i = 0; i < enll; i ++)
    work = enlarger.enlarge2ds(work, enlarger2exds<T, complex<T> >::ENLARGE_Y);
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result(work.rows(), work.cols());
  const int ppratio = (result.rows() + input.rows() - 1) / input.rows();
  ww = work.cols();
  hh = work.rows();
  Eigen::Matrix<T, Eigen::Dynamic, 1> p0(3), p1(3);
  p0[0] = 0;
  p0[1] = work.cols() / 2;
  p0[2] = 0;
  p1[0] = work.rows();
  p1[1] = work.cols() / 2;
  p1[2] = 0;
  cerr << " bump";
  cerr.flush();
  auto lrf(prepareLineAxis(p0, p1, z_max));
  for(int i = 0; i < lrf.rows(); i ++)
    for(int j = 0; j < lrf(i, 0).cols(); j ++) {
      lrf(i, 0)(1, j) = 0;
      lrf(i, 1)(1, j) = 0;
    }
  const int& delta(result.rows());
  for(int i = 0; i < result.cols(); i ++) {
    p0[1] = i;
    p1[1] = i;
    cerr << ".";
    cerr.flush();
    for(int j = 0; j < result.rows(); j ++)
      result(j, i) = - T(1);
    Eigen::Matrix<T, Eigen::Dynamic, 1> zval(result.rows());
    for(int j = 0; j < zval.size(); j ++)
      zval[j] = T(0);
    for(int s = 0; s < delta; s ++) {
      auto pt(indiv(p0, p1, s / T(delta)));
      for(int zz = 0; zz < lrf.rows(); zz ++) {
        Eigen::Matrix<T, Eigen::Dynamic, 1> c(lrf(zz, 0).cols());
        for(int u = 0; u < c.size(); u ++) {
          c[u]  = getImgPt(work, lrf(zz, 0)(0, u) + pt[0],
                                 lrf(zz, 0)(1, u) + pt[1]);
          c[u] -= getImgPt(work, lrf(zz, 1)(0, u) + pt[0],
                                 lrf(zz, 1)(1, u) + pt[1]);
        }
        Eigen::Matrix<T, Eigen::Dynamic, 1> lc(c.size() / 2 + 1);
        Eigen::Matrix<T, Eigen::Dynamic, 1> rc(c.size() / 2 + 1);
        for(int u = 0; u < lc.size(); u ++)
          lc[u] = c[u];
        for(int u = 0; u < rc.size(); u ++)
          rc[u] = c[u - rc.size() + c.size()];
        Eigen::Matrix<T, Eigen::Dynamic, 1> lc2(zband), rc2(zband);
        for(int u = 0; u < lc2.size(); u ++)
          lc2[u] = lc[u - lc2.size() + lc.size()];
        for(int u = 0; u < rc2.size(); u ++)
          rc2[u] = rc[u];
        const Eigen::Matrix<T, Eigen::Dynamic, 1> msl( Dop  * lc);
        const Eigen::Matrix<T, Eigen::Dynamic, 1> msr( Dop  * rc);
        const Eigen::Matrix<T, Eigen::Dynamic, 1> msl2(Dop2 * lc2);
        const Eigen::Matrix<T, Eigen::Dynamic, 1> msr2(Dop2 * rc2);
        auto ms( abs(msl[ msl.size()  - 1] - msr[0]) );
        auto ms2(abs(msl2[msl2.size() - 1] - msr2[0]));
        const T n2(ms2 / (ms + 1.));
        if(isfinite(n2) && gratio < n2 && zval[s] < n2) {
          result(s, i) = zz / T(z_max);
          zval[s]      = n2;
        }
      }
    }
    result.col(i) = complementLine(result.col(i), rband * ppratio);
  }
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> res(input.rows(), input.cols());
  for(int x = 0; x < input.cols(); x ++)
    for(int y0 = 0; y0 < input.rows(); y0 ++)
      res(y0, x) = result(y0 * ppratio, x);
  return res;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::getMosaic(const Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& input, const int& i) {
  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> work(input);
  for(int x = 0; x < (input.cols() + i - 1) / i; x ++) {
    for(int y = 0; y < (input.rows() + i - 1) / i; y ++) {
      std::vector<float> meds;
      for(int xx = x * i;
              xx < std::min(int(work.cols()), (x + 1) * i);
              xx ++)
        for(int yy = y * i;
                yy < std::min(int(work.rows()), (y + 1) * i);
                yy ++) {
          meds.push_back(work(yy, xx));
          work(yy, xx) = - 1.;
        }
      std::sort(meds.begin(), meds.end());
      work(std::min(int(input.rows()) - 1, y * i + i / 2), x * i) = meds[meds.size() / 2];
    }
    work.col(x * i) = complementLine(work.col(x * i), 1);
  }
  for(int y = 0; y < work.rows(); y ++)
    work.row(y) = complementLine(work.row(y), 1);
  return work;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::getPseudoBump(const Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& input, const bool& y_only)
 {
  if(y_only) {
    T ratio = 0;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> out(input.rows(), input.cols());
    for(int i = 0; i < out.rows(); i ++)
      for(int j = 0; j < out.cols(); j ++)
        out(i, j) = 0.;
    for(int i = estp; i < input.rows() / estp; i *= 2) {
      out   += getMosaic(getPseudoBumpSub(getMosaic(input, i)), i);
      ratio += 1.;
    }
    out /= ratio;
    edgedetect<T, complex<T> > edge;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> edges(edge.detect(input, edgedetect<T, complex<T> >::COLLECT_Y));
    std::vector<T> estat;
    for(int i = 0; i < edges.rows(); i ++)
      for(int j = 0; j < edges.cols(); j ++)
        estat.push_back(edges(i, j));
    std::sort(estat.begin(), estat.end());
    for(int i = 0; i < out.rows(); i ++)
      for(int j = 0; j < out.cols(); j ++)
        if(estat[estat.size() * (nedge - 1) / nedge] < edges(i, j))
          out(i, j) = - T(1);
    for(int i = 0; i < out.cols(); i ++)
      out.col(i) = complementLine(out.col(i), rband);
    std::vector<T> stat;
    for(int i = 0; i < out.rows(); i ++)
      for(int j = 0; j < out.cols(); j ++) {
        // XXX fixme:
        // out(i, j) = 1. - out(i, j);
        stat.push_back(out(i, j));
      }
    std::sort(stat.begin(), stat.end());
    const T mm(stat[stat.size() / nlevel]);
    T MM(stat[stat.size() * (nlevel - 1) / nlevel]);
    if(MM == mm)
      MM = mm + 1.;
    for(int i = 0; i < out.rows(); i ++)
      for(int j = 0; j < out.cols(); j ++)
        out(i, j) = (std::max(std::min(out(i, j), MM), mm) - mm) / (MM - mm);
    return out;
  } else {
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> pbx(getPseudoBump(input.transpose(), true).transpose());
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> pby(getPseudoBump(input, true));
    return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>((pbx + pby) / 2.);
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
  auto Ut(svd.matrixU().transpose());
  auto Vt(svd.matrixV());
  auto w(svd.singularValues());
  for(int i = 0; i < w.size(); i ++)
    if(abs(w[i]) > sthresh)
      w[i] = T(1) / w[i];
  Eigen::Matrix<T, Eigen::Dynamic, 1> result(Ut * b);
  for(int i = 0; i < w.size(); i ++)
    result[i] *= w[i];
  result = Vt * result;
  result[0] += avg[0];
  return result;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> PseudoBump<T>::getLineAxis(Eigen::Matrix<T, Eigen::Dynamic, 1> p, Eigen::Matrix<T, Eigen::Dynamic, 1> c, const int& w, const int& h) {
  // suppose straight ray from infinitely distant 1 x 1 size.
  // camera geometry with c, lookup on p[2] z-distance
  // virtual images.
  p[0]  = (p[0] - h / 2.) / (h / 2.);
  p[1]  = (p[1] - w / 2.) / (w / 2.);
  p[2]  = zz(p[2]);
  c[2]  = zz(c[2]);
  // <c + (p - c) * t, [0, 0, 1]> = z
  // fix z = 0, p_z = zz(z_max).
  T t((- c[2]) / (p[2] - c[2]));
  Eigen::Matrix<T, Eigen::Dynamic, 1> work((p - c) * t + c);
  work[0] = (work[0] * h / 2. + h / 2.);
  work[1] = (work[1] * w / 2. + w / 2.);
  return work;
}

template <typename T> Eigen::Matrix<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::prepareLineAxis(const Eigen::Matrix<T, Eigen::Dynamic, 1>& p0, const Eigen::Matrix<T, Eigen::Dynamic, 1>& p1, const int& z0) {
  Eigen::Matrix<T, Eigen::Dynamic, 1> dir(p1 - p0);
  dir /= sqrt(dir.dot(dir));
  
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

template <typename T> T PseudoBump<T>::getImgPt(const Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& img, const T& y, const T& x) {
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

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> PseudoBump<T>::complementLine(const Eigen::Matrix<T, Eigen::Dynamic, 1>& line, const int& rband) {
  std::vector<int> ptsi;
  std::vector<T>   pts;
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
  Eigen::Matrix<T, Eigen::Dynamic, 1> result(line.size());
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
    for(; rng[0] < ptsi.size() - 1; rng[0] ++)
      if(i - rband <= ptsi[rng[0] + 1])
        break;
    for(; rng[1] < ptsi.size(); rng[1] ++)
      if(i - rband <= ptsi[rng[1]] && ptsi[rng[1]] <= i + rband)
        break;
    for(; rng[2] < ptsi.size(); rng[2] ++)
      if(i + rband < ptsi[rng[2]])
        break;
    if(rng[2] < rng[1])
      rng[1] = (rng[2] + rng[0]) / 2;
    if(rng[1] == rng[0] && 0 < rng[0])
      rng[0] --;
    if(rng[1] == rng[2] && rng[2] < ptsi.size() - 1)
      rng[2] ++;
    result[i] = 0.;
    for(int ii = 0; ii < 3; ii ++) {
      T work(1);
      for(int jj = 0; jj < 3; jj ++)
        if(ptsi[rng[ii]] != ptsi[rng[jj]])
          work *= (T(i) - ptsi[rng[jj]]) / T(ptsi[rng[ii]] - ptsi[rng[jj]]);
      result[i] += work * pts[rng[ii]];
    }
  }
  return result;
}

#define _2D3D_PSEUDO_
#endif

