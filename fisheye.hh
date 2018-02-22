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
using std::max;
using std::min;
using std::sqrt;
using std::exp;
using std::pow;
using std::isfinite;
using std::cerr;
using std::flush;
using std::vector;

template <typename T> class PseudoBump {
public:
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1>              Vec;
  typedef Eigen::Matrix<T,   3, 1> Vec3;
  typedef Eigen::Matrix<int, 3, 1> Veci3;
  
  PseudoBump();
  ~PseudoBump();
  void initialize(const int& z_max, const int& stp, const int& renl, const int& nslide);
  void getPseudoVec(const Mat& in, vector<Vec3>& geoms, vector<Veci3>& delaunay, const int& vbox = 4, const T& rz = T(1) / T(3));
  Mat  getPseudoBump(Mat in, const bool& elim0 = true);
  
private:
  Vec  complementLine(const Vec& line, const T& rratio = T(.5), const int& guard = int(1));
  Mat  autoLevel(const Mat& data, int npad = - 1);
  Vec  getPseudoBumpSub(const Vec& work, const vector<Vec>& cf);
  const T& getImgPt(const Vec& img, const T& y);
  Vec  getLineAxis(Vec p, Vec c);
  vector<Vec> prepareLineAxis(const T& rstp);
  
  int z_max;
  int renl;
  
  T   Pi;
  Vec Dop;
  Mat Dop2z;
};

template <typename T> PseudoBump<T>::PseudoBump() {
  initialize(60, 21, 60, 12);
}

template <typename T> PseudoBump<T>::~PseudoBump() {
  ;
}

template <typename T> void PseudoBump<T>::initialize(const int& z_max, const int& stp, const int& renl, const int& nslide) {
  this->z_max = z_max;
  this->renl  = renl;
  assert(1 < z_max);
  assert(0 < stp);
  assert(0 < nslide && nslide < z_max);
  this->Pi      = T(4) * atan2(T(1), T(1));
  {
    Eigen::Matrix<complex<T>, Eigen::Dynamic, Eigen::Dynamic>  DFT(stp, stp);
    Eigen::Matrix<complex<T>, Eigen::Dynamic, Eigen::Dynamic> IDFT(stp, stp);
    for(int i = 0; i < DFT.rows(); i ++)
      for(int j = 0; j < DFT.cols(); j ++) {
        DFT( i, j) = exp(complex<T>(- 2.) * Pi * sqrt(complex<T>(- 1)) * complex<T>(i * j / T(stp)));
        IDFT(i, j) = exp(complex<T>(  2.) * Pi * sqrt(complex<T>(- 1)) * complex<T>(i * j / T(stp))) / T(stp);
      }
    for(int i = 0; i < DFT.rows(); i ++)
      DFT.row(i) *= complex<T>(2.) * Pi * sqrt(complex<T>(- 1)) * T(i) / T(DFT.rows());
    Dop = (IDFT * DFT).real().row(IDFT.rows() / 2);
  }
  {
    Eigen::Matrix<complex<T>, Eigen::Dynamic, Eigen::Dynamic>  DFT(nslide, nslide);
    Eigen::Matrix<complex<T>, Eigen::Dynamic, Eigen::Dynamic> IDFT(nslide, nslide);
    for(int i = 0; i < DFT.rows(); i ++)
      for(int j = 0; j < DFT.cols(); j ++) {
        DFT( i, j) = exp(complex<T>(- 2.) * Pi * sqrt(complex<T>(- 1)) * complex<T>(i * j / T(nslide)));
        IDFT(i, j) = exp(complex<T>(  2.) * Pi * sqrt(complex<T>(- 1)) * complex<T>(i * j / T(nslide))) / T(nslide);
      }
    for(int i = 0; i < DFT.rows(); i ++)
      DFT.row(i) *= pow(complex<T>(2.) * Pi * sqrt(complex<T>(- 1)) * T(i) / T(DFT.rows()), T(2));
    Dop2z = (IDFT * DFT).real();
  }
  return;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> PseudoBump<T>::getPseudoBumpSub(const Vec& work, const vector<Vec>& cf) {
  cerr << "." << flush;
  Vec workv(work.size() * renl);
  for(int i = 0; i < workv.size(); i ++)
    if(i % renl == renl / 2)
      workv[i] = work[i / renl];
    else
      workv[i] = - T(1);
  workv = complementLine(workv);
  Vec result(work.size());
  Vec rsub(work.size());
  for(int j = 0; j < rsub.size(); j ++) {
    result[j] =   T(0);
    rsub[j]   = - T(1);
  }
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int s = 0; s < work.size(); s ++) {
    vector<T> n2;
    n2.resize(cf.size());
    for(int zz = 0; zz < cf.size(); zz ++) {
      // d/dt (local color) / ||local color||:
      Vec c(cf[zz].size());
      for(int u = 0; u < c.size(); u ++)
        c[u] = getImgPt(workv, cf[zz][u] + s * renl + renl / 2);
      n2[zz] = abs(Dop.dot(c) / sqrt(c.dot(c)));
    }
    for(int j = 0; j < z_max - Dop2z.cols(); j ++) {
      Vec zsub(Dop2z.cols());
      for(int k = 0; k < zsub.size(); k ++)
        zsub[k] = n2[k + j];
      zsub = Dop2z * zsub;
      T m(0);
      T score(- z_max * 2);
      for(int zz = 0; zz < zsub.size(); zz ++)
        if(zsub[zz] < m) {
          score = zz / T(cf.size());
          m     = zsub[zz];
        }
      if(m == T(0))
        for(int zz = 0; zz < zsub.size(); zz ++)
          if(m < n2[zz]) {
            score = zz / T(cf.size());
            m     = n2[zz];
          }
      rsub[s] += score;
    }
  }
  return complementLine(rsub);
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::getPseudoBump(Mat in, const bool& elim0) {
  if(elim0) {
    for(int i = 0; i < in.rows(); i ++)
      for(int j = 0; j < in.cols(); j ++)
        if(in(i, j) == T(0))
          in(i, j) = - T(1);
    Mat work(in);
    for(int i = 0; i < in.rows(); i ++)
      in.row(i) = complementLine(in.row(i));
    for(int i = 0; i < in.cols(); i ++)
      in.col(i) = complementLine(in.col(i));
    for(int i = 0; i < work.cols(); i ++)
      work.col(i) = complementLine(work.col(i));
    for(int i = 0; i < work.rows(); i ++)
      work.row(i) = complementLine(work.row(i));
    in = (in + work) / T(2);
  }

  auto cf(prepareLineAxis(sqrt(T(in.rows() * in.cols()))));
  Mat  result(in.rows(), in.cols());
  for(int i = 0; i < result.cols(); i ++)
    result.col(i) = getPseudoBumpSub(in.col(i), cf);
  return autoLevel(result);
}

// get bump with multiple scale and vectored result.
template <typename T> void PseudoBump<T>::getPseudoVec(const Mat& in, vector<Vec3>& geoms, vector<Veci3>& delaunay, const int& vbox, const T& rz) {
  /*
   * Get vector based bumps.
   */
  geoms = vector<Vec3>();
  T aavg(0);
  for(int i = 0; i < in.rows(); i ++)
    for(int j = 0; j < in.cols(); j ++)
      aavg += in(i, j);
  aavg /= in.rows() * in.cols();
  for(int i = 0; i < in.rows() / vbox + 1; i ++)
    for(int j = 0; j < in.cols() / vbox + 1; j ++) {
      if(in.rows() < (i + 1) * vbox ||
         in.cols() < (j + 1) * vbox) {
        if(geoms.size() >= 1) {
          Vec3 gbuf;
          gbuf[0] = i * vbox;
          gbuf[1] = j * vbox;
          gbuf[2] = geoms[geoms.size() - 1][2];
          geoms.push_back(gbuf);
        }
        continue;
      }
      T avg(0);
      for(int ii = i * vbox; ii < (i + 1) * vbox; ii ++)
        for(int jj = j * vbox; jj < (j + 1) * vbox; jj ++)
          avg += in(ii, jj);
      Vec3 work;
      work[0] = i * vbox;
      work[1] = j * vbox;
      const T intens(avg / vbox / vbox - aavg);
      // work[2] = - (intens < T(0) ? - T(1) : T(1)) * pow(abs(intens), T(2)) * sqrt(T(in.rows() * in.cols())) * rz;
      work[2] = - intens * sqrt(T(in.rows() * in.cols())) * rz;
      geoms.push_back(work);
    }
  Vec3 avg;
  avg[0] = avg[1] = avg[2] = T(0);
  for(int i = 0; i < geoms.size(); i ++)
    avg += geoms[i];
  avg /= geoms.size();
  for(int i = 0; i < geoms.size(); i ++)
    geoms[i][2] -= avg[2];
  delaunay = vector<Veci3>();
  for(int i = 1; i < in.rows() / vbox + 1; i ++)
    for(int j = 0; j < in.cols() / vbox; j ++) {
      Veci3 work, work2;
      work[0]  = (i - 1) * (in.cols() / vbox + 1) + j;
      work[1]  =  i      * (in.cols() / vbox + 1) + j;
      work[2]  =  i      * (in.cols() / vbox + 1) + j + 1;
      work2[0] = (i - 1) * (in.cols() / vbox + 1) + j;
      work2[2] = (i - 1) * (in.cols() / vbox + 1) + j + 1;
      work2[1] =  i      * (in.cols() / vbox + 1) + j + 1;
      delaunay.push_back(work);
      delaunay.push_back(work2);
    }
  return;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> PseudoBump<T>::getLineAxis(Vec p, Vec c) {
  // suppose straight ray from infinitely distant 2 x 2 size.
  // camera geometry with c, lookup on p[1] z-distance
  // virtual images.
  // virtually : z = p[1], p[1] = 0, all y axis := 0.
  // <c + (p - c) * t, [0, 0, 1]> = z
  const T t((p[1] - c[1]) / (- c[1]));
  return (p - c) * t + c;
}

template <typename T> vector<Eigen::Matrix<T, Eigen::Dynamic, 1> > PseudoBump<T>::prepareLineAxis(const T& rstp) {
  const int& stp(Dop.size());
  
  // N.B. ray is from infinite far, so same side of these.
  Vec camera(2);
  camera[0] = T(0);
  camera[1] = T(1);
  
  vector<Vec> result;
  result.resize(z_max);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int zi = 0; zi < z_max; zi ++) {
    result[zi] = Vec(stp);
    for(int s = 0; s < stp; s ++) {
      Vec cpoint(2);
      cpoint[0] = (s / (stp - 1.) - 1. / 2.) * stp;
      // XXX checkme:
      cpoint[1] = (zi + 1) / T(z_max + 1);
      result[zi][s] = getLineAxis(cpoint, camera)[0] * rstp;
    }
  }
  return result;
}

template <typename T> const T& PseudoBump<T>::getImgPt(const Vec& img, const T& y) {
  const int& h(img.size());
  const int  yy(abs((int(y + .5) + 3 * h) % (2 * h) - h) % h);
  return img[yy];
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

