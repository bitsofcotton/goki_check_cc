/* BSD 3-Clause License:
 * Copyright (c) 2018, bitsofcotton.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 *    Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *    Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation or other materials provided with the distribution.
 *    Neither the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

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
  typedef complex<T> U;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1>              Vec;
  typedef Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> MatU;
  typedef Eigen::Matrix<T,   3, 1> Vec3;
  typedef Eigen::Matrix<int, 3, 1> Veci3;
  
  PseudoBump();
  ~PseudoBump();
  void initialize(const int& z_max, const int& stp);
  void getPseudoVec(const Mat& in, vector<Vec3>& geoms, vector<Veci3>& delaunay, const int& vbox = 4, const T& rz = T(1) / T(8));
  Mat  getPseudoBump(const Mat& in);
  
private:
  Vec         getPseudoBumpSub(const Vec& work, const vector<Vec>& cf);
  const T&    getImgPt(const Vec& img, const T& y);
  vector<Vec> prepareLineAxis(const T& rstp);
  
  int z_max;
  Vec Dop;
  Vec DDDop;
};

template <typename T> PseudoBump<T>::PseudoBump() {
  initialize(150, 31);
}

template <typename T> PseudoBump<T>::~PseudoBump() {
  ;
}

template <typename T> void PseudoBump<T>::initialize(const int& z_max, const int& stp) {
  this->z_max  = z_max;
  assert(1 < z_max && 0 < stp);
  const T Pi(T(4) * atan2(T(1), T(1)));
  MatU DFT(stp, stp), IDFT(stp, stp);
  for(int i = 0; i < DFT.rows(); i ++)
    for(int j = 0; j < DFT.cols(); j ++) {
      DFT( i, j) = exp(U(- 2.) * Pi * sqrt(U(- 1)) * U(i * j / T(stp)));
      IDFT(i, j) = exp(U(  2.) * Pi * sqrt(U(- 1)) * U(i * j / T(stp))) / T(stp);
    }
  auto DDDFT(DFT);
  for(int i = 0; i < DFT.rows(); i ++) {
    const U rtheta(U(2.) * Pi * sqrt(U(- 1)) * T(i) / T(DFT.rows()));
    DFT.row(i)   *= rtheta;
    DDDFT.row(i) *= rtheta * rtheta * rtheta;
  }
  Dop    = (IDFT *   DFT).row(IDFT.rows() / 2).real().template cast<T>();
  DDDop  = (IDFT * DDDFT).row(IDFT.rows() / 2).real().template cast<T>();
  Dop   /= sqrt(  Dop.dot(  Dop));
  DDDop /= sqrt(DDDop.dot(DDDop));
  return;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> PseudoBump<T>::getPseudoBumpSub(const Vec& work, const vector<Vec>& cf) {
  cerr << "." << flush;
  Vec result(work.size());
  for(int s = 0; s < work.size(); s ++) {
    T mbuf(- 1);
    for(int z = 0; z < cf.size(); z ++) {
      // d/dt (local color) / ||local color||:
      Vec c(cf[z].size());
      for(int u = 0; u < c.size(); u ++)
        c[u] = getImgPt(work, cf[z][u] + s);
      const T score(Dop.dot(c) * DDDop.dot(c) / c.dot(c));
      if(mbuf < score) {
        result[s] = z / T(cf.size());
        mbuf      = score;
      }
    }
    // N.B. if we can't focus, it's from infinitely far.
    if(result[s] == T(0))
      result[s] = T(1);
  }
  return result;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::getPseudoBump(const Mat& in) {
  auto cf(prepareLineAxis(sqrt(T(in.rows() * in.cols()))));
  Mat result(in.rows(), in.cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < in.cols(); i ++)
    result.col(i)  = getPseudoBumpSub(in.col(i), cf);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < in.rows(); i ++)
    result.row(i) += getPseudoBumpSub(in.row(i), cf);
  return result;
}

// get bump with multiple scale and vectorized result.
template <typename T> void PseudoBump<T>::getPseudoVec(const Mat& in, vector<Vec3>& geoms, vector<Veci3>& delaunay, const int& vbox, const T& rz) {
  assert(0 < vbox && T(0) < rz);
  // get vectorize.
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

template <typename T> vector<Eigen::Matrix<T, Eigen::Dynamic, 1> > PseudoBump<T>::prepareLineAxis(const T& rstp) {
  // N.B. ray is from infinite far, so same side of these.
  Vec camera(2);
  camera[0] = T(0);
  camera[1] = T(1);
  
  vector<Vec> result;
  result.resize(z_max, Vec(Dop.size()));
  T absmax(0);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int zi = 0; zi < result.size(); zi ++) {
    for(int s = 0; s < result[zi].size(); s ++) {
      Vec cpoint(2);
      cpoint[0] = (s / T(Dop.size() - 1) - 1 / T(2));
      // [0, 1] not works well, [1, infty[ works better.
      cpoint[1] = 1.5 + zi;
      // x-z plane projection of point p with camera geometry c to z=0.
      // c := camera, p := cpoint.
      // <c + (p - c) * t, [0, 1]> = 0
      const T t(- camera[1] / (cpoint[1] - camera[1]));
      // result increases near to far on differential points.
      result[result.size() - 1 - zi][s] = (camera + (cpoint - camera) * t)[0];
#if defined(_OPENMP)
#pragma omp atomic
#endif
      absmax = max(abs(result[result.size() - 1 - zi][s]), absmax);
    }
  }
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int zi = 0; zi < result.size(); zi ++)
    for(int s = 0; s < result[zi].size(); s ++)
      result[result.size() - 1 - zi][s] *= rstp / absmax;
  return result;
}

template <typename T> const T& PseudoBump<T>::getImgPt(const Vec& img, const T& y) {
  const int& h(img.size());
  const int  yy(abs((int(y + .5) + 3 * h) % (2 * h) - h) % h);
  return img[yy];
}

#define _2D3D_PSEUDO_
#endif

