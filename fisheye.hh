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
  void initialize(const int& stp);
  void getPseudoVec(const Mat& in, vector<Vec3>& geoms, vector<Veci3>& delaunay, const int& vbox = 3, const T& rz = T(1) / T(8));
  Mat  getPseudoBump(const Mat& in);
  
private:
  Vec      getPseudoBumpSub(const Vec& work, const Mat& cf);
  const T& getImgPt(const Vec& img, const T& y);
  Mat      prepareLineAxis(const T& rstp);
  
  Mat Dops;
};

template <typename T> PseudoBump<T>::PseudoBump() {
  initialize(21);
}

template <typename T> PseudoBump<T>::~PseudoBump() {
  ;
}

template <typename T> void PseudoBump<T>::initialize(const int& stp) {
  assert(4 < stp);
  const T Pi(T(4) * atan2(T(1), T(1)));
  MatU DFT(stp / 2, stp / 2), IDFT(stp / 2, stp / 2);
  for(int i = 0; i < DFT.rows(); i ++)
    for(int j = 0; j < DFT.cols(); j ++) {
      DFT( i, j) = exp(U(- 2.) * Pi * sqrt(U(- 1)) * U(i * j / T(DFT.rows())));
      IDFT(i, j) = exp(U(  2.) * Pi * sqrt(U(- 1)) * U(i * j / T(DFT.rows()))) / T(DFT.rows());
    }
  for(int i = 0; i < DFT.rows(); i ++)
    DFT.row(i) *= - U(2.) * Pi * sqrt(U(- 1)) * T(i) / T(DFT.rows());
  const Vec DopL((IDFT.row(IDFT.rows() - 1) * DFT).real().template cast<T>());
  const Vec DopM((IDFT.row(IDFT.rows() / 2) * DFT).real().template cast<T>());
  const Vec DopR((IDFT.row(0)               * DFT).real().template cast<T>());
  Dops = Mat(3, stp);
  for(int i = 0; i < Dops.cols(); i ++)
    Dops(0, i) = Dops(1, i) = Dops(2, i) = T(0);
  for(int i = 0; i < DopL.size(); i ++) {
    Dops(0, i)                                     = DopL[i];
    Dops(1, i - DopM.size() / 2 + Dops.cols() / 2) = DopM[i];
    Dops(2, i - DopR.size()     + Dops.cols()    ) = DopR[i];
  }
  cerr << Dops << endl;
  return;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> PseudoBump<T>::getPseudoBumpSub(const Vec& work, const Mat& cf) {
  cerr << "." << flush;
  Vec result(work.size());
  for(int s = 0; s < work.size(); s ++) {
    T sum(0);
    result[s] = T(0);
    for(int z = 0; z < cf.rows(); z ++) {
      // d/dt (local color) / ||local color||:
      Vec c(cf.cols());
      for(int u = 0; u < c.size(); u ++)
        c[u] = getImgPt(work, cf(z, u) + s);
      const auto buf(Dops * c);
      // N.B. <tangents, const distances> / <tangents, 1> is the result,
      //      this is quite pseudo because this similar to
      //      |tan(theta)| * d / |tan(theta)| / 1..
      const auto score(sqrt(buf.dot(buf) / c.dot(c)));
      if(isfinite(score)) {
        result[s] += score * pow(z / T(cf.rows()), T(2));
        sum       += score;
      }
    }
    if(sum != T(0))
      result[s] /= sum;
  }
  return result;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::getPseudoBump(const Mat& in) {
  const auto cf(prepareLineAxis(sqrt(T(in.rows() * in.cols()))));
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

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PseudoBump<T>::prepareLineAxis(const T& rstp) {
  Vec camera(2);
  camera[0] =   T(0);
  camera[1] = - T(1);
  
  Mat result(int(sqrt(rstp) * T(4)), Dops.cols());
  T   absmax(0);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int zi = 0; zi < result.rows(); zi ++)
    for(int s = 0; s < result.cols(); s ++) {
      Vec cpoint(2);
      cpoint[0] = s / T(result.cols() - 1) - 1 / T(2);
      // [1, rstp] works well and scaling works rapidly better.
      cpoint[1] = (T(1) - pow(zi / T(result.rows()), T(2))) * rstp * T(2);
      // x-z plane projection of point p with camera geometry c to z=0.
      // c := camera, p := cpoint.
      // <c + (p - c) * t, [0, 1]> = 0
      const T t(- camera[1] / (cpoint[1] - camera[1]));
      result(zi, s) = (camera + (cpoint - camera) * t)[0];
#if defined(_OPENMP)
#pragma omp atomic
#endif
      absmax = max(abs(result(zi, s)), absmax);
    }
  return result * rstp / absmax;
}

template <typename T> const T& PseudoBump<T>::getImgPt(const Vec& img, const T& y) {
  const int& h(img.size());
  const int  yy(abs((int(y + .5) + 3 * h) % (2 * h) - h) % h);
  return img[yy];
}

#define _2D3D_PSEUDO_
#endif

