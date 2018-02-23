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

#if !defined(_ENLARGE2X_)

#include <Eigen/Core>
#include <Eigen/LU>

using std::cerr;
using std::flush;
using std::complex;
using std::abs;
using std::sqrt;
using std::exp;
using std::pow;

template <typename T> class enlarger2ex {
public:
  typedef enum {
    ENLARGE_X,
    ENLARGE_Y,
    ENLARGE_FX,
    ENLARGE_FY,
    ENLARGE_BOTH,
    ENLARGE_FBOTH,
    ENLARGE_3BOTH,
    DETECT_X,
    DETECT_Y,
    DETECT_BOTH,
    COLLECT_X,
    COLLECT_Y,
    COLLECT_BOTH,
    IDETECT_X,
    IDETECT_Y,
    IDETECT_BOTH } direction_t;
  typedef complex<T> U;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> MatU;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1>              Vec;
  typedef Eigen::Matrix<U, Eigen::Dynamic, 1>              VecU;
  enlarger2ex();
  Mat compute(const Mat& data, const direction_t& dir);
  
private:
  void initPattern(const int& size);
  U    I;
  T    Pi;
  Mat  D;
  Mat  Dop;
  Mat  Iop;
  Mat  bD;
  Mat  bDop;
  Mat  bIop;
};

template <typename T> enlarger2ex<T>::enlarger2ex() {
  I    = sqrt(U(- 1.));
  Pi   = atan2(T(1), T(1)) * T(4);
  D    = Mat();
  Dop  = Mat();
  Iop  = Mat();
  bD   = Mat();
  bDop = Mat();
  bIop = Mat();
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> enlarger2ex<T>::compute(const Mat& data, const direction_t& dir) {
  Mat result;
  switch(dir) {
  case ENLARGE_BOTH:
    result = (compute(compute(data, ENLARGE_X), ENLARGE_Y) +
              compute(compute(data, ENLARGE_Y), ENLARGE_X)) / T(2);
    break;
  case ENLARGE_FBOTH:
    result = (compute(compute(data, ENLARGE_FX), ENLARGE_FY) +
              compute(compute(data, ENLARGE_FY), ENLARGE_FX)) / T(2);
    break;
  case ENLARGE_3BOTH:
    result = compute(data, ENLARGE_BOTH) +
             compute(data, ENLARGE_FBOTH) / T(2) / T(3);
    break;
  case DETECT_BOTH:
    result = (compute(data, DETECT_X) + compute(data, DETECT_Y)) / 2.;
    break;
  case COLLECT_BOTH:
    result = (compute(data, COLLECT_X) + compute(data, COLLECT_Y)) / 2.;
    break;
  case IDETECT_BOTH:
    result = (compute(data, IDETECT_X) + compute(data, IDETECT_Y)) / 2.;
    break;
  case ENLARGE_X:
    result = compute(data.transpose(), ENLARGE_Y).transpose();
    break;
  case ENLARGE_FX:
    result = compute(data.transpose(), ENLARGE_FY).transpose();
    break;
  case DETECT_X:
    result = compute(data.transpose(), DETECT_Y).transpose();
    break;
  case COLLECT_X:
    result = compute(data.transpose(), COLLECT_Y).transpose();
    break;
  case IDETECT_X:
    result = compute(data.transpose(), IDETECT_Y).transpose();
    break;
  case ENLARGE_Y:
    {
      initPattern(data.rows());
      result = D * data;
      for(int i = 0; i < result.rows(); i ++)
        result.row(i) += data.row(i / 2);
    }
    break;
  case ENLARGE_FY:
    {
      initPattern(data.rows());
      const Mat dcache(Dop * data);
      T ndc(0), nd(0);
      for(int i = 1; i < data.rows(); i ++)
        for(int j = 0; j < data.cols(); j ++)
          nd = max(abs(data(i, j) - data(i - 1, j)), nd);
      for(int i = 0; i < data.rows(); i ++)
        for(int j = 1; j < data.cols(); j ++)
          nd = max(abs(data(i, j) - data(i, j - 1)), nd);
      for(int i = 0; i < dcache.rows(); i ++)
        for(int j = 0; j < dcache.cols(); j ++)
          ndc = max(abs(dcache(i, j)), ndc);
      result = compute(dcache, ENLARGE_Y);
      for(int i = 0; i < result.rows(); i ++)
        result.row(i) += dcache.row(i / 2);
      result *= nd / ndc;
    }
    break;
  case DETECT_Y:
    initPattern(data.rows());
    result = Dop * data;
    break;
  case COLLECT_Y:
    result = compute(data, DETECT_Y);
    for(int i = 0; i < result.rows(); i ++)
      for(int j = 0; j < result.cols(); j ++)
        result(i, j) = abs(result(i, j));
    break;
  case IDETECT_Y:
    {
      initPattern(data.rows());
      Vec avg(data.cols());
      for(int i = 0; i < data.cols(); i ++) {
        avg[i]  = T(0);
        for(int j = 0; j < data.rows(); j ++)
          avg[i] += data(j, i);
        avg[i] /= data.rows();
      }
      result = Iop * data;
      for(int i = 0; i < data.cols(); i ++)
        for(int j = 0; j < data.rows(); j ++)
          result(j, i) += avg[i] * j * j / data.rows() / data.rows();
    }
    break;
  default:
    ;
  }
  return result;
}

template <typename T> void enlarger2ex<T>::initPattern(const int& size) {
  cerr << " initPat" << flush;
  if(Dop.rows() == size)
    return;
  Mat work(bD);
  bD   = D;
  D    = work;
  work = bDop;
  bDop = Dop;
  Dop  = work;
  work = bIop;
  bIop = Iop;
  Iop  = work;
  if(Dop.rows() == size)
    return;
  cerr << " new" << flush;
  MatU DFT( size, size);
  MatU IDFT(size, size);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < DFT.rows(); i ++)
    for(int j = 0; j < DFT.cols(); j ++) {
      DFT( i, j) = exp(U(- 2.) * Pi * I * U(i * j / T(size)));
      IDFT(i, j) = exp(U(  2.) * Pi * I * U(i * j / T(size))) / T(size);
    }
  MatU Dbuf(DFT);
  MatU Ibuf(DFT);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < Dbuf.rows(); i ++)
    Dbuf.row(i) *= U(- 2.) * Pi * I * U(i / T(size));
  Ibuf.row(0) *= T(0);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 1; i < Ibuf.rows(); i ++)
    Ibuf.row(i) /= U(- 2.) * Pi * I * U(i / T(size));
  Dop =   (IDFT * Dbuf).real().template cast<T>();
  Iop = - (IDFT * Ibuf).real().template cast<T>();
  {
    MatU DFTa(DFT), DFTb(DFT);
    for(int i = 0; i < DFT.rows(); i ++) {
      const U phase(T(2) * (i + .5) * Pi / DFT.rows());
      // N.B. refer enlarge.wxm, ideally U(1) but it is flat,
      //      uses with skewness with U(.5) and U(2).
      const T r((sqrt(exp(phase) / (exp(phase) - U(.5)) *
                      exp(phase) / (exp(phase) - U(2.)))).real());
      DFTa.row(i) *= U(       r);
      DFTb.row(i) *= U(T(1) - r);
    }
    const Mat Da((IDFT * DFTa).real().template cast<T>());
    const Mat Db((IDFT * DFTb).real().template cast<T>());
    D = Mat(DFT.rows() * 2, DFT.cols());
    for(int i = 0; i < DFT.rows(); i ++) {
      D.row(i * 2 + 0) = Da.row(i);
      D.row(i * 2 + 1) = Db.row(i);
    }
    D /= T(2) * Pi * T(2);
  }
  return;
}

#define _ENLARGE2X_
#endif

