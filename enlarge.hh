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
using std::tan;
using std::atan2;

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
    IDETECT_BOTH,
    BUMP_X,
    BUMP_Y,
    BUMP_BOTH } direction_t;
  typedef complex<T> U;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> MatU;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1>              Vec;
  typedef Eigen::Matrix<U, Eigen::Dynamic, 1>              VecU;
  enlarger2ex();
  void initialize(const int& stp);
  Mat  compute(const Mat& data, const direction_t& dir);
  
private:
  void initPattern(const int& size);
  int  getImgPt(const T& y, const T& rows);
  void prepareMat(const int& rows, const T& rstp);
  Mat  Dops;
  U    I;
  T    Pi;
  Mat  A;
  Mat  B;
  Mat  D;
  Mat  Dop;
  Mat  Iop;
  Mat  bA;
  Mat  bB;
  Mat  bD;
  Mat  bDop;
  Mat  bIop;
};

template <typename T> enlarger2ex<T>::enlarger2ex() {
  I  = sqrt(U(- 1.));
  Pi = atan2(T(1), T(1)) * T(4);
  initialize(31);
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
    result = (compute(data, DETECT_X)  + compute(data, DETECT_Y)) / 2.;
    break;
  case COLLECT_BOTH:
    result = (compute(data, COLLECT_X) + compute(data, COLLECT_Y)) / 2.;
    break;
  case IDETECT_BOTH:
    result = (compute(data, IDETECT_X) + compute(data, IDETECT_Y)) / 2.;
    break;
  case BUMP_BOTH:
    result = (compute(data, BUMP_X)    + compute(data, BUMP_Y)) / 2.;
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
  case BUMP_X:
    result = compute(data.transpose(), BUMP_Y).transpose();
    break;
  case ENLARGE_FY:
    result = compute(compute(data, DETECT_Y), ENLARGE_Y);
    break;
  case ENLARGE_Y:
    initPattern(data.rows());
    result = D * data;
    break;
  case DETECT_Y:
    initPattern(data.rows());
    result = Dop * data;
    break;
  case COLLECT_Y:
    result = compute(data, DETECT_Y);
    {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          result(i, j) = abs(result(i, j));
    }
    break;
  case IDETECT_Y:
    {
      initPattern(data.rows());
      Vec avg(data.cols());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < data.cols(); i ++) {
        avg[i]  = T(0);
        for(int j = 0; j < data.rows(); j ++)
          avg[i] += data(j, i);
        avg[i] /= data.rows();
      }
      result = Iop * data;
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < data.cols(); i ++)
        for(int j = 0; j < data.rows(); j ++)
          result(j, i) += avg[i] * j * j / data.rows() / data.rows();
    }
    break;
  case BUMP_Y:
    {
      prepareMat(data.rows(), sqrt(T(data.rows() * data.cols())));
      result = Mat(data.rows(), data.cols());
      assert(A.rows() == result.rows() && A.cols() == result.rows() &&
             B.rows() == result.rows() && B.cols() == result.rows() &&
             data.rows() == result.rows() && data.cols() == result.cols());
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          // XXX : is tensor multiplication norm better?
          result(i, j) = abs(A.col(i).dot(data.col(j)) / B.col(i).dot(data.col(j)));
    }
    break;
  default:
    assert(0 && "unknown command in enlarger2ex (should not be reached.)");
  }
  return result;
}

template <typename T> void enlarger2ex<T>::initialize(const int& stp) {
  assert(4 < stp);
  MatU DFT(stp / 2, stp / 2), IDFT(stp / 2, stp / 2);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < DFT.rows(); i ++)
    for(int j = 0; j < DFT.cols(); j ++) {
      DFT( i, j) = exp(U(- 2.) * Pi * sqrt(U(- 1)) * U(i * j / T(DFT.rows())));
      IDFT(i, j) = exp(U(  2.) * Pi * sqrt(U(- 1)) * U(i * j / T(DFT.rows()))) / T(DFT.rows());
    }
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < DFT.rows(); i ++)
    DFT.row(i) *= - U(2.) * Pi * sqrt(U(- 1)) * T(i) / T(DFT.rows());
  Dops = Mat(IDFT.rows(), stp);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < Dops.rows(); i ++)
    for(int j = 0; j < Dops.cols(); j ++)
      Dops(i, j) = T(0);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < IDFT.rows(); i ++) {
    const Vec lDop((IDFT.row(IDFT.rows() - 1 - i) * DFT).real().template cast<T>());
    for(int j = i; j < lDop.size(); j ++)
      Dops(i, j) = lDop[j - i];
  }
  return;
}

template <typename T> void enlarger2ex<T>::initPattern(const int& size) {
  cerr << "." << flush;
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
  cerr << "new" << flush;
  MatU DFT( size, size);
  MatU IDFT(size, size);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < DFT.rows(); i ++)
    for(int j = 0; j < DFT.cols(); j ++) {
      DFT( i, j) = exp(U(- 2.) * Pi * I * U(i * j / T(size)));
      IDFT(i, j) = exp(U(  2.) * Pi * I * U(i * j / T(size))) / T(size);
    }
  MatU Dbuf(DFT); MatU Ibuf(DFT);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < Dbuf.rows(); i ++)
    Dbuf.row(i) *= U(- 2.) * Pi * I * U(i / T(size));
  Ibuf.row(0) *= T(0);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 1; i < Ibuf.rows(); i ++)
    Ibuf.row(i) /= U(- 2.) * Pi * I * U(i / T(size));
  Dop =   (IDFT * Dbuf).real().template cast<T>();
  Iop = - (IDFT * Ibuf).real().template cast<T>();
  {
    MatU DFTa(DFT);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
    for(int i = 0; i < DFTa.rows(); i ++) {
      // N.B. please refer enlarge.wxm, uses each freq.
      //      b(t) -> i * sin(phase) / (cos(phase) - 1) * f(t) for each phase.
      const T phase(Pi * T(i + DFT.rows()) / T(DFT.rows()));
      const U r(sqrt(U(- 1)) * sin(phase) / (cos(phase) - U(1)));
      DFTa.row(i) *= r;
    }
    // N.B. : configured ratio.
    const Mat Da((IDFT * DFTa).real().template cast<T>() / (T(2 * DFT.rows()) * Pi));
    D = Mat(DFT.rows() * 2, DFT.cols());
    for(int i = 0; i < DFT.rows(); i ++) {
      D.row(2 * i + 0) = - Da.row(i);
      D.row(2 * i + 1) =   Da.row(i);
      D(2 * i + 0, i) += T(1);
      D(2 * i + 1, i) += T(1);
    }
    // N.B. shif, so both sides of image are wrong.
    const Mat D0(D * T(2));
    for(int i = 0; i < D.rows(); i ++)
      for(int j = 0; j < D.cols(); j ++) {
        D(i, j)                  += D0((i + 1) % D0.rows(), j);
        D((i + 1) % D.rows(), j) += D0(i, j);
      }
    D /= T(4);
#if defined(_WITH_EXTERNAL_)
    // This works perfectly (from referring https://web.stanford.edu/class/cs448f/lectures/2.1/Sharpening.pdf via reffering Q&A sites.).
    // But I don't know whether this method is open or not.
    MatU DFT2(DFT.rows() * 2, DFT.cols() * 2);
    MatU IDFT2(DFT2.rows(), DFT2.cols());
    for(int i = 0; i < DFT2.rows(); i ++)
      for(int j = 0; j < DFT2.cols(); j ++) {
        DFT2( i, j) = exp(U(- 2.) * Pi * I * U(i * j / T(size * 2)));
        IDFT2(i, j) = exp(U(  2.) * Pi * I * U(i * j / T(size * 2))) / T(size * 2);
      }
    for(int i = 0; i < DFT2.rows(); i ++)
      DFT2.row(i) *= T(1.5) - T(.5) * exp(- pow(T(i) / DFT.rows(), T(2)));
    const Mat sharpen((IDFT2 * DFT2).real().template cast<T>());
    D = sharpen * D;
#endif
  }
  return;
}

template <typename T> void enlarger2ex<T>::prepareMat(const int& rows, const T& rstp) {
  cerr << "." << flush;
  if(A.rows() == rows)
    return;
  Mat work(bA);
  bA   = A;
  A    = work;
  work = bB;
  bB   = B;
  B    = work;
  if(A.rows() == rows)
    return;
  cerr << "new" << flush;
  
  Vec camera(2);
  camera[0] = T(0);
  camera[1] = rstp;
  
  work = Mat(int(sqrt(rstp) * T(4)), Dops.cols());
  T absmax(0);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int zi = 0; zi < work.rows(); zi ++)
    for(int s = 0; s < work.cols(); s ++) {
      Vec cpoint(2);
      cpoint[0] = s / T(work.cols() - 1) - 1 / T(2);
      // [1, rstp] works well.
      cpoint[1] = (zi + T(.5)) / T(work.rows()) * rstp;
      // x-z plane projection of point p with camera geometry c to z=0.
      // c := camera, p := cpoint.
      // <c + (p - c) * t, [0, 1]> = 0
      const T t(- camera[1] / (cpoint[1] - camera[1]));
      work(zi, s) = (camera + (cpoint - camera) * t)[0];
#if defined(_OPENMP)
#pragma omp atomic
#endif
      absmax = max(abs(work(zi, s)), absmax);
    }
  work *= rstp / absmax;
  A = Mat(rows, rows);
  B = Mat(rows, rows);
  for(int i = 0; i < rows; i ++)
    for(int j = 0; j < rows; j ++)
      A(i, j) = B(i, j) = T(0);
  for(int i = 1; i < max(sqrt(sqrt(rstp)), T(1)); i ++)
    for(int j = 0; j < work.rows(); j ++)
      for(int k = 0; k < work.cols(); k ++)
        for(int ii = 0; ii < rows; ii ++)
          for(int jj = 0; jj < Dops.rows(); jj ++) {
            A(getImgPt(ii + work(j, k) * T(i) / sqrt(sqrt(rstp)), rows), ii) += Dops(jj, k) * (jj + T(1)) / Dops.rows();
            B(getImgPt(ii + work(j, k) * T(i) / sqrt(sqrt(rstp)), rows), ii) += T(1);
          }
  return;
}

template <typename T> int enlarger2ex<T>::getImgPt(const T& y, const T& rows) {
  const int& h(rows);
  return abs((int(y + .5) + 3 * h) % (2 * h) - h) % h;
}

#define _ENLARGE2X_
#endif

