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
#include <vector>
#include "redig.hh"

using std::cerr;
using std::flush;
using std::complex;
using std::abs;
using std::sqrt;
using std::exp;
using std::pow;
using std::tan;
using std::atan2;
using std::vector;
using std::sort;
using std::ceil;

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
  Mat  compute(const Mat& data, const direction_t& dir);
  
private:
  void initDop(const int& size);
  void initEop(const int& size);
  void initBump(const int& rows, const T& zmax);
  Vec  minSquare(const Vec& in);
  int  getImgPt(const T& y, const T& h);
  void makeDI(const int& size, Vec& Dop, Vec& Iop);
  MatU seed(const int& size, const bool& idft);
  void xchg(Mat& a, Mat& b);
  Vec  Dop0;
  U    I;
  T    Pi;
  Mat  A;
  Mat  D;
  Mat  Dop;
  Mat  Iop;
  Mat  bA;
  Mat  bD;
  Mat  bDop;
  Mat  bIop;
};

template <typename T> enlarger2ex<T>::enlarger2ex() {
  I  = sqrt(U(- 1.));
  Pi = atan2(T(1), T(1)) * T(4);
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
    initEop(data.rows());
    result = D * data;
    break;
  case DETECT_Y:
    initDop(data.rows());
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
      initDop(data.rows());
      Vec ms[data.cols()];
      Mat work(data);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < data.cols(); i ++) {
        ms[i] = minSquare(data.col(i));
        for(int j = 0; j < data.rows(); j ++)
          work(j, i) -= ms[i][0] + ms[i][1] * j;
      }
      result = Iop * work;
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < data.cols(); i ++)
        for(int j = 0; j < data.rows(); j ++)
          result(j, i) += ms[i][0] * j / data.rows() + ms[i][1] * j * j / 2 / data.rows();
    }
    break;
  case BUMP_Y:
    {
      initBump(data.rows(), sqrt(T(data.rows() * data.cols())));
      result = Mat(data.rows(), data.cols());
      assert(A.rows() == result.rows() && A.cols() == result.rows() &&
             data.rows() == result.rows() && data.cols() == result.cols());
      // N.B. local to global, commutative.
      // (uv)'=u'v+uv', u'v=(uv)'-uv',
      //   u' = average(C'*z_k), v = 1 / average(C').
      //   u  = integrate(average(C'*z_k)) = average(integrate(C')*z_k).
      //   v' = - C'' / average^2(C')
      //        integrate(average(C'*z_k) / average(C')) =
      // integrate(u'v) = average(integrate(C')*z_k) / average(C') +
      //        integrate(average(integrate(C')*z_k) / average^2(C') * C'').
            Mat rdata(data.rows(), data.cols());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < rdata.rows(); i ++)
        rdata.row(i) = data.row(data.rows() - 1 - i);
      const Mat tA( A * compute(data, IDETECT_Y));
            Mat work(   compute(data, DETECT_Y));
      const Mat trA(A * compute(rdata, IDETECT_Y));
            Mat rwork(  compute(rdata, DETECT_Y));
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i <= work.rows() / 2; i ++)
        for(int j = 0; j < work.cols(); j ++) {
          work(i, j)  *= tA(i, j)  / pow(data(i, j) + T(1), T(2));
          rwork(i, j) *= trA(i, j) / pow(data(i, j) + T(1), T(2));
      }
      work  = compute( work, IDETECT_Y);
      rwork = compute(rwork, IDETECT_Y);
      for(int i = 0; i <= result.rows() / 2; i ++)
        for(int j = 0; j < result.cols(); j ++) {
          result(i,                     j) = tA(i, j)  / (data(i, j)  + T(1)) + work(i, j);
          result(result.rows() - 1 - i, j) = trA(i, j) / (rdata(i, j) + T(1)) + rwork(i, j);
        }
      result /= - T(2);
    }
    break;
  default:
    assert(0 && "unknown command in enlarger2ex (should not be reached.)");
  }
  return result;
}

template <typename T> void enlarger2ex<T>::initDop(const int& size) {
  cerr << "." << flush;
  if(Dop.rows() == size)
    return;
  xchg(Dop, bDop);
  xchg(Iop, bIop);
  if(Dop.rows() == size)
    return;
  cerr << "new" << flush;
  Vec vDop;
  Vec vIop;
  makeDI(size, vDop, vIop);
  Dop = Mat(size, size);
  Iop = Mat(size, size);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < Dop.rows(); i ++)
    for(int j = 0; j < Dop.cols(); j ++) {
      Dop(i, j) = vDop[(j - i + Dop.cols() * 3 / 2) % Dop.cols()];
      Iop(i, j) = vIop[(j - i + Dop.cols() * 3 / 2) % Iop.cols()];
    }
  return;
}

template <typename T> void enlarger2ex<T>::initEop(const int& size) {
  cerr << "." << flush;
  if(D.rows() == size * 2)
    return;
  xchg(D, bD);
  if(D.rows() == size * 2)
    return;
  cerr << "new" << flush;
  auto DFTa(seed(size, false));
  auto IDFT(seed(size, true));
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < DFTa.rows(); i ++) {
    // N.B. please refer enlarge.wxm, uses each freq.
    //      b(t) -> i * sin(phase) / (cos(phase) - 1) * f(t) for each phase.
    const T phase(Pi * T(i + DFTa.rows()) / T(DFTa.rows()));
    const U r(sqrt(U(- 1)) * sin(phase) / (cos(phase) - U(1)));
    DFTa.row(i) *= r;
  }
  // N.B. : configured ratio.
  const Mat Da((IDFT * DFTa).real().template cast<T>() / (T(2 * DFTa.rows()) * Pi));
  D = Mat(DFTa.rows() * 2, DFTa.cols());
  for(int i = 0; i < DFTa.rows(); i ++) {
    D.row(2 * i + 0) = - Da.row(i);
    D.row(2 * i + 1) =   Da.row(i);
    D(2 * i + 0, i) += T(1);
    D(2 * i + 1, i) += T(1);
  }
  // N.B. shifts, so both sides of image are wrong.
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
  auto DFT2( seed(DFTa.rows() * 2, false));
  auto IDFT2(seed(DFT2.rows(), true));
  for(int i = 0; i < DFT2.rows(); i ++)
    DFT2.row(i) *= T(1.5) - T(.5) * exp(- pow(T(i) / DFTa.rows(), T(2)));
  const Mat sharpen((IDFT2 * DFT2).real().template cast<T>());
  D = sharpen * D;
#endif
  return;
}

template <typename T> void enlarger2ex<T>::initBump(const int& rows, const T& zmax) {
  cerr << "." << flush;
  if(A.rows() == rows)
    return;
  xchg(A, bA);
  if(A.rows() == rows)
    return;

  cerr << "new" << flush;
  assert(0 < rows && T(0) < zmax);
  A = Mat(rows, rows);
  for(int i = 0; i < rows; i ++)
    for(int j = 0; j < rows; j ++)
      A(i, j) = T(0);
  Vec Dop0;
  Vec Iop0;
  makeDI(rows, Dop0, Iop0);
  // Fixed camera, 0 < t < 1 <=> point_z < camera_z
  //             - 1 < t < 0 <=> point_z in [1, 2] * camera_z
  Vec camera(2);
  camera[0] = T(0);
  camera[1] = zmax;
  for(int zi = 0; zi <= zmax * T(2); zi ++)
    for(int j = 0; j < Dop0.size(); j ++) {
      Vec cpoint(2);
      cpoint[0] = j - T(Dop0.size() - 1) / 2;
      cpoint[1] = zi;
      // x-z plane projection of point p with camera geometry c to z=0.
      // c := camera, p := cpoint.
      // <c + (p - c) * t, [0, 1]> = 0
      const auto t(- camera[1] / (cpoint[1] - camera[1]));
      const auto y0((camera + (cpoint - camera) * t)[0]);
      // N.B. average_k(dC_k / dy * z_k).
      for(int i = 0; i < A.rows(); i ++)
        A(i, getImgPt(y0 + i, rows)) += Dop0[j] * (zi + 1);
    }
  T n2(0);
  for(int i = 0; i < A.rows(); i ++)
    n2 += sqrt(A.row(i).dot(A.row(i)));
  A /= n2 / A.rows();
  return;
}

template <typename T> int enlarger2ex<T>::getImgPt(const T& y, const T& h) {
  return int(abs(int(y - pow(h, int(log(y) / log(h))) + .5 + 2 * h * h + h) % int(2 * h) - h)) % int(h);
}

template <typename T> void enlarger2ex<T>::makeDI(const int& size, Vec& Dop, Vec& Iop) {
  assert(6 < size);
  Dop = Vec(size);
  Iop = Vec(size);
  for(int i = 0; i < Dop.size(); i ++)
    Dop[i] = Iop[i] = T(0);
  assert(Dop.size() == size && Iop.size() == size);
  for(int ss = 3; ss <= size / 2; ss ++) {
          auto DFTD(seed(ss, false));
    const auto IDFT(seed(ss, true));
          auto DFTI(DFTD);
    DFTD.row(0) *= U(0);
    DFTI.row(0) *= U(0);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 1; i < DFTD.rows(); i ++) {
      const U phase(- U(2.) * Pi * sqrt(U(- 1)) * T(i) / T(DFTD.rows()));
      DFTD.row(i) *= phase;
      DFTI.row(i) /= phase;
    }
    for(int i = 0; i < IDFT.rows(); i ++) {
      const Vec lDop((IDFT.row(IDFT.rows() - 1 - i) * DFTD).real().template cast<T>());
      const Vec lIop((IDFT.row(IDFT.rows() - 1 - i) * DFTI).real().template cast<T>());
      const T   norm2D(lDop.dot(lDop));
      const T   norm2I(lIop.dot(lIop));
      for(int j = i; j - i < lDop.size(); j ++) {
        Dop[j + size / 2 - IDFT.rows() + 1] += lDop[j - i] / norm2D;
        Iop[j + size / 2 - IDFT.rows() + 1] += lIop[j - i] / norm2I;
      }
    }
  }
  Dop /= sqrt(Dop.dot(Dop));
  Iop /= sqrt(Iop.dot(Iop));
  return;
}

template <typename T> Eigen::Matrix<complex<T>, Eigen::Dynamic, Eigen::Dynamic> enlarger2ex<T>::seed(const int& size, const bool& idft) {
  MatU result(size, size);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < result.rows(); i ++)
    for(int j = 0; j < result.cols(); j ++)
      result(i, j) = exp(U(- 2. * (idft ? - 1 : 1)) * Pi * I * U(i * j / T(size))) / T(idft ? size : 1);
  return result;
}

template <typename T> void enlarger2ex<T>::xchg(Mat& a, Mat& b) {
  Mat work(b);
  b = a;
  a = work;
  return;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> enlarger2ex<T>::minSquare(const Vec& in) {
  T xsum(0);
  T ysum(0);
  T xdot(0);
  T ydot(0);
  for(int i = 0; i < in.size(); i ++) {
    xsum += T(i);
    ysum += in[i];
    xdot += T(i) * T(i);
    ydot += T(i) * in[i];
  }
  Vec result(2);
  result[0] = (xdot * ysum - ydot * xsum) / (in.size() * xdot - xsum * xsum);
  result[1] = (in.size() * ydot - xsum * ysum) / (in.size() * xdot - xsum * xsum);
  return result;
}

#define _ENLARGE2X_
#endif

