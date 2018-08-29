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

#if defined(_WITHOUT_EIGEN_)
#include "simplelin.hh"
#include <complex>
#else
#include <Eigen/Core>
#include <Eigen/LU>
#endif
#include <vector>

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
#if defined(_WITHOUT_EIGEN_)
  typedef SimpleMatrix<T> Mat;
  typedef SimpleMatrix<U> MatU;
  typedef SimpleVector<T> Vec;
  typedef SimpleVector<U> VecU;
#else
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> MatU;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1>              Vec;
  typedef Eigen::Matrix<U, Eigen::Dynamic, 1>              VecU;
#endif
  enlarger2ex();
  Mat  compute(const Mat& data, const direction_t& dir);
  
private:
  void initDop(const int& size);
  void initBump(const int& rows, const T& zmax);
  Vec  minSquare(const Vec& in);
  int  getImgPt(const T& y, const T& h);
  void makeDI(const int& size, Vec& Dop, Vec& Iop, Vec& Eop);
  MatU seed(const int& size, const bool& idft);
  void xchg(Mat& a, Mat& b);
  U    I;
  T    Pi;
  Mat  A;
  Mat  B;
  Mat  Dop;
  Mat  Eop;
  Mat  Iop;
  Mat  bA;
  Mat  bB;
  Mat  bDop;
  Mat  bEop;
  Mat  bIop;
};

template <typename T> enlarger2ex<T>::enlarger2ex() {
  I  = sqrt(U(- 1.));
  Pi = atan2(T(1), T(1)) * T(4);
}

template <typename T> typename enlarger2ex<T>::Mat enlarger2ex<T>::compute(const Mat& data, const direction_t& dir) {
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
    initDop(data.rows());
    result = Eop * data;
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
      result = data;
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < data.cols(); i ++) {
        ms[i] = minSquare(data.col(i));
        for(int j = 0; j < data.rows(); j ++)
          result(j, i) -= ms[i][0] + ms[i][1] * j;
      }
      result = Iop * result;
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
      assert(A.rows() == data.rows() && A.cols() == data.rows());
      // we assume |average(dC * z_k) / average(dC)| be in a differential space.
      Mat dataA(A * data);
      Mat dataB(B * data);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
      for(int j = 0; j < dataA.cols(); j ++) {
        const T offset(sqrt(data.col(j).dot(data.col(j))) / T(2));
        for(int i = 0; i < dataA.rows(); i ++) {
          dataA(i, j) = abs(dataA(i, j)) + offset;
          dataB(i, j) = abs(dataB(i, j)) + offset;
        }
      }
      const Mat datadA(compute(dataA, DETECT_Y));
      const Mat datadB(compute(dataB, DETECT_Y));
      result = Mat(data.rows(), data.cols());
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < result.cols(); i ++) {
        const T offset(data.col(i).dot(data.col(i)));
        for(int j = 0; j < result.rows(); j ++)
          result(j, i) = (datadA(j, i) * dataB(j, i) - dataA(j, i) * datadB(j, i)) / max(dataB(j, i) * dataB(j, i), offset / T(256) / T(256));
      }
      result = compute(result, IDETECT_Y);
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
  xchg(Eop, bEop);
  if(Dop.rows() == size)
    return;
  cerr << "new" << flush;
  Vec vDop;
  Vec vIop;
  Vec vEop;
  makeDI(size, vDop, vIop, vEop);
  Dop = Mat(size, size);
  Iop = Mat(size, size);
  Eop = Mat(size, size);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < Dop.rows(); i ++)
    for(int j = 0; j < Dop.cols(); j ++) {
      Dop(i, j) =   vDop[(j - i + Dop.cols() * 3 / 2) % Dop.cols()];
      Iop(i, j) = - vIop[(j - i + Dop.cols() * 3 / 2) % Iop.cols()];
      Eop(i, j) =   vEop[(j - i + Dop.cols() * 3 / 2) % Eop.cols()];
    }
  Eop *= T(2);
  Mat newEop(Eop.rows() * 2, Eop.cols());
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < Eop.rows(); i ++) {
    newEop.row(2 * i + 0) =   Eop.row(i);
    newEop.row(2 * i + 1) = - Eop.row(i);
    newEop(2 * i + 0, i) += T(1);
    newEop(2 * i + 1, i) += T(1);
  }
  Eop = newEop * T(2);
  for(int i = 0; i < Eop.rows(); i ++) {
    Eop.row(i) += newEop.row(min(i + 1, int(Eop.rows()) - 1));
    Eop.row(i) += newEop.row(max(i - 1, 0));
  }
  Eop /= T(4);
  return;
}

template <typename T> void enlarger2ex<T>::initBump(const int& rows, const T& zmax) {
  cerr << "." << flush;
  if(A.rows() == rows)
    return;
  xchg(A, bA);
  xchg(B, bB);
  if(A.rows() == rows)
    return;

  cerr << "new" << flush;
  assert(0 < rows && T(0) < zmax);
  A = Mat(rows, rows);
  B = Mat(rows, rows);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < rows; i ++)
    for(int j = 0; j < rows; j ++)
      B(i, j) = A(i, j) = T(0);
  Vec Dop0;
  Vec Iop0;
  Vec Eop0;
  makeDI(rows, Dop0, Iop0, Eop0);
  // Fixed camera, 0 < t < 1 <=> point_z < camera_z
  //             - 1 < t < 0 <=> point_z in [1, 2] * camera_z
  Vec camera(2);
  camera[0] = T(0);
  camera[1] = zmax;
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
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
      for(int i = 0; i < A.rows(); i ++) {
#if defined(_OPENMP)
#pragma omp critical
#endif
        {
          A(i, getImgPt(y0 + i, rows)) += Dop0[j] * (zi + 1);
          B(i, getImgPt(y0 + i, rows)) += Dop0[j];
        }
      }
    }
  T n2(0);
  for(int i = 0; i < A.rows(); i ++)
    n2 += sqrt(A.row(i).dot(A.row(i)));
  A /= n2 / A.rows();
  n2 = T(0);
  for(int i = 0; i < B.rows(); i ++)
    n2 += sqrt(B.row(i).dot(B.row(i)));
  B /= n2 / B.rows();
  return;
}

template <typename T> int enlarger2ex<T>::getImgPt(const T& y, const T& h) {
  return int(abs(int(y - pow(h, int(log(y) / log(h))) + .5 + 2 * h * h + h) % int(2 * h) - h)) % int(h);
}

template <typename T> void enlarger2ex<T>::makeDI(const int& size, Vec& Dop, Vec& Iop, Vec& Eop) {
  assert(6 < size);
  Dop = Vec(size);
  Iop = Vec(size);
  Eop = Vec(size);
  for(int i = 0; i < Dop.size(); i ++)
    Dop[i] = Iop[i] = Eop[i] = T(0);
  assert(Dop.size() == size && Iop.size() == size && Eop.size() == size);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int ss = 3; ss <= size / 2; ss ++) {
          auto DFTD(seed(ss, false));
    const auto IDFT(seed(ss, true));
          auto DFTI(DFTD);
          auto DFTE(DFTD);
    DFTD.row(0) *= U(0);
    DFTI.row(0) *= U(0);
    DFTE.row(0) *= U(0);
    T norme2(0);
    for(int i = 1; i < DFTD.rows(); i ++) {
      const U phase(- U(2.) * Pi * sqrt(U(- 1)) * T(i) / T(DFTD.rows()));
      DFTD.row(i) *= phase;
      DFTI.row(i) /= phase;
      // N.B. please refer enlarge.wxm, uses each freq.
      //      b(t) -> i * sin(phase) / (cos(phase) - 1) * f(t) for each phase.
      const T phase2(Pi * T(i + DFTE.rows()) / T(DFTE.rows()));
      const U r(sqrt(U(- 1)) * sin(phase2) / (cos(phase2) - U(1)));
      DFTE.row(i) *= r;
      norme2      += abs(r) * abs(r);
    }
    DFTE /= sqrt(norme2);
#if defined(_WITHOUT_EIGEN_)
    DFTD = DFTD.transpose();
    DFTI = DFTI.transpose();
    DFTE = DFTE.transpose();
#endif
    for(int i = 0; i < IDFT.rows(); i ++) {
      const int iidx(IDFT.rows() - 1 - i);
#if defined(_WITHOUT_EIGEN_)
      const VecU lDop(DFTD * IDFT.row(iidx));
      const VecU lIop(DFTI * IDFT.row(iidx));
      const VecU lEop(DFTE * IDFT.row(iidx));
#else
      const VecU lDop(IDFT.row(iidx) * DFTD);
      const VecU lIop(IDFT.row(iidx) * DFTI);
      const VecU lEop(IDFT.row(iidx) * DFTE);
#endif
#if defined(_OPENMP)
#pragma omp critical
#endif
      for(int j = i; j - i < lDop.size(); j ++) {
        const int idx(j + size / 2 - IDFT.rows() + 1);
        const int jdx(j - i);
        Dop[idx] += T(lDop[jdx].real()) / lDop.size();
        Iop[idx] += T(lIop[jdx].real()) / lIop.size();
        Eop[idx] += T(lEop[jdx].real()) / lEop.size();
      }
    }
  }
  Dop /= size / 2 - 2;
  Iop /= size / 2 - 2;
  Eop /= size / 2 - 2;
  return;
}

template <typename T> typename enlarger2ex<T>::MatU enlarger2ex<T>::seed(const int& size, const bool& idft) {
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

template <typename T> typename enlarger2ex<T>::Vec enlarger2ex<T>::minSquare(const Vec& in) {
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

