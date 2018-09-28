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
using std::vector;

/*
 * This class is NOT thread safe.
 * When applying this to another sizes of images,
 * please create instances for each, or we get slow results.
 */
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
    DETECT_NOP_Y,
    DETECT_BOTH,
    COLLECT_X,
    COLLECT_Y,
    COLLECT_BOTH,
    IDETECT_X,
    IDETECT_Y,
    IDETECT_BOTH,
    BUMP_X,
    BUMP_Y0,
    BUMP_Y,
    BUMP_BOTH,
    EXTEND_X,
    EXTEND_Y,
    EXTEND_BOTH,
    DIV2_X,
    DIV2_Y,
    DIV2_BOTH,
    REVERSE_X,
    REVERSE_Y,
    REVERSE_BOTH,
    CLIP,
    CLIPPM,
    BCLIP,
    ABS,
    EXPSCALE,
    LOGSCALE,
    NORMALIZE,
    NORMALIZE_CLIP } direction_t;
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
  Mat compute(const Mat& data, const direction_t& dir);
  T   dratio;
  T   offset;
  T   blur;
  
private:
  void initDop(const int& size);
  void initBump(const int& rows, const int& cols);
  Vec  minSquare(const Vec& in);
  int  getImgPt(const T& y, const T& h);
  void makeDI(const int& size, Vec& Dop, Vec& Dhop, Vec& Iop, Vec& Eop);
  MatU seed(const int& size, const bool& idft);
  void xchg(Mat& a, Mat& b);
  Mat  round2y(const Mat& in, const int& h);
  U    I;
  T    Pi;
  Mat  A;
  Mat  B;
  Mat  Dop;
  Mat  Dhop;
  Mat  Eop;
  Mat  Iop;
  Mat  bA;
  Mat  bB;
  Mat  bDop;
  Mat  bDhop;
  Mat  bEop;
  Mat  bIop;
};

template <typename T> enlarger2ex<T>::enlarger2ex() {
  I  = sqrt(U(- 1.));
  Pi = atan2(T(1), T(1)) * T(4);
  dratio = T(.0005);
  offset = T(4) / T(256);
  blur   = T(8);
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
    result = (compute(data, DETECT_X)  + compute(data, DETECT_Y)) / T(2);
    break;
  case COLLECT_BOTH:
    result = (compute(data, COLLECT_X) + compute(data, COLLECT_Y)) / T(2);
    break;
  case IDETECT_BOTH:
    result = (compute(data, IDETECT_X) + compute(data, IDETECT_Y)) / T(2);
    break;
  case BUMP_BOTH:
    result = (compute(data, BUMP_Y) + compute(data, BUMP_X)) / T(2);
    break;
  case EXTEND_BOTH:
    result = (compute(compute(data, EXTEND_X), EXTEND_Y) +
              compute(compute(data, EXTEND_Y), EXTEND_X)) / T(2);
    break;
  case DIV2_BOTH:
    result = compute(compute(data, DIV2_X), DIV2_Y);
    break;
  case REVERSE_BOTH:
    result = compute(compute(data, REVERSE_X), REVERSE_Y);
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
  case EXTEND_X:
    result = compute(data.transpose(), EXTEND_Y).transpose();
    break;
  case DIV2_X:
    result = compute(data.transpose(), DIV2_Y).transpose();
    break;
  case REVERSE_X:
    result = compute(data.transpose(), REVERSE_Y).transpose();
    break;
  case ENLARGE_FY:
    result = compute(compute(data, DETECT_Y), ENLARGE_Y);
    break;
  case ENLARGE_Y:
    initDop(data.rows());
    result = compute(Eop * data, CLIP);
    break;
  case DETECT_Y:
    initDop(data.rows());
    result = Dop * data;
    break;
  case DETECT_NOP_Y:
    initDop(data.rows());
    result = Dhop * data;
    break;
  case COLLECT_Y:
    result = compute(compute(data, DETECT_Y), ABS);
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
          result(j, i) += ms[i][0] * j / data.rows() + ms[i][1] * j * j / 2 / data.rows() / data.rows();
    }
    break;
  case BUMP_Y0:
    {
      result = Mat(data.rows(), data.cols());
      initBump(data.rows(), data.cols());
      assert(A.rows() == data.rows() && A.cols() == data.rows());
      // |average(dC*z_k)/average(dC)| == dataA / dataB.
      const auto dataA(compute(A * data, ABS));
            auto dataB(compute(B * data, ABS));
      const auto datadyA(compute(dataA, DETECT_Y));
      const auto datadyB(compute(dataB, DETECT_Y));
      const auto datadxA(compute(dataA, DETECT_X));
      const auto datadxB(compute(dataB, DETECT_X));
      const auto datad2yxA(compute(datadxA, DETECT_Y));
      const auto datad2xyB(compute(datadyB, DETECT_X));
      dataB = compute(dataB, BCLIP);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          // N.B. simply: d/dy(d/dx(A/B)) =
          result(i, j) = (datad2yxA(i, j) * dataB(i, j) + datadxA(i, j) * datadyB(i, j) - datadyA(i, j) * datadxB(i, j) - dataA(i, j) * datad2xyB(i, j)) / pow(dataB(i, j), T(2)) + (datadxA(i, j) * dataB(i, j) + dataA(i, j) * datadxB(i, j)) / pow(dataB(i, j), T(3)) * (- T(2)) * datadyB(i, j);
      result = compute(compute(result, IDETECT_Y), IDETECT_X);
    }
    break;
  case BUMP_Y:
    result = compute(data, BUMP_Y0) + compute(compute(compute(data, REVERSE_X), BUMP_Y0), REVERSE_X) + compute(compute(compute(data, REVERSE_Y), BUMP_Y0), REVERSE_Y) + compute(compute(compute(data, REVERSE_BOTH), BUMP_Y0), REVERSE_BOTH);
    break;
  case EXTEND_Y:
    {
      initDop(data.rows() + 1);
      result = Mat(data.rows() + 1, data.cols());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < data.rows(); i ++)
        result.row(i) = data.row(i);
#if defined(_OPENMP)
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < data.cols(); i ++)
        result(data.rows(), i) = T(0);
      // N.B. nearest data in differential space.
      const auto ddata(Dop * result);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < result.cols(); i ++) {
        T a(0), b(0), c(0);
        for(int j = 0; j < data.rows(); j ++) {
          const T aa(T(2) * Dop(j, Dop.rows() - 1) - Dop(getImgPt(j - 1, Dop.rows()), Dop.rows() - 1) - Dop(getImgPt(j + 1, Dop.rows()), Dop.rows() - 1));
          const T bb(T(2) * ddata(j, i) - ddata(getImgPt(j - 1, ddata.rows()), i) - ddata(getImgPt(j + 1, ddata.rows()), i));
          a += aa * aa;
          b += aa * bb;
          c += bb * bb;
        }
        if(T(0) < b * b - a * c)
          result(data.rows(), i) = - b / a + sqrt(b * b - a * c) / a;
        else
          result(data.rows(), i) = - b / a;
      }
      // XXX fixme: don't know why:
      result.row(data.rows()) /= T(2);
      for(int i = 0; i < result.cols(); i ++)
        result(data.rows(), i) = max(T(0), min(T(1), result(data.rows(), i)));
    }
    break;
  case DIV2_Y:
    {
      result = Mat(data.rows() / 2 + data.rows() % 2, data.cols());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          result(i, j) = T(0);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < data.rows() / 2; i ++)
        result.row(i / 2) += data.row(i) / T(2);
      if(data.rows() % 2)
        result.row(data.rows() / 2) = data.row(data.rows() - 1);
    }
    break;
  case REVERSE_Y:
    result = Mat(data.rows(), data.cols());
    for(int i = 0; i < data.rows(); i ++)
      result.row(result.rows() - i - 1) = data.row(i);
    break;
  case CLIP:
    {
      result = Mat(data.rows(), data.cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          result(i, j) = min(T(1), max(T(0), data(i, j)));
    }
    break;
  case CLIPPM:
    {
      result = Mat(data.rows(), data.cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          result(i, j) = min(T(1), max(- T(1), data(i, j)));
    }
    break;
  case BCLIP:
    {
      result = Mat(data.rows(), data.cols());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          result(i, j) = (data(i, j) < T(0) ? - T(1) : T(1)) * max(abs(data(i, j)), offset);
    }
    break;
  case ABS:
    {
      result = Mat(data.rows(), data.cols());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          result(i, j) = abs(data(i, j));
    }
    break;
  case EXPSCALE:
    {
      result = Mat(data.rows(), data.cols());
      // N.B. might sigmoid is better.
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          result(i, j) = (data(i, j) < T(0) ? - T(1) : T(1)) * exp(T(1) + abs(data(i, j)));
    }
    break;
  case LOGSCALE:
    {
      result = Mat(data.rows(), data.cols());
      // N.B. might be sigmoid is better.
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          result(i, j) = (data(i, j) < T(0) ? - T(1) : T(1)) * (abs(data(i, j)) < exp(T(1)) ? abs(data(i, j)) / exp(T(1)) : log(abs(data(i, j)) + exp(T(1))));
    }
    break;
  case NORMALIZE:
    {
      result = Mat(data.rows(), data.cols());
      T MM(data(0, 0));
      T mm(data(0, 0));
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++) {
          MM = max(MM, data(i, j));
          mm = min(mm, data(i, j));
        }
      if(mm == MM)
        MM += T(1);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          result(i, j) = (data(i, j) - mm) / (MM - mm);
    }
    break;
  case NORMALIZE_CLIP:
    {
      result = Mat(data.rows(), data.cols());
      vector<T> stat;
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          stat.push_back(data(i, j));
      sort(stat.begin(), stat.end());
      auto mm(stat[data.rows() + data.cols()]);
      auto MM(stat[stat.size() - data.rows() - data.cols() - 1]);
      if(mm <= MM + T(1))
        MM += T(1);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          result(i, j) = max(mm, min(MM, (data(i, j) - mm) / (MM - mm)));
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
  xchg(Dhop, bDhop);
  xchg(Iop, bIop);
  xchg(Eop, bEop);
  if(Dop.rows() == size)
    return;
  cerr << "new" << flush;
  Vec vDop;
  Vec vDhop;
  Vec vIop;
  Vec vEop;
  makeDI(size, vDop, vDhop, vIop, vEop);
  Dop  = Mat(size, size);
  Dhop = Mat(size, size);
  Iop  = Mat(size, size);
  Eop  = Mat(size, size);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < Dop.rows(); i ++)
    for(int j = 0; j < Dop.cols(); j ++) {
      Dop(i, j)  =   vDop[ getImgPt(j - i - size / 2, Dop.cols())];
      Dhop(i, j) =   vDhop[getImgPt(j - i - size / 2, Dop.cols())];
      Iop(i, j)  = - vIop[ getImgPt(j - i - size / 2, Dop.cols())];
      Eop(i, j)  =   vEop[ getImgPt(j - i - size / 2, Dop.cols())];
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

template <typename T> void enlarger2ex<T>::initBump(const int& rows, const int& cols) {
  cerr << "." << flush;
  if(A.rows() == rows)
    return;
  xchg(A, bA);
  xchg(B, bB);
  if(A.rows() == rows)
    return;

  cerr << "new" << flush;
  assert(0 < rows);
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
  Vec Dhop0;
  Vec Iop0;
  Vec Eop0;
  makeDI(int(sqrt(T(rows))), Dop0, Dhop0, Iop0, Eop0);
  Vec camera(2);
  camera[0] = T(0);
  camera[1] = T(1);
  T MM(0);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int zi = 0; MM < rows; zi ++)
    for(int j = 0; j < Dop0.size(); j ++) {
      Vec cpoint(2);
      cpoint[0] = j - T(Dop0.size() - 1) / 2;
      cpoint[1] = zi * dratio;
      // x-z plane projection of point p with camera geometry c to z=0.
      // c := camera, p := cpoint.
      // <c + (p - c) * t, [0, 1]> = 0
      const auto t(- camera[1] / (cpoint[1] - camera[1]));
      const auto y0((camera + (cpoint - camera) * t)[0]);
      MM = max(y0, MM);
      // N.B. average_k(dC_k / dy * z_k).
      for(int i = 0; i < A.rows(); i ++) {
#if defined(_OPENMP)
#pragma omp critical
#endif
        {
          A(i, getImgPt(i + y0, rows)) += Dop0[j] * (zi + 1);
          B(i, getImgPt(i + y0, rows)) += Dop0[j];
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

template <typename T> void enlarger2ex<T>::makeDI(const int& size, Vec& Dop, Vec& Dhop, Vec& Iop, Vec& Eop) {
  assert(6 < size);
  Dop  = Vec(size);
  Dhop = Vec(size);
  Iop  = Vec(size);
  Eop  = Vec(size);
  for(int i = 0; i < Dop.size(); i ++)
    Dop[i] = Dhop[i] = Iop[i] = Eop[i] = T(0);
  assert(Dop.size() == size && Dhop.size() == size && Iop.size() == size && Eop.size() == size);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int ss = 3; ss <= size / 2; ss ++) {
          auto DFTD(seed(ss, false));
    const auto IDFT(seed(ss, true));
          auto DFTH(DFTD);
          auto DFTI(DFTD);
          auto DFTE(DFTD);
    DFTD.row(0) *= U(0);
    DFTH.row(0) *= U(0);
    DFTI.row(0) *= U(0);
    DFTE.row(0) *= U(0);
    T norme2(0);
    for(int i = 1; i < DFTD.rows(); i ++) {
      const U phase(- U(2.) * Pi * sqrt(U(- 1)) * T(i) / T(DFTD.rows()));
      // N.B. d/dy, integrate.
      DFTD.row(i) *= phase;
      DFTI.row(i) /= phase;
      // N.B. (d^(log(h))/dy^(log(h)) f, lim h -> 1. : nop.
      DFTH.row(i) *= log(phase);
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
    DFTH = DFTH.transpose();
    DFTI = DFTI.transpose();
    DFTE = DFTE.transpose();
#endif
    for(int i = 0; i < IDFT.rows(); i ++) {
      const int iidx(IDFT.rows() - 1 - i);
#if defined(_WITHOUT_EIGEN_)
      const VecU lDop( DFTD * IDFT.row(iidx));
      const VecU lDhop(DFTH * IDFT.row(iidx));
      const VecU lIop( DFTI * IDFT.row(iidx));
      const VecU lEop( DFTE * IDFT.row(iidx));
#else
      const VecU lDop( IDFT.row(iidx) * DFTD);
      const VecU lDhop(IDFT.row(iidx) * DFTH);
      const VecU lIop( IDFT.row(iidx) * DFTI);
      const VecU lEop( IDFT.row(iidx) * DFTE);
#endif
#if defined(_OPENMP)
#pragma omp critical
#endif
      for(int j = i; j - i < lDop.size(); j ++) {
        const int idx(j + size / 2 - IDFT.rows() + 1);
        const int jdx(j - i);
        Dop[idx]  += T(lDop[ jdx].real()) / lDop.size();
        Dhop[idx] += T(lDhop[jdx].real()) / lDhop.size();
        Iop[idx]  += T(lIop[ jdx].real()) / lIop.size();
        Eop[idx]  += T(lEop[ jdx].real()) / lEop.size();
      }
    }
  }
  Dop  /= size / 2 - 2;
  Dhop /= size / 2 - 2;
  Iop  /= size / 2 - 2;
  Eop  /= size / 2 - 2;
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

template <typename T> typename enlarger2ex<T>::Mat enlarger2ex<T>::round2y(const Mat& in, const int& h) {
  Mat result(min(h, int(in.rows())), in.cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < result.rows(); i ++)
    for(int j = 0; j < result.cols(); j ++)
      result(i, j) = in(i, j);
  return result;
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

