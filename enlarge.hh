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
using std::max;
using std::min;

/*
 * This class is NOT thread safe.
 * When applying this to another sizes of images,
 * please create instances for each, or we get slow results.
 */
template <typename T> class enlarger2ex {
public:
  typedef enum {
    ENLARGE_X,
    ENLARGE_XQ,
    ENLARGE_XQS,
    ENLARGE_Y,
    ENLARGE_YQ,
    ENLARGE_YQS,
    ENLARGE_FX,
    ENLARGE_FY,
    ENLARGE_BOTH,
    ENLARGE_BOTHQ,
    ENLARGE_BOTHQS,
    ENLARGE_FBOTH,
    ENLARGE_3BOTH,
    DETECT_X,
    DETECT_XQ,
    DETECT_XQS,
    DETECT_Y,
    DETECT_YQ,
    DETECT_YQS,
    DETECT_NOP_Y,
    DETECT_BOTH,
    DETECT_BOTHQ,
    DETECT_BOTHQS,
    COLLECT_X,
    COLLECT_XQ,
    COLLECT_XQS,
    COLLECT_Y,
    COLLECT_YQ,
    COLLECT_YQS,
    COLLECT_BOTH,
    COLLECT_BOTHQ,
    COLLECT_BOTHQS,
    IDETECT_X,
    IDETECT_XQ,
    IDETECT_XQS,
    IDETECT_X_BOTH,
    IDETECT_Y,
    IDETECT_YQ,
    IDETECT_YQS,
    IDETECT_Y_BOTH,
    IDETECT_BOTH,
    IDETECT_BOTHQ,
    IDETECT_BOTHQS,
    IDETECT_QUAD,
    BLUR_X,
    BLUR_Y,
    BLUR_BOTH,
    BUMP_X,
    BUMP_XQ,
    BUMP_XQS,
    BUMP_Y00,
    BUMP_Y0,
    BUMP_Y1,
    BUMP_Y2,
    BUMP_Y,
    BUMP_YP,
    BUMP_YQ,
    BUMP_YQS,
    BUMP_BOTH,
    BUMP_BOTHQ,
    BUMP_BOTHQS,
    EXTEND_X,
    EXTEND_XQ,
    EXTEND_XQS,
    EXTEND_Y0,
    EXTEND_Y,
    EXTEND_YQ,
    EXTEND_YQS,
    EXTEND_BOTH,
    EXTEND_BOTHQ,
    EXTEND_BOTHQS,
    DIV2_X,
    DIV2_Y,
    DIV2_BOTH,
    REVERSE_X,
    REVERSE_Y,
    REVERSE_BOTH,
    CLIP,
    CLIPPM,
    CLIPPLUS,
    BCLIP,
    ABS,
    SQRTSCALE,
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
  Mat  compute(const Mat& data, const direction_t& dir);
  MatU seed(const int& size, const bool& idft);
  T    dratio;
  T    offset;
  T    blur;
  int  sz_cell;
  int  st_cell;
  
private:
  void initDop(const int& size);
  void initBump(const int& rows, const int& cols);
  Vec  minSquare(const Vec& in);
  int  getImgPt(const T& y, const T& h);
  void makeDI(const int& size, Vec& Dop, Vec& Dhop, Vec& Iop, Vec& Eop);
  void xchg(Mat& a, Mat& b);
  Mat  round2y(const Mat& in, const int& h);
  Mat  recursive(const Mat& data, const direction_t& dir, const direction_t& dir0);
  Mat  recursiveSumup(const Mat& data, const direction_t& dir, const direction_t& dir0);
  U    I;
  T    Pi;
  Mat  A;
  Mat  dataA;
  Mat  B;
  Mat  C;
  Mat  Dop;
  Mat  Dhop;
  Mat  Eop;
  Mat  Iop;
  Mat  bA;
  Mat  bB;
  Mat  bC;
  Mat  bDop;
  Mat  bDhop;
  Mat  bEop;
  Mat  bIop;
};

template <typename T> enlarger2ex<T>::enlarger2ex() {
  I  = sqrt(U(- 1.));
  Pi = atan2(T(1), T(1)) * T(4);
  dratio  = T(.0005);
  offset  = T(4) / T(256);
  blur    = T(8);
  sz_cell = 96;
  st_cell = 12;
}

template <typename T> typename enlarger2ex<T>::Mat enlarger2ex<T>::compute(const Mat& data, const direction_t& dir) {
  Mat result;
  switch(dir) {
  case ENLARGE_BOTH:
    result = (compute(compute(data, ENLARGE_X), ENLARGE_Y) +
              compute(compute(data, ENLARGE_Y), ENLARGE_X)) / T(2);
    break;
  case ENLARGE_BOTHQ:
    result = (compute(compute(data, ENLARGE_XQ), ENLARGE_YQ) +
              compute(compute(data, ENLARGE_YQ), ENLARGE_XQ)) / T(2);
    break;
  case ENLARGE_BOTHQS:
    result = (compute(compute(data, ENLARGE_XQS), ENLARGE_YQS) +
              compute(compute(data, ENLARGE_YQS), ENLARGE_XQS)) / T(2);
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
  case DETECT_BOTHQ:
    result = (compute(data, DETECT_XQ) + compute(data, DETECT_YQ)) / T(2);
    break;
  case DETECT_BOTHQS:
    result = (compute(data, DETECT_XQS) + compute(data, DETECT_YQS)) / T(2);
    break;
  case COLLECT_BOTH:
    result = (compute(data, COLLECT_X) + compute(data, COLLECT_Y)) / T(2);
    break;
  case COLLECT_BOTHQ:
    result = (compute(data, COLLECT_XQ) + compute(data, COLLECT_YQ)) / T(2);
    break;
  case COLLECT_BOTHQS:
    result = (compute(data, COLLECT_XQS) + compute(data, COLLECT_YQS)) / T(2);
    break;
  case IDETECT_BOTH:
    result = (compute(data, IDETECT_X) + compute(data, IDETECT_Y)) / T(2);
    break;
  case IDETECT_BOTHQ:
    result = (compute(data, IDETECT_XQ) + compute(data, IDETECT_YQ)) / T(2);
    break;
  case IDETECT_BOTHQS:
    result = (compute(data, IDETECT_XQS) + compute(data, IDETECT_YQS)) / T(2);
    break;
  case IDETECT_QUAD:
    result = (compute(data, IDETECT_Y_BOTH) + compute(data, IDETECT_X_BOTH)) / T(2);
    break;
  case BLUR_BOTH:
    result = (compute(data, BLUR_X)    + compute(data, BLUR_Y)) / T(2);
    break;
  case BUMP_BOTH:
    result = (compute(data, BUMP_Y)    + compute(data, BUMP_X)) / T(2);
    break;
  case EXTEND_BOTH:
    result = (compute(compute(data, EXTEND_X), EXTEND_Y) +
              compute(compute(data, EXTEND_Y), EXTEND_X)) / T(2);
    break;
  case EXTEND_BOTHQ:
    result = (compute(compute(data, EXTEND_XQ), EXTEND_YQ) +
              compute(compute(data, EXTEND_YQ), EXTEND_XQ)) / T(2);
    break;
  case EXTEND_BOTHQS:
    result = (compute(compute(data, EXTEND_XQS), EXTEND_YQS) +
              compute(compute(data, EXTEND_YQS), EXTEND_XQS)) / T(2);
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
  case ENLARGE_XQ:
    result = compute(data.transpose(), ENLARGE_YQ).transpose();
    break;
  case ENLARGE_XQS:
    result = compute(data.transpose(), ENLARGE_YQS).transpose();
    break;
  case ENLARGE_FX:
    result = compute(data.transpose(), ENLARGE_FY).transpose();
    break;
  case DETECT_X:
    result = compute(data.transpose(), DETECT_Y).transpose();
    break;
  case DETECT_XQ:
    result = compute(data.transpose(), DETECT_YQ).transpose();
    break;
  case DETECT_XQS:
    result = compute(data.transpose(), DETECT_YQS).transpose();
    break;
  case COLLECT_X:
    result = compute(data.transpose(), COLLECT_Y).transpose();
    break;
  case COLLECT_XQ:
    result = compute(data.transpose(), COLLECT_YQ).transpose();
    break;
  case COLLECT_XQS:
    result = compute(data.transpose(), COLLECT_YQS).transpose();
    break;
  case IDETECT_X:
    result = compute(data.transpose(), IDETECT_Y).transpose();
    break;
  case IDETECT_XQ:
    result = compute(data.transpose(), IDETECT_YQ).transpose();
    break;
  case IDETECT_XQS:
    result = compute(data.transpose(), IDETECT_YQS).transpose();
    break;
  case IDETECT_X_BOTH:
    result = compute(data.transpose(), IDETECT_Y_BOTH).transpose();
    break;
  case BLUR_X:
    result = compute(data.transpose(), BLUR_Y).transpose();
    break;
  case BUMP_X:
    result = compute(data.transpose(), BUMP_Y).transpose();
    break;
  case BUMP_XQ:
    result = compute(data.transpose(), BUMP_YQ).transpose();
    break;
  case BUMP_XQS:
    result = compute(data.transpose(), BUMP_YQS).transpose();
    break;
  case EXTEND_X:
    result = compute(data.transpose(), EXTEND_Y).transpose();
    break;
  case EXTEND_XQ:
    result = compute(data.transpose(), EXTEND_YQ).transpose();
    break;
  case EXTEND_XQS:
    result = compute(data.transpose(), EXTEND_YQS).transpose();
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
  case ENLARGE_YQ:
    result = recursive(data, ENLARGE_YQ, ENLARGE_Y);
    break;
  case ENLARGE_YQS:
    result = recursiveSumup(data, ENLARGE_YQ, ENLARGE_Y);
    break;
  case DETECT_Y:
    initDop(data.rows());
    result = Dop * data;
    break;
  case DETECT_YQ:
    result = recursive(data, DETECT_YQ, DETECT_Y);
    break;
  case DETECT_YQS:
    result = recursiveSumup(data, DETECT_YQ, DETECT_Y);
    break;
  case DETECT_NOP_Y:
    initDop(data.rows());
    result = Dhop * data;
    break;
  case BLUR_Y:
    initDop(data.rows());
    result = C * data;
    break;
  case COLLECT_Y:
    result = compute(compute(data, DETECT_Y), ABS);
    break;
  case COLLECT_YQ:
    result = compute(compute(data, DETECT_YQ), ABS);
    break;
  case COLLECT_YQS:
    result = compute(compute(data, DETECT_YQS), ABS);
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
  case IDETECT_YQ:
    result = recursive(data, IDETECT_YQ, IDETECT_Y);
    break;
  case IDETECT_YQS:
    result = recursiveSumup(data, IDETECT_YQS, IDETECT_Y);
    break;
  case IDETECT_Y_BOTH:
    result = compute(data, IDETECT_Y) + compute(compute(compute(data, REVERSE_Y), IDETECT_Y), REVERSE_Y);
    break;
  case BUMP_Y00:
    {
      assert(A.rows() == data.rows() && A.cols() == data.rows());
      assert(dataA.rows() == data.rows() && dataA.cols() == data.cols());
      // |average(dC*z_k)/average(dC)| == dataA / dataB.
      const auto dataB(compute(B * data, ABS));
      const auto dataBc(compute(dataB, BCLIP));
      const auto datadyA(compute(dataA, DETECT_Y));
      const auto datadyB(compute(dataB, DETECT_Y));
      const auto datadxA(compute(dataA, DETECT_X));
      const auto datadxB(compute(dataB, DETECT_X));
      const auto datad2yyA(compute(datadyA, DETECT_Y));
      const auto datad2yyB(compute(datadyB, DETECT_Y));
      const auto datad2xyA(compute(datadyA, DETECT_X));
      const auto datad2xyB(compute(datadyB, DETECT_X));
      const auto datad2xxA(compute(datadxA, DETECT_X));
      const auto datad2xxB(compute(datadxB, DETECT_X));
      const auto datad2yxA(compute(datadxA, DETECT_Y));
      const auto datad2yxB(compute(datadxB, DETECT_Y));
      Mat d2xx(data.rows(), data.cols());
      Mat d2xy(data.rows(), data.cols());
      Mat d2yx(data.rows(), data.cols());
      Mat d2yy(data.rows(), data.cols());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < data.rows(); i ++)
        for(int j = 0; j < data.cols(); j ++) {
          d2xx(i, j) = (datad2xxA(i, j) * dataB(i, j) + datadxA(i, j) * datadxB(i, j) - datadxA(i, j) * datadxB(i, j) - dataA(i, j) * datad2xxB(i, j)) / pow(dataBc(i, j), T(2)) + (datadxA(i, j) * dataB(i, j) + dataA(i, j) * datadxB(i, j)) / pow(dataBc(i, j), T(3)) * (- T(2)) * datadxB(i, j);
          d2xy(i, j) = (datad2xyA(i, j) * dataB(i, j) + datadyA(i, j) * datadxB(i, j) - datadxA(i, j) * datadyB(i, j) - dataA(i, j) * datad2xyB(i, j)) / pow(dataBc(i, j), T(2)) + (datadyA(i, j) * dataB(i, j) + dataA(i, j) * datadyB(i, j)) / pow(dataBc(i, j), T(3)) * (- T(2)) * datadxB(i, j);
          d2yx(i, j) = (datad2yxA(i, j) * dataB(i, j) + datadxA(i, j) * datadyB(i, j) - datadyA(i, j) * datadxB(i, j) - dataA(i, j) * datad2yxB(i, j)) / pow(dataBc(i, j), T(2)) + (datadxA(i, j) * dataB(i, j) + dataA(i, j) * datadxB(i, j)) / pow(dataBc(i, j), T(3)) * (- T(2)) * datadyB(i, j);
          d2yy(i, j) = (datad2yyA(i, j) * dataB(i, j) + datadyA(i, j) * datadyB(i, j) - datadyA(i, j) * datadyB(i, j) - dataA(i, j) * datad2yyB(i, j)) / pow(dataBc(i, j), T(2)) + (datadyA(i, j) * dataB(i, j) + dataA(i, j) * datadyB(i, j)) / pow(dataBc(i, j), T(3)) * (- T(2)) * datadyB(i, j);
        }
      // N.B. (d/dt)^2 f = d/dt(df/dx dx/dt + df/dy dy/dt)
      // = d^2f/dx^2 dx/dt + df/dx d^2x/dt^2 +
      //   d^2f/dydx dx/dt + df/dx dxdy/dt^2 +
      //   d^2f/dxdy dy/dt + df/dy dydx/dt^2 +
      //   d^2f/dy^2 dy/dt + df/dy d^2y/dt^2
      // dx/dt == dy/dt == const. ->
      // result(i, j) = d2dxdx + d2dxdy + d2dydx + d2dydy;
      // N.B. in fact, in this case, so function don't have rotation related
      //      information, rotation + d/dt sumup is needed, but now, not so.
      result = compute(compute(compute(d2xx, IDETECT_X), IDETECT_X) + compute(compute(d2xy, IDETECT_Y), IDETECT_X) + compute(compute(d2yx, IDETECT_X), IDETECT_Y) + compute(compute(d2yy, IDETECT_Y), IDETECT_Y), LOGSCALE);
    }
    break;
  case BUMP_Y0:
    initBump(data.rows(), data.cols());
    dataA  = compute(A * data, ABS);
    result = compute(data, BUMP_Y00);
    break;
  case BUMP_Y1:
    initBump(data.rows(), data.cols());
    dataA  = - compute(A * data, ABS);
    result = - compute(data, BUMP_Y00);
    break;
  case BUMP_Y2:
    result = compute(data, BUMP_Y0) + compute(data, BUMP_Y1);
    break;
  case BUMP_Y:
    result = compute(data, BUMP_Y2) + compute(compute(compute(data, REVERSE_X), BUMP_Y2), REVERSE_X) + compute(compute(compute(data, REVERSE_Y), BUMP_Y2), REVERSE_Y) + compute(compute(compute(data, REVERSE_BOTH), BUMP_Y2), REVERSE_BOTH);
    break;
  case BUMP_YQ:
    result = recursive(data, BUMP_YQ, BUMP_Y);
    break;
  case BUMP_YQS:
    result = recursiveSumup(data, BUMP_YQ, BUMP_Y);
    break;
  case EXTEND_Y0:
    {
      result = Mat(data.rows() + 1, data.cols());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < data.rows(); i ++)
        result.row(i) = data.row(i);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < data.cols(); i ++)
        result(data.rows(), i) = T(0);
      // N.B. nearest data in differential space.
      const auto d0data(compute(data, DETECT_Y));
      const auto ddata(compute(result, DETECT_Y));
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < result.cols(); i ++) {
        T a(0), b(0), c(0);
        for(int j = 0; j < d0data.rows(); j ++) {
          const T aa(- Dop(j, Dop.rows() - 1));
          const T bb(ddata(j, i) - d0data(j, i));
          a += aa * aa;
          b += aa * bb;
          c += bb * bb;
        }
        if(T(0) < b * b - a * c)
          result(data.rows(), i) = - b / a + sqrt(b * b - a * c) / a;
        else
          result(data.rows(), i) = - b / a;
      }
      for(int i = 0; i < result.cols(); i ++)
        result(data.rows(), i) = max(T(0), min(T(1), result(data.rows(), i)));
    }
    break;
  case EXTEND_Y:
    {
      const auto buf0(compute(data, EXTEND_Y0));
      const auto buf1(compute(compute(compute(data, REVERSE_Y), EXTEND_Y0), REVERSE_Y));
      result = Mat(data.rows() + 2, data.cols());
      result.row(0) = buf1.row(0);
      for(int i = 0; i < data.rows(); i ++)
        result.row(i + 1) = data.row(i);
      result.row(data.rows() + 1) = buf0.row(data.rows() - 1);
    }
    break;
  case EXTEND_YQ:
    result = recursive(data, EXTEND_YQ, EXTEND_Y);
    break;
  case EXTEND_YQS:
    result = recursiveSumup(data, EXTEND_YQ, EXTEND_Y);
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
  case CLIPPLUS:
    {
      result = Mat(data.rows(), data.cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          result(i, j) = max(T(0), data(i, j));
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
  case SQRTSCALE:
    {
      result = Mat(data.rows(), data.cols());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          result(i, j) = (data(i, j) < T(0) ? - T(1) : T(1)) * sqrt(abs(data(i, j)));
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
  xchg(C, bC);
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
  C = Mat(Dop.rows(), Dop.cols());
  for(int i = 0; i < C.rows(); i ++)
    for(int j = 0; j < C.cols(); j ++)
      C(i, j) = T(1) / (T(1 + abs(i - j)) / C.rows() * blur);
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
  T sumup(0);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int zi = 0; zi < T(1) / dratio; zi ++) {
    for(int j = 0; j < Dop0.size(); j ++) {
      Vec cpoint(2);
      cpoint[0] = (j - T(Dop0.size() - 1) / 2);
      cpoint[1] = T(zi + 1) * dratio;
      // x-z plane projection of point p with camera geometry c to z=0.
      // c := camera, p := cpoint.
      // <c + (p - c) * t, [0, 1]> = 0
      const auto t(- camera[1] / (cpoint[1] - camera[1]));
      const auto y0((camera + (cpoint - camera) * t)[0]);
      // N.B. average_k(dC_k / dy * z_k).
#if defined(_OPENMP)
#pragma omp critical
#endif
      {
        for(int i = 0; i < A.rows(); i ++) {
          // N.B. exp exp then log log to get better results.
          //      this is because in fact we need to do this with multiple
          //      sets of A, B and get maximum abs index as a result.
          A(i, getImgPt(i + y0, rows)) += Dop0[j] * exp(T(zi + 1) * dratio);
          B(i, getImgPt(i + y0, rows)) += Dop0[j];
        }
      }
    }
#if defined(_OPENMP)
#pragma omp atomic
#endif
    sumup += exp(T(zi + 1) * dratio);
  }
  const T n2(T(1) / dratio);
  A /= n2 * sumup;
  B /= n2;
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
      {
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

template <typename T> typename enlarger2ex<T>::Mat enlarger2ex<T>::recursive(const Mat& data, const direction_t& dir, const direction_t& dir0) {
  assert(96 <= sz_cell && 0 < st_cell && st_cell < sz_cell);
  Mat result(data.rows(), data.cols());
  if(sz_cell < data.rows()) {
    Mat former(data.rows() / 2, data.cols());
    Mat latter(data.rows() - data.rows() / 2, data.cols());
    Mat shrink(compute(data, DIV2_Y));
    for(int i = 0; i < former.rows(); i ++)
      former.row(i) = data.row(i);
    for(int i = 0; i < latter.rows(); i ++)
      latter.row(i) = data.row(former.rows() + i);
    // omit parallelize.
    former = compute(former, dir);
    latter = compute(latter, dir);
    shrink = compute(shrink, dir);
    for(int i = 0; i < former.rows(); i ++)
      result.row(i) = former.row(i);
    for(int i = 0; i < latter.rows(); i ++)
      result.row(former.rows() + i) = latter.row(i);
    for(int i = 0; i < data.rows(); i ++)
      result.row(i) += shrink.row(i / 2);
  } else
    result = compute(data, dir0);
  return result;
}

template <typename T> typename enlarger2ex<T>::Mat enlarger2ex<T>::recursiveSumup(const Mat& data, const direction_t& dir, const direction_t& dir0) {
  assert(96 <= sz_cell && 0 < st_cell && st_cell < sz_cell);
  Mat result(recursive(data, dir, dir0));
  if(sz_cell < result.rows())
    for(int ii = 1; ii < sz_cell; ii += st_cell) {
      Mat data2(data.rows() + ii, data.cols());
      for(int i = 0; i < ii; i ++)
        data2.row(i) = data.row(ii - i);
      for(int i = 0; i < data.rows(); i ++)
        data2.row(i + ii) = data.row(i);
      data2 = recursive(data, dir, dir0);
      for(int i = 0; i < result.rows(); i ++)
        result.row(i) += data2.row(i + ii);
    }
  return result;
}

#define _ENLARGE2X_
#endif

