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
    ENLARGE_Y,
    ENLARGE_FX,
    ENLARGE_FY,
    ENLARGE_BOTH,
    ENLARGE_BOTHQ,
    ENLARGE_BOTHQS,
    ENLARGE_FBOTH,
    ENLARGE_3BOTH,
    DETECT_X,
    DETECT_Y,
    DETECT_NOP_Y,
    DETECT_BOTH,
    DETECT_BOTHQ,
    DETECT_BOTHQS,
    COLLECT_X,
    COLLECT_Y,
    COLLECT_BOTH,
    COLLECT_BOTHQ,
    COLLECT_BOTHQS,
    IDETECT_X,
    IDETECT_Y,
    IDETECT_BOTH,
    IDETECT_BOTHQ,
    IDETECT_BOTHQS,
    BLUR_X,
    BLUR_Y,
    BLUR_BOTH,
    BUMP_X,
    BUMP_Y,
    BUMP_BOTH,
    BUMP_BOTHQ,
    BUMP_BOTHQS,
    EXTEND_X,
    EXTEND_XQ,
    EXTEND_XQS,
    EXTEND_Y0,
    EXTEND_Y1,
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
  void initBump(const int& size);
  Vec  minSquare(const Vec& in);
  int  getImgPt(const T& y, const T& h);
  void makeDI(const int& size, Mat& Dop, Mat& Dhop, Mat& Iop, Mat& Eop);
  Mat  recursive(const Mat& data, const direction_t& dir, const direction_t& dir0);
  Mat  recursiveSumup(const Mat& data, const direction_t& dir, const direction_t& dir0);
  U    I;
  T    Pi;
  T    recursive_ratio;
  vector<Mat> A;
  vector<Mat> B;
  vector<Mat> C;
  vector<Mat> Dop;
  vector<Mat> Dhop;
  vector<Mat> Eop;
  vector<Mat> Iop;
  int idx_d;
  int idx_b;
};

template <typename T> enlarger2ex<T>::enlarger2ex() {
  I  = sqrt(U(- 1.));
  Pi = atan2(T(1), T(1)) * T(4);
  dratio  = T(.0005);
  offset  = T(4) / T(256);
  blur    = T(8);
  sz_cell = 4;
  st_cell = 8;
  idx_d   = - 1;
  idx_b   = - 1;
  recursive_ratio = T(1.);
}

template <typename T> typename enlarger2ex<T>::Mat enlarger2ex<T>::compute(const Mat& data, const direction_t& dir) {
  Mat result;
  switch(dir) {
  case ENLARGE_BOTH:
    result = (compute(compute(data, ENLARGE_X), ENLARGE_Y) +
              compute(compute(data, ENLARGE_Y), ENLARGE_X)) / T(2);
    break;
  case ENLARGE_BOTHQ:
    result = recursive(data, ENLARGE_BOTHQ, ENLARGE_BOTH);
    break;
  case ENLARGE_BOTHQS:
    result = recursiveSumup(data, ENLARGE_BOTHQ, ENLARGE_BOTH);
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
    result = recursive(data, DETECT_BOTHQ, DETECT_BOTH);
    break;
  case DETECT_BOTHQS:
    result = recursiveSumup(data, DETECT_BOTHQ, DETECT_BOTH);
    break;
  case COLLECT_BOTH:
    result = (compute(data, COLLECT_X) + compute(data, COLLECT_Y)) / T(2);
    break;
  case COLLECT_BOTHQ:
    result = compute(compute(data, COLLECT_BOTHQ), ABS);
    break;
  case COLLECT_BOTHQS:
    result = compute(compute(data, COLLECT_BOTHQS), ABS);
    break;
  case IDETECT_BOTH:
    result = (compute(data, IDETECT_X) + compute(data, IDETECT_Y)) / T(2);
    break;
  case IDETECT_BOTHQ:
    result = recursive(data, IDETECT_BOTHQ, IDETECT_BOTH);
    break;
  case IDETECT_BOTHQS:
    result = recursiveSumup(data, IDETECT_BOTHQ, IDETECT_BOTH);
    break;
  case BLUR_BOTH:
    result = (compute(data, BLUR_X)    + compute(data, BLUR_Y)) / T(2);
    break;
  case BUMP_BOTH:
    result = (compute(data, BUMP_Y)    + compute(data, BUMP_X)) / T(2);
    break;
  case BUMP_BOTHQ:
    result = recursive(data, BUMP_BOTHQ, BUMP_BOTH);
    break;
  case BUMP_BOTHQS:
    result = recursiveSumup(data, BUMP_BOTHQ, BUMP_BOTH);
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
  case BLUR_X:
    result = compute(data.transpose(), BLUR_Y).transpose();
    break;
  case BUMP_X:
    result = compute(data.transpose(), BUMP_Y).transpose();
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
    result = compute(Eop[idx_d] * data, CLIP);
    break;
  case DETECT_Y:
    initDop(data.rows());
    result = Dop[idx_d] * data;
    break;
  case DETECT_NOP_Y:
    initDop(data.rows());
    result = Dhop[idx_d] * data;
    break;
  case BLUR_Y:
    initDop(data.rows());
    result = C[idx_d] * data;
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
      result = Iop[idx_d] * result;
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < data.cols(); i ++)
        for(int j = 0; j < data.rows(); j ++)
          result(j, i) += ms[i][0] * j / data.rows() + ms[i][1] * j * j / 2 / data.rows() / data.rows();
    }
    break;
  case BUMP_Y:
    {
      // |average(dC*z_k)/average(dC)| == dataA / dataB.
      initBump(data.rows());
      const auto dataA(compute(A[idx_b] * data, ABS));
      const auto dataB(compute(B[idx_b] * data, ABS));
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
      result = Mat(data.rows(), data.cols());
      Mat ddx(data.rows(), data.cols());
      Mat ddy(data.rows(), data.cols());
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
          result(i, j) = dataA(i, j) / dataBc(i, j);
          ddx(i, j)  = (datadxA(i, j) * dataB(i, j) - dataA(i, j) * datadxB(i, j)) / dataBc(i, j) / dataBc(i, j);
          ddy(i, j)  = (datadyA(i, j) * dataB(i, j) - dataA(i, j) * datadyB(i, j)) / dataBc(i, j) / dataBc(i, j);
          d2xx(i, j) = (datad2xxA(i, j) * dataB(i, j) + datadxA(i, j) * datadxB(i, j) - datadxA(i, j) * datadxB(i, j) - dataA(i, j) * datad2xxB(i, j)) / pow(dataBc(i, j), T(2)) + (datadxA(i, j) * dataB(i, j) + dataA(i, j) * datadxB(i, j)) / pow(dataBc(i, j), T(3)) * (- T(2)) * datadxB(i, j);
          d2xy(i, j) = (datad2xyA(i, j) * dataB(i, j) + datadyA(i, j) * datadxB(i, j) - datadxA(i, j) * datadyB(i, j) - dataA(i, j) * datad2xyB(i, j)) / pow(dataBc(i, j), T(2)) + (datadyA(i, j) * dataB(i, j) + dataA(i, j) * datadyB(i, j)) / pow(dataBc(i, j), T(3)) * (- T(2)) * datadxB(i, j);
          d2yx(i, j) = (datad2yxA(i, j) * dataB(i, j) + datadxA(i, j) * datadyB(i, j) - datadyA(i, j) * datadxB(i, j) - dataA(i, j) * datad2yxB(i, j)) / pow(dataBc(i, j), T(2)) + (datadxA(i, j) * dataB(i, j) + dataA(i, j) * datadxB(i, j)) / pow(dataBc(i, j), T(3)) * (- T(2)) * datadyB(i, j);
          d2yy(i, j) = (datad2yyA(i, j) * dataB(i, j) + datadyA(i, j) * datadyB(i, j) - datadyA(i, j) * datadyB(i, j) - dataA(i, j) * datad2yyB(i, j)) / pow(dataBc(i, j), T(2)) + (datadyA(i, j) * dataB(i, j) + dataA(i, j) * datadyB(i, j)) / pow(dataBc(i, j), T(3)) * (- T(2)) * datadyB(i, j);
        }
      // N.B. f ~ f0 + f0' * dy + f0'' * dy^2 / 2! + ... + f0^(n) * dy^n / n! + ...
      // N.B. (d/dt)^2 f = d/dt(df/dx dx/dt + df/dy dy/dt)
      // = d^2f/dx^2 dx/dt + df/dx d^2x/dt^2 +
      //   d^2f/dydx dx/dt + df/dx dxdy/dt^2 +
      //   d^2f/dxdy dy/dt + df/dy dydx/dt^2 +
      //   d^2f/dy^2 dy/dt + df/dy d^2y/dt^2
      // dx/dt == dy/dt == const. ->
      // d^2/dt^2(result(i, j)) = d2dxdx + d2dxdy + d2dydx + d2dydy;
      // N.B. in fact, in this case, so function don't have rotation related
      //      information, rotation + d/dt sumup is needed, but now, not so.
      result += (compute(ddx, IDETECT_X) + compute(ddy, IDETECT_Y)) / T(2);
      result += (compute(compute(d2xx, IDETECT_X) +
                         compute(d2yx, IDETECT_Y), IDETECT_X) +
                 compute(compute(d2xy, IDETECT_X) +
                         compute(d2yy, IDETECT_Y), IDETECT_Y))
                / T(4) / T(2);
      // N.B. logscale is artificial,
      //      so this should not be needed but works well.
      result  = compute(result, LOGSCALE);
    }
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
          const T aa(- Dop[idx_d](j, Dop[idx_d].rows() - 1));
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
  case EXTEND_Y1:
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
  case EXTEND_Y:
    result = compute(data, EXTEND_Y1);
    break;
  case EXTEND_YQ:
    result = recursive(data, EXTEND_YQ, EXTEND_Y1);
    break;
  case EXTEND_YQS:
    result = recursiveSumup(data, EXTEND_YQ, EXTEND_Y1);
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
          result(i, j) = (data(i, j) < T(0) ? - T(1) : T(1)) * (exp(T(1) + abs(data(i, j))) - exp(T(1)));
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
          result(i, j) = (data(i, j) < T(0) ? - T(1) : T(1)) * (log(exp(T(1)) + abs(data(i, j))) - T(1));
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
  for(int i = 0; i < Dop.size(); i ++)
    if(Dop[i].rows() == size) {
      idx_d = i;
      return;
    }
  cerr << "n" << flush;
  idx_d = Dop.size();
  Dop.push_back(Mat());
  Dhop.push_back(Mat());
  Iop.push_back(Mat());
  Eop.push_back(Mat());
  C.push_back(Mat());
  Mat vDop;
  Mat vDhop;
  Mat vIop;
  Mat vEop;
  assert(2 <= size);
  makeDI(size, vDop, vDhop, vIop, vEop);
  Dop[idx_d]  = Mat(size, size);
  Dhop[idx_d] = Mat(size, size);
  Iop[idx_d]  = Mat(size, size);
  Eop[idx_d]  = Mat(size, size);
  C[idx_d]    = Mat(size, size);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < Dop[idx_d].rows(); i ++)
    for(int j = 0; j < Dop[idx_d].cols(); j ++)
      Dop[idx_d](i, j) = Dhop[idx_d](i, j) = Iop[idx_d](i, j) = Eop[idx_d](i, j) = T(0);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < Dop[idx_d].rows(); i ++)
    for(int j = 0; j < Dop[idx_d].cols() / 2; j ++) {
      Dop[idx_d]( i, i / 2 + j) =   vDop( i / 2, j);
      Dhop[idx_d](i, i / 2 + j) =   vDhop(i / 2, j);
      Iop[idx_d]( i, i / 2 + j) = - vIop( i / 2, j);
      Eop[idx_d]( i, i / 2 + j) =   vEop( i / 2, j);
    }
  Eop[idx_d] *= T(2);
  Mat newEop(Eop[idx_d].rows() * 2, Eop[idx_d].cols());
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < Eop[idx_d].rows(); i ++) {
    newEop.row(2 * i + 0) =   Eop[idx_d].row(i);
    newEop.row(2 * i + 1) = - Eop[idx_d].row(i);
    newEop(2 * i + 0, i) += T(1);
    newEop(2 * i + 1, i) += T(1);
  }
  Eop[idx_d] = newEop * T(2);
  for(int i = 0; i < Eop[idx_d].rows(); i ++) {
    Eop[idx_d].row(i) += newEop.row(min(i + 1, int(Eop[idx_d].rows()) - 1));
    Eop[idx_d].row(i) += newEop.row(max(i - 1, 0));
  }
  Eop[idx_d] /= T(4);
  for(int i = 0; i < C[idx_d].rows(); i ++)
    for(int j = 0; j < C[idx_d].cols(); j ++)
      C[idx_d](i, j) = T(1) / (T(1 + abs(i - j)) / C[idx_d].rows() * blur);
  return;
}

template <typename T> void enlarger2ex<T>::initBump(const int& size) {
  for(int i = 0; i < A.size(); i ++)
    if(A[i].rows() == size) {
      idx_b = i;
      return;
    }
  cerr << "n" << flush;
  idx_b = A.size();
  A.push_back(Mat());
  B.push_back(Mat());
  assert(2 <= size);
  A[idx_b] = Mat(size, size);
  B[idx_b] = Mat(size, size);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < size; i ++)
    for(int j = 0; j < size; j ++)
      B[idx_b](i, j) = A[idx_b](i, j) = T(0);
  Mat Dop0;
  Mat Dhop0;
  Mat Iop0;
  Mat Eop0;
  makeDI(max(3, int(T(2) * sqrt(T(size)))), Dop0, Dhop0, Iop0, Eop0);
  Vec camera(2);
  camera[0] = T(0);
  camera[1] = T(1);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int zi = 0; zi < T(1) / dratio; zi ++) {
    for(int j = 0; j < Dop0.rows() / 2; j ++) {
      Vec cpoint(2);
      cpoint[0] = (j - T(Dop0.rows() / 2 - 1) / 2);
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
        if(Dop0.rows() % 2 == 1 || Dop0.rows() <= 3)
          for(int i = 0; i < A[idx_b].rows(); i ++) {
            A[idx_b](i, getImgPt(i + y0, size)) += Dop0(Dop0.rows() / 2, j) * T(zi + 1) * dratio;
            B[idx_b](i, getImgPt(i + y0, size)) += Dop0(Dop0.rows() / 2, j);
          }
        else
          for(int i = 0; i < A[idx_b].rows(); i ++) {
            A[idx_b](i, getImgPt(i + y0, size)) += (Dop0(Dop0.rows() / 2, j) + Dop0(Dop0.rows() / 2 + 1, j)) / T(2) * T(zi + 1) * dratio;
            B[idx_b](i, getImgPt(i + y0, size)) += (Dop0(Dop0.rows() / 2, j) + Dop0(Dop0.rows() / 2 + 1, j)) / T(2);
          }
      }
    }
  }
  A[idx_b] *= dratio;
  B[idx_b] *= dratio;
  return;
}

template <typename T> int enlarger2ex<T>::getImgPt(const T& y, const T& h) {
  return int(abs(int(y - pow(h, int(log(y) / log(h))) + .5 + 2 * h * h + h) % int(2 * h) - h)) % int(h);
}

template <typename T> void enlarger2ex<T>::makeDI(const int& size, Mat& Dop, Mat& Dhop, Mat& Iop, Mat& Eop) {
  assert(2 <= size);
  const auto ss((size + 1) / 2);
        auto DFTD(seed(ss, false));
  DFTD.row(0) *= U(0);
  const auto IDFT(seed(ss, true));
        auto DFTH(DFTD);
        auto DFTI(DFTD);
        auto DFTE(DFTD);
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
  }
  DFTD /= T(2) * Pi;
  DFTH /= T(2) * Pi;
  DFTI /= T(2) * Pi;
  DFTE /= T(2) * Pi;
#if defined(_WITHOUT_EIGEN_)
  Dop  = (IDFT * DFTD).template real<T>();
  Dhop = (IDFT * DFTH).template real<T>();
  Iop  = (IDFT * DFTI).template real<T>();
  Eop  = (IDFT * DFTE).template real<T>();
#else
  Dop  = (IDFT * DFTD).real().template cast<T>();
  Dhop = (IDFT * DFTH).real().template cast<T>();
  Iop  = (IDFT * DFTI).real().template cast<T>();
  Eop  = (IDFT * DFTE).real().template cast<T>();
#endif
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

// over all O(const. * lg(y) * (compute(f(y), dir0) + 2. * compute(f(y/2), dir0) + ...)
template <typename T> typename enlarger2ex<T>::Mat enlarger2ex<T>::recursive(const Mat& data, const direction_t& dir, const direction_t& dir0) {
  assert(4 <= sz_cell);
  assert(T(0) <= recursive_ratio);
  Mat result;
  // N.B. might have better result with data.rows() - 1
  //      other than data.rows() / 2, but it's extremely slow.
  if(sz_cell < data.rows() && sz_cell < data.cols()) {
    Mat a00((data.rows() + 1) / 2, (data.cols() + 1) / 2);
    Mat a01(a00);
    Mat a10(a00);
    Mat a11(a00);
    for(int i = 0; i < a00.rows(); i ++)
      for(int j = 0; j < a00.cols(); j ++) {
        a00(i, j) = data(i, j);
        a01(i, j) = data(i, j - a00.cols() + data.cols());
        a10(i, j) = data(i - a00.rows() + data.rows(), j);
        a11(i, j) = data(i - a00.rows() + data.rows(), j - a00.cols() + data.cols());
      }
    a00 = compute(a00, dir);
    a01 = compute(a01, dir);
    a10 = compute(a10, dir);
    a11 = compute(a11, dir);
    result = compute(data, dir0) * recursive_ratio;
    switch(dir0) {
    case EXTEND_Y1:
      for(int i = 0; i < result.cols(); i ++) {
        result(0, i) += i < a00.cols() ? a00(0, i) : a01(0, i - result.cols() + a01.cols());
        result(data.rows() + 1, i) += i < a00.cols() ? a10(a10.rows() - 1, i) : a11(a11.rows() - 1, i - result.cols() + a11.cols());
      }
      result.row(0) /= T(1) + recursive_ratio;
      result.row(result.rows() - 1) /= T(1) + recursive_ratio;
      break;
    default:
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++) {
          result(i, j) += (i < a00.rows() ?
                  (j < a00.cols() ? a00(i, j) :
                                    a01(i, j - result.cols() + a01.cols())) :
                  (j < a10.cols() ? a10(i - result.rows() + a10.rows(), j) :
                                    a11(i - result.rows() + a11.rows(),
                                        j - result.cols() + a11.cols())));
        }
      result /= T(1) + recursive_ratio;
    }
  } else if(sz_cell < data.rows()) { 
    Mat former((data.rows() + 1) / 2, data.cols());
    Mat latter((data.rows() + 1) / 2, data.cols());
    for(int i = 0; i < former.rows(); i ++) {
      former.row(i) = data.row(i);
      latter.row(i) = data.row(i - former.rows() + data.rows());
    }
    former = compute(former, dir);
    latter = compute(latter, dir);
    result = compute(data, dir0) * recursive_ratio;
    switch(dir0) {
    case EXTEND_Y1:
      result.row(0) += former.row(0);
      result.row(result.rows() - 1) += latter.row(latter.rows() - 1);
      result.row(0) /= T(1) + recursive_ratio;
      result.row(result.rows() - 1) /= T(1) + recursive_ratio;
      break;
    default:
      for(int i = 0; i < former.rows(); i ++)
        result.row(i) += former.row(i);
      for(int i = former.rows(); i < result.rows(); i ++)
        result.row(i) += latter.row(i - result.rows() + latter.rows());
      result /= T(1) + recursive_ratio;
    }
  } else if(sz_cell < data.cols()) {
    Mat former(data.rows(), data.cols() / 2);
    Mat latter(data.rows(), data.cols() - former.cols());
#if defined(_WITHOUT_EIGEN_)
    for(int i = 0; i < former.cols(); i ++)
      former.setCol(i, data.col(i));
    for(int i = former.cols(); i < data.cols(); i ++)
      latter.setCol(i - data.cols() + latter.cols(), data.col(i));
#else
    for(int i = 0; i < former.cols(); i ++)
      former.col(i) = data.col(i);
    for(int i = former.cols(); i < data.cols(); i ++)
      latter.col(i - data.cols() + latter.cols()) = data.col(i);
#endif
    former = compute(former, dir);
    latter = compute(latter, dir);
    result = compute(data, dir0) * recursive_ratio;
    switch(dir0) {
    case EXTEND_Y1:
      for(int i = 0; i < result.cols(); i ++) {
        result(0, i) += i < former.cols() ? former(0, i) : latter(0, i - result.cols() + latter.cols());
        result(result.rows() - 1, i) += i < former.cols() ? former(former.rows() - 1, i) : latter(latter.rows() - 1, i - result.cols() + latter.cols());
      }
      result.row(0) /= T(1) + recursive_ratio;
      result.row(result.rows() - 1) /= T(1) + recursive_ratio;
      break;
    default:
#if defined(_WITHOUT_EIGEN_)
      for(int i = 0; i < former.cols(); i ++)
        result.setCol(i, result.col(i) + former.col(i));
      for(int i = former.cols(); i < result.cols(); i ++)
        result.setCol(i, result.col(i) + latter.col(i - result.cols() + latter.cols()));
#else
      for(int i = 0; i < former.cols(); i ++)
        result.col(i) += former.col(i);
      for(int i = former.cols(); i < result.cols(); i ++)
        result.col(i) += latter.col(i - result.cols() + latter.cols());
#endif
      result /= T(1) + recursive_ratio;
    }
  } else
    result = compute(data, dir0);
  return result;
}

template <typename T> typename enlarger2ex<T>::Mat enlarger2ex<T>::recursiveSumup(const Mat& data, const direction_t& dir, const direction_t& dir0) {
  assert(0 < st_cell && 0 < sz_cell);
  Mat result(recursive(data, dir, dir0));
  if(sz_cell < result.rows())
    for(int ii = 1; ii < int(data.rows() / 2); ii += max(int(data.rows() / 2) / st_cell, 1)) {
      Mat data2(data.rows(), data.cols());
      for(int i = ii; i < data.rows(); i++)
        data2.row(i) = data.row(i - ii);
      for(int i = 0; i < ii; i ++)
        data2.row(i - ii + data.rows()) = data.row(i);
      data2 = recursive(data2, dir, dir0);
      for(int i = 0; i < result.rows(); i ++)
        result.row(i) += data2.row((i + ii) % data.rows());
    }
  else
    return result;
  return result / st_cell;
}

#define _ENLARGE2X_
#endif

