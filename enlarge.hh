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
    IDETECT_X_BOTH,
    IDETECT_Y,
    IDETECT_Y_BOTH,
    IDETECT_BOTH,
    IDETECT_QUAD,
    BUMP_X,
    BUMP_Y0,
    BUMP_Y1,
    BUMP_Y2,
    BUMP_Y3,
    BUMP_Y,
    BUMP_YP,
    BUMP_BOTH,
    EXTEND_X,
    EXTEND_Y0,
    EXTEND_Y1,
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
  enlarger2ex(const enlarger2ex<T>& src);
  enlarger2ex<T>& operator = (const enlarger2ex<T>& src);
  Mat  compute(const Mat& data, const direction_t& dir);
  MatU seed(const int& size, const bool& idft);
  T    dratio;
  T    offset;
  bool di_mode[5];
  bool di_bump_mode[5];
  
private:
  void initDop(const int& size);
  void initBump(const int& size);
  Vec  minSquare(const Vec& in);
  int  getImgPt(const T& y, const T& h);
  void makeDI(const int& size, Vec& Dop, Vec& Dhop, Vec& Iop, Vec& Eop);
  void xchg(Mat& a, Mat& b);
  Mat  round2y(const Mat& in, const int& h);
  U    I;
  T    Pi;
  Mat  A[4];
  Mat  B[4];
  Mat  Dop[4];
  Mat  Dhop[4];
  Mat  Eop[4];
  Mat  Iop[4];
};

template <typename T> enlarger2ex<T>::enlarger2ex() {
  I  = sqrt(U(- 1.));
  Pi = atan2(T(1), T(1)) * T(4);
  dratio  = T(1e-3);
  offset  = T(4) / T(256);
  for(int i = 0; i < 5; i ++)
    di_bump_mode[i] = di_mode[i] = false;
}

template <typename T> enlarger2ex<T>::enlarger2ex(const enlarger2ex<T>& src) {
  *this = src;
}

template <typename T> enlarger2ex<T>& enlarger2ex<T>::operator = (const enlarger2ex<T>& src) {
  I  = src.I;
  Pi = src.Pi;
  dratio  = src.dratio;
  offset  = src.offset;
  for(int i = 0; i < 5; i ++) {
    di_mode[i]      = src.di_mode[i];
    di_bump_mode[i] = src.di_bump_mode[i];
  }
  for(int i = 0; i < 4; i ++) {
    A[i]     = src.A[i];
    B[i]     = src.B[i];
    Dop[i]   = src.Dop[i];
    Dhop[i]  = src.Dhop[i];
    Eop[i]   = src.Eop[i];
    Iop[i]   = src.Iop[i];
  }
  return *this;
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
  case IDETECT_QUAD:
    result = (compute(data, IDETECT_Y_BOTH) + compute(data, IDETECT_X_BOTH)) / T(2);
    break;
  case BUMP_BOTH:
    result = (compute(data, BUMP_Y)    + compute(data, BUMP_X)) / T(2);
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
  case IDETECT_X_BOTH:
    result = compute(data.transpose(), IDETECT_Y_BOTH).transpose();
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
    di_mode[0] = false;
    initDop(data.rows());
    result  = compute(Eop[0] * data, CLIP);
    di_mode[0] = true;
    initDop(data.rows());
    result += compute(Eop[0] * data, CLIP);
    break;
    break;
  case DETECT_Y:
    di_mode[0] = false;
    initDop(data.rows());
    result  = Dop[0] * data;
    di_mode[0] = true;
    initDop(data.rows());
    result += Dop[0] * data;
    break;
  case DETECT_NOP_Y:
    di_mode[0] = false;
    initDop(data.rows());
    result  = Dhop[0] * data;
    di_mode[0] = true;
    initDop(data.rows());
    result += Dhop[0] * data;
    break;
  case COLLECT_Y:
    result = compute(compute(data, DETECT_Y), ABS);
    break;
  case IDETECT_Y:
    {
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
      di_mode[0] = false;
      initDop(data.rows());
      result  = Iop[0] * result;
      di_mode[0] = true;
      initDop(data.rows());
      result += Iop[0] * result;
      result /= T(2);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < data.cols(); i ++)
        for(int j = 0; j < data.rows(); j ++)
          result(j, i) += ms[i][0] * j / data.rows() + ms[i][1] * j * j / 2 / data.rows() / data.rows();
    }
    break;
  case IDETECT_Y_BOTH:
    result = compute(data, IDETECT_Y) + compute(compute(compute(data, REVERSE_Y), IDETECT_Y), REVERSE_Y);
    break;
  case BUMP_Y0:
    {
      initBump(data.rows());
      // |average(dC*z_k)/average(dC)| == dataA / dataB.
      const auto dataA(compute(A[0] * data, ABS));
      const auto dataB(compute(B[0] * data, ABS));
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
      result = compute(compute(d2xx, IDETECT_X), IDETECT_X) +
               compute(compute(d2xy, IDETECT_Y), IDETECT_X) +
               compute(compute(d2yx, IDETECT_X), IDETECT_Y) +
               compute(compute(d2yy, IDETECT_Y), IDETECT_Y);
    }
    break;
  case BUMP_Y:
    di_bump_mode[0] = false;
    result  = compute(data, BUMP_Y0) + compute(compute(compute(data, REVERSE_X), BUMP_Y0), REVERSE_X) + compute(compute(compute(data, REVERSE_Y), BUMP_Y0), REVERSE_Y) + compute(compute(compute(data, REVERSE_BOTH), BUMP_Y0), REVERSE_BOTH);
    di_bump_mode[0] = true;
    result += compute(data, BUMP_Y0) + compute(compute(compute(data, REVERSE_X), BUMP_Y0), REVERSE_X) + compute(compute(compute(data, REVERSE_Y), BUMP_Y0), REVERSE_Y) + compute(compute(compute(data, REVERSE_BOTH), BUMP_Y0), REVERSE_BOTH);
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
          const T aa(- Dop[0](j, Dop[0].rows() - 1));
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
  for(int i = 0; i < 4; i ++)
    if(Dop[i].rows() == size && di_mode[0] == di_mode[i + 1]) {
      if(i) {
        xchg(Dop[0],  Dop[i]);
        xchg(Dhop[0], Dhop[i]);
        xchg(Iop[0],  Iop[i]);
        xchg(Eop[0],  Eop[i]);
      }
      di_mode[i + 1] = di_mode[1];
      di_mode[1]     = di_mode[0];
      return;
    }
  cerr << "new" << flush;
  Vec vDop,  vDop0;
  Vec vDhop;
  Vec vIop;
  Vec vEop;
  di_bump_mode[1] = di_bump_mode[0];
  makeDI(size, vDop, vDhop, vIop, vEop);
  vEop *= T(2);
  Dop[0]  = Mat(size, size);
  Dhop[0] = Mat(size, size);
  Iop[0]  = Mat(size, size);
  Eop[0]  = Mat(size, size);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < Dop[0].rows(); i ++)
    for(int j = 0; j < Dop[0].cols(); j ++) {
      Dop[0](i, j)  =   vDop[ getImgPt(j - i - size / 2, Dop[0].cols())];
      Dhop[0](i, j) =   vDhop[getImgPt(j - i - size / 2, Dop[0].cols())];
      Iop[0](i, j)  = - vIop[ getImgPt(j - i - size / 2, Dop[0].cols())];
      Eop[0](i, j)  =   vEop[ getImgPt(j - i - size / 2, Dop[0].cols())];
    }
  Mat newEop(Eop[0].rows() * 2, Eop[0].cols());
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < Eop[0].rows(); i ++) {
    newEop.row(2 * i + 0) =   Eop[0].row(i);
    newEop.row(2 * i + 1) = - Eop[0].row(i);
    newEop(2 * i + 0, i) += T(1);
    newEop(2 * i + 1, i) += T(1);
  }
  Eop[0] = newEop * T(2);
  for(int i = 0; i < Eop[0].rows(); i ++) {
    Eop[0].row(i) += newEop.row(min(i + 1, int(Eop[0].rows()) - 1));
    Eop[0].row(i) += newEop.row(max(i - 1, 0));
  }
  Eop[0] /= T(4);
  return;
}

template <typename T> void enlarger2ex<T>::initBump(const int& size) {
  cerr << "." << flush;
  for(int i = 0; i < 4; i ++)
    if(A[i].rows() == size && di_bump_mode[0] == di_bump_mode[i + 1]) {
      if(i) {
        xchg(A[0], A[i]);
        xchg(B[0], B[i]);
      }
      di_bump_mode[i + 1] = di_bump_mode[1];
      di_bump_mode[1]     = di_bump_mode[0];
      return;
    }
  cerr << "new" << flush;
  assert(0 < size);
  A[0] = Mat(size, size);
  B[0] = Mat(size, size);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < size; i ++)
    for(int j = 0; j < size; j ++)
      B[0](i, j) = A[0](i, j) = T(0);
  Vec Dop0;
  Vec Dhop0;
  Vec Iop0;
  Vec Eop0;
  di_bump_mode[1] = di_bump_mode[0];
  makeDI(min(max(7, int(sqrt(T(size)))), int(size)), Dop0, Dhop0, Iop0, Eop0);
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
        for(int i = 0; i < A[0].rows(); i ++) {
          A[0](i, getImgPt(i + y0, size)) += Dop0[j] * T(zi + 1) * dratio;
          B[0](i, getImgPt(i + y0, size)) += Dop0[j];
        }
      }
    }
#if defined(_OPENMP)
#pragma omp atomic
#endif
    sumup += T(zi + 1) * dratio;
  }
  const T n2(T(1) / dratio);
  A[0] /= n2 * sumup;
  B[0] /= n2;
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
  T sumup(0);
  assert(Dop.size() == size && Dhop.size() == size && Iop.size() == size && Eop.size() == size);
#if defined(_RECURSIVE_RECURSIVE_)
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int ss = 3; ss <= size / 2; ss ++) {
#else
  for(int ss = size / 2; ss <= size / 2; ss ++) {
#endif
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
      VecU lDop( DFTD * IDFT.row(iidx));
      VecU lDhop(DFTH * IDFT.row(iidx));
      VecU lIop( DFTI * IDFT.row(iidx));
      VecU lEop( DFTE * IDFT.row(iidx));
#else
      VecU lDop( IDFT.row(iidx) * DFTD);
      VecU lDhop(IDFT.row(iidx) * DFTH);
      VecU lIop( IDFT.row(iidx) * DFTI);
      VecU lEop( IDFT.row(iidx) * DFTE);
#endif
      // N.B. averate effect for each diffs.
      const T ratio(di_bump_mode[0] ? T(ss) : T(1) / ss);
      sumup += ratio;
      lDop  *= ratio;
      lDhop *= ratio;
      lIop  *= ratio;
      lEop  *= ratio;
#if defined(_OPENMP)
#pragma omp critical
#endif
      {
        for(int j = i; j - i < lDop.size(); j ++) {
          const int idx(j + size / 2 - IDFT.rows() + 1);
          const int jdx(j - i);
          Dop[idx]  += T(lDop[ jdx].real());
          Dhop[idx] += T(lDhop[jdx].real());
          Iop[idx]  += T(lIop[ jdx].real());
          Eop[idx]  += T(lEop[ jdx].real());
        }
      }
    }
  }
  Dop  /= sumup;
  Dhop /= sumup;
  Iop  /= sumup;
  Eop  /= sumup;
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

