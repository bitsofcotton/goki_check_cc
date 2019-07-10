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
using std::sin;
using std::cos;
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
 * Please re-initialize when parameters changed.
 */
template <typename T> class Filter {
public:
  typedef enum {
    ENLARGE_X,
    ENLARGE_Y,
    ENLARGE_BOTH,
    DETECT_X,
    DETECT_Y,
    DETECT_BOTH,
    COLLECT_X,
    COLLECT_Y,
    COLLECT_BOTH,
    BUMP_X,
    BUMP_Y,
    BUMP_BOTH,
    EXTEND_X,
    EXTEND_Y0,
    EXTEND_Y,
    EXTEND_BOTH,
    BCLIP,
    CLIP,
    ABS,
    EXPSCALE,
    LOGSCALE } direction_t;
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
  Filter();
  Mat  compute(const Mat& data, const direction_t& dir);
  MatU seed(const int& size, const bool& idft);
  Mat  gmean(const Mat& a, const Mat& b);
  void reinit();
  T    dratio;
  T    offset;
  
private:
  void initDop(const int& size);
  int  getImgPt(const T& y, const T& h);
  U    I;
  T    Pi;
  vector<Mat> Dop;
  vector<Mat> Eop;
  int  idx;
};

template <typename T> Filter<T>::Filter() {
  I  = sqrt(U(- 1.));
  Pi = atan2(T(1), T(1)) * T(4);
  // N.B. from accuracy reason, low depth.
  dratio  = T(.05);
  offset  = T(1) / T(64);
  idx     = - 1;
}

template <typename T> void Filter<T>::reinit() {
  Dop = Eop = vector<Mat>();
  idx = - 1;
}

template <typename T> typename Filter<T>::Mat Filter<T>::compute(const Mat& data, const direction_t& dir) {
  Mat result;
  switch(dir) {
  case ENLARGE_BOTH:
    // N.B. commutative.
    result = compute(compute(data, ENLARGE_X), ENLARGE_Y);
    break;
  case DETECT_BOTH:
    result = gmean(compute(data, DETECT_X), compute(data, DETECT_Y));
    break;
  case COLLECT_BOTH:
    result = gmean(compute(data, COLLECT_X), compute(data, COLLECT_Y));
    break;
  case BUMP_BOTH:
    result = gmean(compute(data, BUMP_X), compute(data, BUMP_Y));
    break;
  case EXTEND_BOTH:
    result = gmean(compute(compute(data, EXTEND_X), EXTEND_Y),
                   compute(compute(data, EXTEND_Y), EXTEND_X));
    break;
  case ENLARGE_X:
    result = compute(data.transpose(), ENLARGE_Y).transpose();
    break;
  case DETECT_X:
    result = compute(data.transpose(), DETECT_Y).transpose();
    break;
  case COLLECT_X:
    result = compute(data.transpose(), COLLECT_Y).transpose();
    break;
  case BUMP_X:
    result = compute(data.transpose(), BUMP_Y).transpose();
    break;
  case EXTEND_X:
    result = compute(data.transpose(), EXTEND_Y).transpose();
    break;
  case ENLARGE_Y:
    {
      initDop(data.rows());
      const Mat  diff(Dop[idx] * data);
            auto delta(compute(Eop[idx] * data, ABS));
      result = Mat(data.rows() * 2, data.cols());
      for(int i = 0; i < data.rows(); i ++) {
        result.row(i * 2 + 0) = data.row(i);
        result.row(i * 2 + 1) = data.row(i);
        for(int j = 0; j < data.cols(); j ++)
          if(diff(i, j) < T(0)) {
            result(i * 2 + 0, j) += delta(i, j);
            result(i * 2 + 1, j) -= delta(i, j);
          } else {
            result(i * 2 + 0, j) -= delta(i, j);
            result(i * 2 + 1, j) += delta(i, j);
          }
      }
    }
    break;
  case DETECT_Y:
    initDop(data.rows());
    result = Dop[idx] * data;
    break;
  case COLLECT_Y:
    result = compute(compute(data, DETECT_Y), ABS);
    break;
  case BUMP_Y:
    {
      result = Mat(data.rows(), data.cols());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          result(i, j) = T(0);
      initDop(max(3, min(int(data.rows()) / 16, int(T(1) / dratio / dratio))));
      Vec camera(2);
      camera[0] = T(0);
      camera[1] = T(1);
      assert(0 < dratio);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
      for(int zi = 0; zi < T(1) / dratio; zi ++) {
        Mat A(data.rows(), data.rows());
        for(int i = 0; i < A.rows(); i ++)
          for(int j = 0; j < A.cols(); j ++)
            A(i, j) = T(0);
        const Vec Dop0(Dop[idx].row(Dop[idx].rows() / 2) * exp(T(zi) * sqrt(dratio)));
        for(int j = 0; j < Dop0.size(); j ++) {
          Vec cpoint(2);
          cpoint[0] = j - T(Dop0.size() - 1) / 2;
          cpoint[1] = T(zi) * dratio;
          // x-z plane projection of point p with camera geometry c to z=0.
          // c := camera, p := cpoint.
          // <c + (p - c) * t, [0, 1]> = 0
          const auto t(- camera[1] / (cpoint[1] - camera[1]));
          const auto y0((camera + (cpoint - camera) * t)[0]);
          // N.B. average_k(dC_k / dy * z_k).
          for(int i = 0; i < A.rows(); i ++)
            A(i, getImgPt(i + y0, data.rows())) += Dop0[j];
        }
#if defined(_OPENMP)
#pragma omp atomic
#endif
        result += compute(A * data, ABS);
      }
      // N.B.
      // From hypothesis, it is correct to use:
      //   log(|average(d_k C * exp(z_k))| / |dC| == result.
      // But it's not here, because of some calculation experiment
      // causes some of false detection. So to avoid that, we use:
      //   log(|average(d_k C * exp(z_k * sqrt(dratio)))| / |dC|) == result.
      const auto dC(compute(compute(data, COLLECT_Y), BCLIP));
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          result(i, j) /= dC(i, j);
      result = compute(result, LOGSCALE) * dratio;
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
      initDop(data.rows());
      // N.B. nearest data in differential space.
      const Mat   d0data(Dop[idx] * data);
      initDop(result.rows());
      const Mat   ddata(Dop[idx] * result);
      const auto& Dop0(Dop[idx]);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < result.cols(); i ++) {
        T a(0), b(0);
        for(int j = 0; j < d0data.rows(); j ++) {
          // <d0data, (ddata + Diff * t * e_k)> == <d0data, d0data> in cut.
          a += d0data(j, i) * Dop0(j, Dop0.cols() - 1);
          b += d0data(j, i) * (d0data(j, i) - ddata(j, i));
        }
        if(a == 0) {
          cerr << "Error in EXTEND_Y0" << endl;
          result(data.rows(), i) = data(data.rows() - 1, i);
        } else
          result(data.rows(), i) = b / a;
      }
      result = compute(result, CLIP);
    }
    break;
  case EXTEND_Y:
    {
      result = Mat(data.rows() + 2, data.cols());
      Mat revdata(data.rows(), data.cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int i = 0; i < data.rows(); i ++)
        revdata.row(i) = data.row(data.rows() - i - 1);
      result.row(0) = compute(revdata, EXTEND_Y0).row(data.rows());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int i = 0; i < data.rows(); i ++)
        result.row(i + 1) = data.row(i);
      result.row(data.rows() + 1) = compute(data, EXTEND_Y0).row(data.rows());
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
  case CLIP:
    {
      result = Mat(data.rows(), data.cols());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          result(i, j) = max(T(0), min(T(1), data(i, j)));
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
  default:
    assert(0 && "unknown command in Filter (should not be reached.)");
  }
  return result;
}

template <typename T> void Filter<T>::initDop(const int& size) {
  for(int i = 0; i < Dop.size(); i ++)
    if(Dop[i].rows() == size) {
      idx = i;
      return;
    }
  cerr << "n" << flush;
  assert(2 <= size);
        auto DFTD(seed(size, false));
  DFTD.row(0) *= U(0);
  const auto IDFT(seed(size, true));
        auto DFTE(DFTD);
  T ni(0);
  T nd(0);
  for(int i = 1; i < DFTD.rows(); i ++) {
    const U phase(- U(2.) * Pi * sqrt(U(- 1)) * T(i) / T(DFTD.rows()));
    const U phase2( U(2.) * Pi * sqrt(U(- 1)) * T(i) / T(DFTD.rows()));
    // N.B. d/dy.
    DFTD.row(i) *= phase;
    nd += abs(phase)  * abs(phase);
    ni += abs(phase2) * abs(phase2);
    // N.B. integrate.
    // DFTI.row(i) /= phase;
    // N.B. (d^(log(h))/dy^(log(h)) f, lim h -> 1. : nop.
    // DFTH.row(i) *= log(phase);
    // N.B. please refer enlarge.wxm, half freq space refer and uses each.
    DFTE.row(i) /= exp(sqrt(U(- 1)) * Pi / T(2 * DFTE.rows())) - U(T(1));
  }
  // N.B. similar to det(Dop * Iop) == 1
  DFTD /= sqrt(nd * ni);
  idx   = Dop.size();
#if defined(_WITHOUT_EIGEN_)
  Dop.push_back((IDFT * DFTD).template real<T>());
  Eop.push_back((IDFT * DFTE).template real<T>());
#else
  Dop.push_back((IDFT * DFTD).real().template cast<T>());
  Eop.push_back((IDFT * DFTE).real().template cast<T>());
#endif
  return;
}

template <typename T> int Filter<T>::getImgPt(const T& y, const T& h) {
  int yy(int(y) % int(2 * h));
  if(yy < 0)
    yy = - yy;
  if(yy >= int(h))
    yy = int(h) - (yy - int(h));
  return yy % int(h);
}

template <typename T> typename Filter<T>::MatU Filter<T>::seed(const int& size, const bool& idft) {
  MatU result(size, size);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < result.rows(); i ++)
    for(int j = 0; j < result.cols(); j ++)
      result(i, j) = exp(I * U(- 2. * Pi * T(idft ? - 1 : 1) * i * j / T(size))) / T(idft ? - size : 1);
  return result;
}

template <typename T> typename Filter<T>::Mat Filter<T>::gmean(const Mat& a, const Mat& b) {
  assert(a.rows() == b.rows() && a.cols() == b.cols());
  Mat res(a.rows(), a.cols());
  for(int i = 0; i < a.rows(); i ++)
    for(int j = 0; j < a.cols(); j ++) {
      const auto lval(a(i, j) * b(i, j));
      res(i, j) = (lval < T(0) ? - T(1) : T(1)) * sqrt(abs(lval));
    }
  return res;
}

#define _ENLARGE2X_
#endif

