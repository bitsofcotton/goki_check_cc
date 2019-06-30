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
    DEDGE,
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
  void reinit();
  T    dratio;
  T    offset;
  T    thedge;
  int  sq;
  Mat  gmean(const Mat& a, const Mat& b);
  
private:
  void initDop(const int& size);
  void initBump(const int& size);
  int  getImgPt(const T& y, const T& h);
  void makeDI(const int& size, Mat& Dop, Mat& Eop, const bool& recursive = false);
  U    I;
  T    Pi;
  vector<Mat> A;
  vector<Mat> Dop;
  vector<Mat> Eop;
  int idx_d;
  int idx_b;
};

template <typename T> Filter<T>::Filter() {
  I  = sqrt(U(- 1.));
  Pi = atan2(T(1), T(1)) * T(4);
  dratio  = T(.005);
  offset  = T(1) / T(64);
  thedge  = T(.05);
  sq      = 4;
  idx_d   = - 1;
  idx_b   = - 1;
}

template <typename T> void Filter<T>::reinit() {
  A     = Dop   = Eop = vector<Mat>();
  idx_d = idx_b = - 1;
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
    result = (compute(compute(data, EXTEND_X), EXTEND_Y) +
              compute(compute(data, EXTEND_Y), EXTEND_X)) / T(2);
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
      const Mat  diff(Dop[idx_d] * data);
      const auto delta(compute(Eop[idx_d] * data, ABS));
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
    result = Dop[idx_d] * data;
    break;
  case COLLECT_Y:
    result = compute(compute(data, DETECT_Y), ABS);
    break;
  case BUMP_Y:
    {
      initBump(data.rows());
      // log(|average(d_k C * exp(z_k))| / |dC|) == result.
      result = compute(A[idx_b] * data, ABS);
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
      Mat Dop0, Eop0;
      makeDI(data.rows() * 2 - 1, Dop0, Eop0);
      // N.B. nearest data in differential space.
      const Mat d0data(Dop0 * data / data.rows());
      makeDI(result.rows() * 2 - 1, Dop0, Eop0);
      const Mat ddata(Dop0 * result / result.rows());
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < result.cols(); i ++) {
        T a(0), b(0);
        for(int j = 0; j < d0data.rows(); j ++) {
          // <d0data, (ddata + Diff * t * e_k)> == <d0data, d0data> in cut.
          a += Dop0(j, Dop0.cols() - 1);
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
  case DEDGE:
    {
      assert(T(0) <= thedge && 0 < sq);
      vector<T> stat;
      for(int i = 0; i < data.rows(); i ++)
        for(int j = 0; j < data.cols(); j ++)
          stat.push_back(data(i, j));
      sort(stat.begin(), stat.end());
      result = data;
      for(int i = 0; i < data.rows(); i ++)
        for(int j = 0; j < data.cols(); j ++)
          if(data(i, j) <= stat[int((stat.size() - 1) * thedge)])
            result(i, j) = - T(1);
      bool flag(true);
      while(flag) {
        flag = false;
        bool flag2(false);
        for(int i = 0; i < result.rows(); i ++)
          for(int j = 0; j < result.cols(); j ++)
            if(result(i, j) < T(0)) {
              int cnt(0);
              T   sum(0);
              for(int ii = max(0, i - sq + 1); ii < min(i + sq, int(result.rows())); ii ++)
                for(int jj = max(0, j - sq + 1); jj < min(j + sq, int(result.cols())); jj ++)
                  if(T(0) <= result(ii, jj)) {
                    sum += result(ii, jj);
                    cnt ++;
                  }
              if(cnt) {
                result(i, j) = sum / cnt;
                flag2 = true;
              } else
                flag  = true;
            }
        flag = flag && flag2;
      }
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
      idx_d = i;
      return;
    }
  cerr << "n" << flush;
  idx_d = Dop.size();
  Mat vDop;
  Mat vEop;
  assert(2 <= size);
  makeDI(size, vDop, vEop);
  Dop.push_back(Mat(size, size));
  Eop.push_back(Mat(size, size));
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < Dop[idx_d].rows(); i ++)
    for(int j = 0; j < Dop[idx_d].cols(); j ++)
      Dop[idx_d](i, j) = Eop[idx_d](i, j) = T(0);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < Dop[idx_d].rows(); i ++)
    for(int j = 0; j < Dop[idx_d].cols() / 2; j ++) {
      Dop[idx_d](i, i / 2 + j) = vDop(i / 2, j);
      Eop[idx_d](i, i / 2 + j) = vEop(i / 2, j);
    }
  return;
}

template <typename T> void Filter<T>::initBump(const int& size) {
  for(int i = 0; i < A.size(); i ++)
    if(A[i].rows() == size) {
      idx_b = i;
      return;
    }
  cerr << "n" << flush;
  idx_b = A.size();
  assert(2 <= size);
  A.push_back(Mat(size, size));
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < size; i ++)
    for(int j = 0; j < size; j ++)
      A[idx_b](i, j) = T(0);
  Mat Dop0;
  Mat Eop0;
  makeDI(max(3, size / 16), Dop0, Eop0);
  Vec camera(2);
  camera[0] = T(0);
  camera[1] = T(1);
  assert(0 < dratio);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int zi = 0; zi < T(1) / dratio; zi ++) {
    for(int j = 0; j < Dop0.rows() / 2; j ++) {
      Vec cpoint(2);
      cpoint[0] = (j - T(Dop0.rows() / 2 - 1) / 2);
      cpoint[1] = T(zi) * dratio;
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
        for(int i = 0; i < A[idx_b].rows(); i ++)
          A[idx_b](i, getImgPt(i + y0, size)) += Dop0(Dop0.rows() / 2, j) * exp(T(int(T(1) / dratio) - zi) * sqrt(dratio));
      }
    }
  }
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

template <typename T> void Filter<T>::makeDI(const int& size, Mat& Dop, Mat& Eop, const bool& recursive) {
  assert(2 <= size);
  const auto ss((size + 1) / 2);
  Dop = Mat(ss, ss);
  for(int i = 0; i < Dop.rows(); i ++)
    for(int j = 0; j < Dop.cols(); j ++)
      Dop(i, j) = T(0);
  Eop = Mat(Dop);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int s = (recursive ? 2 : ss); s <= ss; s ++) {
          auto DFTD(seed(s, false));
    DFTD.row(0) *= U(0);
    const auto IDFT(seed(s, true));
          auto DFTE(DFTD);
    T ni(0);
    T nd(0);
    for(int i = 1; i < DFTD.rows(); i ++) {
      const U phase(- U(2.) * Pi * sqrt(U(- 1)) * T(i) / T(DFTD.rows()));
      // N.B. d/dy.
      DFTD.row(i) *= phase;
      // N.B. integrate.
      // DFTI.row(i) /= phase;
      nd += abs(phase) * abs(phase);
      ni += T(1) / (abs(phase) * abs(phase));
      // N.B. (d^(log(h))/dy^(log(h)) f, lim h -> 1. : nop.
      // DFTH.row(i) *= log(phase);
      // N.B. please refer enlarge.wxm, half freq space refer and uses each.
      DFTE.row(i) /= exp(sqrt(U(- 1)) * Pi / T(2 * DFTE.rows())) - U(T(1));
    }
    // similar to det(Iop * Dop) == 1.
    DFTD /= sqrt(ni * nd) * DFTD.rows();
    DFTE /= T(DFTE.rows() - 1);
#if defined(_WITHOUT_EIGEN_)
    Mat lDop((IDFT * DFTD).template real<T>());
    Mat lEop((IDFT * DFTE).template real<T>());
#else
    Mat lDop((IDFT * DFTD).real().template cast<T>());
    Mat lEop((IDFT * DFTE).real().template cast<T>());
#endif
#if defined(_OPENMP)
#pragma omp critical
#endif
    {
      for(int i = 0; i < Dop.rows(); i ++)
        for(int j = 0; j < lDop.cols(); j ++) {
          Dop(i, i * (Dop.cols() - lDop.cols()) / Dop.rows() + j) +=
            lDop(i * lDop.rows() / Dop.rows(), j);
          Eop(i, i * (Eop.cols() - lEop.cols()) / Eop.rows() + j) +=
            lEop(i * lEop.rows() / Eop.rows(), j);
        }
    }
  }
  if(recursive) {
    Dop /= ss - 1;
    Eop /= ss - 1;
  }
  return;
}

template <typename T> typename Filter<T>::MatU Filter<T>::seed(const int& size, const bool& idft) {
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

