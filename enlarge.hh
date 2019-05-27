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
 */
template <typename T> class enlarger2ex {
public:
  typedef enum {
    ENLARGE_X,
    ENLARGE_Y,
    ENLARGE_Y1,
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
    REVERSE_X,
    REVERSE_Y,
    REVERSE_BOTH,
    LANCZOS_X,
    LANCZOS_Y,
    LANCZOS_BOTH,
    D2Y,
    DEDGE,
    BCLIP,
    CLIP,
    ABS,
    EXPSCALE,
    LOGSCALE} direction_t;
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
  T    thedge;
  T    lanczos;
  T    sharpen;
  int  sq;
  // true for photo, false for illust.
  bool photo_illust;
  
private:
  void initDop(const int& size);
  void initBump(const int& size);
  int  getImgPt(const T& y, const T& h);
  void makeDI(const int& size, Mat& Dop, Mat& Eop, const bool& recursive = true);
  U    I;
  T    Pi;
  vector<Mat> A;
  vector<Mat> Dop;
  vector<Mat> Eop;
#if defined(_WITH_EXTERNAL_)
  vector<Mat> Hop;
#endif
  int idx_d;
  int idx_b;
};

template <typename T> enlarger2ex<T>::enlarger2ex() {
  I  = sqrt(U(- 1.));
  Pi = atan2(T(1), T(1)) * T(4);
  dratio  = T(.005);
  offset  = T(1) / T(256);
  thedge  = T(.05);
  sq      = 4;
  lanczos = T(1);
  sharpen = T(1);
  photo_illust = ! true;
  idx_d   = - 1;
  idx_b   = - 1;
}

template <typename T> typename enlarger2ex<T>::Mat enlarger2ex<T>::compute(const Mat& data, const direction_t& dir) {
  Mat result;
  switch(dir) {
  case ENLARGE_BOTH:
    result = (compute(compute(data, ENLARGE_X), ENLARGE_Y) +
              compute(compute(data, ENLARGE_Y), ENLARGE_X)) / T(2);
    break;
  case DETECT_BOTH:
    result = (compute(data, DETECT_X)  + compute(data, DETECT_Y)) / T(2);
    break;
  case COLLECT_BOTH:
    result = (compute(data, COLLECT_X) + compute(data, COLLECT_Y)) / T(2);
    break;
  case BUMP_BOTH:
    result = (compute(data, BUMP_X)    + compute(data, BUMP_Y) +
              compute(compute(compute(data, REVERSE_X), BUMP_X), REVERSE_X) +
              compute(compute(compute(data, REVERSE_Y), BUMP_Y), REVERSE_Y))
             / T(4);
    break;
  case LANCZOS_BOTH:
    result = (compute(data, LANCZOS_X) + compute(data, LANCZOS_Y)) / T(2);
    break;
  case EXTEND_BOTH:
    result = (compute(compute(data, EXTEND_X), EXTEND_Y) +
              compute(compute(data, EXTEND_Y), EXTEND_X)) / T(2);
    break;
  case REVERSE_BOTH:
    result = compute(compute(data, REVERSE_X), REVERSE_Y);
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
  case LANCZOS_X:
    result = compute(data.transpose(), LANCZOS_Y).transpose();
    break;
  case EXTEND_X:
    result = compute(data.transpose(), EXTEND_Y).transpose();
    break;
  case REVERSE_X:
    result = compute(data.transpose(), REVERSE_Y).transpose();
    break;
  case ENLARGE_Y1:
    {
      initDop(data.rows());
      const Mat dc(compute(Eop[idx_d] * data, ABS));
      const Mat delta(compute(data, DETECT_Y));
      result = Mat(data.rows() * 2, data.cols());
      for(int i = 0; i < data.rows(); i ++)
        for(int j = 0; j < data.cols(); j ++)
          if(delta(i, j) < 0) {
            result(2 * i,     j) = data(i, j) + dc(i, j);
            result(2 * i + 1, j) = data(i, j) - dc(i, j);
          } else {
            result(2 * i,     j) = data(i, j) - dc(i, j);
            result(2 * i + 1, j) = data(i, j) + dc(i, j);
          }
      result = compute(result, LANCZOS_Y);
    }
    break;
  case ENLARGE_Y:
    result = compute(data, ENLARGE_Y1);
#if defined(_WITH_EXTERNAL_)
    initDop(data.rows());
    result = Hop[idx_d] * result;
#endif
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
      // |average((dC)/C*z_k)| == dataA / dataB.
      result = Mat(data.rows(), data.cols());
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          result(i, j) = T(0);
      // N.B. there's better method with low freq, but that has a patent matter.
      initBump(data.rows());
      Mat work(data);
      for(int i0 = 1; 4 < work.rows(); i0 ++) {
        Mat w0(work);
        for( ; w0.rows() < data.rows(); )
          w0 = compute(w0, ENLARGE_Y1);
        Mat w1(data.rows(), data.cols());
        for(int i = 0; i < w1.rows(); i ++)
          w1.row(i) = w0.row(i);
        const auto dataA(compute(compute(A[idx_b] * w1, ABS), BCLIP));
        const auto dataB(compute(compute(w1, ABS), BCLIP));
        for(int i = 0; i < result.rows(); i ++)
          for(int j = 0; j < result.cols(); j ++)
            result(i, j) += dataA(i, j) / dataB(i, j);
        work = compute(work, D2Y);
      }
      result = - compute(result, LOGSCALE) * dratio;
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
      const Mat d0data(Dop0 * data);
      makeDI(result.rows() * 2 - 1, Dop0, Eop0);
      const Mat ddata(Dop0 * result);
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
      result.row(0) = compute(compute(data, REVERSE_Y), EXTEND_Y0).row(data.rows());
      for(int i = 0; i < data.rows(); i ++)
        result.row(i + 1) = data.row(i);
      result.row(data.rows() + 1) = compute(data, EXTEND_Y0).row(data.rows());
    }
    break;
  case LANCZOS_Y:
    {
      result = data;
      for(int i = 0; i < result.rows(); i ++)
        for(int j = - int(lanczos); j <= int(lanczos); j ++)
          if(j)
            result.row(i) += data.row(max(min(j + i, int(data.rows() - 1)), 0)) * sin(Pi * j / T(2)) * sin(Pi * j / T(2) / lanczos) / (Pi * Pi * j * j / T(2) / T(2) / lanczos);
      T sum(1);
      for(int j = - int(lanczos); j <= int(lanczos); j ++)
        if(j)
          sum += sin(Pi * j / T(2)) * sin(Pi * j / T(2) / lanczos) / (Pi * Pi * j * j / T(2) / T(2) / lanczos);
      result /= sum;
    }
    break;
  case REVERSE_Y:
    {
      result = Mat(data.rows(), data.cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int i = 0; i < data.rows(); i ++)
        result.row(result.rows() - i - 1) = data.row(i);
    }
    break;
  case D2Y:
    result = Mat((data.rows() + 1) / 2, data.cols());
    for(int i = 0; i < result.rows(); i ++) {
      result.row(i) = data.row(i * 2);
      if(i * 2 + 1 < data.rows()) {
        result.row(i) += data.row(i * 2 + 1);
        result.row(i) /= T(2);
      }
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
  Eop.push_back(Mat());
  Mat vDop;
  Mat vEop;
  assert(2 <= size);
  makeDI(size, vDop, vEop);
  Dop[idx_d]  = Mat(size, size);
  Eop[idx_d]  = Mat(size, size);
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
#if defined(_WITH_EXTERNAL_)
  // This works perfectly (from referring https://web.stanford.edu/class/cs448f/lectures/2.1/Sharpening.pdf via reffering Q&A sites.).
  // But I don't know whether this method is open or not.
  auto DFT2(seed(Eop[idx_d].rows() * 2, false));
  for(int i = 0; i < DFT2.rows(); i ++)
    DFT2.row(i) *= sharpen + T(1) - sharpen * exp(- pow(T(i) / DFT2.rows(), T(2)));
  Hop.push_back(Mat());
#if defined(_WITHOUT_EIGEN_)
  Hop[idx_d] = (seed(DFT2.rows(), true) * DFT2).template real<T>();
#else
  Hop[idx_d] = (seed(DFT2.rows(), true) * DFT2).real().template cast<T>();
#endif
#endif
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
  assert(2 <= size);
  A[idx_b] = Mat(size, size);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < size; i ++)
    for(int j = 0; j < size; j ++)
      A[idx_b](i, j) = T(0);
  Mat Dop0;
  Mat Eop0;
  makeDI(max(3, int(T(2) * sqrt(T(size)))), Dop0, Eop0);
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
            if(photo_illust)
              A[idx_b](i, getImgPt(i + y0, size)) += Dop0(Dop0.rows() / 2, j) * (exp(T(1) / dratio) - exp(T(zi + 1)));
            else
              A[idx_b](i, getImgPt(i + y0, size)) += Dop0(Dop0.rows() / 2, j) * (- exp(T(zi + 1)));
          }
        else
          for(int i = 0; i < A[idx_b].rows(); i ++) {
            if(photo_illust)
              A[idx_b](i, getImgPt(i + y0, size)) += (Dop0(Dop0.rows() / 2, j) + Dop0(Dop0.rows() / 2 + 1, j)) / T(2) * (exp(T(1) / dratio) - exp(T(zi + 1)));
            else
              A[idx_b](i, getImgPt(i + y0, size)) += (Dop0(Dop0.rows() / 2, j) + Dop0(Dop0.rows() / 2 + 1, j)) / T(2) * (- exp(T(zi + 1)));
          }
      }
    }
  }
  return;
}

template <typename T> int enlarger2ex<T>::getImgPt(const T& y, const T& h) {
  int yy(int(y) % int(2 * h));
  if(yy < 0)
    yy = - yy;
  if(yy >= int(h))
    yy = int(h) - (yy - int(h));
  return yy % int(h);
}

template <typename T> void enlarger2ex<T>::makeDI(const int& size, Mat& Dop, Mat& Eop, const bool& recursive) {
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
    T nd(0), ni(0);
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
      // N.B. please refer enlarge.wxm, uses each freq.
      //      b(t) -> i * sin(phase) / (cos(phase) - 1) * f(t) for each phase.
      const T phase2(Pi * T(i) / T(DFTE.rows()));
      const U r(sqrt(U(- 1)) * sin(phase2) / (cos(phase2) - U(1)));
      DFTE.row(i) *= r;
    }
    // N.B. similar to DFTI * DFTD == id.
    const T ratio(sqrt(nd * ni));
    DFTD /= ratio;
    DFTE /= T(DFTE.rows()) * T(DFTE.rows());
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

#define _ENLARGE2X_
#endif

