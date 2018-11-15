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
 * When applying this to another sizes of images,
 * please create instances for each, or we get slow results.
 */
template <typename T> class enlarger2ex {
public:
  typedef enum {
    ENLARGE_X,
    ENLARGE_Y,
    ENLARGE_BOTH,
    ENLARGE_FBOTH,
    DETECT_X,
    DETECT_Y,
    DETECT_BOTH,
    DETECT_NOP_X,
    DETECT_NOP_Y,
    DETECT_NOP_BOTH,
    COLLECT_X,
    COLLECT_Y,
    COLLECT_BOTH,
    IDETECT_X,
    IDETECT_Y,
    IDETECT_BOTH,
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
    BCLIP,
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
  int  rec_tayl;
  
private:
  void initDop(const int& size);
  void initBump(const int& size);
  Vec  minSquare(const Vec& in);
  int  getImgPt(const T& y, const T& h);
  void makeDI(const int& size, Mat& Dop, Mat& Dhop, Mat& Iop, Mat& Eop);
  Mat  recursivePTayl(const Mat& A, const Mat& B, const Mat& ddxB, const Mat& ddyB, const Mat& B0, const int count, const T& dt = T(.25), const int count2 = 1);
  Mat  recursiveETayl(const Mat& A, const int count, const T& dt = T(.25), const int count2 = 1);
  Mat  recursiveNopInt(const Mat& A, const int count, const T& dt = T(.25), const int count2 = 1);
  U    I;
  T    Pi;
  vector<Mat> A;
  vector<Mat> B;
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
  dratio   = T(.0005);
  offset   = T(2) / T(256);
  rec_tayl = 4;
  idx_d    = - 1;
  idx_b    = - 1;
}

template <typename T> typename enlarger2ex<T>::Mat enlarger2ex<T>::compute(const Mat& data, const direction_t& dir) {
  Mat result;
  switch(dir) {
  case ENLARGE_BOTH:
    result = (compute(compute(data, ENLARGE_X), ENLARGE_Y) +
              compute(compute(data, ENLARGE_Y), ENLARGE_X)) / T(2);
    break;
  case ENLARGE_FBOTH:
    result = recursiveETayl(data, rec_tayl);
    break;
  case DETECT_BOTH:
    result = (compute(data, DETECT_X)  + compute(data, DETECT_Y)) / T(2);
    break;
  case DETECT_NOP_BOTH:
    result = (compute(data, DETECT_NOP_X) + compute(data, DETECT_NOP_Y)) / T(2);
    break;
  case COLLECT_BOTH:
    result = (compute(data, COLLECT_X) + compute(data, COLLECT_Y)) / T(2);
    break;
  case IDETECT_BOTH:
    result = (compute(data, IDETECT_X) + compute(data, IDETECT_Y)) / T(2);
    break;
  case BUMP_BOTH:
    result = (compute(data, BUMP_Y) +
              compute(data, BUMP_X) +
              compute(compute(compute(data, REVERSE_Y), BUMP_Y), REVERSE_Y) +
              compute(compute(compute(data, REVERSE_X), BUMP_X), REVERSE_X)) /
              T(4);
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
  case DETECT_NOP_X:
    result = compute(data.transpose(), DETECT_NOP_Y).transpose();
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
  case REVERSE_X:
    result = compute(data.transpose(), REVERSE_Y).transpose();
    break;
  case ENLARGE_Y:
    initDop(data.rows());
    result = Eop[idx_d] * data;
    break;
  case DETECT_Y:
    initDop(data.rows());
    result = Dop[idx_d] * data;
    break;
  case DETECT_NOP_Y:
    initDop(data.rows());
    result = Dhop[idx_d] * data;
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
      // N.B. similar to f(x) ~ f(x0) + f'(x0) * (x - x0) + f''(x0) * (x - x0)^ 2 / 2! + ...
      result = recursivePTayl(dataA, dataB,
                              compute(dataB, DETECT_X),
                              compute(dataB, DETECT_Y),
                              dataB, rec_tayl);
      // N.B. artificial.
      result = compute(recursiveNopInt(result, rec_tayl, T(2)), LOGSCALE);
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
      Mat Dop0, Dhop0, Iop0, Eop0;
      makeDI(data.rows() * 2 - 1, Dop0, Dhop0, Iop0, Eop0);
      // N.B. nearest data in differential space.
      const Mat d0data(Dop0 * data);
      makeDI(result.rows() * 2 - 1, Dop0, Dhop0, Iop0, Eop0);
      const Mat ddata(Dop0 * result);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < result.cols(); i ++) {
        T d0d(0), d0t(0), d0n(0);
        for(int j = 0; j < d0data.rows(); j ++) {
          // <d0data, ddata + Diff * t * e_k> == ||d0data||^2
          // <d0data, ddata> + <d0data, Diff * t * e_k> == ||d0data||^2
          d0d += d0data(j, i) * ddata(j, i);
          d0t += d0data(j, i) * Dop0(j, Dop0.cols() - 1);
          d0n += d0data(j, i) * d0data(j, i);
        }
        if(d0t == T(0)) {
          cerr << "NG in EXTEND_Y0" << endl;
          result(data.rows(), i) = data(data.rows() - 1, i);
        } else
          result(data.rows(), i) = (d0n - d0d) / d0t;
      }
      for(int i = 0; i < result.cols(); i ++)
        result(data.rows(), i) = max(T(0), min(T(1), result(data.rows(), i)));
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
  case REVERSE_Y:
    result = Mat(data.rows(), data.cols());
    for(int i = 0; i < data.rows(); i ++)
      result.row(result.rows() - i - 1) = data.row(i);
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
  Dhop.push_back(Mat());
  Iop.push_back(Mat());
  Eop.push_back(Mat());
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
#if defined(_WITH_EXTERNAL_)
  // reverted.
  // This works perfectly (from referring https://web.stanford.edu/class/cs448f/lectures/2.1/Sharpening.pdf via reffering Q&A sites.).
  // But I don't know whether this method is open or not.
  auto DFT2( seed(Eop[idx_d].rows(), false));
  auto IDFT2(seed(DFT2.rows(), true));
  for(int i = 0; i < DFT2.rows(); i ++)
    DFT2.row(i) *= T(1.5) - T(.5) * exp(- pow(T(i) / DFT2.rows(), T(2)));
#if defined(_WITHOUT_EIGEN_)
  Eop[idx_d] = (IDFT2 * DFT2).template real<T>() * Eop[idx_d];
#else
  Eop[idx_d] = (IDFT2 * DFT2).real().template cast<T>() * Eop[idx_d];
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
  B[idx_b] *= T(.5);
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
  DFTE /= T(2) * Pi * sqrt(T(DFTE.rows()));
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

template <typename T> typename enlarger2ex<T>::Mat enlarger2ex<T>::recursivePTayl(const Mat& A, const Mat& B, const Mat& ddxB, const Mat& ddyB, const Mat& B0, const int count, const T& dt, const int count2) {
  assert(0 <= count);
  const auto datadxA(compute(A, DETECT_X));
  const auto datadyA(compute(A, DETECT_Y));
        auto BB(B);
  const auto Bc(compute(B, BCLIP));
        Mat  datax(A.rows(), A.cols());
        Mat  datay(A.rows(), A.cols());
        Mat  datai(A.rows(), A.cols());
  for(int i = 0; i < BB.rows(); i ++)
    for(int j = 0; j < BB.cols(); j ++) {
      BB(i, j)   *= B0(i, j);
      // d/dt (A / (B0^n)) = ((d/dt A) * B0 - n * A * (d/dt B0)) / (B0^(n+1))
      datax(i, j) = datadxA(i, j) * B0(i, j) - count2 * A(i, j) * ddxB(i, j);
      datay(i, j) = datadyA(i, j) * B0(i, j) - count2 * A(i, j) * ddyB(i, j);
      datai(i, j) = A(i, j) / Bc(i, j);
    }
  if(count)
    // res = (A / B * 2 + integrate(d/dx (A / B), dx) + same for y) / n / 4.
    return (datai * T(2) +
            (compute(recursivePTayl(datax, BB, ddxB, ddyB, B0,
                                    count - 1, dt, count2 + 1), IDETECT_X) +
             compute(recursivePTayl(datay, BB, ddxB, ddyB, B0,
                                    count - 1, dt, count2 + 1), IDETECT_Y))
            * dt) / T(4) / count2;
  return datai;
}

template <typename T> typename enlarger2ex<T>::Mat enlarger2ex<T>::recursiveETayl(const Mat& A, const int count, const T& dt, const int count2) {
  assert(0 <= count);
  if(count)
    return (compute(A, ENLARGE_BOTH) * T(2) +
            (recursiveETayl(compute(A, DETECT_Y), count - 1, dt, count2 + 1) +
             recursiveETayl(compute(A, DETECT_X), count - 1, dt, count2 + 1)) *
           dt) / T(4) / count2;
  return compute(A, ENLARGE_BOTH);
}

template <typename T> typename enlarger2ex<T>::Mat enlarger2ex<T>::recursiveNopInt(const Mat& A, const int count, const T& dt, const int count2) {
  assert(0 <= count);
  if(count)
    return (A * T(2) +
            (recursiveNopInt(compute(A, DETECT_NOP_Y), count - 1, dt, count2 + 1) +
             recursiveNopInt(compute(A, DETECT_NOP_X), count - 1, dt, count2 + 1)) * dt) / T(4) / count2;
  return A;
}

#define _ENLARGE2X_
#endif

