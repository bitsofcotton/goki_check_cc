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
using std::vector;
using std::sort;
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
    BUMP_Y_SHORT,
    BUMP_Y0,
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
  Mat  bump2(const Mat& data0, const Mat& data1, const T& pixels = T(1));
  MatU seed(const int& size, const bool& idft);
  Mat  gmean(const Mat& a, const Mat& b);
  void reinit();
  T    dratio;
  T    offset;
  int  pstart;
  int  pend;
  
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
  I  = sqrt(U(- T(1)));
  Pi = atan2(T(1), T(1)) * T(4);
  // N.B. from accuracy reason, low depth.
  dratio  = T(005) / T(100);
  offset  = T(1) / T(64);
  pstart  = 1;
  pend    = 8;
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
    result = gmean(compute(data, BUMP_X ), compute(data, BUMP_Y ));
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
    initDop(data.rows());
    result = Eop[idx] * data;
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
      result = compute(data, BUMP_Y0);
      Mat shrink(data);
      for(int i0 = 1; 16 <= shrink.rows(); i0 ++) {
        Mat work((shrink.rows() + 1) / 2, shrink.cols());
        for(int i = 0; i < work.rows(); i ++)
          if(i * 2 + 1 < shrink.rows())
            work.row(i) = (shrink.row(i * 2) + shrink.row(i * 2 + 1)) / T(2);
          else
            work.row(i) =  shrink.row(i * 2);
        shrink = work;
        while(work.rows() < data.rows())
          work = compute(work, ENLARGE_Y);
        Mat work2(data.rows(), data.cols());
        for(int i = 0; i < work2.rows(); i ++)
          work2.row(i) = work.row(i);
        result += compute(work2, BUMP_Y0);
      }
      result = - compute(result, LOGSCALE) * sqrt(dratio);
    }
    break;
  case BUMP_Y_SHORT:
    result = - compute(compute(data, BUMP_Y0), LOGSCALE) * sqrt(dratio);
    break;
  case BUMP_Y0:
    {
      result = Mat(data.rows(), data.cols());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          result(i, j) = T(0);
      // XXX:
      // initDop(max(3, min(int(data.rows()) / 16, int(T(1) / dratio / dratio))));
      initDop(60);
      Vec camera(2);
      camera[0] = T(0);
      camera[1] = T(1);
      assert(T(0) < dratio);
      const auto rxy(sqrt(T(data.rows()) * T(data.cols())));
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
      for(int zi = 0; T(zi) < T(1) / dratio; zi ++) {
        Mat A(data.rows(), data.rows());
        for(int i = 0; i < A.rows(); i ++)
          for(int j = 0; j < A.cols(); j ++)
            A(i, j) = T(0);
        const Vec Dop0(Dop[idx].row(Dop[idx].rows() / 2) * exp(T(zi)));
        for(int j = 0; j < Dop0.size(); j ++) {
          Vec cpoint(2);
          cpoint[0] = (T(j) - T(Dop0.size() - 1) / T(2)) / rxy;
          cpoint[1] = T(zi) * dratio;
          // x-z plane projection of point p with camera geometry c to z=0.
          // c := camera, p := cpoint.
          // <c + (p - c) * t, [0, 1]> = 0
          const auto t(- camera[1] / (cpoint[1] - camera[1]));
          const auto y0((camera + (cpoint - camera) * t)[0] * rxy);
          // N.B. average_k(dC_k / dy * z_k).
          for(int i = 0; i < A.rows(); i ++)
            A(i, getImgPt(i + int(y0), data.rows())) += Dop0[j];
        }
#if defined(_OPENMP)
#pragma omp atomic
#endif
        result += compute(A * data, ABS);
      }
      // N.B. log(|average(d_k C * exp(z_k))| / |dC|) == result.
      const auto dC(compute(compute(data, COLLECT_Y), BCLIP));
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          result(i, j) /= dC(i, j);
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
      assert(0 < pstart && pstart < pend);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < data.cols(); i ++) {
        for(int k = pstart; k < min(pend, int(data.rows()) / 8); k ++) {
          // out of range prediction causes low frequency.
          P0<T, complex<T> > p(data.rows() * 2 / k);
          for(int j = (data.rows() % k + k - 1) % k; j < data.rows() - k; j += k)
            p.nextNoreturn(data(j, i));
          result(data.rows(), i) += p.next(data(data.rows() - 1, i));
        }
        result(data.rows(), i) /= T(min(pend, int(data.rows()) / 8) - pstart);
      }
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
      result = compute(result, CLIP);
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

template <typename T> typename Filter<T>::Mat Filter<T>::bump2(const Mat& data0, const Mat& data1, const T& pixels) {
  assert(data0.rows() == data1.rows() && data0.cols() == data0.cols());
  Mat result(data0.rows(), data0.cols());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < result.rows(); i ++)
    for(int j = 0; j < result.cols(); j ++)
      result(i, j) = T(1);
  assert(T(0) < dratio);
  initDop(max(3, min(int(data0.rows()) / 16, int(T(1) / dratio / dratio))));
  const auto dC(compute(compute(data0, COLLECT_Y) + compute(data1, COLLECT_Y), ABS));
  const auto rxy(sqrt(T(data0.rows()) * T(data0.cols())));
  Vec camera0(2);
  Vec camera1(2);
  camera0[0] =   pixels / rxy;
  camera0[1] =   T(1);
  camera1[0] = - pixels / rxy;
  camera1[1] =   T(1);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int zi = 0; T(zi) < T(1) / dratio; zi ++) {
    Mat A(data0.rows(), data0.rows());
    for(int i = 0; i < A.rows(); i ++)
      for(int j = 0; j < A.cols(); j ++)
        A(i, j) = T(0);
    Mat B(A);
    const Vec Dop0(Dop[idx].row(Dop[idx].rows() / 2) * exp(T(zi)));
    for(int i = 0; i < A.rows(); i ++)
      for(int j = 0; j < Dop0.size(); j ++) {
        Vec cpoint(2);
        cpoint[0] = (T(j) - T(Dop0.size() - 1) / T(2)) / rxy;
        cpoint[1] = T(zi) * dratio;
        // x-z plane projection of point p with camera geometry c to z=0.
        // c := camera, p := cpoint.
        // <c + (p - c) * t, [0, 1]> = 0
        const auto t0(- camera0[1] / (cpoint[1] - camera0[1]));
        const auto t1(- camera1[1] / (cpoint[1] - camera1[1]));
        const auto y0((camera0 + (cpoint - camera0) * t0)[0] * rxy);
        const auto y1((camera1 + (cpoint - camera1) * t1)[0] * rxy);
        // N.B. average_k(dC_k / dy * z_k).
        A(i, getImgPt(i + int(y0), data0.rows())) += Dop0[j];
        B(i, getImgPt(i + int(y1), data0.rows())) += Dop0[j];
      }
    const auto work(compute(compute(A * data0 - B * data1, ABS), BCLIP));
    for(int i = 0; i < result.rows(); i ++)
      for(int j = 0; j < result.cols(); j ++)
        result(i, j) += dC(i, j) / work(i, j);
  }
  return compute(result, LOGSCALE) * dratio;
}

template <typename T> void Filter<T>::initDop(const int& size) {
  for(int i = 0; i < Dop.size(); i ++)
    if(Dop[i].rows() == size) {
      idx = i;
      return;
    }
  cerr << "n" << flush;
  idx = Dop.size();
  Dop.push_back(Mat(size, size));
  int cnt(0);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < size; i ++)
    for(int j = 0; j < size; j ++)
      Dop[idx](i, j) = T(0);
  Eop.push_back(Dop[idx]);
  assert(2 <= size);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
#if defined(_RECURSIVE_)
  for(int lsize = 2; lsize <= size; lsize *= 2) {
#else
  for(int lsize = size; lsize <= size; lsize *= 2) {
#endif
          auto DFTD(seed(lsize, false));
    const auto IDFT(seed(lsize, true));
    DFTD.row(0) *= U(T(0));
          auto DFTE(DFTD);
    T ni(0);
    T nd(0);
    for(int i = 1; i < DFTD.rows(); i ++) {
      const auto phase(- U(T(2)) * Pi * sqrt(U(- T(1))) * T(i) / T(DFTD.rows()));
      const auto phase2(U(T(1)) / phase);
      // N.B. d/dy.
      DFTD.row(i) *= phase;
      nd += abs(phase)  * abs(phase);
      ni += abs(phase2) * abs(phase2);
      // N.B. integrate.
      // DFTI.row(i) /= phase;
      // N.B. (d^(log(h))/dy^(log(h)) f, lim h -> 1. : nop.
      // DFTH.row(i) *= log(phase);
      // N.B. please refer enlarge.wxm, half freq space refer and uses each.
      DFTE.row(i) /= exp(sqrt(U(- T(1))) * Pi * U(T(i)) / T(DFTE.rows())) - U(T(1));
    }
    // N.B. similar to det(Dop * Iop) == det(Dop) * det(Iop) == 1,
    //      but Dop * Iop == I in ideal (Iop.row(0) == NaN) case.
    //      in matrix-matrix operation:
    //      ||Dop * Iop * x|| / ||x|| == sum((d_k*i_k*x_k)^2)/sum(x_k^2)
    //                                == sum(d_k^2*i_k^2)*cos theta cos psi
    //                                == sqrt(n - 1) * cos psi
    //      in matrix-vector operation:
    //      ||Dop * Iop * x|| / ||x|| == sum((d_k*i_k*x_k)^2)/sum(x_k^2)
    //                                == sum(d_k^2)*cos theta'*sum(i_k^2)cos phi
    //                                == ||Dop|| ||Iop|| cos theta' cos phi
    //      so we choose matrix-vector operation with matrix-matrix style,
    //      because of cosine range, we choose:
    //        Dop' := Dop sqrt(n - 1) / sqrt(||Dop|| ||Iop||).
    //      (sqrt instead of sqrt(sqrt(...)) is because of
    //        the ratio is applied to differential operator itself.)
    //      And, if we change coefficients ratio on differential operator,
    //        and its inverse of integrate operator, it causes invalid
    //        on the meaning of DFT core, but in experiment,
    //          if ||Dop|| ||Iop|| == 1 in that meaning,
    //          exists r in R, Dop * x == r * x results gains.
    DFTD *= sqrt(T(DFTD.rows() - 1) / sqrt(nd * ni));
    DFTE /= T(DFTE.rows()) - T(1);
#if defined(_WITHOUT_EIGEN_)
    const Mat lDop((IDFT * DFTD).template real<T>());
    const Mat lEop((IDFT * DFTE).template real<T>());
#else
    const Mat lDop((IDFT * DFTD).real().template cast<T>());
    const Mat lEop((IDFT * DFTE).real().template cast<T>());
#endif
#if defined(_OPENMP)
#pragma omp critical
#endif
    for(int i = 0; i < Dop[idx].rows(); i ++)
      for(int j = 0; j < lDop.cols(); j ++) {
        Dop[idx](i, i * (Dop[idx].cols() - lDop.cols()) / Dop[idx].cols() + j) +=
          lDop(i * lDop.rows() / Dop[idx].rows(), j);
        Eop[idx](i, i * (Eop[idx].cols() - lEop.cols()) / Eop[idx].cols() + j) +=
          lEop(i * lEop.rows() / Eop[idx].rows(), j);
      }
    cnt ++;
  }
  Dop[idx] /= T(cnt);
  Eop[idx] /= T(cnt);
  Mat Eop2x(Eop[idx].rows() * 2, Eop[idx].cols());
  for(int i = 0; i < Eop[idx].rows(); i ++) {
    Eop2x.row(2 * i + 0) = - Eop[idx].row(i);
    Eop2x.row(2 * i + 1) =   Eop[idx].row(i);
    Eop2x(2 * i + 0, i) += T(1);
    Eop2x(2 * i + 1, i) += T(1);
  }
  Eop[idx] = Eop2x;
  return;
}

template <typename T> int Filter<T>::getImgPt(const T& y, const T& h) {
  int yy(int(y) % int(T(2) * h));
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
      result(i, j) = exp(I * U(- T(2) * Pi * T(idft ? - 1 : 1) * T(i) * T(j) / T(size))) / T(idft ? size : 1);
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

