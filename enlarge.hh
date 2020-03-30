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
using std::sort;
using std::vector;
using std::max;
using std::min;
template <typename T> class reDig;

// This class is NOT thread safe.
template <typename T> class Filter {
public:
  typedef enum {
    SHARPEN_X,
    SHARPEN_Y,
    SHARPEN_BOTH,
    DETECT_X,
    DETECT_Y,
    DETECT_BOTH,
    COLLECT_X,
    COLLECT_Y,
    COLLECT_BOTH,
    BUMP_X,
    BUMP_Y_SHALLOW,
    BUMP_Y,
    BUMP_BOTH,
    EXTEND_X,
    EXTEND_Y,
    EXTEND_BOTH,
    BCLIP,
    CLIP,
    ABS,
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
  ~Filter();
  Mat  compute(const Mat& data, const direction_t& dir);
  MatU seed(const int& size, const bool& idft);
  Mat  gmean(const Mat& a, const Mat& b);
  T    dratio;
  T    offset;
  T    dbratio;
  int  plen;
  int  lrecur;
  int  bumpd;

private:
  void initDop(const int& size);
  int  getImgPt(const int& y, const int& h);
  T    Pi;
  vector<Mat> Dop;
  vector<Mat> Sop;
  int  idx;
};

template <typename T> Filter<T>::Filter() {
  Pi = atan2(T(1), T(1)) * T(4);
  // N.B. from accuracy reason, low depth.
  dratio  = T(01) / T(10);
  dbratio = T(01) / T(10);
  offset  = T(1)  / T(64);
  plen    = 1;
  lrecur  = 8;
  bumpd   = 65;
  idx     = - 1;
}

template <typename T> Filter<T>::~Filter() {
  ;
}

template <typename T> typename Filter<T>::Mat Filter<T>::compute(const Mat& data, const direction_t& dir) {
  Mat result;
  switch(dir) {
  case SHARPEN_BOTH:
    result = compute(compute(data, SHARPEN_X), SHARPEN_Y);
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
    result = compute(compute(data, EXTEND_X), EXTEND_Y);
    break;
  case SHARPEN_X:
    result = compute(data.transpose(), SHARPEN_Y).transpose();
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
  case SHARPEN_Y:
    initDop(data.rows());
    result = Sop[idx] * data;
    break;
  case DETECT_Y:
    initDop(data.rows());
    result = Dop[idx] * data;
    break;
  case COLLECT_Y:
    result = compute(compute(data, DETECT_Y), ABS);
    break;
  case BUMP_Y_SHALLOW:
    {
      result = Mat(data.rows(), data.cols());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          result(i, j) = T(0);
      initDop(bumpd);
      Vec camera(2);
      camera[0] = T(0);
      camera[1] = T(1);
      assert(T(0) < dratio);
      const auto rxy(sqrt(T(data.rows()) * T(data.cols())));
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
      for(int zi = 0; T(zi) < T(1) / dratio; zi ++) {
        Mat A(data.rows(), data.cols());
        for(int i = 0; i < A.rows(); i ++)
          for(int j = 0; j < A.cols(); j ++)
            A(i, j) = T(0);
        const Vec Dop0(Dop[idx].row(Dop[idx].rows() / 2) * exp(T(zi) / sqrt(dratio)));
        for(int j = 0; j < Dop0.size(); j ++) {
          Vec cpoint(2);
          cpoint[0] = (T(j) - T(Dop0.size() - 1) / T(2)) / rxy;
          cpoint[1] = T(zi) * dratio / T(2);
          // x-z plane projection of point p with camera geometry c to z=0.
          // c := camera, p := cpoint.
          // <c + (p - c) * t, [0, 1]> = 0
          const auto t(- camera[1] / (cpoint[1] - camera[1]));
          const auto y0((camera + (cpoint - camera) * t)[0] * rxy);
          //  N.B. average_k(dC_k / dy * z_k).
          for(int i = 0; i < A.rows(); i ++)
            A.row(i) += data.row(getImgPt(i + int(y0), data.rows())) * Dop0[j];
        }
#if defined(_OPENMP)
#pragma omp atomic
#endif
        result += compute(A, ABS);
      }
      // N.B. log(|average(d_k C * exp(z_k))| / |dC|) == result.
      const auto dC(compute(data, COLLECT_Y));
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          result(i, j) /= dC(i, j);
      result = - compute(compute(result, BCLIP), LOGSCALE) * sqrt(dratio);
    }
    break;
  case BUMP_Y:
    {
      reDig<T> redig;
      const auto bump(compute(data, BUMP_Y_SHALLOW));
      const auto t0(redig.tilt(redig.makeRefMatrix(data, 1), bump,
                               redig.tiltprep(bump, 0, 4, dbratio)));
      const auto data0(redig.pullRefMatrix(t0, 1, data));
      const auto data1(redig.tilt(data, bump, redig.tiltprep(bump, 2, 4, dbratio)) );
      result = Mat(data.rows(), data.cols());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          result(i, j) = T(0);
      assert(T(0) < dratio);
      initDop(bumpd);
      const auto rxy(sqrt(T(data0.rows()) * T(data0.cols())));
      Vec camera(2);
      camera[0] = T(0);
      camera[1] = T(20);
      auto work(result);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
      for(int zi = 0; T(zi) < T(1) / dratio; zi ++) {
        Mat A(data0.rows(), data0.cols());
        for(int i = 0; i < A.rows(); i ++)
          for(int j = 0; j < A.cols(); j ++)
            A(i, j) = T(0);
        Mat C(A);
        const Vec Dop0(Dop[idx].row(Dop[idx].rows() / 2) * exp(T(zi) / sqrt(dratio)));
        for(int i = 0; i < A.rows(); i ++) {
          for(int j = 0; j < Dop0.size(); j ++) {
            Vec cpoint(2);
            cpoint[0] = (T(j) - T(Dop0.size() - 1) / T(2) + T(i) - T(data0.rows() - 1) / T(2)) / rxy;
            cpoint[1] = T(zi) * dratio / T(2);
            // x-z plane projection of point p with camera geometry c to z=0.
            // c := camera, p := cpoint.
            // <c + (p - c) * t, [0, 1]> = 0
            const auto t(- camera[1] / (cpoint[1] - camera[1]));
            const auto y((camera + (cpoint - camera) * t)[0] * rxy + T(data0.rows() - 1) / T(2));
            // N.B. average_k(dC_k / dy * z_k).
            const auto idx(getImgPt(int(y), data0.rows()));
            if(j <= Dop0.size() / 2)
              A.row(i) += data0.row(idx) * Dop0[j];
            if(Dop0.size() / 2 <= j)
              A.row(i) += data1.row(idx) * Dop0[j];
            if(j == Dop0.size() / 2)
              C.row(i) += data0.row(idx) + data1.row(idx);
          }
        }
        auto lres(compute(A, ABS));
        C = compute(compute(C, ABS), BCLIP);
        for(int i = 0; i < C.rows(); i ++)
          for(int j = 0; j < C.cols(); j ++)
            lres(i, j) /= C(i, j);
#if defined(_OPENMP)
#pragma omp atomic
#endif
        result += lres;
      }
      result = compute(compute(result, BCLIP), LOGSCALE) * sqrt(dratio);
      int ii(0);
      for(int i = 0; i < t0.rows(); i ++) {
        bool flg(true);
        for(int j = 0; j < t0.cols(); j ++)
          flg = flg &&
            t0(i, j) != T(0) && t0(t0.rows() - i - 1, j) != T(0);
        if(flg) {
          ii = i;
          break;
        }
      }
      MatU b0(t0.rows() - ii * 2, t0.cols());
      for(int i = 0; i < b0.rows(); i ++)
        for(int j = 0; j < b0.cols(); j ++)
          b0(i, j) = U(result(i + ii, j));
      b0 = seed(b0.rows(), false) * b0 / sqrt(T(b0.rows()));
      MatU b1(t0.rows(), t0.cols());
      for(int i = 0; i < b0.rows(); i ++)
        b1.row(i) = b0.row(i);
      for(int i = b0.rows(); i < b1.rows(); i ++)
        b1.row(i) *= U(T(0));
#if defined(_WITHOUT_EIGEN_)
      result = (seed(t0.rows(), true) * b1).template real<T>() * sqrt(T(t0.rows()));
#else
      result = (seed(t0.rows(), true) * b1).real().template cast<T>() * sqrt(T(t0.rows()));
#endif
    }
    break;
  case EXTEND_Y:
    {
      assert(0 < plen);
      result = Mat(data.rows() + plen * 2, data.cols());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < data.rows(); i ++)
        result.row(i + plen) = data.row(i);
      Mat rdata(data.rows(), data.cols());
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
      for(int i = 0; i < data.rows(); i ++)
        rdata.row(data.rows() - 1 - i) = data.row(i);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
      for(int j = 0; j < plen; j ++) {
        P0<T> p(data.rows(), 2, j + 1);
        for(int i = 0; i < data.cols(); i ++) {
          result(data.rows() + j + plen, i) = p.next(data.col(i));
          result(plen - j - 1, i) = p.next(rdata.col(i));
        }
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
  case LOGSCALE:
    {
      result = Mat(data.rows(), data.cols());
      // N.B. before to use, please BCLIP.
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          result(i, j) = (data(i, j) < T(0) ? - T(1) : T(1)) * log(abs(data(i, j)) / offset) * offset;
    }
    break;
  default:
    assert(0 && "unknown command in Filter (should not be reached.)");
  }
  return result;
}

template <typename T> void Filter<T>::initDop(const int& size) {
  assert(2 <= size);
  for(int i = 0; i < Dop.size(); i ++)
    if(Dop[i].rows() == size) {
      idx = i;
      return;
    }
  cerr << "n" << flush;
  idx = Dop.size();
  Dop.push_back(Mat(size, size));
  for(int i = 0; i < Dop[idx].rows(); i ++)
    for(int j = 0; j < Dop[idx].cols(); j ++)
      Dop[idx](i, j) = T(0);
  Sop.push_back(Dop[idx]);
  int cnt(1);
  for(int ss = 2; ss <= size; ss *= 2, cnt ++) {
          auto DFTD(seed(ss, false));
          auto DFTL(seed(ss * 2, false));
    const auto IDFT(seed(ss, true));
    DFTD.row(0) *= U(T(0));
          auto DFTE(DFTD);
    T ni(0);
    T nd(0);
    for(int i = 1; i < DFTD.rows(); i ++) {
      const auto phase(- U(T(2)) * Pi * U(T(0), T(1)) * T(i) / T(DFTD.rows()));
      const auto phase2(U(T(1)) / phase);
      // N.B. d/dy with sampling theorem.
      if(i < DFTD.rows() / 2) {
        DFTD.row(i) *= phase;
        nd += abs(phase)  * abs(phase);
        ni += abs(phase2) * abs(phase2);
      } else
        DFTD.row(i) *= U(T(0));
      // N.B. integrate.
      // DFTI.row(i) /= phase;
      // N.B. (d^(log(h))/dy^(log(h)) f, lim h -> 1. : nop.
      // DFTH.row(i) *= log(phase);
      // N.B. please refer enlarge.wxm, half freq space refer and uses each.
      DFTE.row(i) /= exp(U(T(0), T(1)) * Pi * U(T(i)) / T(DFTE.rows())) - U(T(1));
    }
    // N.B. similar to det(Dop * Iop) == det(Dop) * det(Iop) == 1,
    //      but Dop * Iop == I in ideal (Iop.row(0) == NaN) case.
    //      with omitting IDFT and DFT on both side, ||.||_2^2,
    //      in matrix-matrix operation:
    //      ||Dop * Iop * x|| / ||x|| == sum(sum(d_k*i_l*x_m)^2)/sum(x_k^2)
    //                                == sum(d_k^2)*sum(i_l^2)*cos theta cos psi
    //                                == (n - 1) * cos psi
    //      in matrix-vector operation:
    //      ||Dop * Iop * x|| / ||x|| == sum(sum(d_k*i_l*x_m)^2)/sum(x_k^2)
    //                                == sum(d_k^2)*cos theta'*sum(i_l^2)cos phi
    //                                == ||Dop|| ||Iop|| cos theta' cos phi
    //      so we choose matrix-vector operation with matrix-matrix style,
    //      because of cosine range, with normal ||.||_2, we choose:
    //        Dop' := Dop sqrt(sqrt(n - 1)) / sqrt(||Dop|| ||Iop||).
    //        Iop' := Iop sqrt(sqrt(n - 1)) / sqrt(||Dop|| ||Iop||).
    //      then we get:
    //        Iop' * Dop' == Iop * Dop * sqrt(n - 1) / (||Dop|| * ||Iop||)
    //      And, if we change coefficients ratio on differential operator,
    //        and its inverse of integrate operator, it causes invalid
    //        on the meaning of DFT core, but in experiment,
    //          if ||Dop|| ||Iop|| == 1 in that meaning,
    //          exists r in R, Dop * x == r * x results gains.
    DFTD *= nd * ni == T(0) ? T(0) : sqrt(sqrt(T(DFTD.rows() - 1) / (nd * ni)));
    DFTE /= T(DFTE.rows()  - 1);
#if defined(_WITHOUT_EIGEN_)
    const Mat lDop((IDFT * DFTD).template real<T>());
    const Mat EE((IDFT * DFTE).template real<T>());
    const Mat LL((seed(ss * 2, true) * DFTL).template real<T>());
#else
    const Mat lDop((IDFT * DFTD).real().template cast<T>());
    const Mat EE((IDFT * DFTE).real().template cast<T>());
    const Mat LL((seed(ss * 2, true) * DFTL).real().template cast<T>());
#endif
    Mat Eop(ss * 2, ss);
    Mat Lop(ss, ss * 2);
    for(int i = 0; i < EE.rows(); i ++) {
      Eop.row(i * 2 + 0) = - EE.row(i);
      Eop.row(i * 2 + 1) =   EE.row(i);
      Eop(i * 2 + 0, i) += T(1);
      Eop(i * 2 + 1, i) += T(1);
      Lop.row(i)         = LL.row(i * 2);
    }
    const Mat lLop(Lop * Eop);
    for(int i = 0; i < Dop[idx].rows(); i ++)
      for(int j = 0; j < lDop.rows(); j ++) {
        int ij(i - lDop.rows() / 2 + j);
        int jj(lDop.rows() / 2);
        if(i < lDop.rows() / 2) {
          ij = j;
          jj = i;
        } else if(Dop[idx].rows() - i < (lDop.rows() + 1) / 2) {
          ij = j - lDop.rows() + Dop[idx].rows();
          jj = i - Dop[idx].rows() + lDop.rows();
        }
        Dop[idx](i, ij) += lDop(jj, j);
        Sop[idx](i, ij) += lLop(jj, j);
      }
  }
  Dop[idx] /= T(cnt);
  const Mat SS(Sop[idx]);
  for(int i = 0; i < SS.rows(); i ++)
     for(int j = 0; j < SS.cols(); j ++)
      Sop[idx](SS.rows() - i - 1, SS.cols() - j - 1) += SS(i, j);
  T mnorm(0);
  for(int i = 0; i < Sop[idx].rows(); i ++)
    mnorm = max(mnorm, Sop[idx].row(i).dot(Sop[idx].row(i)));
  Sop[idx] /= sqrt(mnorm);
  for(int i = 0; i < lrecur; i ++)
    Sop[idx] = Sop[idx] * Sop[idx];
  mnorm = T(0);
  for(int i = 0; i < Sop[idx].rows(); i ++)
    mnorm = max(mnorm, Sop[idx].row(i).dot(Sop[idx].row(i)));
  Sop[idx] /= sqrt(mnorm);
  mnorm = T(0);
  for(int i = 0; i < Dop[idx].rows(); i ++)
    mnorm = max(mnorm, Dop[idx].row(i).dot(Dop[idx].row(i)));
  Dop[idx] /= sqrt(mnorm);
  return;
}

template <typename T> int Filter<T>::getImgPt(const int& y, const int& h) {
  int yy(y % (2 * h));
  if(yy < 0)
    yy = - yy;
  if(yy >= h)
    yy = h - (yy - h);
  return yy % h;
}

template <typename T> typename Filter<T>::MatU Filter<T>::seed(const int& size, const bool& idft) {
  MatU result(size, size);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < result.rows(); i ++)
    for(int j = 0; j < result.cols(); j ++) {
      const auto phase(- T(2) * Pi * T(idft ? - 1 : 1) * T(i) * T(j) / T(size));
      result(i, j) = U(cos(phase), sin(phase)) / T(idft ? size : 1);
    }
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

