/* BSD 3-Clause License:
 * Copyright (c) 2018-2020, bitsofcotton.
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
using std::abs;

// This class is NOT thread safe.
template <typename T> class Filter {
public:
  typedef enum {
    SHARPEN_X,
    SHARPEN_Y,
    SHARPEN_BOTH,
    ENLARGE_X,
    ENLARGE_Y,
    ENLARGE_BOTH,
    INTEG_X,
    INTEG_Y,
    RINTEG_X,
    RINTEG_Y,
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
    EXTEND_Y,
    EXTEND_BOTH,
    CLIP,
    ABS } direction_t;
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
  int  dist;
  int  dratio;
  int  plen;
  int  lrecur;
  int  bumpd;

private:
  void initDop(const int& size);
  int  getImgPt(const int& y, const int& h);
  T    Pi;
  vector<Mat> Dop;
  vector<Mat> Iop;
  vector<Mat> RIop;
  vector<Mat> Eop;
  vector<Mat> Sop;
  int  idx;
};

template <typename T> Filter<T>::Filter() {
  Pi = atan2(T(1), T(1)) * T(4);
  dist   = 1;
  dratio = 128;
  plen   = 1;
  lrecur = 8;
  bumpd  = 65;
  idx    = - 1;
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
  case ENLARGE_BOTH:
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
    result = compute(compute(data, EXTEND_X), EXTEND_Y);
    break;
  case SHARPEN_X:
    result = compute(data.transpose(), SHARPEN_Y).transpose();
    break;
  case ENLARGE_X:
    result = compute(data.transpose(), ENLARGE_Y).transpose();
    break;
  case INTEG_X:
    result = compute(data.transpose(), INTEG_Y).transpose();
    break;
  case RINTEG_X:
    result = compute(data.transpose(), RINTEG_Y).transpose();
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
  case ENLARGE_Y:
    initDop(data.rows());
    result = Eop[idx] * data;
    break;
  case INTEG_Y:
    initDop(data.rows());
    result = Iop[idx] * data;
    break;
  case RINTEG_Y:
    initDop(data.rows());
    result = RIop[idx] * data;
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
      Mat zscore(result);
      initDop(bumpd);
      Vec camera(2);
      camera[0] = T(0);
      camera[1] = T(1);
      assert(T(0) < dratio);
      const auto rxy(sqrt(T(data.rows()) * T(data.cols())));
      const auto Dop0((Dop[idx] * Dop[idx]).row(Dop[idx].rows() / 2));
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
      for(int zi = 0; zi <= dratio; zi ++) {
        Mat A(data.rows(), data.cols());
        for(int i = 0; i < A.rows(); i ++)
          for(int j = 0; j < A.cols(); j ++)
            A(i, j) = T(0);
        for(int j = 0; j < Dop0.size(); j ++) {
          Vec cpoint(2);
          cpoint[0] = (T(j) - T(Dop0.size() - 1) / T(2)) / rxy;
          cpoint[1] = T(zi) / T(dratio) * T(dist);
          // x-z plane projection of point p with camera geometry c to z=0.
          // c := camera, p := cpoint.
          // <c + (p - c) * t, [0, 1]> = 0
          const auto t(- camera[1] / (cpoint[1] - camera[1]));
          const auto y0((camera + (cpoint - camera) * t)[0] * rxy);
          //  N.B. d^2C_k/dy^2 on zi.
          for(int i = 0; i < A.rows(); i ++)
            A.row(i) += data.row(getImgPt(i + int(y0), data.rows())) * Dop0[j];
        }
#if defined(_OPENMP)
#pragma omp critical
#endif
        {
          for(int i = 0; i < A.rows(); i ++)
            for(int j = 0; j < A.cols(); j ++)
              if(zi == 0 || abs(A(i, j)) <= zscore(i, j)) {
                result(i, j) = T(1) - T(zi) / T(dratio);
                zscore(i, j) = abs(A(i, j));
              }
        }
      }
    }
    result = gmean(compute(result, INTEG_Y), compute(result, RINTEG_Y));
    result = gmean(compute(result, INTEG_Y), compute(result, RINTEG_Y));
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
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
      for(int j = 0; j < plen; j ++) {
        P0<T> p(data.rows() / (j + 1));
        Mat fdata(data.rows() / (j + 1), data.cols());
        Mat rdata(data.rows() / (j + 1), data.cols());
        for(int i = 0; i < fdata.rows(); i ++) {
          for(int k = 0; k < fdata.cols(); k ++)
            fdata(fdata.rows() - i - 1, k) = rdata(rdata.rows() - i - 1, k) = T(0);
          for(int k = 0; k < j + 1; k ++) {
            fdata.row(fdata.rows() - i - 1) += data.row(data.rows() - 1 - i * (j + 1) - k);
            rdata.row(rdata.rows() - i - 1) += data.row(i * (j + 1) + k);
          }
        }
        for(int i = 0; i < data.cols(); i ++) {
          result(data.rows() + j + plen, i) = p.next(fdata.col(i), true);
          result(plen - j - 1, i) = p.next(rdata.col(i), true);
        }
        for(int i = 0; i < j; i ++) {
          result.row(data.rows() + j + plen) -= result.row(data.rows() + i + plen);
          result.row(plen - j - 1) -= result.row(plen - i - 1);
        }
      }
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
  Iop.push_back(Dop[idx]);
  RIop.push_back(Dop[idx]);
  Eop.push_back(Mat(Dop[idx].rows() * 2, Dop[idx].cols()));
  Sop.push_back(Dop[idx]);
  int cnt(1);
  for(int ss = 2; ss <= size; ss *= 2, cnt ++) {
          auto DFTD(seed(ss, false));
          auto DFTL(seed(ss * 2, false));
    const auto IDFT(seed(ss, true));
          auto DFTI(DFTD);
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
        DFTI.row(i) /= phase;
        nd += abs(phase)  * abs(phase);
        ni += abs(phase2) * abs(phase2);
      } else {
        DFTD.row(i) *= U(T(0));
        DFTI.row(i) *= U(T(0));
      }
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
    DFTI *= nd * ni == T(0) ? T(0) : sqrt(sqrt(T(DFTD.rows() - 1) / (nd * ni)));
    DFTE /= T(DFTE.rows()  - 1);
#if defined(_WITHOUT_EIGEN_)
    const Mat lDop((IDFT * DFTD).template real<T>());
    const Mat lIop((IDFT * DFTI).template real<T>());
    const Mat EE((IDFT * DFTE).template real<T>());
    const Mat LL((seed(ss * 2, true) * DFTL).template real<T>());
#else
    const Mat lDop((IDFT * DFTD).real().template cast<T>());
    const Mat lIop((IDFT * DFTI).real().template cast<T>());
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
        Iop[idx](i, ij) += lIop(jj, j);
        Sop[idx](i, ij) += lLop(jj, j);
      }
  }
  Dop[idx] /= T(cnt);
  Iop[idx] /= T(cnt);
  for(int i = 0; i < Sop[idx].rows(); i ++) {
    Eop[idx].row(i * 2 + 0) = Sop[idx].row(i);
    for(int j = 0; j < Sop[idx].cols(); j ++)
      Eop[idx](i * 2 + 1, j) = Sop[idx](Sop[idx].rows() - i - 1, Sop[idx].cols() - j - 1);
  }
  const Mat SS(Sop[idx]);
  for(int i = 0; i < SS.rows(); i ++)
    for(int j = 0; j < SS.cols(); j ++)
     Sop[idx](SS.rows() - i - 1, SS.cols() - j - 1) += SS(i, j);
  const Mat II(Iop[idx]);
  for(int i = 0; i < II.rows(); i ++)
    for(int j = 0; j < II.cols(); j ++)
      Iop[idx](II.rows() - i - 1, II.cols() - j - 1) += II(i, j);
  Iop[idx] /= T(2);
  for(int i = 0; i < Iop[idx].rows(); i ++)
    for(int j = 0; j < Iop[idx].cols(); j ++)
      RIop[idx](i, j) = Iop[idx](Iop[idx].rows() - i - 1, Iop[idx].cols() - j - 1);
  T mnorm(0);
  for(int i = 0; i < Sop[idx].rows(); i ++)
    mnorm = max(mnorm, Sop[idx].row(i).dot(Sop[idx].row(i)));
  Sop[idx] /= sqrt(mnorm);
  for(int i = 0; i < lrecur; i ++)
    Sop[idx] = Sop[idx] * Sop[idx];
  Sop[idx] *= sqrt(mnorm);
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

