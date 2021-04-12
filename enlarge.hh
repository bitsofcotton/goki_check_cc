/* BSD 3-Clause License:
 * Copyright (c) 2018-2021, bitsofcotton.
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

template <typename T> class reDig;

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
    FLARGE_X,
    FLARGE_Y,
    FLARGE_BOTH,
    DETECT_X,
    DETECT_Y,
    DETECT_BOTH,
    INTEG_X,
    INTEG_Y,
    INTEG_BOTH,
    COLLECT_X,
    COLLECT_Y,
    COLLECT_BOTH,
    BUMP_X,
    BUMP_Y,
    BUMP_BOTH,
    EXTEND_X,
    EXTEND_Y,
    EXTEND_BOTH,
    EXTEND2_X,
    EXTEND2_Y,
    EXTEND2_BOTH,
    BLINK_X,
    BLINK_Y,
    BLINK_BOTH,
    REPRESENT,
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
  Filter(const int& recur = 2);
  ~Filter();
  Mat compute(const Mat& data, const direction_t& dir, const int& n = 0);
  inline int getImgPt(const int& y, const int& h);
  inline Mat gmean(const Mat& a, const Mat& b);

private:
  T   Pi;
  int recur;
};

template <typename T> Filter<T>::Filter(const int& recur) {
  Pi  = atan2(T(1), T(1)) * T(4);
  this->recur = recur;
}

template <typename T> Filter<T>::~Filter() {
  ;
}

template <typename T> typename Filter<T>::Mat Filter<T>::compute(const Mat& data, const direction_t& dir, const int& n) {
  assert(0 <= n);
  if(n <= 1 || dir == REPRESENT) {
    static P0<T> p;
    switch(dir) {
    case SHARPEN_BOTH:
      return compute(compute(data, SHARPEN_X), SHARPEN_Y);
    case ENLARGE_BOTH:
      return compute(compute(data, ENLARGE_X), ENLARGE_Y);
    case FLARGE_BOTH:
      return compute(compute(data, FLARGE_X), FLARGE_Y);
    case DETECT_BOTH:
      return (compute(data, DETECT_X) + compute(data, DETECT_Y)) / T(2);
    case INTEG_BOTH:
      return (compute(data, INTEG_X) + compute(data, INTEG_Y)) / T(2);
    case COLLECT_BOTH:
      return gmean(compute(data, COLLECT_X), compute(data, COLLECT_Y));
    case BUMP_BOTH:
      // eigen sum on curvature.
      return (compute(data, BUMP_X) + compute(data, BUMP_Y)) / T(2);
    case EXTEND_BOTH:
      return compute(compute(data, EXTEND_X), EXTEND_Y);
    case EXTEND2_BOTH:
      return compute(compute(data, EXTEND2_X), EXTEND2_Y);
    case BLINK_BOTH:
      return compute(compute(data, BLINK_X), BLINK_Y);
    case SHARPEN_X:
      return compute(data.transpose(), SHARPEN_Y).transpose();
    case ENLARGE_X:
      return compute(data.transpose(), ENLARGE_Y).transpose();
    case FLARGE_X:
      return compute(data.transpose(), FLARGE_Y).transpose();
    case DETECT_X:
      return compute(data.transpose(), DETECT_Y).transpose();
    case INTEG_X:
      return compute(data.transpose(), INTEG_Y).transpose();
    case COLLECT_X:
      return compute(data.transpose(), COLLECT_Y).transpose();
    case BUMP_X:
      return compute(data.transpose(), BUMP_Y).transpose();
    case EXTEND_X:
      return compute(data.transpose(), EXTEND_Y).transpose();
    case EXTEND2_X:
      return compute(data.transpose(), EXTEND2_Y).transpose();
    case BLINK_X:
      return compute(data.transpose(), BLINK_Y).transpose();
    case DETECT_Y:
      return p.diff(  data.rows()) * data;
    case INTEG_Y:
      return p.diff(- data.rows()) * data;
    case COLLECT_Y:
      return compute(compute(data, DETECT_Y), ABS);
    case SHARPEN_Y:
      {
        assert(2 <= data.rows());
        static vector<vector<Mat> > Sop;
        int idx(0);
        for(int i = 0; i < Sop.size(); i ++)
          if(Sop[i][0].rows() == data.rows()) {
            idx = i;
            goto sopi;
          }
        cerr << "s" << flush;
        idx = Sop.size();
        Sop.emplace_back(vector<Mat>());
        Sop[idx].emplace_back(Mat(data.rows(), data.rows()));
        {
          auto& sop(Sop[idx][0]);
          for(int i = 0; i < sop.rows(); i ++)
            for(int j = 0; j < sop.cols(); j ++)
              sop(i, j) = T(0);
          int cnt(1);
          for(int ss = 2; ss <= data.rows(); ss *= 2, cnt ++) {
            auto DFTS(p.seed(ss));
            DFTS.row(0) *= U(T(0));
            for(int i = 1; i < DFTS.rows(); i ++) {
              // N.B. d/dt((d^(t)/dy^(t)) f), differential-integral space tilt on f.
              // DFTH.row(i) *= log(phase);
              // N.B. please refer enlarge.wxm, half freq space refer and uses each.
              //   -> This is sharpen operation at all because this is same as original
              //      picture when {x0 + x0.5, x0.5 + x1, x1 + x1.5, x1.5 + x2, ...}
              //      series, and both picture of dft is same, them, pick {x0, x1, ...}.
              DFTS.row(i) /= exp(U(T(0), T(1)) * Pi * U(T(i)) / T(DFTS.rows())) - U(T(1));
            }
            DFTS /= T(DFTS.rows() - 1);
#if defined(_WITHOUT_EIGEN_)
            Mat lSop(- (p.seed(- ss) * DFTS).template real<T>());
#else
            Mat lSop(- (p.seed(- ss) * DFTS).template real<T>());
#endif
            for(int i = 0; i < lSop.rows(); i ++)
              lSop(i, i) += T(1);
            for(int i = 0; i < sop.rows(); i ++)
              for(int j = 0; j < lSop.rows(); j ++) {
                int ij(i - lSop.rows() / 2 + j);
                int jj(lSop.rows() / 2);
                if(i < lSop.rows() / 2) {
                  ij = j;
                  jj = i;
                } else if(sop.rows() - i < (lSop.rows() + 1) / 2) {
                  ij = j - lSop.rows() + sop.rows();
                  jj = i - sop.rows() + lSop.rows();
                }
                sop(i, ij) += lSop(jj, j);
              }
          }
          sop /= T(cnt);
          const Mat SS(sop);
          for(int i = 0; i < SS.rows(); i ++)
            for(int j = 0; j < SS.cols(); j ++)
              sop(SS.rows() - i - 1, SS.cols() - j - 1) += SS(i, j);
          sop /= T(2);
        }
       sopi:
        // N.B. insufficient:
        for( ; Sop[idx].size() <= recur; )
          Sop[idx].emplace_back(Sop[idx][Sop[idx].size() - 1] * Sop[idx][0]);
        return Sop[idx][recur] * data;
      }
      break;
    case ENLARGE_Y:
      {
        assert(2 <= data.rows());
        static vector<vector<Mat> > Eop;
        const auto& size(data.rows());
        if(Eop.size() <= size)
          Eop.resize(size + 1, vector<Mat>());
        else if(recur < Eop[size].size())
          goto eopi;
        if(Eop[size].size() <= recur)
          Eop[size].resize(recur + 1, Mat());
        {
          auto& eop(Eop[size][recur]);
          if(eop.cols() == size)
            goto eopi;
          cerr << "e" << flush;
          eop.resize((size - 1) * recur + 1, size);
          for(int j = 0; j < eop.rows(); j ++)
            eop.row(j) = p.taylor(eop.cols(), T(j) / T(eop.rows() - 1) * T(eop.cols() - 1));
        }
       eopi:
        return Eop[size][recur] * data;
      }
      break;
    case FLARGE_Y:
      {
        Mat work(data);
        for(int i = 0; i < work.rows(); i ++)
          for(int j = 0; j < work.cols(); j ++)
            work(i, j) += T(1) / T(256);
        Decompose<T> e(work.rows());
        Mat result(work.rows() * recur, work.cols());
        for(int i = 0; i < work.cols(); i ++)
          result.setCol(i, e.enlarge(work.col(i), recur));
        return result;
      }
      break;
    case BUMP_Y:
      {
        Mat result(data.rows(), data.cols());
        Mat zscore(data.rows(), data.cols());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
        for(int i = 0; i < result.rows(); i ++)
          for(int j = 0; j < result.cols(); j ++) {
            result(i, j) = T(0);
            zscore(i, j) = - T(1);
          }
        const auto rxy(sqrt(T(data.rows()) * T(data.cols())));
        const int  dratio(sqrt(sqrt(rxy)));
              Vec  camera(2);
              Vec  cpoint(2);
        camera[0] = T(0);
        camera[1] = T(1);
        cpoint[0] = T(1) / T(2 * dratio);
        for(int zi = 0; zi < dratio; zi ++) {
          Mat A(data.rows(), data.cols());
          for(int i = 0; i < A.rows(); i ++)
            for(int j = 0; j < A.cols(); j ++)
              A(i, j) = T(0);
          // N.B. projection scale is linear.
          cpoint[1] = T(zi) / T(dratio);
          // x-z plane projection of point p with camera geometry c to z=0.
          // c := camera, p := cpoint.
          // <c + (p - c) * t, [0, 1]> = 0
          const auto t(- camera[1] / (cpoint[1] - camera[1]));
          const auto y0((camera + (cpoint - camera) * t)[0] * rxy);
          if(abs(int(y0 * T(2))) < 3 || data.rows() / 2 < int(y0 * T(2)))
            continue;
          assert(int(y0) * 2 <= data.rows());
          const auto& Dop(p.diff(abs(int(y0 * T(2))) & ~ int(1)));
          const auto& Dop0(Dop.row(Dop.rows() / 2));
          //  N.B. dC_k/dy on zi.
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
          for(int i = 0; i < A.rows(); i ++) {
            for(int j = 0; j < Dop0.size(); j ++)
              A.row(i) += data.row(getImgPt(i + j - Dop0.size() / 2, data.rows())) * Dop0[j];
          }
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
          for(int i = 0; i < A.rows(); i ++) {
            for(int j = 0; j < A.cols(); j ++)
              if(zscore(i, j) < abs(A(i, j))) {
                result(i, j) = T(zi + 1);
                zscore(i, j) = abs(A(i, j));
              }
          }
        }
        return result;
      }
      break;
    case EXTEND_Y:
      {
        MatU result(data.rows() + 2 * recur, data.cols());
        auto DFT(data.template cast<complex<T> >() * p.seed(data.cols()));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
        for(int i = 0; i < data.rows(); i ++) {
          result.row(i + recur) = std::move(DFT.row(i));
        }
        for(int i = 0; i < recur; i ++) {
          const auto& next(p.next(int(data.rows()) / (i + 1)));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
          for(int j = 0; j < data.cols(); j ++) {
            result(data.rows() + recur + i, j) =
              result(recur - i - 1, j) = complex<T>(T(0));
            for(int k = 0; k < next.size(); k ++) {
              result(data.rows() + recur + i, j) +=
                result(recur + (k - next.size() + 1) * (i + 1) + data.rows() - 1, j) * next[k];
              result(recur - i - 1, j) +=
                result(recur - (k - next.size() + 1) * (i + 1), j) * next[k];
            }
          }
        }
        return (result * p.seed(- result.cols())).template real<T>();
      }
      break;
    case EXTEND2_Y:
      {
        Mat result(data.rows() + 2 * recur, data.cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
        for(int i = 0; i < data.rows(); i ++)
          result.row(recur + i) = data.row(i);
        for(int i = 0; i < recur; i ++)
          for(int j = 0; j < data.cols(); j ++) {
            SimpleVector<T> bnext(data.rows() / (i + 1));
            SimpleVector<T> bback(bnext.size());
            for(int k = 0; k < bnext.size(); k ++) {
              bnext[k] = data(  (k - bnext.size() + 1) * (i + 1) + data.rows() - 1, j);
              bback[k] = data(- (k - bback.size() + 1) * (i + 1), j);
            }
            const auto varlen(min(bnext.size(), bnext.size() / 3));
            const auto nbn(sqrt(bnext.dot(bnext)));
            const auto nbb(sqrt(bback.dot(bback)));
            const auto next(invariantP1(bnext / nbn, varlen, T(n)));
            const auto back(invariantP1(bback / nbb, varlen, T(n)));
                  auto bbnext(next);
                  auto bbback(back);
            for(int k = 1; k < varlen; k ++) {
              bbnext[k - 1] = bnext[k - varlen + bnext.size()] / nbn;
              bbback[k - 1] = bback[k - varlen + bnext.size()] / nbb;
            }
            bbnext[varlen - 1] = bbnext[varlen - 2];
            bbback[varlen - 1] = bbback[varlen - 2];
            bbnext[varlen + 1] = bbnext[varlen] =
              bbback[varlen + 1] = bbback[varlen] =
                T(1) / sqrt(T((bnext.size() - varlen + 1) * 2) * T(varlen + 2));
            result(data.rows() + recur + i, j) = (next.dot(bbnext) - next[varlen - 1] * bbnext[varlen - 1]) / next[varlen - 1] * nbn;
            result(recur - i - 1, j) = (back.dot(bbback) - back[varlen - 1] * bbback[varlen - 1]) / back[varlen - 1] * nbb;
          }
        return result;
      }
      break;
    case BLINK_Y:
      {
        auto diff(p.seed(data.rows()) * data.template cast<complex<T> >());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
        for(int i = 1; i < data.rows(); i ++) {
          diff.row(i) *= - complex<T>(T(0), T(2)) * T(i) / T(data.rows());
        }
        diff = p.seed(- data.rows()) * diff;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
        for(int i = 1; i < recur - 1; i ++) {
          diff.row(i) += (diff.row(i - 1) + diff.row(i + 1)) * complex<T>(T(recur) / T(256));
        }
        diff = p.seed(data.rows()) * diff;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
        for(int i = 1; i < data.rows(); i ++) {
          diff.row(i) /= - complex<T>(T(0), T(2)) * T(i) / T(data.rows());
        }
        return (p.seed(- data.rows()) * diff).template real<T>();
      }
      break;
    case REPRESENT:
      return Decompose<T>(recur).decompose(data, n);
    case CLIP:
      {
        Mat result(data.rows(), data.cols());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
        for(int i = 0; i < result.rows(); i ++)
          for(int j = 0; j < result.cols(); j ++)
            result(i, j) = max(T(0), min(T(1), data(i, j)));
        return result;
      }
      break;
    case ABS:
      {
        Mat result(data.rows(), data.cols());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
        for(int i = 0; i < result.rows(); i ++)
          for(int j = 0; j < result.cols(); j ++)
            result(i, j) = abs(data(i, j));
        return result;
      }
      break;
    }
    assert(0 && "unknown command in Filter (should not be reached.)");
    return Mat();
  }
  if(dir == EXTEND_Y || dir == EXTEND_X || dir == EXTEND_BOTH || dir == EXTEND2_Y || dir == EXTEND2_X || dir == EXTEND2_BOTH || dir == REPRESENT)
    return compute(data, dir, 0);
  vector<Mat> res;
  res.reserve(n);
  res.emplace_back(compute(data, dir));
  for(int i = 1; i < n; i ++) {
    cerr << "r" << flush;
          Mat  lres(res[0]);
    const auto theta(T(i) * Pi / T(2 * n));
    const auto c(cos(theta));
    const auto s(sin(theta));
    Mat work(abs(int(c * T(data.rows()) + s * T(data.cols()))),
             abs(int(s * T(data.rows()) + c * T(data.cols()))));
    for(int j = 0; j < work.rows(); j ++)
      for(int k = 0; k < work.cols(); k ++)
        work(j, k) = T(0);
    for(int j = - (work.rows() + work.cols());
            j <   (work.rows() + work.cols()) * 2; j ++)
      for(int k = - (work.rows() + work.cols());
              k <   (work.rows() + work.cols()) * 2; k ++) {
        const int yy(c * T(j) - s * T(k) + T(1) / T(2));
        const int xx(s * T(j) + c * T(k) + T(1) / T(2));
        if(0 <= yy && yy < work.rows() &&
           0 <= xx && xx < work.cols()) {
          const auto dyy(((j % data.rows()) + data.rows()) % data.rows());
          const auto dxx(((k % data.cols()) + data.cols()) % data.cols());
          work(yy, xx) = work(min(yy + 1, int(work.rows()) - 1), xx) =
            work(yy, min(xx + 1, int(work.cols()) - 1)) =
            work(min(yy + 1, int(work.rows()) - 1),
                 min(xx + 1, int(work.cols()) - 1)) =
              data(dyy, dxx);
        }
      }
    work = compute(work, dir);
    // XXX inefficient:
    for(int j = - (work.rows() + work.cols());
            j <   (work.rows() + work.cols()) * 2; j ++)
      for(int k = - (work.rows() + work.cols());
              k <   (work.rows() + work.cols()) * 2; k ++) {
        const int yy(c * T(j) - s * T(k) + T(1) / T(2));
        const int xx(s * T(j) + c * T(k) + T(1) / T(2));
        if(0 <= yy && yy < work.rows() &&
           0 <= xx && xx < work.cols()) {
          lres(((j % lres.rows()) + lres.rows()) % lres.rows(),
               ((k % lres.cols()) + lres.cols()) % lres.cols()) = work(yy, xx);
        }
      }
    res.emplace_back(lres);
  }
  for(int i = 1; i < res.size(); i ++)
    res[0] += res[i];
  return res[0] /= T(res.size());
}

template <typename T> inline int Filter<T>::getImgPt(const int& y, const int& h) {
  int yy(y % (2 * h));
  if(yy < 0)
    yy = - yy;
  if(yy >= h)
    yy = h - (yy - h);
  return yy % h;
}

template <typename T> inline typename Filter<T>::Mat Filter<T>::gmean(const Mat& a, const Mat& b) {
  assert(a.rows() == b.rows() && a.cols() == b.cols());
  Mat res(a.rows(), a.cols());
  for(int i = 0; i < a.rows(); i ++)
    for(int j = 0; j < a.cols(); j ++) {
      const auto lval(a(i, j) * b(i, j));
      res(i, j) = (abs(a(i, j)) < abs(b(i, j)) ? (b(i, j) < T(0) ? - T(1) : T(1)) : (a(i, j) < T(0) ? - T(1) : T(1))) * sqrt(abs(lval));
    }
  return res;
}

#define _ENLARGE2X_
#endif

