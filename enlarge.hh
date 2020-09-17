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
  Filter(const int& recur = 2);
  ~Filter();
  Mat rotcompute(const Mat& data, const direction_t& dir, const int& n);
  Mat compute(const Mat& data, const direction_t& dir);
  inline int getImgPt(const int& y, const int& h);
  inline Mat gmean(const Mat& a, const Mat& b);

private:
  T   Pi;
  int recur;
};

template <typename T> Filter<T>::Filter(const int& recur) {
  Pi  = atan2(T(1), T(1)) * T(4);
  this->recur = recur;
  assert(0 < recur);
}

template <typename T> Filter<T>::~Filter() {
  ;
}

template <typename T> typename Filter<T>::Mat Filter<T>::rotcompute(const Mat& data, const direction_t& dir, const int& n) {
  if(dir == EXTEND_X)
    return rotcompute(data.transpose(), EXTEND_Y, n).transpose();
  if(dir == EXTEND_BOTH)
    return rotcompute(rotcompute(data, EXTEND_Y, n), EXTEND_X, n);
  vector<Mat> res;
  res.emplace_back(compute(data, dir));
  for(int i = 1; i < n; i ++) {
    cerr << "r" << flush;
          auto lres(res[0]);
    const auto theta(T(i) * Pi / T(2 * n));
    if(dir == EXTEND_Y) {
      const auto c(cos(theta * T(2)));
      const auto s(sin(theta * T(2)));
      const auto rows(abs(int(c * T(data.rows()) + s * T(data.cols()))));
      const auto cols(abs(int(s * T(data.rows()) + c * T(data.cols()))));
      if(rows / (recur + 1) - 1 < 8 || cols < 8) continue;
      Mat workt(rows, cols);
      Mat workb(rows, cols);
      for(int j = - (workt.rows() + workt.cols());
              j <   (workt.rows() + workt.cols()) * 2; j ++)
        for(int k = - (workt.rows() + workt.cols());
                k <   (workt.rows() + workt.cols()) * 2; k ++) {
          const int yy(c * T(j));
          const int xx(s * T(j) + c * T(k) + T(1) / T(2));
          if(0 <= yy && yy < workt.rows() &&
             0 <= xx && xx < workt.cols()) {
            const auto dyy(((j % data.rows()) + data.rows()) % data.rows());
            const auto dxx(((k % data.cols()) + data.cols()) % data.cols());
            workt(yy, xx) =
              workt(min(yy + 1, int(workt.rows()) - 1), xx) =
              workt(yy, min(xx + 1, int(workt.cols()) - 1)) =
              workt(min(yy + 1, int(workt.rows()) - 1),
                     min(xx + 1, int(workt.cols()) - 1)) =
                data(dyy, dxx);
            workb(yy, xx) =
              workb(min(yy + 1, int(workb.rows()) - 1), xx) =
              workb(yy, min(xx + 1, int(workb.cols()) - 1)) =
              workb(min(yy + 1, int(workb.rows()) - 1),
                     min(xx + 1, int(workb.cols()) - 1)) =
                data(data.rows() - dyy - 1, data.cols() - dxx - 1);
          }
        }
      workt = compute(workt, dir);
      workb = compute(workb, dir);
      for(int i = 0; i < data.rows(); i ++)
        for(int j = 0; j < res[0].cols(); j ++)
          for(int k = 1; k <= recur; k ++) {
            const int bx(j - s * T(k - 1));
            const int xx(j - s * T(k));
            if(0 <= xx && bx < workt.cols()) {
              const auto yy(recur - k);
              const auto yyy(k - recur - 1 + res[0].rows());
              if(0 <= yy && yy + 1 <= workt.rows()) {
                const auto lst(abs(workt(yy, xx) - workt(yy + 1, bx)));
                const auto lsb(abs(workb(yy, xx) - workb(yy + 1, bx)));
                if(0 <= yy  && yy  < lres.rows())
                  lres(yy,  j) = workt(yy, xx);
                if(0 <= yyy && yyy < lres.rows())
                  lres(yyy, j) = workb(yy, xx);
              }
            }
          }
    } else {
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
          const int jj(((j % lres.rows()) + lres.rows()) % lres.rows());
          const int kk(((k % lres.cols()) + lres.cols()) % lres.cols());
          if(0 <= yy && yy < work.rows() &&
             0 <= xx && xx < work.cols()) {
            lres(((j % lres.rows()) + lres.rows()) % lres.rows(),
                 ((k % lres.cols()) + lres.cols()) % lres.cols()) = work(yy, xx);
          }
        }
    }
    res.emplace_back(lres);
  }
  Mat rres(res[0].rows(), res[0].cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < rres.rows(); i ++) {
    for(int j = 0; j < rres.cols(); j ++) {
      vector<T> sr;
      sr.reserve(res.size());
      for(int k = 0; k < res.size(); k ++)
        sr.emplace_back(res[k](i, j));
      sort(sr.begin(), sr.end());
      rres(i, j) = sr[sr.size() / 2];
    }
  }
  return rres;
}

template <typename T> typename Filter<T>::Mat Filter<T>::compute(const Mat& data, const direction_t& dir) {
  static P0<T> p;
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
      result = Sop[idx][recur] * data;
    }
    break;
  case ENLARGE_Y:
    {
      assert(2 <= data.rows());
      static vector<pair<vector<Mat>, int> > Eop;
      int idx(- 1);
      for(int i = 0; i < Eop.size(); i ++)
        if(Eop[i].second == data.rows()) {
          idx = i;
          break;
        }
      if(idx < 0) {
        idx = Eop.size();
        Eop.emplace_back(make_pair(vector<Mat>(), data.rows()));
      }
      cerr << "e" << flush;
      if(Eop[idx].first.size() <= recur)
        Eop[idx].first.resize(recur + 1, Mat());
      auto& eop(Eop[idx].first[recur]);
      if(eop.cols() == data.rows())
        goto eopi;
      eop.resize((data.rows() - 1) * recur + 1, data.rows());
      for(int i = 0; i < eop.rows(); i ++)
        eop.row(i) = p.taylor(eop.cols(), T(i) / T(recur));
      {
        const Mat EEop(eop);
        for(int i = 0; i < EEop.rows(); i ++)
          for(int j = 0; j < EEop.cols(); j ++)
            eop(i, j) += EEop(EEop.rows() - i - 1,
                              EEop.cols() - j - 1);
        eop /= T(2);
      }
     eopi:
      result = Eop[idx].first[recur] * data;
    }
    break;
  case DETECT_Y:
    result = p.diff(data.rows()) * data;
    break;
  case COLLECT_Y:
    result = compute(compute(data, DETECT_Y), ABS);
    break;
  case BUMP_Y:
    {
      result = Mat(data.rows(), data.cols());
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
      const int  dratio(sqrt(rxy));
            Vec  camera(2);
            Vec  cpoint(2);
      camera[0] = T(0);
      camera[1] = T(1);
      cpoint[0] = T(1) / T(2 * dratio);
      cerr << dratio << "depth ";
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
        const auto& Dop(p.diff(abs(int(y0 * T(2)))));
        const auto& Dop0(Dop.row(Dop.rows() / 2));
        //const auto  Dop0((Dop * Dop).row(Dop.rows() / 2));
        //  N.B. d^2C_k/dy^2 on zi.
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
            if(zscore(i, j) < T(0) || abs(A(i, j)) <= zscore(i, j)) {
              result(i, j) = T(zi) / T(dratio);
              zscore(i, j) = abs(A(i, j));
            }
        }
      }
    }
    break;
  case EXTEND_Y:
    {
      result = Mat(data.rows() + 2 * recur, data.cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int i = 0; i < data.rows(); i ++)
        result.row(i + recur) = data.row(i);
      for(int i = 0; i < recur; i ++) {
        const auto  size(min(120, int(data.rows()) / (i + 1) - 1));
        const auto& comp(p.next(size - 1));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
        for(int j = 0; j < data.cols(); j ++) {
          result(data.rows() + recur + i, j) =
            result(recur - i - 1, j) = T(0);
          for(int k = 1; k < size; k ++) {
            result(data.rows() + recur + i, j) +=
              result(data.rows() + recur - 1 + (k - size) * (i + 1), j) *
                comp[k - 1];
            result(recur - i - 1, j) +=
              result(recur + (size - k) * (i + 1), j) * comp[k - 1];
          }
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

