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

#if !defined(_FILTER_)

using std::cerr;
using std::flush;
using std::sort;
using std::vector;
using std::max;
using std::min;
using std::abs;

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
  BLINK_X,
  BLINK_Y,
  BLINK_BOTH,
  REPRESENT,
  CLIP,
  ABS } direction_t;

template <typename T> static inline T getImgPt(const T& y, const T& h) {
  auto yy(y % (2 * h));
  if(yy < 0)
    yy = - yy;
  if(yy >= h)
    yy = h - (yy - h);
  return yy % h;
}

template <typename T> SimpleMatrix<T> rotate(const SimpleMatrix<T>& d, const T& theta) {
  assert(abs(theta) < atan(T(1)));
  const auto c(cos(theta));
  const auto s(sin(theta));
  const auto h0(abs(int(c * T(d.rows()) - s * T(d.cols()))));
  const auto h1(h0 + abs(int(s * T(d.cols()))) * 2);
  const auto w0(abs(int(s * T(d.rows()) + c * T(d.cols()))));
  const auto w1(w0 + abs(int(s * T(d.rows()))) * 2);
  SimpleMatrix<T> res(h0 < d.rows() ? h1 : h0,
                      w0 < d.cols() ? w1 : w0);
  const T offy(h0 < d.rows() ? abs(int(s * T(d.cols()))) : 0);
  const T offx(w0 < d.cols() ? abs(int(s * T(d.rows()))) : 0);
  res.O();
  const auto diag(int(sqrt(res.rows() * res.rows() +
                           res.cols() * res.cols())) + 1);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int j = - diag; j < diag; j ++)
    for(int k = - diag; k < diag; k ++) {
      const int yy(c * T(j) - s * T(k) + offy);
      const int xx(s * T(j) + c * T(k) + offx);
      if(0 <= yy && yy < res.rows() &&
         0 <= xx && xx < res.cols()) {
        const auto dyy(getImgPt<int>(j, d.rows()));
        const auto dxx(getImgPt<int>(k, d.cols()));
#if defined(_OPENMP)
#pragma omp critical
#endif
        {
          res(yy, xx) = res(min(yy + 1, int(res.rows()) - 1), xx) =
            res(yy, min(xx + 1, int(res.cols()) - 1)) =
            res(min(yy + 1, int(res.rows()) - 1),
                min(xx + 1, int(res.cols()) - 1)) =
              d(dyy, dxx);
        }
      }
    }
  return res;
}

template <typename T> static inline SimpleMatrix<T> center(const SimpleMatrix<T>& dr, const SimpleMatrix<T>& d) {
  SimpleMatrix<T> res(d.rows(), d.cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < res.rows(); i ++)
    for(int j = 0; j < res.cols(); j ++)
      res(i, j) = dr(max(0, min(i + (dr.rows() - d.rows()) / 2, dr.rows() - 1)),
                     max(0, min(j + (dr.cols() - d.cols()) / 2, dr.cols() - 1)));
  return res;
}

template <typename T> SimpleMatrix<T> sharpen(const int& size) {
  assert(0 < size);
  SimpleMatrix<T> s;
  const auto file(std::string("./.cache/lieonn/sharpen-") + std::to_string(size) +
#if defined(_FLOAT_BITS_)
    std:string("-") + std::to_string(_FLOAT_BITS_)
#else
    std::string("-ld")
#endif
  );
  ifstream cache(file.c_str());
  if(cache.is_open()) {
    cache >> s;
    cache.close();
  } else {
    if(2 < size) {
      const auto s0(sharpen<T>(size - 1) * T(size - 1));
      s = SimpleMatrix<T>(size, size).O().setMatrix(0, 0, s0);
      s.setMatrix(1, 1, s.subMatrix(1, 1, size - 1, size - 1) + s0);
      s.row(0) *= T(2);
      s.row(s.rows() - 1) *= T(2);
      s /= T(2);
    } else
      s  = SimpleMatrix<T>(size, size).O();
    auto dfts(dft<T>(size));
    dfts.row(0) *= complex<T>(T(0));
    static const auto Pi(atan(T(1)) * T(4));
    for(int i = 1; i < dfts.rows(); i ++) {
      // N.B. d/dt((d^(t)/dy^(t)) f), differential-integral space tilt on f.
      // DFTH.row(i) *= log(phase);
      // N.B. please refer enlarge.wxm, half freq space refer and uses each.
      //   -> This is sharpen operation at all because this is same as original
      //      picture when {x0 + x0.5, x0.5 + x1, x1 + x1.5, x1.5 + x2, ...}
      //      series, and both picture of dft is same, them, pick {x0, x1, ...}.
      dfts.row(i) /= exp(complex<T>(T(0), T(1)) * Pi * complex<T>(T(i)) / T(dfts.rows())) - complex<T>(T(1));
    }
    dfts /= T(dfts.rows() - 1);
    s += (dft<T>(- size) * dfts).template real<T>();
    s /= T(size);
    ofstream ocache(file.c_str());
    ocache << s;
    ocache.close();
    cerr << "." << flush;
  }
  return s;
}

// N.B. this function is NOT thread safe.
template <typename T> SimpleMatrix<T> filter(const SimpleMatrix<T>& data, const direction_t& dir, const int& n = 0, const int& recur = 2) {
  assert(0 <= n && 0 < recur);
  static const auto Pi(atan2(T(1), T(1)) * T(4));
  if(n <= 1 || dir == REPRESENT || dir == EXTEND_Y || dir == EXTEND_X || dir == EXTEND_BOTH || dir == ABS) {
    switch(dir) {
    case SHARPEN_BOTH:
      return filter<T>(filter<T>(data, SHARPEN_X, n, recur), SHARPEN_Y, n, recur);
    case ENLARGE_BOTH:
      return filter<T>(filter<T>(data, ENLARGE_X, n, recur), ENLARGE_Y, n, recur);
    case FLARGE_BOTH:
      return filter<T>(filter<T>(data, FLARGE_X, n, recur), FLARGE_Y, n, recur);
    case DETECT_BOTH:
      return (filter<T>(data, DETECT_X, n, recur) + filter<T>(data, DETECT_Y, n, recur)) / T(2);
    case INTEG_BOTH:
      return (filter<T>(data, INTEG_X, n, recur) + filter<T>(data, INTEG_Y, n, recur)) / T(2);
    case COLLECT_BOTH:
      return (filter<T>(data, COLLECT_X, n, recur) + filter<T>(data, COLLECT_Y, n, recur)) / T(2);
    case BUMP_BOTH:
      // eigen sum on curvature.
      return (filter<T>(data, BUMP_X, n, recur) + filter<T>(data, BUMP_Y, n, recur)) / T(2);
    case EXTEND_BOTH:
      return filter<T>(filter<T>(filter<T>(filter<T>(data, EXTEND_X, n, recur), CLIP, n, recur), EXTEND_Y, n, recur), CLIP, n, recur);
    case BLINK_BOTH:
      return filter<T>(filter<T>(data, BLINK_X, n, recur), BLINK_Y, n, recur);
    case SHARPEN_X:
      return filter<T>(data.transpose(), SHARPEN_Y, n, recur).transpose();
    case ENLARGE_X:
      return filter<T>(data.transpose(), ENLARGE_Y, n, recur).transpose();
    case FLARGE_X:
      return filter<T>(data.transpose(), FLARGE_Y, n, recur).transpose();
    case DETECT_X:
      return filter<T>(data.transpose(), DETECT_Y, n, recur).transpose();
    case INTEG_X:
      return filter<T>(data.transpose(), INTEG_Y, n, recur).transpose();
    case COLLECT_X:
      return filter<T>(data.transpose(), COLLECT_Y, n, recur).transpose();
    case BUMP_X:
      return filter<T>(data.transpose(), BUMP_Y, n, recur).transpose();
    case EXTEND_X:
      return filter<T>(data.transpose(), EXTEND_Y, n, recur).transpose();
    case BLINK_X:
      return filter<T>(data.transpose(), BLINK_Y, n, recur).transpose();
    case DETECT_Y:
      return diff<T>(  data.rows()) * data;
    case INTEG_Y:
      // return diff<T>(- data.rows()) * data;
      {
        // instead of integrate, we can use normalize:
        // N.B. d^exp(t)/dx^exp(t) f(x) == f(x + t dx), t != 0.
        //   so d^exp(- t)/dx^exp(- t)
        //   == d^(- exp(t))/dx^(- exp(t)) == f(x - t dx).
        //   we normalize with f(x - dx) + f(x + dx) in weak differential meaning.
        // Cor: d/dt d^t/d(x^t) f(x)
        //   == d/dt f(x + log(t) dx) == f'(x + log(t) dx) / g(x, t).
        // N.B. d^0/dx^0 f(x) == f(x)
        //   == d^exp(- inf)/dx^exp(- inf) f(x) == f(x - inf dx)
        //   == d^(- exp(- inf))/dx^(- exp(- inf dx))
        //   == d^(exp(inf))/dx^exp(inf) == f(x + inf dx)
        auto normalize(dft<T>(data.rows()) * data.template cast<complex<T> >());
        for(int i = 0; i < normalize.rows(); i ++) {
          const auto n(complex<T>(T(0), - T(2) * Pi * T(i) / T(normalize.rows())));
          normalize.row(i) *= exp(n) + exp(- n);
        }
        return (dft<T>(- data.rows()) * normalize).template real<T>();
      }
      break;
    case COLLECT_Y:
      return filter<T>(filter<T>(data, DETECT_Y, n, recur), ABS, n, recur);
    case SHARPEN_Y:
      {
        auto shp(sharpen<T>(int(data.rows())));
        for(int i = 0; i < recur; i ++)
          shp = shp * shp;
        return shp * data;
      }
      break;
    case ENLARGE_Y:
      {
        assert(2 <= data.rows());
        static vector<vector<SimpleMatrix<T> > > Eop;
        const auto& size(data.rows());
        if(Eop.size() <= size)
          Eop.resize(size + 1, vector<SimpleMatrix<T> >());
        else if(recur < Eop[size].size())
          goto eopi;
        if(Eop[size].size() <= recur)
          Eop[size].resize(recur + 1, SimpleMatrix<T>());
        {
          auto& eop(Eop[size][recur]);
          if(eop.cols() == size)
            goto eopi;
          cerr << "e" << flush;
          eop.resize((size - 1) * recur + 1, size);
          for(int j = 0; j < eop.rows(); j ++)
            eop.row(j) = taylor<T>(eop.cols(), T(j) / T(eop.rows() - 1) * T(eop.cols() - 1));
        }
       eopi:
        return Eop[size][recur] * data;
      }
      break;
    case FLARGE_Y:
      {
        SimpleMatrix<T> work(data);
        for(int i = 0; i < work.rows(); i ++)
          for(int j = 0; j < work.cols(); j ++)
            work(i, j) += T(1) / T(256);
        Decompose<T> e(work.rows());
        SimpleMatrix<T> result(work.rows() * recur, work.cols());
        for(int i = 0; i < work.cols(); i ++)
          result.setCol(i, e.enlarge(work.col(i), recur));
        return result;
      }
      break;
    case BUMP_Y:
      {
        SimpleMatrix<T> result(data.rows(), data.cols());
        SimpleMatrix<T> zscore(data.rows(), data.cols());
        result.O();
        zscore.I(- T(1));
        const auto rxy(T(min(data.rows(), data.cols())));
        const int  dratio(sqrt(sqrt(rxy)));
              SimpleVector<T> camera(2);
              SimpleVector<T> cpoint(2);
        camera[0] = T(0);
        camera[1] = T(1);
        cpoint[0] = T(1) / T(2 * dratio);
        for(int zi = 0; zi < dratio; zi ++) {
          SimpleMatrix<T> A(data.rows(), data.cols());
          A.O();
          // N.B. projection scale is linear.
          cpoint[1] = T(zi) / T(dratio);
          // x-z plane projection of point p with camera geometry c to z=0.
          // c := camera, p := cpoint.
          // <c + (p - c) * t, [0, 1]> = 0
          const auto t(- camera[1] / (cpoint[1] - camera[1]));
          const auto y0((camera + (cpoint - camera) * t)[0] * rxy);
          if(abs(int(y0)) < 3 || rxy < abs(y0) * T(2)) continue;
          const auto Dop(diff<T>(abs(int(y0)) & ~ int(1)));
          const auto Dop0(Dop.row(Dop.rows() / 2) + Dop.row(Dop.rows() / 2 + 1));
          //  N.B. dC_k/dy on zi.
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
          for(int i = 0; i < A.rows(); i ++) {
            for(int j = 0; j < Dop0.size(); j ++)
              A.row(i) += data.row(getImgPt<int>(i + j - Dop0.size() / 2, data.rows())) * Dop0[j];
          }
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
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
        SimpleMatrix<T> result(data.rows() + 2 * recur, data.cols());
        SimpleMatrix<complex<T> > ddft(data.rows(), data.cols());
        const auto cdft(dft<T>(   data.cols()));
        const auto cidft(dft<T>(- data.cols()));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
        for(int i = 0; i < data.rows(); i ++) {
          ddft.row(i) = cdft * data.row(i).template cast<complex<T> >();
          result.row(i + recur) = data.row(i);
        }
        for(int i = 0; i < recur; i ++) {
          const auto next(nextP0<T>(int(ddft.rows()) / (i + 1)));
          SimpleVector<complex<T> > bwd(data.cols());
          for(int j = 0; j < bwd.size(); j ++)
            bwd[j] = complex<T>(T(0));
          auto fwd(bwd);
          for(int j = 0; j < next.size(); j ++) {
            bwd += ddft.row(- (j - next.size() + 1) * (i + 1)) * next[j];
            fwd += ddft.row(  (j - next.size() + 1) * (i + 1) + ddft.rows() - 1) * next[j];
          }
          result.row(recur - i - 1) = (cidft * bwd).template real<T>();
          result.row(data.rows() + recur + i) = (cidft * fwd).template real<T>();
        }
        return result;
      }
      break;
    case BLINK_Y:
      {
        auto dif(dft<T>(data.rows()) * data.template cast<complex<T> >());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
        for(int i = 1; i < data.rows(); i ++) {
          dif.row(i) *= - complex<T>(T(0), T(2)) * T(i) / T(data.rows());
        }
        dif = dft<T>(- data.rows()) * dif;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
        for(int i = 1; i < recur - 1; i ++) {
          dif.row(i) += (dif.row(i - 1) + dif.row(i + 1)) * complex<T>(T(recur) / T(256));
        }
        dif = dft<T>(data.rows()) * dif;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
        for(int i = 1; i < data.rows(); i ++) {
          dif.row(i) /= - complex<T>(T(0), T(2)) * T(i) / T(data.rows());
        }
        return (dft<T>(- data.rows()) * dif).template real<T>();
      }
      break;
    case REPRESENT:
      return Decompose<T>(recur).represent(data, n);
    case CLIP:
      {
        SimpleMatrix<T> result(data.rows(), data.cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
        for(int i = 0; i < result.rows(); i ++)
          for(int j = 0; j < result.cols(); j ++)
            result(i, j) = max(T(0), min(T(1), data(i, j)));
        return result;
      }
      break;
    case ABS:
      {
        SimpleMatrix<T> result(data.rows(), data.cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
        for(int i = 0; i < result.rows(); i ++)
          for(int j = 0; j < result.cols(); j ++)
            result(i, j) = abs(data(i, j));
        return result;
      }
      break;
    }
    assert(0 && "unknown command in filter (should not be reached.)");
    return SimpleMatrix<T>();
  }
  auto res(filter<T>(data, dir, 0, recur));
  for(int i = 0; i < n; i ++) {
    cerr << "r" << flush;
    const auto theta((T(i) - T(n - 1) / T(2)) * atan(T(1)) / (T(n) / T(2)));
    res += center<T>(rotate<T>(filter<T>(rotate<T>(data, theta), dir, 0, recur), - theta), res);
  }
  return res /= T(n + 1);
}

#define _FILTER_
#endif

