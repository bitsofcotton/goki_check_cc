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
  BLUR_X,
  BLUR_Y,
  BLUR_BOTH,
  INTEG_BOTH,
  DETECT_BOTH,
  COLLECT_BOTH,
  BUMP_BOTH,
  EXTEND_X,
  EXTEND_Y,
  EXTEND_BOTH,
  BLINK_X,
  BLINK_Y,
  BLINK_BOTH,
  LPF_X,
  LPF_Y,
  LPF_BOTH,
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
  for(int j = - diag; j < diag; j ++)
    for(int k = - diag; k < diag; k ++) {
      const int yy(c * T(j) - s * T(k) + offy);
      const int xx(s * T(j) + c * T(k) + offx);
      if(0 <= yy && yy < res.rows() &&
         0 <= xx && xx < res.cols()) {
        const auto dyy(getImgPt<int>(j, d.rows()));
        const auto dxx(getImgPt<int>(k, d.cols()));
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
    std::string("-") + std::to_string(_FLOAT_BITS_)
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
    static const auto Pi(atan(T(1)) * T(4));
    dfts.row(0) *= complex<T>(T(0));
    for(int i = 1; i < dfts.rows(); i ++) {
      // N.B. d/dt((d^(t)/dy^(t)) f), differential-integral space tilt on f.
      // DFTH.row(i) *= log(phase);
      // N.B. please refer enlarge.wxm, half freq space refer and uses each.
      //   -> This is sharpen operation at all because this is same as original
      //      picture when {x0 + x0.5, x0.5 + x1, x1 + x1.5, x1.5 + x2, ...}
      //      series, and both picture of dft is same, them, pick {x0, x1, ...}.
      dfts.row(i) /= exp(complex<T>(T(0), Pi * T(i) / T(dfts.rows()))) - complex<T>(T(1));
    }
    s += (dft<T>(- size) * dfts).template real<T>() / T(size - 1);
    if(2 < size)
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
  assert(0 <= n && (dir == BLINK_Y || dir == BLINK_X || dir == BLINK_BOTH || 0 < recur));
  static const auto Pi(atan2(T(1), T(1)) * T(4));
  if(n <= 1 || dir == REPRESENT || dir == EXTEND_Y || dir == EXTEND_X || dir == EXTEND_BOTH || dir == CLIP || dir == ABS) {
    switch(dir) {
    case SHARPEN_BOTH:
      return filter<T>(filter<T>(data, SHARPEN_X, n, recur), SHARPEN_Y, n, recur);
    case ENLARGE_BOTH:
      return filter<T>(filter<T>(data, ENLARGE_X, n, recur), ENLARGE_Y, n, recur);
    case FLARGE_BOTH:
      return filter<T>(filter<T>(data, FLARGE_X, n, recur), FLARGE_Y, n, recur);
    case BLUR_BOTH:
      return filter<T>(filter<T>(data, BLUR_Y, n, recur), BLUR_X, n, recur);
    case EXTEND_BOTH:
      return filter<T>(filter<T>(filter<T>(filter<T>(data, EXTEND_X, n, recur), CLIP, n, recur), EXTEND_Y, n, recur), CLIP, n, recur);
    case BLINK_BOTH:
      return filter<T>(filter<T>(data, BLINK_X, n, recur), BLINK_Y, n, recur);
    case LPF_BOTH:
      return filter<T>(filter<T>(data, LPF_X, n, recur), LPF_Y, n, recur);
    case SHARPEN_X:
      return filter<T>(data.transpose(), SHARPEN_Y, n, recur).transpose();
    case ENLARGE_X:
      return filter<T>(data.transpose(), ENLARGE_Y, n, recur).transpose();
    case FLARGE_X:
      return filter<T>(data.transpose(), FLARGE_Y, n, recur).transpose();
    case BLUR_X:
      return filter<T>(data.transpose(), BLUR_Y, n, recur).transpose();
    case EXTEND_X:
      return filter<T>(data.transpose(), EXTEND_Y, n, recur).transpose();
    case BLINK_X:
      return filter<T>(data.transpose(), BLINK_Y, n, recur).transpose();
    case LPF_X:
      return filter<T>(data.transpose(), LPF_Y, n, recur).transpose();
    case BLUR_Y:
      {
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
    case SHARPEN_Y:
      {
        const auto      shp(sharpen<T>(int(data.rows())) * data);
        SimpleMatrix<T> res(data.rows() * 2, data.cols());
        for(int i = 0; i < data.rows(); i ++) {
          res.row(2 * i)     = data.row(i) - shp.row(i);
          res.row(2 * i + 1) = data.row(i) + shp.row(i);
        }
        return res;
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
    case DETECT_BOTH:
      {
        const auto zy(diffRecur<T>(data.rows()) * data);
        const auto zx(data * diffRecur<T>(data.cols()).transpose());
              auto res(data);
        res.O();
        for(int i = 0; i < res.rows(); i ++)
          for(int j = 0; j < res.cols(); j ++) {
            const auto E(T(1) + zy(i, j) * zy(i, j));
            const auto F(       zy(i, j) * zx(i, j));
            const auto G(T(1) + zx(i, j) * zx(i, j));
            res(i, j) = E * G - F * F;
          }
        return res;
      }
      break;
    case INTEG_BOTH:
      {
        const auto zx(diffRecur<T>(- data.rows()) * data);
        const auto zy(data * diffRecur<T>(- data.cols()).transpose());
              auto res(data);
        res.O();
        for(int i = 0; i < res.rows(); i ++)
          for(int j = 0; j < res.cols(); j ++) {
            const auto E(T(1) + zy(i, j) * zy(i, j));
            const auto F(       zy(i, j) * zx(i, j));
            const auto G(T(1) + zx(i, j) * zx(i, j));
            res(i, j) = E * G - F * F;
          }
        return res;
      }
      break;
    case COLLECT_BOTH:
      return filter<T>(filter<T>(data, DETECT_BOTH, n, recur), ABS, n, recur);
    case BUMP_BOTH:
      {
        SimpleMatrix<T> result(data.rows(), data.cols());
        SimpleMatrix<T> zscore(data.rows(), data.cols());
        result.O();
        zscore.O(- T(1));
        const auto rxy(T(min(data.rows(), data.cols())));
        const int  dratio(sqrt(sqrt(rxy)));
              SimpleVector<T> camera(2);
              SimpleVector<T> cpoint(2);
        camera[0] = T(0);
        camera[1] = T(1);
        cpoint[0] = T(1) / T(2 * dratio);
        for(int zi = 0; zi < dratio; zi ++) {
          // N.B. projection scale is linear.
          cpoint[1] = T(zi) / T(dratio);
          // x-z plane projection of point p with camera geometry c to z=0.
          // c := camera, p := cpoint.
          // <c + (p - c) * t, [0, 1]> = 0
          const auto t(- camera[1] / (cpoint[1] - camera[1]));
          const auto y0((camera + (cpoint - camera) * t)[0] * rxy);
          if(abs(int(y0)) < 3 || rxy < abs(y0) * T(2)) continue;
          const auto Dop(diff<T>(abs(int(y0) & ~ int(1))));
          const auto Dop0((Dop.row(y0 / 2) + Dop.row(y0 / 2 + 1)) / T(2));
          // N.B. curvature matrix det == EG - F^2, we see only \< relation.
          for(int i = 0; i < data.rows(); i ++) {
            for(int j = 0; j < data.cols(); j ++) {
              T zy(0), zx(0);
              for(int kk = 0; kk < Dop0.size(); kk ++) {
                zy += data(getImgPt<int>(i + kk - Dop0.size() / 2,
                  data.rows()), j) * Dop0[kk];
                zx += data(i, getImgPt<int>(j + kk - Dop0.size() / 2,
                  data.cols()) ) * Dop0[kk];
              }
              const auto E(T(1) + zy * zy);
              const auto F(       zy * zx);
              const auto G(T(1) + zx * zx);
              const auto lscore(abs(E * G - F * F));
              if(zscore(i, j) < lscore) {
                result(i, j) = T(zi + 1);
                zscore(i, j) = lscore;
              }
            }
          }
        }
        return result;
      }
      break;
    case EXTEND_Y:
      {
        SimpleMatrix<T> result(data.rows() + 2 * recur, data.cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
        for(int i = 0; i < data.rows(); i ++)
          result.row(i + recur) = data.row(i);
        for(int i = 0; i < recur; i ++) {
          const auto next(pnext<T>(int(data.rows()) / (i + 1)));
          result.row(recur - i - 1).O();
          result.row(data.rows() + recur + i).O();
          for(int j = 0; j < next.size(); j ++) {
            result.row(recur - i - 1) += data.row(- (j - next.size() + 1) * (i + 1)) * next[j];
            result.row(data.rows() + recur + i) += data.row(  (j - next.size() + 1) * (i + 1) + data.rows() - 1) * next[j];
          }
        }
        return result;
      }
      break;
    case BLINK_Y:
      {
        auto dif(dft<T>(data.rows()) * data.template cast<complex<T> >());
        for(int i = 1; i < data.rows(); i ++) {
          dif.row(i) *= - complex<T>(T(0), T(2)) * T(i) / T(data.rows());
        }
        dif = dft<T>(- data.rows()) * dif;
        for(int i = 1; i < dif.rows() - 1; i ++) {
          dif.row(i) += (dif.row(i - 1) + dif.row(i + 1)) * complex<T>(T(recur) / T(256));
        }
        dif = dft<T>(data.rows()) * dif;
        for(int i = 1; i < data.rows(); i ++) {
          dif.row(i) /= - complex<T>(T(0), T(2)) * T(i) / T(data.rows());
        }
        return (dft<T>(- data.rows()) * dif).template real<T>();
      }
      break;
    case LPF_Y:
      return (dft<T>(- data.rows()).subMatrix(0, 0, data.rows(), recur) * (dft<T>(data.rows()).subMatrix(0, 0, recur, data.rows()) * data.template cast<complex<T> >())).template real<T>();
    case REPRESENT:
      return Decompose<T>(recur).represent(data, n);
    case CLIP:
      {
        SimpleMatrix<T> result(data.rows(), data.cols());
        for(int i = 0; i < result.rows(); i ++)
          for(int j = 0; j < result.cols(); j ++)
            result(i, j) = max(T(0), min(T(1), data(i, j)));
        return result;
      }
      break;
    case ABS:
      {
        SimpleMatrix<T> result(data.rows(), data.cols());
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
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < n; i ++) {
    cerr << "r" << flush;
    const auto theta((T(i) - T(n - 1) / T(2)) * atan(T(1)) / (T(n) / T(2)));
          auto work(center<T>(rotate<T>(filter<T>(rotate<T>(data, theta), dir, 0, recur), - theta), res));
#if defined(_OPENMP)
#pragma omp critical
#endif
    {
      res += move(work);
    }
  }
  return res /= T(n + 1);
}

#define _FILTER_
#endif

