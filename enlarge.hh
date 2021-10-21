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
  BUMP_SIDE,
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
template <typename T> SimpleMatrix<T> filter(const SimpleMatrix<T>& data, const direction_t& dir, const int& recur = 2) {
  assert(0 < recur);
  static const auto Pi(atan2(T(1), T(1)) * T(4));
  switch(dir) {
  case SHARPEN_BOTH:
    return filter<T>(filter<T>(data, SHARPEN_X, recur), SHARPEN_Y, recur);
  case ENLARGE_BOTH:
    return filter<T>(filter<T>(data, ENLARGE_X, recur), ENLARGE_Y, recur);
  case FLARGE_BOTH:
    return filter<T>(filter<T>(data, FLARGE_X, recur), FLARGE_Y, recur);
  case BLUR_BOTH:
    return filter<T>(filter<T>(data, BLUR_Y, recur), BLUR_X, recur);
  case EXTEND_BOTH:
    return filter<T>(filter<T>(filter<T>(filter<T>(data, EXTEND_X, recur), CLIP, recur), EXTEND_Y, recur), CLIP, recur);
  case BLINK_BOTH:
    return filter<T>(filter<T>(data, BLINK_X, recur), BLINK_Y, recur);
  case LPF_BOTH:
    return filter<T>(filter<T>(data, LPF_X, recur), LPF_Y, recur);
  case SHARPEN_X:
    return filter<T>(data.transpose(), SHARPEN_Y, recur).transpose();
  case ENLARGE_X:
    return filter<T>(data.transpose(), ENLARGE_Y, recur).transpose();
  case FLARGE_X:
    return filter<T>(data.transpose(), FLARGE_Y, recur).transpose();
  case BLUR_X:
    return filter<T>(data.transpose(), BLUR_Y, recur).transpose();
  case EXTEND_X:
    return filter<T>(data.transpose(), EXTEND_Y, recur).transpose();
  case BLINK_X:
    return filter<T>(data.transpose(), BLINK_Y, recur).transpose();
  case LPF_X:
    return filter<T>(data.transpose(), LPF_Y, recur).transpose();
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
      const auto zxx(diff<T>(data.rows()) * diff<T>(data.rows()) * data);
      const auto zxy(diff<T>(data.rows()) * data * diff<T>(data.cols()).transpose());
      const auto zyy(data * diff<T>(data.cols()).transpose() * diff<T>(data.cols()).transpose());
            auto res(data);
      res.O();
      for(int i = 0; i < res.rows(); i ++)
        for(int j = 0; j < res.cols(); j ++)
          res(i, j) = zxx(i, j) * zyy(i, j) - zxy(i, j) * zxy(i, j);
      return res;
    }
    break;
  case INTEG_BOTH:
    {
      const auto  zyy(diff<T>(- data.rows()) * data * diff<T>(  data.cols()).transpose());
      const auto& zxy(data);
      const auto  zxx(diff<T>(  data.rows()) * data * diff<T>(- data.cols()).transpose());
            auto res(data);
      res.O();
      for(int i = 0; i < res.rows(); i ++)
        for(int j = 0; j < res.cols(); j ++)
          res(i, j) = zxx(i, j) * zyy(i, j) - zxy(i, j) * zxy(i, j);
      return res;
    }
    break;
  case COLLECT_BOTH:
    return filter<T>(filter<T>(data, DETECT_BOTH, recur), ABS, recur);
  case BUMP_BOTH:
  case BUMP_SIDE:
    {
      assert(dir == BUMP_BOTH || data.cols() == 2);
      SimpleMatrix<T> result(data.rows(), data.cols());
      SimpleMatrix<T> zscore(data.rows(), data.cols());
      if(dir == BUMP_BOTH) {
        result.resize(data.rows(), data.cols());
        zscore.resize(data.rows(), data.cols());
      } else {
        result.resize(1, 2);
        zscore.resize(1, 2);
      }
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
        const auto Dop0((Dop.row(int(y0) / 2) + Dop.row(int(y0) / 2 + 1)) / T(2));
        const auto DDop(Dop * Dop);
        const auto DDop0((DDop.row(int(y0) / 2) + DDop.row(int(y0) / 2 + 1)) / T(2));
        // N.B. curvature matrix det == EG - F^2, we see only \< relation.
        if(dir == BUMP_BOTH)
          for(int i = 0; i < data.rows(); i ++)
            for(int j = 0; j < data.cols(); j ++) {
              T L(0), M(0), N(0);
              for(int kk = 0; kk < Dop0.size(); kk ++) {
                L += data(getImgPt<int>(i + kk - Dop0.size() / 2,
                  data.rows()), j) * DDop0[kk];
                N += data(i, getImgPt<int>(j + kk - Dop0.size() / 2,
                  data.cols()) ) * DDop0[kk];
                for(int ll = 0; ll < Dop0.size(); ll ++)
                  M += Dop0[kk] * Dop0[ll] *
                    data(getImgPt<int>(i + kk - Dop0.size() / 2, data.rows()),
                         getImgPt<int>(j + ll - Dop0.size() / 2, data.cols()));
              }
              const auto lscore(abs(L * N - M * M));
              if(zscore(i, j) < lscore) {
                result(i, j) = T(zi + 1);
                zscore(i, j) = lscore;
              }
            }
        else {
          SimpleVector<T> work(2);
          work.O();
          for(int i = 0; i < Dop0.size(); i ++)
            work += data.row(getImgPt<int>(data.rows() / 2 + i - Dop0.size() / 2, data.rows())) * Dop0[i];
          for(int i = 0; i < work.size(); i ++) {
            const auto lscore(abs(work[i]));
            if(zscore(0, i) < lscore) {
              result(0, i) = T(zi + 1);
              zscore(0, i) = lscore;
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
      for(int i = 1; i < data.rows(); i ++)
        dif.row(i) *= - complex<T>(T(0), T(2)) * T(i) / T(data.rows());
      dif = dft<T>(- data.rows()) * dif;
      for(int i = 1; i < dif.rows() - 1; i ++)
        dif.row(i) += (dif.row(i - 1) + dif.row(i + 1)) * complex<T>(T(recur) / T(256));
      dif = dft<T>(data.rows()) * dif;
      for(int i = 1; i < data.rows(); i ++)
        dif.row(i) /= - complex<T>(T(0), T(2)) * T(i) / T(data.rows());
      return (dft<T>(- data.rows()) * dif).template real<T>();
    }
    break;
  case LPF_Y:
    return (dft<T>(- data.rows()).subMatrix(0, 0, data.rows(), recur) * (dft<T>(data.rows()).subMatrix(0, 0, recur, data.rows()) * data.template cast<complex<T> >())).template real<T>();
  case REPRESENT:
    return Decompose<T>(recur).represent(data, 2);
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

#define _FILTER_
#endif

