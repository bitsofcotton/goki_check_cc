/*
 BSD 3-Clause License

Copyright (c) 2019-2022, bitsofcotton (kazunobu watatsu)
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#if !defined(_P0_)

using std::ifstream;
using std::ofstream;
using std::cerr;
using std::flush;
using std::endl;

template <typename T> SimpleVector<T> pnext(const int& size, const int& step = 1) {
  SimpleVector<T> p;
  if(size <= 1) {
    p.resize(1);
    p[0] = T(1);
  } else if(size <= 2) {
    p.resize(2);
    p[0] = T(1) / T(2) + T(step < 0 ? 1 : 0);
    p[1] = T(1) / T(2) + T(step < 0 ? 0 : 1);
    p   /= T(2);
  } else {
    const auto file(std::string("./.cache/lieonn/next-") +
      std::to_string(size) + std::string("-") + std::to_string(step) +
#if defined(_FLOAT_BITS_)
      std::string("-") + std::to_string(_FLOAT_BITS_)
#else
      std::string("-ld")
#endif
    );
    ifstream cache(file.c_str());
    if(cache.is_open()) {
      cache >> p;
      cache.close();
    } else {
      p = taylor(size, T(step < 0 ? step : size + step - 1));
      cerr << "." << flush;
      if(abs(step) * 2 < size) {
        const auto pp(pnext<T>(size - 1, step));
        for(int j = 0; j < pp.size(); j ++)
          p[step < 0 ? j : j - pp.size() + p.size()] += pp[j] * T(size - 1);
        p /= T(size);
        ofstream ocache(file.c_str());
        ocache << p << endl;
        ocache.close();
      }
    }
  }
  assert(p.size() == size);
  return p;
}

template <typename T> const SimpleVector<T>& pnextcache(const int& size, const int& step) {
  assert(0 < size && 0 <= step);
  static vector<vector<SimpleVector<T> > > cp;
  if(cp.size() <= size) cp.resize(size + 1, vector<SimpleVector<T> >());
  if(cp[size].size() <= step) cp[size].resize(step + 1, SimpleVector<T>());
  if(cp[size][step].size()) return cp[size][step];
  return cp[size][step] = (pnext<T>(size, step) + pnext<T>(size, step + 1)) / T(int(2));
}

template <typename T> class Pnull {
public:
  inline Pnull() { ; }
  inline ~Pnull() { ; };
  inline T next(const T& in) { return T(int(0)); }
};

template <typename T, typename feeder> class P0 {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleMatrix<T> Mat;
  inline P0() { ; }
  inline P0(const int& size, const int& step = 1) {
    f = feeder(size);
    this->step = step;
  }
  inline ~P0() { ; };
  inline T next(const T& in) {
    const auto& ff(f.next(in));
    return f.full ? pnextcache<T>(ff.size(), step).dot(ff) : T(int(0));
  }
  int step;
  feeder f;
};

template <typename T> const SimpleMatrix<complex<T> >& dftcache(const int& size) {
  assert(size != 0);
  static vector<SimpleMatrix<complex<T> > > cdft;
  static vector<SimpleMatrix<complex<T> > > cidft;
  if(0 < size) {
    if(cdft.size() <= size) cdft.resize(size + 1, SimpleMatrix<complex<T> >());
    if(cdft[size].rows() && cdft[size].cols()) return cdft[size];
    return cdft[size] = dft<T>(size);
  }
  if(cidft.size() <= abs(size)) cidft.resize(abs(size) + 1, SimpleMatrix<complex<T> >());
  if(cidft[abs(size)].rows() && cidft[abs(size)].cols()) return cidft[abs(size)];
  return cidft[abs(size)] = dft<T>(size);
}

template <typename T, typename P, typename feeder> class P0DFT {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleMatrix<T> Mat;
  inline P0DFT() { ; }
  inline P0DFT(P&& p, const int& size) {
    f = feeder(size);
    (this->p).reserve(size);
    q.reserve(size);
    (this->p).emplace_back(p);
    q.emplace_back((this->p)[0]);
    for(int i = 1; i < size; i ++) {
      (this->p).emplace_back((this->p)[0]);
      q.emplace_back((this->p)[0]);
    }
  }
  inline ~P0DFT() { ; };
  inline T next(const T& in) {
    const auto& fn(f.next(in));
    if(! f.full) return T(int(0));
    auto ff(dftcache<T>(fn.size()) * fn.template cast<complex<T> >());
    assert(ff.size() == p.size() && p.size() == q.size());
    for(int i = 0; i < ff.size(); i ++)
      if(! (ff[i].real() == T(int(0)) && ff[i].imag() == T(int(0)) ) )
        ff[i] = abs(p[i].next(abs(ff[i]))) * exp(complex<T>(T(int(0)), q[i].next(arg(ff[i]))));
    return dftcache<T>(- fn.size()).row(fn.size() - 1).dot(ff).real();
  }
  vector<P> p;
  vector<P> q;
  feeder f;
};

template <typename T, typename P> class northPole {
public:
  inline northPole() { ; }
  inline northPole(P&& p) { this->p = p; }
  inline ~northPole() { ; }
  inline T next(const T& in) {
    static const T zero(int(0));
    static const T one(int(1));
    static const T M(atan(one / sqrt(SimpleMatrix<T>().epsilon)));
    if(! isfinite(in) || in == zero) return in;
    auto ain(atan(in));
    //assert(- M < ain && ain < M);
    auto bin(atan(one / move(ain)));
    //assert(- M < bin && bin < M);
    auto pn(p.next(move(bin)));
    if(! isfinite(pn) || pn == zero) return in;
    auto res(tan(max(- M, min(M, one / tan(max(- M, min(M, move(pn))))))));
    if(isfinite(res)) return move(res);
    return in;
  }
  P p;
};

template <typename T, typename P, bool avg = false> class sumChain {
public:
  inline sumChain() { S = T(t ^= t); }
  inline sumChain(P&& p) { this->p = p; S = T(t ^= t); }
  inline ~sumChain() { ; }
  inline T next(const T& in) {
    S += in;
    if(! avg) return p.next(S) - S;
    const auto A(S / T(++ t));
    return p.next(in - A) + A;
  }
  T S;
  P p;
  myuint t;
};

template <typename T, typename P> class logChain {
public:
  inline logChain() { ; }
  inline logChain(P&& p) { this->p = p; S = T(int(0)); }
  inline ~logChain() { ; }
  inline T next(const T& in) {
    static const T zero(int(0));
    static const T one(int(1));
    static const auto epsilon(sqrt(sqrt(SimpleMatrix<T>().epsilon)));
    const auto bS(S);
    S += in;
    if(bS == zero) return zero;
    const auto dd(S / bS - one);
    if(! isfinite(dd)) return zero;
    return p.next(dd / epsilon) * S * epsilon;
  }
  T S;
  P p;
};

template <typename T, typename P> class P0Expect {
public:
  inline P0Expect() { ; }
  inline P0Expect(P&& p, const int& nyquist = 2, const int& offset = 0) {
    Mx = M = d = T(t ^= t);
    t -= offset;
    tM = nyquist;
    assert(0 < tM);
    this->p = p;
  }
  inline ~P0Expect() { ; }
  inline const T& next(const T& in) {
    if(0 <= t) d += in;
    if(++ t < tM) return M;
    Mx = max(Mx, abs(d /= T(tM * tM)) * T(int(2)));
    M  = max(- Mx, min(Mx, p.next(d)));
    d  = T(t ^= t);
    return M;
  }
  int t;
  int tM;
  T d;
  T M;
  T Mx;
  P p;
};

template <typename T> class P0maxRank {
public:
  inline P0maxRank() { ; }
  inline P0maxRank(const int& status, const int& var) {
    assert(0 < status && 0 < var);
    p0 = p0_st(p0_s6t(p0_s5t(p0_s4t(p0_s3t(p0_s2t(p0_1t(p0_0t(status, var), var) )) ) )) );
    rr = qq = q = r = p = p0;
    bt = - (this->status = status);
    t  = (btt ^= btt);
    M  = SimpleMatrix<T>(4, max(4, status)).O();
  }
  inline ~P0maxRank() { ; }
  inline T next(const T& in) {
    static const T epsilon(sqrt(sqrt(SimpleMatrix<T>().epsilon)));
    static const T zero(int(0));
    auto res(zero);
    for(int i = 0; i < M.cols() - 1; i ++)
      M.setCol(i, M.col(i + 1));
    M(0, M.cols() - 1) = p.next(in);
    M(1, M.cols() - 1) = q.next(in * (
      T(status) + T(t - bt) / T(status) ));
    M(2, M.cols() - 1) = r.next(in * (
      T(status) - T(t - bt) / T(status) ));
    M(3, M.cols() - 1) = - avg.next(in);
    qq.next(in * (T(status) + T(t - btt) / T(status)));
    rr.next(in * (T(status) - T(t - btt) / T(status)));
    bvg.next(in);
    auto MM(M);
    for(int i = 0; i < MM.rows(); i ++) {
      const auto norm2(MM.row(i).dot(MM.row(i)));
      if(norm2 != zero) MM.row(i) /= sqrt(norm2);
    }
    const auto lsvd(MM.SVD());
    const auto svd(lsvd * MM);
    vector<T> stat;
    stat.reserve(svd.rows());
    for(int i = 0; i < svd.rows(); i ++)
      stat.emplace_back(sqrt(svd.row(i).dot(svd.row(i))));
    auto sstat(stat);
    sort(sstat.begin(), sstat.end());
    for(int i = 0; i < svd.rows(); i ++)
      if(sstat[sstat.size() - 1] * epsilon < stat[i]) {
        auto sum(lsvd(i, 0));
        for(int j = 1; j < lsvd.cols(); j ++)
          sum += lsvd(i, j);
        res += svd(i, svd.cols() - 1) / stat[i] * sgn<T>(sum);
      }
    if(! isfinite(res)) res = zero;
    if(! ((++ t) % status)) {
      q    = qq;
      r    = rr;
      qq   = p0;
      rr   = p0;
      bt   = btt;
      btt += status;
      avg  = bvg;
      bvg  = cvg;
      M.O();
    }
    return res;
  }
  // N.B. on existing taylor series.
  //      if the sampling frequency is not enough, middle range of the original
  //      function frequency (enough large bands) will effect prediction fail.
  //      this is because we only observes highest and lowest frequency on
  //      sampling points, so omitted part exists.
  //      even if the parameter on P0 is large, situation unchange.
  //      so we should use sectional measurement for them.
  typedef P0<T, idFeeder<T> > p0_0t;
  // N.B. sectional measurement, also expected value.
  typedef shrinkMatrix<T, p0_0t> p0_1t;
/*
  // N.B. make information-rich not to associative/commutative.
  //      2 dimension semi-order causes (x, status) from input as sedenion.
  typedef P0DFT<T, p0_1t, idFeeder<T> > p0_2t;
  typedef P0DFT<T, p0_2t, idFeeder<T> > p0_3t;
  typedef P0DFT<T, p0_3t, idFeeder<T> > p0_4t;
  typedef P0DFT<T, p0_4t, idFeeder<T> > p0_5t;
  // N.B. on any R to R into reasonable taylor.
  typedef northPole<T, p0_5t> p0_6t;
  typedef northPole<T, p0_6t> p0_7t;
  // N.B. we make the prediction on (delta) summation.
  typedef sumChain<T, p0_7t>  p0_8t;
  // N.B. we treat periodical part as non aligned complex arg part.
  typedef logChain<T, p0_8t>  p0_9t;
  typedef logChain<T, p0_9t>  p0_10t;
  // N.B. we take average as origin of input.
  typedef sumChain<T, p0_10t, true> p0_t;
  // N.B. if original sample lebesgue integrate is not enough continuous,
  //      imitate original function by some of sample points,
  //      but move origin point to average one, so a little better
  //      original function estimation.
  // N.B. frequency space *= 2 causes nyquist frequency ok.
  // N.B. but this is equivalent from jammer on PRNG, and probe on some
  //      measurable phenomenon.
  //typedef P0Expect<T, p0_11t> p0_t;
  // N.B. this needs huge memory to run.
*/
  // N.B. plain complex form.
  typedef northPole<T, p0_1t>  p0_s2t;
  typedef northPole<T, p0_s2t> p0_s3t;
  typedef sumChain<T, p0_s3t>  p0_s4t;
  typedef logChain<T, p0_s4t>  p0_s5t;
  typedef logChain<T, p0_s5t>  p0_s6t;
  typedef sumChain<T, p0_s6t, true> p0_st;
  int t;
  int bt;
  int btt;
  int status;
  p0_st p0;
  p0_st p;
  p0_st q;
  p0_st r;
  p0_st qq;
  p0_st rr;
  sumChain<T, Pnull<T>, true> avg;
  sumChain<T, Pnull<T>, true> bvg;
  sumChain<T, Pnull<T>, true> cvg;
  SimpleMatrix<T> M;
};

template <typename T, typename P> class P0recur {
public:
  inline P0recur() { ; }
  inline P0recur(const int& status) {
    // N.B. parameters are not optimal but we use this.
    for(int i = status; i > 4; i = int(max(T(int(2)), ceil(exp(sqrt(log(T(i)))))))) {
      p.emplace_back(P(i, int(max(T(int(3)), ceil(exp(sqrt(log(T(i)))))))));
      std::cerr << i << ", ";
    }
    if(! p.size()) { p.emplace_back(P(status, 1)); std::cerr << status; }
    std::cerr << std::endl;
    M.resize(p.size(), T(int(0)));
  }
  inline ~P0recur() { ; }
  inline T next(const T& in) {
    auto d(M);
    T    res(int(1));
    for(int i = 0; i < d.size(); i ++) d[i] *= in;
    for(int i = 0; i < p.size(); i ++)
      res *= (M[i] = p[i].next(i ? d[i - 1] : in));
    return res;
  }
  vector<T> M;
  vector<P> p;
};

template <typename T, typename P> class P0ContRand {
public:
  inline P0ContRand() { ; }
  inline P0ContRand(P&& p, const int& para) {
    (this->p).reserve(para);
    (this->p).emplace_back(p);
    for(int i = 1; i < para; i ++)
      (this->p).emplace_back((this->p)[0]);
    r.resize(para, T(t ^= t));
    br.resize(para, T(t));
  }
  inline ~P0ContRand() { ; }
  inline T next(const T& in) {
    t ++;
    T res(0);
    for(int i = 0; i < p.size(); i ++) {
      const auto rr(t & 1 ? r[i] + br[i] : r[i] + r[i]);
      res += p[i].next(in * rr);
      if(! (t & 1)) {
        br[i] = r[i];
        r[i]  = T((random() & 0x7ffffff) + 1) / T(int(0x8000000));
      }
    }
    return res /= T(int(p.size()));
  }
  vector<P> p;
  vector<T> r;
  vector<T> br;
  int t;
};

template <typename T, typename P> class P0Binary01 {
public:
  inline P0Binary01() { ; }
  inline P0Binary01(P&& p, const int& depth = 1) { (this->p).resize(depth, p); }
  inline ~P0Binary01() { ; }
  inline T next(const T& in) {
    T pw(int(1));
    T res(int(0));
    for(int i = 0; i < p.size(); i ++) {
      res += pw * (p[i].next(T(int(in / pw) & 1) - T(int(1)) / T(int(2))) + T(int(1)) / T(int(2)) );
      pw /= T(int(2));
    }
    return res;
  }
  vector<P> p;
};

#define _P0_
#endif

