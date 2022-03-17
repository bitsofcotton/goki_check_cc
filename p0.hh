/*
 BSD 3-Clause License

Copyright (c) 2019-2021, bitsofcotton (kazunobu watatsu)
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
          p[j - pp.size() + p.size()] += pp[j] * T(size - 1);
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

template <typename T> const SimpleVector<T> pnextcache(const int& size, const int& step) {
  assert(0 < size && 0 <= step);
  vector<vector<SimpleVector<T> > > cp;
  if(cp.size() <= size) cp.resize(size + 1, vector<SimpleVector<T> >());
  if(cp[size].size() <= step) cp[size].resize(step + 1, SimpleVector<T>());
  if(cp[size][step].size()) return cp[size][step];
  return cp[size][step] = pnext<T>(size, step);
}

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
    const auto ff(f.next(in));
    return f.full ? pnextcache<T>(ff.size(), step).dot(ff) : T(int(0));
  }
  int step;
  feeder f;
};

template <typename T> const SimpleMatrix<complex<T> > dftcache(const int& size) {
  assert(size != 0);
  vector<SimpleMatrix<complex<T> > > cdft;
  vector<SimpleMatrix<complex<T> > > cidft;
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
    (this->p).resize(size, p);
    q = this->p;
  }
  inline ~P0DFT() { ; };
  inline T next(const T& in) {
    const auto& fn(f.next(in));
    if(! f.full) return T(int(0));
    auto ff(dftcache<T>(fn.size()) * fn.template cast<complex<T> >());
    assert(ff.size() == p.size() && p.size() == q.size());
    for(int i = 0; i < ff.size(); i ++)
      ff[i] = complex<T>(p[i].next(ff[i].real()), q[i].next(ff[i].imag()));
    return (dftcache<T>(- fn.size()) * ff)[ff.size() - 1].real();
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
  inline sumChain() { ; }
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

#define _P0_
#endif

