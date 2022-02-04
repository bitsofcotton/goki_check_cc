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
    p[0] = T(1) / T(2);
    p[1] = T(1) / T(2) + T(1);
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
      p = taylor(size, T(size + step - 1));
      cerr << "." << flush;
      const auto pp(pnext<T>(size - 1, step));
      for(int j = 0; j < pp.size(); j ++)
        p[j - pp.size() + p.size()] += pp[j] * T(size - 1);
      p /= T(size);
      ofstream ocache(file.c_str());
      ocache << p << endl;
      ocache.close();
    }
  }
  assert(p.size() == size);
  return p;
}

template <typename T, typename feeder> class P0 {
public:
  typedef SimpleVector<T> Vec;
  inline P0() { ; }
  inline P0(const int& size, const int& step = 1) {
    f = feeder(size);
    p = pnext<T>(size, step);
  }
  inline ~P0() { ; };
  inline T next(const T& in) {
    const auto& ff(f.next(in));
    static const T zero(int(0));
    return f.full ? p.dot(ff) : zero;
  }
  Vec p;
  feeder f;
};

template <typename T, typename P> class P0D {
public:
  typedef SimpleVector<T> Vec;
  inline P0D() { ; }
  inline P0D(const int& size, const int& step = 1) {
    p = P(size, step);
    q = P(size, step);
  }
  inline ~P0D() { ; };
  inline T next(const T& in) {
    const static T zero(int(0));
    const static T one(int(1));
    const static T two(int(2));
    if(in == zero) return zero;
    const auto qn(q.next(one / in));
    if(qn == zero) return p.next(in);
    return (p.next(in) + one / qn) / two;
  }
  P p;
  P q;
};

template <typename T, typename P> class PC {
public:
  inline PC() { ; }
  inline PC(const P& p0, const int& sz = 20) {
    p = P(p0);
    q = P(p0);
    bb.resize(sz, M = T(int(0)));
    assert(! (bb.size() & 1));
  }
  inline ~PC() { ; };
  inline T next(const T& in) {
    for(int i = 1; i < bb.size() - 2; i ++)
      bb[i - 1] = move(bb[i]);
    bb[bb.size() - 2]  = bb[bb.size() - 1];
    bb[bb.size() - 1] += in;
    vector<T> H;
    vector<T> h;
    vector<T> g;
    H.reserve(bb.size() / 2);
    h.reserve(bb.size() / 2);
    g.reserve(bb.size() / 2);
    for(int i = 0; i < bb.size() / 2; i ++) {
      g.emplace_back(bb[i * 2 + 0] + bb[i * 2 + 0]);
      H.emplace_back(bb[i * 2 + 0] + bb[i * 2 + 1]);
      h.emplace_back(bb[i * 2 + 1] + bb[i * 2 + 1]);
    }
    const auto MM(q.next(M * (H[H.size() - 1] - h[h.size() - 1] + g[g.size() - 2] - H[H.size() - 2])));
    auto pp(p);
    for(int i = 1; i < H.size(); i ++)
      M = pp.next(H[i] - h[i] + g[g.size() - 1] - H[H.size() - 1]);
    return - M * MM;
  }
  P p;
  P q;
  T M;
  vector<T> bb;
};

template <typename T, typename P> class pC {
public:
  inline pC() { ; }
  inline pC(const P& p0, const T& lower, const T& upper, const int& bit) {
    assert(T(int(0)) < lower && lower < upper);
    p.resize(bit, p0);
    r = exp(log(upper - lower) / T(int(bit)));
    this->lower = lower;
    bres = T(int(0));
  }
  inline ~pC() { ; };
  inline T next(const T& in) {
    auto res(bres * T(int(2)));
    bres = T(int(0));
    for(int i = 0; i < p.size(); i ++) {
      if(p[i].next(T(int(myint(in / lower / pow(r, T(i)))) & 1) - T(int(1)) / T(int(2))) > T(int(0)))
        bres += lower * pow(r, T(i));
    }
    return res -= in;
  }
  vector<P> p;
  T lower;
  T r;
  T bres;
};

#define _P0_
#endif

