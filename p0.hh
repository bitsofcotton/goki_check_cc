/*
 BSD 3-Clause License

Copyright (c) 2019-2020, bitsofcotton (kazunobu watatsu)
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

template <typename T, bool recur = true> class P0 {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleMatrix<T> Mat;
  typedef SimpleMatrix<complex<T> > MatU;
  inline P0();
  inline ~P0();
  const MatU& seed(const int& size0);
  const Mat&  diff(const int& size);
  inline Vec  taylor(const int& size, const T& step);
  const Vec&  next(const int& size);
  const T&    Pi() const;
  const complex<T>& J() const;
};

template <typename T, bool recur> inline P0<T,recur>::P0() {
  ;
}

template <typename T, bool recur> inline P0<T,recur>::~P0() {
  ;
}

template <typename T, bool recur> const T& P0<T,recur>::Pi() const {
  const static auto pi(atan2(T(1), T(1)) * T(4));
  return pi;
}

template <typename T, bool recur> const complex<T>& P0<T,recur>::J() const {
  const static auto i(complex<T>(T(0), T(1)));
  return i;
}

template <typename T, bool recur> const typename P0<T,recur>::MatU& P0<T,recur>::seed(const int& size0) {
  const auto size(abs(size0));
  assert(size);
  static vector<MatU> dft;
  static vector<MatU> idft;
  if(dft.size() <= size)
    dft.resize(size + 1, MatU());
  if(idft.size() <= size)
    idft.resize(size + 1, MatU());
  auto& edft( dft[size]);
  auto& eidft(idft[size]);
  if(edft.rows() != size || edft.cols() != size) {
    edft.resize(size, size);
    eidft.resize(size, size);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < edft.rows(); i ++) {
      for(int j = 0; j < edft.cols(); j ++) {
        const auto theta(- T(2) * Pi() * T(i) * T(j) / T(edft.rows()));
        edft(i, j)  = complex<T>(cos(  theta), sin(  theta));
        eidft(i, j) = complex<T>(cos(- theta), sin(- theta)) / complex<T>(T(size));
      }
    }
  }
  return size0 < 0 ? eidft : edft;
}

template <typename T, bool recur> const typename P0<T,recur>::Mat& P0<T,recur>::diff(const int& size0) {
  assert(size0);
  const auto size(abs(size0));
  static vector<Mat> D;
  static vector<Mat> I;
  if(D.size() <= size)
    D.resize(size + 1, Mat());
  if(I.size() <= size)
    I.resize(size + 1, Mat());
  auto& dd(D[size]);
  auto& ii(I[size]);
  if(dd.rows() != size || dd.cols() != size) {
    auto DD(seed(size));
    auto II(seed(size));
    for(int i = 0; i < DD.rows(); i ++)
      DD.row(i) *= J() * T(2) * Pi() * T(i) / T(DD.rows());
    for(int i = 1; i < DD.rows(); i ++)
      II.row(i) /= J() * T(2) * Pi() * T(i) / T(DD.rows()) / Pi();
    dd = (seed(- size) * DD).template real<T>() / Pi();
    ii = (seed(- size) * II).template real<T>();
  }
  return size0 < 0 ? ii : dd;
}

template <typename T, bool recur> inline typename P0<T,recur>::Vec P0<T,recur>::taylor(const int& size, const T& step) {
  const int  step00(max(0, min(size - 1, int(floor(step)))));
  const auto residue0(step - T(step00));
  const auto step0(step00 == size - 1 || abs(residue0) <= T(1) / T(2) ? step00 : step00 + 1);
  const auto residue(step - T(step0));
        Vec  res(size);
  for(int i = 0; i < size; i ++)
    res[i] = i == step0 ? T(1) : T(0);
  if(residue == T(0))
    return res;
  // N.B.
  // if we deal with (D *= r, residue /= r), it is identical with (D, residue)
  // So ||D^n * residue^n|| / T(n!) < 1 case, this loop converges.
  // but with n^n v.s. n!, differential of n! is faster than n^n.
  // (n! < n^n but a^n < n! somewhere).
  // And, we treat D * residue as a block, so Readme.md's condition 1/x^k needs
  // to be in the series in this.
  const auto& D(diff(size));
        auto  dt(D.col(step0) * residue);
  for(int i = 2; ; i ++) {
    const auto last(res);
    res += dt;
    if(last == res) break;
    dt   = D * dt * residue / T(i);
  }
  return res;
}

template <typename T, bool recur> const typename P0<T,recur>::Vec& P0<T,recur>::next(const int& size) {
  assert(0 < size);
  static vector<Vec> P;
  if(P.size() <= size)
    P.resize(size + 1, Vec());
  auto& p(P[size]);
  if(p.size() != size) {
    p = taylor(size, T(size));
    if(recur) {
      std::cerr << "p" << std::flush;
      if(1 < size) {
        const auto& back(next(size - 1));
        for(int i = 0; i < back.size(); i ++)
          p[i - back.size() + p.size()] += back[i] * T(size - 1);
        p /= T(size);
      }
    }
  }
  return p;
}


template <typename T, bool recur = true> class P0B {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleVector<complex<T> > VecU;
  inline P0B();
  inline P0B(const int& size);
  inline ~P0B();
  inline T next(const T& in);
private:
  Vec buf;
};

template <typename T, bool recur> inline P0B<T,recur>::P0B() {
  ;
}

template <typename T, bool recur> inline P0B<T,recur>::P0B(const int& size) {
  buf.resize(size);
  for(int i = 0; i < buf.size(); i ++)
    buf[i] = T(0);
}

template <typename T, bool recur> inline P0B<T,recur>::~P0B() {
  ;
}

template <typename T, bool recur> inline T P0B<T,recur>::next(const T& in) {
  static P0<T,recur> p;
  for(int i = 0; i < buf.size() - 1; i ++)
    buf[i] = buf[i + 1];
  buf[buf.size() - 1] = in;
  return p.next(buf.size()).dot(buf);
}

#define _P0_
#endif

