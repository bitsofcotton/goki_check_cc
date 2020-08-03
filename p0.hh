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

template <typename T> class P0 {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleMatrix<T> Mat;
  typedef SimpleMatrix<complex<T> > MatU;
  inline P0();
  inline ~P0();
  inline T    next(const Vec& in, const T& err = T(1) / T(8000));
  inline T    next0(const Vec& in, const T& err = T(1) / T(8000));
  inline Vec  taylor(const int& size, const T& step);
  const MatU& seed(const int& size0);
  const Mat&  diff(const int& size);
  inline Mat  diffinv(const int& size, const int& k);
private:
  const Mat&  lpf(const int& size0);
  const Vec&  nextP(const int& size);
  const Vec&  minSq(const int& size);
  const T&    Pi() const;
  const complex<T>& J() const;
};

template <typename T> inline P0<T>::P0() {
  ;
}

template <typename T> inline P0<T>::~P0() {
  ;
}

template <typename T> inline T P0<T>::next(const Vec& in, const T& err) {
  assert(in.size());
  auto res(next0(in, err));
  auto hpf(lpf(- in.size()) * in);
  while(in.dot(in) * err < hpf.dot(hpf)) {
    auto work(hpf);
    for(int i = 0; i < work.size(); i ++)
      work[i] *= (work.size() + i) & 1 ? - T(1) : T(1);
    res += next0(work, err);
    hpf  = lpf(- work.size()) * work;
  }
  return res;
}

template <typename T> inline T P0<T>::next0(const Vec& in, const T& err) {
  assert(in.size());
  Vec   work(in.size() + 1);
  for(int i = 0; i < in.size(); i ++)
    work[i] = in[i];
  auto& res(work[work.size() - 1]);
  res = nextP(in.size()).dot(in);
  const auto normin(sqrt(in.dot(in)) / T(in.size()));
  auto  tilt(normin);
  while(in.dot(in) != T(0) && err * normin < abs(tilt)) {
    tilt  = minSq(work.size()).dot(work);
    Vec buf(in.size());
    for(int i = 0; i < buf.size(); i ++)
      buf[i] = in[i] - tilt * T(i);
    res   = nextP(buf.size()).dot(buf) + tilt * T(buf.size());
    tilt -= minSq(work.size()).dot(work);
  }
  return res;
}

template <typename T> const T& P0<T>::Pi() const {
  const static auto pi(atan2(T(1), T(1)) * T(4));
  return pi;
}

template <typename T> const complex<T>& P0<T>::J() const {
  const static auto i(complex<T>(T(0), T(1)));
  return i;
}

template <typename T> const typename P0<T>::MatU& P0<T>::seed(const int& size0) {
  const auto size(abs(size0));
  assert(0 < size);
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

template <typename T> const typename P0<T>::Mat& P0<T>::diff(const int& size) {
  assert(0 < size);
  static vector<Mat> D;
  if(D.size() <= size)
    D.resize(size + 1, Mat());
  auto& dd(D[size]);
  if(dd.rows() != size || dd.cols() != size) {
    auto DD(seed(size));
    DD.row(0) *= complex<T>(T(0));
    for(int i = 1; i < DD.rows(); i ++)
      DD.row(i) *= J() * T(2) * Pi() * T(i) / T(DD.rows());
    dd = (seed(- size) * DD).template real<T>();
    Vec calibrate(dd.rows());
    for(int i = 0; i < calibrate.size(); i ++)
      calibrate[i] = sin(T(i) / T(calibrate.size()) * T(2) * Pi());
    dd /= - dd.row(dd.rows() / 2).dot(calibrate);
  }
  return dd;
}

template <typename T> inline typename P0<T>::Mat P0<T>::diffinv(const int& size, const int& k) {
  auto DD(seed(size));
  DD.row(0) *= complex<T>(T(0));
  for(int i = 1; i < DD.rows(); i ++)
    DD.row(i) *= pow(J() * T(2) * Pi() * T(i) / T(DD.rows()), T(1) / T(k));
  return (seed(- size) * DD).template real<T>();
}

template <typename T> const typename P0<T>::Mat& P0<T>::lpf(const int& size0) {
  const auto size(abs(size0));
  static vector<Mat> L;
  static vector<Mat> H;
  if(L.size() <= size)
    L.resize(size + 1, Mat());
  if(H.size() <= size)
    H.resize(size + 1, Mat());
  auto& l(L[size]);
  auto& h(H[size]);
  if(l.rows() != size || l.cols() != size) {
    auto LL(seed(size));
    auto HH(LL);
    for(int i = 0; i < LL.rows(); i ++)
      if(i < LL.rows() / 2)
        HH.row(i) *= complex<T>(T(0));
      else
        LL.row(i) *= complex<T>(T(0));
    l = (seed(- size) * LL).template real<T>();
    h = (seed(- size) * HH).template real<T>();
  }
  return size0 < 0 ? h : l;
}

template <typename T> inline typename P0<T>::Vec P0<T>::taylor(const int& size, const T& step) {
  const int  step0(max(0, min(size - 1, int(floor(step)))));
  const auto residue(step - T(step0));
        Vec  res(size);
  for(int i = 0; i < size; i ++)
    res[i] = i == step0 ? T(1) : T(0);
  if(residue != T(0)) {
    const auto& D(diff(size));
          auto  dt(D * residue);
    for(int i = 2; ; i ++) {
      const auto last(res);
      res += dt.col(step0);
      if(last == res) break;
      dt   = D * dt * residue / T(i);
    }
  }
  return res;
}

template <typename T> const typename P0<T>::Vec& P0<T>::nextP(const int& size) {
  assert(1 < size);
  static vector<Vec> P;
  if(P.size() <= size)
    P.resize(size + 1, Vec());
  auto& p(P[size]);
  if(p.size() != size) {
    const auto reverse(lpf(size).transpose() * taylor(size, - T(1)));
    p = lpf(size).transpose() * taylor(size, T(size));
    for(int i = 0; i < reverse.size(); i ++)
      p[i] += reverse[reverse.size() - i - 1];
    p /= T(2);
  }
  return p;
}

template <typename T> const typename P0<T>::Vec& P0<T>::minSq(const int& size) {
  assert(1 < size);
  static vector<Vec> S;
  if(S.size() <= size)
    S.resize(size + 1, Vec());
  auto& s(S[size]);
  if(s.size() != size) {
    s.resize(size);
    const T    xsum(size * (size - 1) / 2);
    const T    xdot(size * (size - 1) * (2 * size - 1) / 6);
    const auto denom(xdot * T(size) - xsum * xsum);
    for(int i = 0; i < s.size(); i ++)
      s[i] = (T(i) * T(size) - xsum) / denom;
  }
  return s;
}


template <typename T> class P0B {
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

template <typename T> inline P0B<T>::P0B() {
  ;
}

template <typename T> inline P0B<T>::P0B(const int& size) {
  buf.resize(size);
  for(int i = 0; i < buf.size(); i ++)
    buf[i] = T(0);
}

template <typename T> inline P0B<T>::~P0B() {
  ;
}

template <typename T> inline T P0B<T>::next(const T& in) {
  static P0<T> p;
  for(int i = 0; i < buf.size() - 1; i ++)
    buf[i] = buf[i + 1];
  buf[buf.size() - 1] = in;
  return p.next(buf);
}

#define _P0_
#endif

