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
  inline Vec  taylor(const int& size, const T& step);
  const MatU& seed(const int& size0);
  const Mat&  diff(const int& size0);
  const Mat&  diffCalibrate(const int& size);
  const Vec&  nextP(const int& size);
  const Vec&  nextQ(const int& size);
  const Vec&  nextR(const int& size);
  const Vec&  nextS(const int& size);
  inline const Vec& next(const int& size);
  const Vec&  minSq(const int& size);
  const T&    Pi() const;
  inline T    dot1(const Vec& x);
  const complex<T>& J() const;
};

template <typename T> inline P0<T>::P0() {
  ;
}

template <typename T> inline P0<T>::~P0() {
  ;
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

template <typename T> const typename P0<T>::Mat& P0<T>::diff(const int& size0) {
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
    auto II(DD);
    DD.row(0) *= complex<T>(T(0));
    for(int i = 1; i < DD.rows(); i ++) {
      DD.row(i) *= J() * T(2) * Pi() * T(i) / T(DD.rows());
      II.row(i) /= J() * T(2) * Pi() * T(i) / T(DD.rows());
    }
    dd = (seed(- size) * DD).template real<T>();
    ii = (seed(- size) * II).template real<T>();
  }
  return size0 < 0 ? ii : dd;
}

template <typename T> const typename P0<T>::Mat& P0<T>::diffCalibrate(const int& size) {
  assert(0 < size);
  static vector<Mat> D;
  if(D.size() <= size)
    D.resize(size + 1, Mat());
  auto& dd(D[size]);
  dd = diff(size);
  Vec calibrate(dd.rows());
  for(int i = 0; i < calibrate.size(); i ++)
    calibrate[i] = sin(T(i) / T(calibrate.size()) * T(2) * Pi());
  return dd *= - T(2) * Pi() / T(dd.rows()) / dd.row(dd.rows() / 2).dot(calibrate);
}

template <typename T> inline typename P0<T>::Vec P0<T>::taylor(const int& size, const T& step) {
  const int  step00(max(0, min(size - 1, int(floor(step)))));
  const auto residue0(step - T(step00));
  const auto step0(step00 == size - 1 || abs(residue0) <= T(1) / T(2) ? step00 : step00 + 1);
  const auto residue(step - T(step0));
        Vec  res(size);
  for(int i = 0; i < size; i ++)
    res[i] = i == step0 ? T(1) : T(0);
  if(residue == T(0))
    return res;
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

template <typename T> const typename P0<T>::Vec& P0<T>::nextP(const int& size) {
  assert(1 < size);
  static vector<Vec> P;
  if(P.size() <= size)
    P.resize(size + 1, Vec());
  auto& p(P[size]);
  if(p.size() != size) {
    const auto reverse(taylor(size, - T(1)));
    p = taylor(size, T(size));
    for(int i = 0; i < reverse.size(); i ++)
      p[i] += reverse[reverse.size() - i - 1];
    p /= dot1(p);
  }
  return p;
}

template <typename T> const typename P0<T>::Vec& P0<T>::nextQ(const int& size) {
  assert(1 < size);
  static vector<Vec> P;
  if(P.size() <= size)
    P.resize(size + 1, Vec());
  auto& p(P[size]);
  if(p.size() != size) {
    const auto& pp(nextP(size));
    p.resize(size);
    p[0] = T(1);
    for(int i = 1; i < p.size(); i ++)
      p[p.size() - i] = - pp[i];
    p /= pp[0];
    p += pp;
    p /= dot1(p);
    std::cerr << "q" << std::flush;
  }
  return p;
}

template <typename T> const typename P0<T>::Vec& P0<T>::nextR(const int& size) {
  assert(1 < size);
  static vector<Vec> P;
  if(P.size() <= size)
    P.resize(size + 1, Vec());
  auto& p(P[size]);
  if(p.size() != size) {
    p = nextQ(size);
    for(int i = 3; i < size; i ++)
      for(int j = 0; j < i; j ++)
        p[j - i + p.size()] += nextQ(i)[j];
    p /= dot1(p);
  }
  return p;
}

template <typename T> const typename P0<T>::Vec& P0<T>::nextS(const int& size) {
  assert(1 < size);
  static vector<Vec> P;
  if(P.size() <= size)
    P.resize(size + 1, Vec());
  auto& p(P[size]);
  if(p.size() != size) {
    Mat minusMinSq(size + 1, size + 1);
    Mat predict(size + 1, size + 1);
    for(int i = 0; i < minusMinSq.rows(); i ++)
      for(int j = 0; j < minusMinSq.cols(); j ++)
        minusMinSq(i, j) = predict(i, j) = (i == j ? T(1) : T(0));
    for(int i = 1; i < minusMinSq.rows(); i ++)
      minusMinSq.row(i) -= minSq(size + 1) * T(i);
    predict(predict.rows() - 1, predict.cols() - 1) = T(0);
    const auto& pp(nextR(size));
    for(int i = 0; i < pp.size(); i ++)
      predict(predict.rows() - 1, i) = pp[i];
    auto retry(minusMinSq * predict);
    while(true) {
      const auto  btry(retry.row(retry.rows() - 1));
      retry = retry * retry;
      const auto& rtry(retry.row(retry.rows() - 1));
      if(abs(rtry.dot(btry) / sqrt(rtry.dot(rtry) * btry.dot(btry)) - T(1)) <= T(0)) break;
    }
    p.resize(size);
    for(int i = 0; i < p.size(); i ++)
      p[i] = retry(retry.rows() - 1, i) - (retry(1, i) - (i == 1 ? T(1) : T(0))) * T(p.size());
    p /= dot1(p);
  }
  return p;
}

template <typename T> inline const typename P0<T>::Vec& P0<T>::next(const int& size) {
  return nextS(size);
}

template <typename T> inline T P0<T>::dot1(const Vec& x) {
  auto sum(x[0]);
  for(int i = 1; i < x.size(); i ++)
    sum += x[i];
  return sum;
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
  return p.next(buf.size()).dot(buf);
}

#define _P0_
#endif

