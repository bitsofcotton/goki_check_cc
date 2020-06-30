/*
 BSD 3-Clause License

Copyright (c) 2019-2020, bitsofcotton
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

template <typename T, int ratio = 2> class P0 {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleMatrix<T> Mat;
  typedef SimpleMatrix<complex<T> > MatU;
  inline P0();
  inline ~P0();
  inline T    next(const Vec& in, const T& err = T(1) / T(80));
  inline Vec  taylor(const int& size, const T& step);
  const Mat&  diff(const int& size);
  const Vec&  minSq(const int& size);
private:
  const MatU& seed(const int& size, const bool& idft);
  const Vec&  nextP(const int& size);
  const Vec&  nextDeepP(const int& size);
  const Vec&  nextReverseDeepP(const int& size);
  const Vec&  nextBothP(const int& size);
  const T&    Pi() const;
  const complex<T>& J() const;
};

template <typename T, int ratio> inline P0<T,ratio>::P0() {
  ;
}

template <typename T, int ratio> inline P0<T,ratio>::~P0() {
  ;
}

template <typename T, int ratio> inline T P0<T,ratio>::next(const Vec& in, const T& err) {
  assert(in.size());
  Vec   work(in.size() + 1);
  for(int i = 0; i < in.size(); i ++)
    work[i] = in[i];
  auto& res(work[work.size() - 1]);
  res = nextBothP(in.size()).dot(in);
  const auto normin(sqrt(in.dot(in)));
  auto  tilt(normin);
  while(err * normin < abs(tilt)) {
    tilt  = minSq(work.size()).dot(work);
    Vec buf(in.size());
    for(int i = 0; i < buf.size(); i ++)
      buf[i] = in[i] - tilt * T(i);
    res   = nextBothP(buf.size()).dot(buf) + tilt * T(buf.size());
    tilt -= minSq(work.size()).dot(work);
  }
  return res;
}

template <typename T, int ratio> const T& P0<T,ratio>::Pi() const {
  const static auto pi(atan2(T(1), T(1)) * T(4));
  return pi;
}

template <typename T, int ratio> const complex<T>& P0<T,ratio>::J() const {
  const static auto i(complex<T>(T(0), T(1)));
  return i;
}

template <typename T, int ratio> const typename P0<T,ratio>::MatU& P0<T,ratio>::seed(const int& size, const bool& f_idft) {
  assert(0 < size);
  static vector<MatU> dft;
  static vector<MatU> idft;
  if(dft.size() <= size)
    dft.resize(size + 1, MatU());
  if(idft.size() <= size)
    idft.resize(size + 1, MatU());
  if((!f_idft) &&  dft[size].rows() == size &&  dft[size].cols() == size)
    return dft[size];
  if(  f_idft  && idft[size].rows() == size && idft[size].cols() == size)
    return idft[size];
  auto& edft( dft[size]);
  auto& eidft(idft[size]);
  edft.resize(size, size);
  eidft.resize(size, size);
  for(int i = 0; i < edft.rows(); i ++) {
    for(int j = 0; j < edft.cols(); j ++) {
      const auto theta(- T(2) * Pi() * T(i) * T(j) / T(edft.rows()));
      edft(i, j)  = complex<T>(cos(  theta), sin(  theta));
      eidft(i, j) = complex<T>(cos(- theta), sin(- theta)) / complex<T>(T(size));
    }
  }
  if(f_idft)
    return eidft;
  return edft;
}

template <typename T, int ratio> const typename P0<T,ratio>::Mat& P0<T,ratio>::diff(const int& size) {
  assert(0 < size);
  static vector<Mat> D;
  if(D.size() <= size)
    D.resize(size + 1, Mat());
  if(D[size].rows() == size && D[size].cols() == size)
    return D[size];
  auto& d(D[size]);
  d.resize(size, size);
  auto  DD(seed(size, false));
  DD.row(0) *= complex<T>(T(0));
  T nd(0);
  T ni(0);
  for(int i = 1; i < DD.rows(); i ++) {
    const auto phase(J() * T(2) * Pi() * T(i) / T(DD.rows()));
    const auto phase2(complex<T>(T(1)) / phase);
    DD.row(i) *= phase;
    nd += abs(phase)  * abs(phase);
    ni += abs(phase2) * abs(phase2);
  }
  return d = (seed(size, true) * DD).template real<T>() * sqrt(sqrt(T(DD.rows() - 1) / (nd * ni)));
}

template <typename T, int ratio> inline typename P0<T,ratio>::Vec P0<T,ratio>::taylor(const int& size, const T& step) {
  const auto  spt(min(size - 1, max(0, int(ceil(step)))));
  const auto  residue(step - T(spt));
        Vec   tayl(size);
  const auto& D(diff(size));
        auto  dt(D * residue);
  for(int i = 0; i < tayl.size(); i ++)
    tayl[i] = T(0);
  tayl[spt] = T(1);
  if(residue != T(0))
    for(int i = 2; 0 <= i; i ++) {
      const auto bt(tayl);
      tayl += dt.row(spt);
      if(bt == tayl)
        break;
      dt = (D * dt) * (residue / T(i));
    }
  return tayl;
}

template <typename T, int ratio> const typename P0<T,ratio>::Vec& P0<T,ratio>::nextP(const int& size) {
  assert(1 < size);
  static vector<Vec> P;
  if(P.size() <= size)
    P.resize(size + 1, Vec());
  if(P[size].size() == size)
    return P[size];
  auto& p(P[size]);
  Mat   extends(size * ratio - ratio + 1, size);
  for(int i = 0; i < extends.rows(); i ++)
    extends.row(i) = taylor(size, T(i) / T(ratio));
  const auto reverse(taylor(extends.rows(), - T(1)));
  p = taylor(extends.rows(), T(extends.rows()));
  for(int i = 0; i < reverse.size(); i ++)
    p[i] += reverse[reverse.size() - i - 1];
  p /= T(2);
  std::vector<Vec> pa;
  pa.reserve(ratio);
  pa.emplace_back(p);
  for(int i = 0; i < ratio - 1; i ++) {
    for(int j = 0; j < i + 1; j ++)
      p[j] = T(0);
    for(int j = i + 1; j < pa[0].size(); j ++)
      p[j] = pa[0][j - i - 1];
    for(int j = 0; j <= i; j ++)
      p   += pa[j] * pa[0][j - i + p.size() - 1];
    pa.emplace_back(p);
  }
  return p = extends.transpose() * p;
}

template <typename T, int ratio> const typename P0<T,ratio>::Vec& P0<T,ratio>::nextDeepP(const int& size) {
  assert(1 < size);
  static vector<Vec> P;
  if(P.size() <= size)
    P.resize(size + 1, Vec());
  if(P[size].size() == size)
    return P[size];
  auto& p(P[size]);
  if(size <= 3)
    return p = nextP(size);
  const auto pp(nextDeepP(size - 1) * T(size - 3));
  p = nextP(size);
  for(int i = 1; i < p.size(); i ++)
    p[i] += pp[i - 1];
  return p /= T(size - 2);
}

template <typename T, int ratio> const typename P0<T,ratio>::Vec& P0<T,ratio>::nextReverseDeepP(const int& size) {
  assert(1 < size);
  static vector<Vec> P;
  if(P.size() <= size)
    P.resize(size + 1, Vec());
  if(P[size].size() == size)
    return P[size];
  auto& p(P[size]);
  p.resize(size);
  if(size <= 3) {
    const auto pq(nextP(size));
    p[0] = T(1);
    for(int i = 1; i < pq.size(); i ++)
      p[i] = - pq[pq.size() - i];
    return p /= pq[0];
  }
  const auto pp(nextDeepP(size - 1) * T(size - 3));
  const auto pq(nextP(size));
  p[0] = T(1);
  for(int i = 1; i < pq.size(); i ++)
    p[i] = - pq[pq.size() - i];
  p /= pq[0];
  for(int i = 1; i < p.size(); i ++)
    p[i] += pp[i - 1];
  return p /= T(size - 2);
}

template <typename T, int ratio> const typename P0<T,ratio>::Vec& P0<T,ratio>::nextBothP(const int& size) {
  assert(1 < size);
  static vector<Vec> P;
  if(P.size() <= size)
    P.resize(size + 1, Vec());
  if(P[size].size() == size)
    return P[size];
  return P[size] = (nextDeepP(size) + nextReverseDeepP(size)) / T(2);
}

template <typename T, int ratio> const typename P0<T,ratio>::Vec& P0<T,ratio>::minSq(const int& size) {
  assert(1 < size);
  static vector<Vec> S;
  if(S.size() <= size)
    S.resize(size + 1, Vec());
  if(S[size].size() == size)
    return S[size];
  auto& s(S[size]);
  s.resize(size);
  const T    xsum(size * (size - 1) / 2);
  const T    xdot(size * (size - 1) * (2 * size - 1) / 6);
  const auto denom(xdot * T(size) - xsum * xsum);
  for(int i = 0; i < s.size(); i ++)
    s[i] = (T(i) * T(size) - xsum) / denom;
  return s;
}


template <typename T> class P0B {
public:
  typedef SimpleVector<T> Vec;
  inline P0B();
  inline P0B(const int& size);
  inline ~P0B();
  inline T next(const T& in);
private:
  Vec buf;
  Vec buf2;
  T   M;
  T   Mb;
};

template <typename T> inline P0B<T>::P0B() {
  ;
}

template <typename T> inline P0B<T>::P0B(const int& size) {
  buf.resize(size);
  for(int i = 0; i < buf.size(); i ++)
    buf[i] = T(0);
  buf2 = buf;
  M = Mb = T(0);
}

template <typename T> inline P0B<T>::~P0B() {
  ;
}

template <typename T> inline T P0B<T>::next(const T& in) {
  static P0<T> p;
  for(int i = 0; i < buf.size() - 1; i ++) {
    buf[ i] = buf[i + 1];
    buf2[i] = buf2[i + 1];
  }
  buf[buf.size() - 1] = in;
  buf2[buf2.size() - 1] += (Mb < num_t(0) ? M - in : in - M);
  const auto ms(p.minSq(buf2.size()).dot(buf2));
  Mb = in - M;
  M  = p.next(buf);
  return (Mb < num_t(0) ? - (M - in + ms) : M - in + ms) + in;
}

#define _P0_
#endif

