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

template <typename T> class P0 {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleVector<complex<T> > VecU;
  typedef SimpleMatrix<T> Mat;
  typedef SimpleMatrix<complex<T> > MatU;
  inline P0();
  inline P0(const int& range);
  inline ~P0();
  inline T next(const Vec& in, const bool& deep = true);
  const Mat&  lpf(const int& size);
private:
  Vec pred;
  const MatU& seed(const int& size, const bool& idft);
  const Mat&  diff(const int& size);
  const Vec&  nextTaylor(const int& size, const int& step);
  const T&    Pi() const;
  const complex<T>& J() const;
};

template <typename T> inline P0<T>::P0() {
  ;
}

template <typename T> inline P0<T>::P0(const int& range) {
  assert(1 < range);
  const auto  look(1);
  // with convolution meaning, exchange diff and integrate.
  const auto& reverse(nextTaylor(range, - look));
  pred.resize(reverse.size());
  for(int i = 0; i < reverse.size(); i ++)
    pred[i] = reverse[reverse.size() - 1 - i];
  pred = lpf(range).transpose() * ((pred + nextTaylor(range, look)) / T(2));
}

template <typename T> inline P0<T>::~P0() {
  ;
}

template <typename T> inline T P0<T>::next(const Vec& in, const bool& deep) {
  assert(pred.size() == in.size());
  auto res(pred.dot(in));
  if(! deep)
    return res;
  auto residue(in - lpf(in.size()) * in);
  while(T(1) / T(1000) / T(1000) < residue.dot(residue)) {
    for(int i = 0; i < residue.size(); i ++)
      if(i % 2 == (residue.size() + 1) % 2)
        residue[i] = - residue[i];
    res     += pred.dot(residue);
    residue -= lpf(residue.size()) * residue;
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

template <typename T> const typename P0<T>::MatU& P0<T>::seed(const int& size, const bool& f_idft) {
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
  for(int i = 0; i < edft.rows(); i ++)
    for(int j = 0; j < edft.cols(); j ++) {
      const auto theta(- T(2) * Pi() * T(i) * T(j) / T(edft.rows()));
      edft(i, j)  = complex<T>(cos(  theta), sin(  theta));
      eidft(i, j) = complex<T>(cos(- theta), sin(- theta)) / complex<T>(T(size));
    }
  if(f_idft)
    return eidft;
  return edft;
}

template <typename T> const typename P0<T>::Mat& P0<T>::diff(const int& size) {
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
    const auto phase(- J() * T(2) * Pi() * T(i) / T(DD.rows()));
    const auto phase2(complex<T>(T(1)) / phase);
    DD.row(i) *= phase;
    nd += abs(phase)  * abs(phase);
    ni += abs(phase2) * abs(phase2);
  }
  return d = (seed(size, true) * DD).template real<T>() * sqrt(sqrt(T(DD.rows() - 1) / (nd * ni)));
}

template <typename T> const typename P0<T>::Mat& P0<T>::lpf(const int& size) {
  assert(0 < size);
  static vector<Mat> L;
  if(L.size() <= size)
    L.resize(size + 1, Mat());
  if(L[size].rows() == size && L[size].cols() == size)
    return L[size];
  auto& l(L[size]);
  auto  ll(seed(size, false));
  for(int i = size / 2; i < size; i ++)
    ll.row(i) *= complex<T>(T(0));
  return l = (seed(size, true) * ll).template real<T>();
}

template <typename T> const typename P0<T>::Vec& P0<T>::nextTaylor(const int& size, const int& step) {
  assert(0 < size && (step == 1 || step == - 1));
  static vector<Vec> P;
  static vector<Vec> M;
  if(P.size() <= size)
    P.resize(size + 1, Vec());
  if(M.size() <= size)
    M.resize(size + 1, Vec());
  if(0 < step && P[size].size() == size)
    return P[size];
  if(step < 0 && M[size].size() == size)
    return M[size];
        auto& p(P[size]);
        auto& m(M[size]);
  const auto& D(diff(size));
        auto  ddp(  D);
        auto  ddm(- D);
  p.resize(size);
  for(int i = 0; i < p.size(); i ++)
    p[i] = T(0);
  m = p;
  p[p.size() - 1] = T(1);
  m[0] = T(1);
  for(int i = 2; 0 <= i; i ++) {
    const auto bp(p);
    const auto bm(m);
    p += ddp.row(ddp.rows() - 1);
    m += ddm.row(0);
    if(bm == m && bp == p)
      break;
    ddp = (D * ddp) * (  T(1) / T(i));
    ddm = (D * ddm) * (- T(1) / T(i));
  }
  if(0 < step)
    return p;
  return m;
}

#define _P0_
#endif

