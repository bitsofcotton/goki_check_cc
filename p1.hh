/*
BSD 3-Clause License

Copyright (c) 2020, kazunobu watatsu
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
#if !defined(_P1_)

// data depend prediction:
//   in taylor series meaning,
//     f(x)=sum a_k*x_k^k, if such with relation, f(x)=sum a_k*x_k*(x_k\pm b_k).
//                    so with better conversion it'll be f(x)=sum a_k*x_k^k.
//   <a,x>==x+ series, if we suppose Diffmat and Integrate mat,
//                     we get for all k in Z, x^k.
//                     and with some start point x_0, we get <a,x>==x+ + const.
//   if x is the same method, <a,Bx> == <B^ta,x> == <a',x>.
//   in x, y, exp(k*log(x)+l*log(y)) series, there exists x^ky^l
//     if <x,y> == 0, and, <x, y> != 0 case:
//       series x + series x - <x,y>x/<x,x> converges to it.
//   if we predict with d_k := d_{k-1} + sum(d0_k - d0_k-n),
//     predict with A a, <a, x>.
//   if we predict with copycat, it's in trivial in this <a+ - a, x>.
//   if we prepare new function g(x) from f(x) with such prediction,
//      it's the form g_k(x) = (b_k-<a_k,x>/||x||)*f_k(x) ..., recursive,
//        so it's in such prediction.
//      it causes g(x)_k:=<a'_k,x>, with lagrange polynomial, g(x):=<b,x>
//   if original x series contains x^{{k}} term partially,
//      k in R, there exists x^{{k}}^{{-k}} == x term, and, it can concat.
//   if we predict with last saved value error,
//      it's x+ + <a,x>/x[0], <=> x+ * x[0] + <a,x> with certain ratio in taylor
//      x is one of A.row(sigma), so each <a,x> is to be A_sigma * x.
//      And, A_sigma is below g_k(x), it's in the form A with x'.
//   if there's lpf, hpf in prediction, it's trivial in this,
//      a_k in lpf, shift A a, then, <a, x>.
//   if there exists so to call triangular series like DFT,
//      its series are trivial in this.
template <typename T> class P1 {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleMatrix<T> Mat;
  inline P1();
  inline P1(const int& statlen, const int& varlen, const bool& addconst = false);
  inline ~P1();
  const Vec& next(const Vec& in, const int& step = 1, const bool& vanish = false);
  T    lasterr;
private:
  Vec  fvec;
  Mat  A;
  Vec  b;
  int  statlen;
  int  varlen;
  bool addconst;
  T    threshold_feas;
  T    threshold_p0;
  T    threshold_inner;
  Vec  one;
  SimpleVector<bool> checked;
  SimpleVector<bool> fix;
  Mat  Pverb;
  Vec  norm;
  Vec  orth;
  Mat  F;
  Vec  f;
  Mat  Pt;
};

template <typename T> inline P1<T>::P1() {
  statlen  = varlen = 0;
  addconst = false;
  threshold_feas = threshold_p0 = threshold_inner = lasterr = T(0);
}

template <typename T> inline P1<T>::P1(const int& statlen, const int& varlen, const bool& addconst) {
  assert(1 < varlen && varlen < statlen);
  this->statlen  = statlen;
  this->varlen   = varlen;
  this->addconst = addconst;
  b.resize(statlen * 2 + 1);
  A.resize(b.size(), varlen + (addconst ? 1 : 0));
  lasterr = T(0);
  const auto epsilon(std::numeric_limits<T>::epsilon());
  // const auto epsilon(T(1) >> short(62));
  threshold_feas    = pow(epsilon, T(5) / T(6));
  threshold_p0      = pow(epsilon, T(4) / T(6));
  threshold_inner   = pow(epsilon, T(2) / T(6));
  one.resize(b.size());
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int i = 0; i < one.size(); i ++)
    one[i] = T(1);
  norm.resize(one.size());
  checked.resize(one.size());
  fix.resize(one.size());
  Pt.resize(A.cols(), A.rows());
  F.resize(A.cols(), A.cols());
  f.resize(A.cols());
  fvec.resize(A.cols());
}

template <typename T> inline P1<T>::~P1() {
  ;
}

template <typename T> const typename P1<T>::Vec& P1<T>::next(const Vec& in, const int& step, const bool& vanish) {
  assert(in.size() == statlen + varlen * step);
  T MM(0);
  for(int i = 0; i < statlen; i ++) {
    for(int j = 0; j < varlen; j ++)
      MM = max(MM, abs(A(i, j) = in[i + (varlen - j - 1) * step + (vanish ? step : 0)]));
    MM = max(MM, abs(b[i] = vanish ? T(0) : in[i + varlen * step]));
  }
  for(int i = 0; i < statlen; i ++) {
    if(addconst) A(i, varlen) = MM;
    A.row(statlen + i) = - A.row(i);
    b[statlen + i]     = - b[i];
  }
  A /= MM;
  b /= MM;
  for(int i = 1; i < A.cols(); i ++)
    A(statlen * 2, i) = T(0);
  A(statlen * 2, 0) = b[statlen * 2] = - T(1);
  for(int i = 0; i < fvec.size(); i ++)
    fvec[i] = T(0);
  lasterr = A.rows() + A.cols();
  for(int i = 0; i < Pt.rows(); i ++)
    for(int j = 0; j < Pt.cols(); j ++)
      Pt(i, j) = T(0);
  for(int i = 0; i < A.cols(); i ++) {
    const auto Atrowi(A.col(i));
    const auto work(Atrowi - Pt.projectionPt(Atrowi));
    Pt.row(i) = work / sqrt(work.dot(work));
  }
  const auto R(Pt * A);
  for(auto ratio0(lasterr / T(2));
           threshold_inner <= ratio0;
           ratio0 /= T(2)) {
    const auto ratio(lasterr - ratio0);
    int n_fixed;
    T   ratiob;
    T   normb0;
    Vec rvec;
    Vec on;
    Vec deltab;
    Vec mbb;
    Vec bb;
    if(A.cols() == A.rows()) {
      rvec = Pt * b;
      goto pnext;
    }
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int i = 0; i < one.size(); i ++)
      fix[i]  = false;
    bb = b - Pt.projectionPt(b);
    if(sqrt(bb.dot(bb)) <= threshold_feas * sqrt(b.dot(b))) {
      for(int i = 0; i < bb.size(); i ++)
        bb[i] = sqrt(Pt.col(i).dot(Pt.col(i)));
      const auto bbb(bb - Pt.projectionPt(bb));
      if(sqrt(bbb.dot(bbb)) <= threshold_feas * sqrt(bb.dot(bb))) {
        rvec  = Pt * (b - bb - bbb);
        goto pnext;
      }
      bb = bbb;
    }
    mbb    = - bb;
    normb0 = sqrt(mbb.dot(mbb));
    Pverb  = Pt;
    for(n_fixed = 0 ; n_fixed < Pverb.rows(); n_fixed ++) {
#if defined(_OPENMP)
#pragma omp simd
#endif
      for(int j = 0; j < Pverb.cols(); j ++) {
        norm[j]    = sqrt(Pverb.col(j).dot(Pverb.col(j)));
        checked[j] = fix[j] || norm[j] <= threshold_p0;
      }
      auto mb(mbb + norm * normb0 * ratio);
      mb -= (deltab = Pverb.projectionPt(mb));
      mb /= (ratiob = sqrt(mb.dot(mb)));
      on  = Pverb.projectionPt(- one) + mb * mb.dot(- one);
      int fidx(0);
      for( ; fidx < on.size(); fidx ++)
        if(!checked[fidx])
          break;
      for(int j = fidx + 1; j < on.size(); j ++)
        if(!checked[j] && on[fidx] / norm[fidx] < on[j] / norm[j])
          fidx = j;
      if(fidx >= one.size() || on[fidx] / norm[fidx] <= threshold_inner) {
        if(fidx < one.size())
          on /= sqrt(norm.dot(norm));
        if(one.size() <= fidx || n_fixed < Pverb.rows() - 1) {
          n_fixed --;
          break;
        }
      }
      orth = Pverb.col(fidx);
      const auto norm2orth(orth.dot(orth));
      const auto mbb0(mbb[fidx]);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int j = 0; j < Pverb.cols(); j ++) {
        const auto work(Pverb.col(j).dot(orth) / norm2orth);
        Pverb.setCol(j, Pverb.col(j) - orth * work);
        mbb[j] -= mbb0 * work;
      }
      mbb[fidx] = T(0);
      fix[fidx] = true;
    }
    if(n_fixed == Pt.rows()) {
      int j(0);
      for(int i = 0; i < Pt.cols() && j < f.size(); i ++)
        if(fix[i]) {
          const auto lratio(sqrt(Pt.col(i).dot(Pt.col(i)) + b[i] * b[i]));
          F.row(j) = Pt.col(i) / lratio;
          f[j]     = b[i]      / lratio + ratio;
          j ++;
        }
      assert(j == f.size());
      try {
        rvec = F.solve(f);
      } catch (const char* e) {
        continue;
      }
    } else
      rvec = Pt * (on * ratiob + deltab + b);
   pnext:
    SimpleVector<T> err0(Pt.cols());
    for(int i = 0; i < err0.size(); i ++)
      err0[i] = Pt.col(i).dot(rvec);
    auto err(err0 - b - one * ratio);
    for(int i = 0; i < b.size(); i ++)
      if(err[i] <= T(0)) err[i] = T(0);
    if(sqrt(err.dot(err)) <= sqrt(threshold_inner * err0.dot(err0))) {
      try {
        fvec = R.solve(rvec);
      } catch(const char* e) {
        continue;
      }
      lasterr -= ratio0;
    }
  }
  return fvec;
}


template <typename T> class P1B {
public:
  typedef SimpleVector<T> Vec;
  inline P1B();
  inline P1B(const int& stat, const int& var);
  inline ~P1B();
  inline T next(const T& in);
private:
  P1<T> p;
  Vec   buf;
};

template <typename T> inline P1B<T>::P1B() {
  ;
}

template <typename T> inline P1B<T>::P1B(const int& stat, const int& var) {
  buf.resize(stat + var);
  for(int i = 0; i < buf.size(); i ++)
    buf[i] = T(0);
  p = P1<T>(stat, var);
}

template <typename T> inline P1B<T>::~P1B() {
  ;
}

template <typename T> inline T P1B<T>::next(const T& in) {
  for(int i = 0; i < buf.size() - 1; i ++)
    buf[i] = buf[i + 1];
  buf[buf.size() - 1] = in;
  const auto& fvec(p.next(buf));
  T res(0);
  for(int i = 0; i < fvec.size(); i ++)
    res += fvec[i] * buf[buf.size() - 1 - i];
  if(! isfinite(res) || isnan(res))
    res  = in;
  return res;
}


#define _P1_
#endif

