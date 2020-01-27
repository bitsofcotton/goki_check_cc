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
  inline P1(const int& statlen, const int& varlen);
  inline ~P1();
  const T& next(const Vec& in);
  T    lasterr;
  Vec  fvec;
private:
  Vec  a;
  Mat  A;
  Vec  b;
  T    M;
  int  statlen;
  int  varlen;
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
  statlen = varlen = 0;
  M = threshold_feas = threshold_p0 = threshold_inner = lasterr = T(0);
}

template <typename T> inline P1<T>::P1(const int& statlen, const int& varlen) {
  assert(1 < varlen && varlen < statlen);
  this->statlen = statlen;
  this->varlen  = varlen;
  b.resize(statlen * 2);
  A.resize(b.size(), varlen + 1);
  a.resize(A.cols());
  M = lasterr = T(0);
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

template <typename T> const T& P1<T>::next(const Vec& in) {
  assert(in.size() == statlen + varlen);
  T MM(0);
  for(int i = 0; i < in.size(); i ++)
    MM = max(MM, abs(in[i]));
  for(int i = 0; i < statlen; i ++) {
    for(int j = 0; j < varlen; j ++)
      A(i, j) = in[i + varlen - j - 1];
    b[i]         = in[i + varlen];
    A(i, varlen) = MM;
  }
  for(int i = 0; i < statlen; i ++) {
    A.row(statlen + i) = - A.row(i);
    b[statlen + i]     = - b[i];
  }
  for(int i = 0; i < varlen; i ++)
    a[i] = in[statlen + varlen - i - 1];
  a[varlen] = MM;
  try {
    for(int i = 0; i < fvec.size(); i ++)
      fvec[i] = T(0);
    lasterr = MM * T(4);
    for(auto ratio0(MM * T(2)); MM * threshold_inner <= ratio0; ratio0 /= T(2)) {
      const auto ratio(lasterr - ratio0);
      int n_fixed;
      T   ratiob;
      T   normb0;
      Vec rvec;
      Vec on;
      Vec deltab;
      Vec mbb;
      Vec bb;
      for(int i = 0; i < Pt.rows(); i ++)
        for(int j = 0; j < Pt.cols(); j ++)
          Pt(i, j) = T(0);
      for(int i = 0; i < A.cols(); i ++) {
        const auto Atrowi(A.col(i));
        const auto work(Atrowi - Pt.projectionPt(Atrowi));
        Pt.row(i) = work / sqrt(work.dot(work));
      }
      const auto R(Pt * A);
      if(A.cols() == A.rows()) {
        rvec = Pt * (one * ratio + b);
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
          if(j < statlen)
            norm[j] /= T(2);
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
            f[j]     = b[i]      / lratio + (i < statlen ? ratio / T(2) : ratio);
            j ++;
          }
        assert(j == f.size());
        rvec = F.solve(f);
      } else
        rvec = Pt * (on * ratiob + deltab + b);
     pnext:
      const auto err(Pt.transpose() * rvec - b - one * ratio);
            T    errorM(0);
      for(int i = 0; i < b.size(); i ++) if(errorM < err[i])
        errorM = err[i];
      if(isfinite(errorM) && errorM <= sqrt(threshold_inner) * normb0) {
        lasterr -= ratio0;
        fvec     = R.solve(rvec);
      }
    }
    M = a.dot(fvec);
  } catch (const char* e) {
    M = T(0);
  }
  if(isfinite(M) && ! isinf(M))
    return M;
  return M = T(0);
}

#define _P1_
#endif

