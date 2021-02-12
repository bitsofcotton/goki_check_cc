/*
BSD 3-Clause License

Copyright (c) 2020-2021, kazunobu watatsu
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

// Get invariant structure that
// \[0, &alpha;\[ register computer with deterministic calculation
// that whole in/out is shown in the varlen.
// Please refer bitsofcotton/randtools .
// [[A I], [I A^-1]] [x y], if rank A = dim x = dim y,
// rank [[A I], [I A^-1]] == rank A.
// so with this, U [[B], [O]] [x y] == 0.
// this returns one of u_k that concerns 0.
// skip attemps maximum skip status number if original data has a noised ones.
template <typename T> SimpleVector<T> invariantP1(const SimpleVector<T>& in, const int& varlen, const int& skip = 0, const int& guard = 0) {
  assert(varlen * 2 < in.size());
#if defined(_FLOAT_BITS_)
  const auto epsilon(T(1) >> int64_t(mybits - 2));
#else
  const auto epsilon(std::numeric_limits<T>::epsilon());
#endif
  const auto threshold_inner(sqrt(epsilon));
  SimpleMatrix<T> A((in.size() - varlen - guard) * 2 + 1, varlen + 1);
  SimpleVector<T> fvec(A.cols());
  SimpleVector<T> one(A.rows());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < in.size() - varlen - guard; i ++) {
    for(int j = 0; j < varlen; j ++)
      A(i, j) = in[i + j + guard];
    A(i, varlen) = T(0);
    const auto norm(sqrt(A.row(i).dot(A.row(i))));
    if(norm != T(0))
      A.row(i) /= norm;
    A(i, varlen) = T(1) / sqrt(T(varlen));
    A.row(i + A.rows() / 2) = - A.row(i);
  }
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < A.rows(); i ++)
    one[i] = T(1);
  one /= sqrt(one.dot(one));
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < fvec.size() - 1; i ++)
    A(A.rows() - 1, i) = fvec[i] = T(0);
  A(A.rows() - 1, A.cols() - 1) = T(1);
  fvec[fvec.size() - 1] = T(0);
  SimpleMatrix<T> Pt(A.cols(), A.rows());
  for(int i = 0; i < Pt.rows(); i ++)
    for(int j = 0; j < Pt.cols(); j ++)
      Pt(i, j) = T(0);
  for(int i = 0; i < A.cols(); i ++) {
    const auto Atrowi(A.col(i));
    const auto work(Atrowi - Pt.projectionPt(Atrowi));
    Pt.row(i) = work / sqrt(work.dot(work));
    const auto npt(Pt.row(i).dot(Pt.row(i)));
    if(! (isfinite(npt) && ! isnan(npt)))
      return fvec;
  }
  const auto Pt0(Pt);
  const auto R(Pt * A);
  SimpleVector<T> on;
  SimpleVector<T> orth;
  for(int n_fixed = 0 ; n_fixed < Pt.rows() - 1; n_fixed ++) {
    on = Pt.projectionPt(- one);
    std::vector<pair<T, int> > onM;
    onM.reserve(on.size());
    for(int j = 0; j < on.size(); j ++)
      onM.emplace_back(std::make_pair(on[j], j));
    std::sort(onM.begin(), onM.end());
    auto fidx(- 1);
    for(int i = onM.size() - 1 - skip; i < onM.size(); i ++)
      if(threshold_inner < onM[i].first) {
        fidx = onM[i].second;
        break;
      }
    if(fidx < 0)
      break;
    orth = Pt.col(fidx);
    const auto norm2orth(orth.dot(orth));
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
    for(int j = 0; j < Pt.cols(); j ++)
      Pt.setCol(j, Pt.col(j) - orth * Pt.col(j).dot(orth) / norm2orth);
  }
  fvec = R.solve(Pt0 * (- one) - Pt * on);
  const auto nfv(fvec.dot(fvec));
  return (isfinite(nfv) ? (nfv == T(0) ? fvec : fvec / sqrt(nfv))
    : fvec * T(0));
}


template <typename T> class P1Istatus {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleMatrix<T> Mat;
  inline P1Istatus();
  inline P1Istatus(const int& stat, const int& var, const int& instat);
  inline ~P1Istatus();
         bool invariantNext(const T& in, const int& skip = 0, const int& guard = 0);
  inline T    next(const T& in, const int& skip = 0, const int& guard = 0);
  Vec invariant_abs;
  Vec invariant_sgn;
  Vec status_abs;
  Vec status_sgn;
  int s_abs;
  int s_sgn;
  T   rabs;
  T   rsgn;
private:
  inline const T& sgn(const T& x) const;
  Vec buf;
  Vec sbuf;
  Vec projinv_abs;
  Vec projinv_sgn;
  Vec v_abs;
  Vec v_sgn;
  Mat m_abs;
  Mat m_sgn;
  Vec wabs;
  Vec wsgn;
  int varlen;
  int t;
  int si;
  bool once;
};

template <typename T> inline P1Istatus<T>::P1Istatus() {
  rabs = rsgn = T(t = s_abs = s_sgn = si = 0);
  once = true;
}

template <typename T> inline P1Istatus<T>::P1Istatus(const int& stat, const int& var, const int& instat) {
  buf.resize(stat + (varlen = var) + instat - 1);
  for(int i = 0; i < buf.size(); i ++)
    buf[i] = T(0);
  sbuf = buf;
  status_abs.resize(instat);
  for(int i = 0; i < status_abs.size(); i ++)
    status_abs[i] = T(0);
  status_sgn = status_abs;
  v_abs.resize(var + stat);
  for(int i = 0; i < v_abs.size(); i ++)
    v_abs[i] = T(0);
  v_sgn = v_abs;
  m_abs.resize(status_abs.size(), v_abs.size());
  for(int i = 0; i < m_abs.rows(); i ++)
    for(int j = 0; j < m_abs.cols(); j ++)
      m_abs(i, j) = T(0);
  m_sgn = m_abs;
  rabs = rsgn = T(t = s_abs = s_sgn = si = 0);
  once = true;
}

template <typename T> inline P1Istatus<T>::~P1Istatus() {
  ;
}

template <typename T> inline const T& P1Istatus<T>::sgn(const T& x) const {
  static const T zero(0);
  static const T one(1);
  static const T mone(- one);
  return x != zero ? (x < zero ? mone : one) : zero;
}

// Get invariant with inner status.
// Get <a, [x0 y0 x1 y1]> = 0 case, which is <a', [x0 x1]> = <a'', [y0 y1]> .
// With non invariant case on them, we get X a' == b == B [y0 y1].
template <typename T> bool P1Istatus<T>::invariantNext(const T& in, const int& skip, const int& guard) {
  assert(buf.size() == sbuf.size());
  for(int i = 0; i < buf.size() - 1; i ++) {
    buf[i]  = buf[i + 1];
    sbuf[i] = sbuf[i + 1];
  }
  buf[ buf.size()  - 1] = abs(in);
  sbuf[sbuf.size() - 1] = sgn(in) + T(2);
  if(t ++ < buf.size()) return false;
  else if(! projinv_abs.size() || ! projinv_sgn.size()) {
    const auto normbuf(sqrt(buf.dot(buf)));
    const auto normsbuf(sqrt(sbuf.dot(sbuf)));
    projinv_abs = invariantP1<T>(buf / normbuf, varlen, skip, guard);
    projinv_sgn = invariantP1<T>(sbuf / normsbuf, varlen, skip, guard);
  } else {
    assert(projinv_abs.size() == projinv_sgn.size());
    const auto normbuf(sqrt(buf.dot(buf)));
    const auto normsbuf(sqrt(sbuf.dot(sbuf)));
    for(int i = 0; i < v_abs.size() - 1; i ++) {
      v_abs[i] = v_abs[i + 1];
      v_sgn[i] = v_sgn[i + 1];
    }
    auto& vva(v_abs[v_abs.size() - 1]);
    auto& vvs(v_sgn[v_sgn.size() - 1]);
    vva = projinv_abs[projinv_abs.size() - 1] / sqrt(T(varlen));
    vvs = projinv_sgn[projinv_abs.size() - 1] / sqrt(T(varlen));
    for(int i = 0; i < projinv_abs.size() - 1; i ++) {
      vva +=  buf[i + 1] * projinv_abs[i];
      vvs += sbuf[i + 1] * projinv_sgn[i];
    }
    si ++;
    if(si < status_abs.size()) return false;
    if(s_abs < status_abs.size() || s_sgn < status_sgn.size()) {
      if(s_abs < status_abs.size()) {
        const auto p_abs(v_abs - m_abs.projectionPt(v_abs));
        const auto n2p_abs(p_abs.dot(p_abs));
        if(T(1) / T(varlen) < n2p_abs)
          m_abs.row(s_abs ++) = p_abs / sqrt(n2p_abs);
      }
      if(s_sgn < status_sgn.size()) {
        const auto p_sgn(v_sgn - m_sgn.projectionPt(v_sgn));
        const auto n2p_sgn(p_sgn.dot(p_sgn));
        if(T(1) / T(varlen) < n2p_sgn)
          m_sgn.row(s_sgn ++) = p_sgn / sqrt(n2p_sgn);
      }
      return false;
    }
    if(once) {
      invariant_abs.resize(projinv_abs.size() + status_abs.size());
      invariant_sgn.resize(projinv_sgn.size() + status_sgn.size());
      for(int i = 0; i < projinv_abs.size(); i ++) {
        invariant_abs[i] = projinv_abs[i];
        invariant_sgn[i] = projinv_sgn[i];
      }
      for(int i = 0; i < status_abs.size(); i ++) {
        invariant_abs[i + projinv_abs.size()] = m_abs(i, m_abs.cols() - 1);
        invariant_sgn[i + projinv_sgn.size()] = m_sgn(i, m_sgn.cols() - 1);
      }
      once = false;
    }
    status_abs = m_abs * v_abs;
    status_sgn = m_sgn * v_sgn;
    assert(status_abs.size() == status_sgn.size());
    wabs.resize(invariant_abs.size());
    wsgn.resize(invariant_sgn.size());
    for(int i = 0; i < varlen - 1; i ++) {
      wabs[i] =  buf[i - varlen +  buf.size() + 1 - guard];
      wsgn[i] = sbuf[i - varlen + sbuf.size() + 1 - guard];
    }
    wabs[varlen - 1] = T(0);
    wsgn[varlen - 1] = T(0);
    wabs[varlen] = T(1) / sqrt(T(varlen));
    wsgn[varlen] = T(1) / sqrt(T(varlen));
    for(int i = 0; i < status_abs.size(); i ++) {
      wabs[varlen + 1 + i] = status_abs[i];
      wsgn[varlen + 1 + i] = status_sgn[i];
    }
    rabs = abs(invariant_abs.dot(wabs) / invariant_abs[varlen - 1]);
    rsgn = sgn(abs(invariant_sgn.dot(wsgn) / invariant_sgn[varlen - 1]) - T(2));
  }
  return true;
}

template <typename T> inline T P1Istatus<T>::next(const T& in, const int& skip, const int& guard) {
  invariantNext(in, skip, guard);
  return rabs * rsgn;
}

#define _P1_
#endif

