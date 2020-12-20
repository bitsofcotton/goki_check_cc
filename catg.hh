/*
 BSD 3-Clause License

Copyright (c) 2020, bitsofcotton
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
#if !defined(_CATG_)

template <typename T> class Catg {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleMatrix<T> Mat;
  inline Catg();
  inline Catg(const int& size);
  inline ~Catg();
  inline void inq(const Vec& in);
  inline void compute();
  Mat Left;
  Vec lambda;
private:
  Mat roughQR(const Mat& At) const;
  Mat AAt;
};

template <typename T> inline Catg<T>::Catg() {
  ;
}

template <typename T> inline Catg<T>::Catg(const int& size) {
  AAt.resize(size, size);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < AAt.rows(); i ++)
    for(int j = 0; j < AAt.cols(); j ++)
      AAt(i, j) = T(0);
}

template <typename T> inline Catg<T>::~Catg() {
  ;
}

template <typename T> inline void Catg<T>::inq(const Vec& in) {
  assert(AAt.rows() == in.size() && AAt.cols() == in.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < in.size(); i ++)
    AAt.row(i) += in * in[i];
}

template <typename T> inline void Catg<T>::compute() {
  Left  = roughQR(AAt);
  const auto Right(Left.transpose() * AAt);
  lambda.resize(Right.rows());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < lambda.size(); i ++)
    lambda[i] = sqrt(Right.row(i).dot(Right.row(i)));
  return;
}

template <typename T> inline typename Catg<T>::Mat Catg<T>::roughQR(const Mat& At) const {
  Mat Q(At.rows(), At.cols());
  for(int i = 0; i < Q.rows(); i ++)
    for(int j = 0; j < Q.cols(); j ++)
      Q(i, j) = T(0);
  for(int i = 0; i < At.rows(); i ++) {
    const auto work(At.row(i) - Q.projectionPt(At.row(i)));
    // generally, assert norm > error is needed.
    // in this case, not.
    Q.row(i) = work / sqrt(work.dot(work));
  }
  return Q;
}


template <typename T> class CatG {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleMatrix<T> Mat;
  inline CatG();
  inline CatG(const int& size);
  inline ~CatG();
  inline void inq(const Vec& in);
  inline void compute();
  Vec cut;
  T   distance;
  T   origin;
  Catg<T> catg;
  std::vector<Vec> cache;
private:
  T threshold_feas;
  T threshold_p0;
  T threshold_inner;
};

template <typename T> inline CatG<T>::CatG() {
#if defined(_FLOAT_BITS_)
  const auto epsilon(T(1) >> int64_t(mybits - 2));
#else
  const auto epsilon(std::numeric_limits<T>::epsilon());
#endif
  threshold_feas  = pow(epsilon, T(5) / T(6));
  threshold_p0    = pow(epsilon, T(4) / T(6));
  threshold_inner = pow(epsilon, T(2) / T(6));
}

template <typename T> inline CatG<T>::CatG(const int& size) {
#if defined(_FLOAT_BITS_)
  const auto epsilon(T(1) >> int64_t(mybits - 2));
#else
  const auto epsilon(std::numeric_limits<T>::epsilon());
#endif
  threshold_feas  = pow(epsilon, T(5) / T(6));
  threshold_p0    = pow(epsilon, T(4) / T(6));
  threshold_inner = pow(epsilon, T(2) / T(6));
  catg = Catg<T>(size);;
}

template <typename T> inline CatG<T>::~CatG() {
  ;
}

template <typename T> inline void CatG<T>::inq(const Vec& in) {
  if(cache.size())
    assert(cache[0].size() == in.size());
  cache.push_back(in);;
  catg.inq(in);
}

template <typename T> inline void CatG<T>::compute() {
  Mat A(cache[0].size(), cache.size());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < A.cols(); i ++)
    A.setCol(i, cache[i]);
  catg.compute();
  Mat Pt(A.rows(), A.cols() * 2);
  Vec b(Pt.cols());
  Vec one(b.size());
  SimpleVector<bool> fix(b.size());
  {
    auto Pt0(catg.Left.transpose() * A);
    for(int i = 0; i < Pt0.rows(); i ++)
      Pt0.row(i) /= sqrt(Pt0.row(i).dot(Pt0.row(i))) * T(2);
    for(int i = 0; i < Pt0.cols(); i ++) {
      Pt.setCol(2 * i,       Pt0.col(i));
      Pt.setCol(2 * i + 1, - Pt0.col(i));
      b[2 * i]       = T(0);
      b[2 * i + 1]   = T(0);
      one[2 * i]     = T(1);
      one[2 * i + 1] = T(1);
      fix[2 * i]     = 0;
      fix[2 * i + 1] = 0;
    }
  }
  auto checked(fix);
  auto norm(b);
  auto norm2(b);
  Mat  F(Pt.rows(), Pt.rows());
  Vec  f(Pt.rows());
  Mat  Pverb;
  Vec  orth;
  T    lasterr(Pt.rows() + Pt.cols());
  distance = T(0);
  cut      = Vec();
  // from bitsofcotton/p1/p1.hh
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
    if(Pt.cols() == Pt.rows()) {
      rvec = Pt * b;
      goto pnext;
    }
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
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
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
      for(int j = 0; j < Pverb.cols(); j ++) {
        norm[j]    = sqrt(Pverb.col(j).dot(Pverb.col(j)));
        norm2[j]   = (j & 1 ? - norm[j] : norm[j]) + ratio;
        checked[j] = fix[j] || norm[j] <= threshold_p0;
      }
      auto mb(mbb + norm2 * normb0 * ratio);
      mb -= (deltab = Pverb.projectionPt(mb));
      mb /= (ratiob = sqrt(mb.dot(mb)));
      on  = Pverb.projectionPt(- one) + mb * mb.dot(- one);
      int fidx(0);
      for( ; fidx < on.size(); fidx ++)
        if(!checked[fidx])
          break;
      for(int j = (fidx + 1) & ~1; j < on.size() / 2; j ++)
        if(!checked[j] && !checked[j + 1] &&
           (on[fidx] / norm[fidx] < on[j] / norm[j] ||
            on[fidx] / norm[fidx] < on[j + 1] / norm[j + 1]) &&
           T(0) <= on[j] && T(0) <= on[j + 1])
          fidx = j;
      if(fidx >= one.size())
        break;
      on /= abs(mb.dot(on));
      if(on[fidx] * sqrt(norm.dot(norm)) / norm[fidx] <= threshold_inner)
        break;
      orth = Pverb.col(fidx);
      const auto norm2orth(orth.dot(orth));
      const auto mbb0(mbb[fidx]);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
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
          f[j]     = b[i]      / lratio + ratio * ratio;
          j ++;
        }
      assert(j == f.size());
      try {
        rvec = F.solve(f);
      } catch (const char* e) {
        std::cerr << e << std::endl;
        continue;
      }
    } else
      rvec = Pt * (on * ratiob + deltab + b);
   pnext:
    std::vector<T> s;
    T newdist(0);
    T neworigin(0);
    rvec /= sqrt(rvec.dot(rvec));
    for(int i = 0; i < rvec.size(); i ++)
      if(! isfinite(rvec[i]))
        goto next;
      else
        rvec[i] *= catg.lambda[i];
    rvec  = catg.Left * rvec;
    rvec /= sqrt(rvec.dot(rvec));
    s.reserve(cache.size());
    for(int i = 0; i < cache.size(); i ++)
      s.emplace_back(cache[i].dot(rvec));
    std::sort(s.begin(), s.end());
    for(int i = 0; i < s.size() - 1; i ++)
      if(newdist < s[i + 1] - s[i]) {
        newdist   =  s[i + 1] - s[i];
        neworigin = (s[i + 1] + s[i]) / T(2);
      }
    if(distance < newdist) {
      cut      = rvec;
      distance = newdist;
      origin   = neworigin;
    }
   next:
    ;
  }
  return;
}

#define _CATG_
#endif

