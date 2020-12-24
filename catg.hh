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
  Mat R;
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
  Left = roughQR(AAt);
  R    = Left.transpose() * AAt;
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
  inline void inqRecur(const Vec& in);
  inline void compute(const bool& recur = false);
  inline void computeRecur();
  inline T    lmrS(const Vec& in);
  inline int  lmr(const Vec& in);
  inline std::pair<int, int> lmrRecur(const Vec& in);
  Vec cut;
  T   distance;
  T   origin;
  Catg<T> catg;
  std::vector<Vec> cache;
private:
  const vector<Vec>& tayl(const int& in);
  T threshold_feas;
  T threshold_p0;
  T threshold_inner;
  int size;
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
  size = 0;
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
  catg = Catg<T>(this->size = size);;
}

template <typename T> inline CatG<T>::~CatG() {
  ;
}

template <typename T> const vector<typename CatG<T>::Vec>& CatG<T>::tayl(const int& in) {
  static vector<vector<Vec> > t;
  static P0<T> p;
  if(in < t.size() && t[in].size())
    return t[in];
  t.resize(in + 1, vector<Vec>());
  t[in].reserve(size);
  for(int i = 0; i < size; i ++)
    t[in].emplace_back(p.taylor(in, T(i) * T(in) / T(size)));
  return t[in];
}

template <typename T> inline void CatG<T>::inq(const Vec& in) {
  if(in.size() == size) {
    cache.push_back(in);
    catg.inq(in);
  } else {
    const auto& t(tayl(in.size()));
    Vec work(size);
    for(int i = 0; i < work.size(); i ++)
      work[i] = t[i].dot(in);
    cache.push_back(work);
    catg.inq(work);
  }
  return;
}

template <typename T> inline void CatG<T>::inqRecur(const Vec& in) {
  auto work(in);
  for(int i = 0; i < in.size(); i ++) {
    inq(work);
    if(i == in.size() - 1) break;
    auto tmp(work[0]);
    for(int j = 0; j < work.size() - 1; j ++)
      work[j] = work[j + 1];
    work[work.size() - 1] = tmp;
  }
  return;
}

template <typename T> inline void CatG<T>::computeRecur() {
  return compute(true);
}

template <typename T> inline void CatG<T>::compute(const bool& recur) {
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
    for(int i = 0; i < Pt0.cols(); i ++) {
      const auto pp(catg.R.solve(Pt0.col(i)));
      Pt.setCol(2 * i,       pp);
      Pt.setCol(2 * i + 1, - pp);
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
  const auto block(recur ? size * 2 : 2);
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
        norm2[j]   = (j & 1 ? - ratio : ratio) + norm[j];
        checked[j] = fix[j] || norm[j] <= threshold_p0;
      }
      auto mb(mbb + norm2 * normb0 * ratio);
      mb -= (deltab = Pverb.projectionPt(mb));
      mb /= (ratiob = sqrt(mb.dot(mb)));
      on  = Pverb.projectionPt(- one) + mb * mb.dot(- one);
      int fidx(- 1);
      for(int i = 0; i < Pverb.cols() / block; i ++) {
        std::pair<T, int> mm;
        mm = std::make_pair(T(0), - 1);
        for(int j = 0; j < block; j ++) {
          const auto jj(i * block + j);
          if(fix[jj]) {
            // no matter other dimensions if one of them is fixed.
            mm.second = - 1;
            break;
          }
          if(checked[jj])
            continue;
          const auto score(on[jj] / norm[jj]);
          if(score < mm.first || mm.second < 0)
            mm = std::make_pair(score, jj);
        }
        if(0 <= mm.second && (fidx < 0 ||
            (on[fidx] / norm[fidx] < on[mm.second] / norm[mm.second] &&
             T(0) <= on[mm.second])))
          fidx = mm.second;
      }
      on /= abs(mb.dot(on));
      if(fidx < 0 ||
         on[fidx] * sqrt(norm.dot(norm)) / norm[fidx] <= threshold_inner)
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
          f[j]     = (b[i] + (i & 1 ? - ratio : ratio)) / lratio * ratio;
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
    SimpleVector<T> err0(Pt.cols());
    for(int i = 0; i < err0.size(); i ++)
      err0[i] = Pt.col(i).dot(rvec);
    auto err(err0 - b - one * ratio * ratio);
    for(int i = 0; i < b.size() / block; i ++) {
      auto work(err[i * block] + ratio * ratio);
      for(int j = 1; j < block; j ++) {
        const auto jj(i * block + j);
        work = min(work, err[jj] + (jj & 1 ? - ratio * ratio : ratio * ratio));
      }
      for(int j = 0; j < block; j ++) {
        const auto jj(i * block + j);
        if(work <= T(0) || err[jj] + ratio * ratio < work)
          err[jj] = T(0);
      }
    }
    if(! (sqrt(err.dot(err)) <= sqrt(threshold_inner * err0.dot(err0)) && T(0) < rvec.dot(rvec)))
      continue;
    std::vector<T> s;
    T newdist(0);
    T neworigin(0);
    rvec  = catg.R * rvec;
    rvec /= sqrt(rvec.dot(rvec));
    for(int i = 0; i < rvec.size(); i ++)
      if(! isfinite(rvec[i]))
        goto next;
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
    lasterr -= ratio0;
   next:
    ;
  }
  return;
}

template <typename T> inline T CatG<T>::lmrS(const Vec& in) {
  Vec work(size);
  if(in.size() == size)
    work = in;
  else {
    const auto& t(tayl(in.size()));
    for(int i = 0; i < work.size(); i ++)
      work[i] = t[i].dot(in);
  }
  return work.dot(cut) - origin;
}

template <typename T> inline int CatG<T>::lmr(const Vec& in) {
  return lmrS(in) < T(0) ? - 1 : 1;
}

template <typename T> inline std::pair<int, int> CatG<T>::lmrRecur(const Vec& in) {
  std::pair<int, int> res(make_pair(0, 0));
  T    dM(0);
  auto work(in);
  for(int i = 0; i < in.size(); i ++) {
    const auto score(lmrS(work));
    if(abs(dM) < abs(score)) {
      res = make_pair(score < 0 ? - 1 : 1, i);
      dM  = abs(score);
    }
    auto tmp(work[0]);
    for(int j = 0; j < work.size() - 1; j ++)
      work[j] = work[j + 1];
    work[work.size() - 1] = tmp;
  }
  return res;
}


template <typename T> std::vector<std::pair<std::vector<std::pair<SimpleVector<T>, int> >, Catg<T> > > crush(const std::vector<SimpleVector<T> >& v, const int& cs) {
  std::vector<std::pair<std::vector<std::pair<SimpleVector<T>, int> >, Catg<T> > > result;
  if(! v.size()) return result;
  int t(0);
  T   Mdist(0);
  result.emplace_back(std::pair<std::vector<std::pair<SimpleVector<T>, int> >, Catg<T> >());
  result[0].first.reserve(v.size());
  for(int i = 0; i < v.size(); i ++)
    result[0].first.emplace_back(std::make_pair(v[i], i));
  while(t < result.size()) {
    if(! result[t].first.size()) {
      t ++;
      continue;
    }
    CatG<T> cat(cs);
    for(int i = 0; i < result[t].first.size(); i ++)
      cat.inq(result[t].first[i].first);
    std::cerr << "." << std::flush;
    cat.compute();
    std::cerr << cat.distance << std::flush;
    if(! t && Mdist == T(0))
      Mdist = cat.distance;
    if(Mdist <= cat.distance && cat.cut.size()) {
      std::vector<std::pair<SimpleVector<T>, int> > left;
      std::vector<std::pair<SimpleVector<T>, int> > right;
      for(int i = 0; i < result[t].first.size(); i ++)
        (cat.lmr(result[t].first[i].first) < 0 ? left : right).emplace_back(result[t].first[i]);
      if(left.size() && right.size()) {
        CatG<T> lC(cs);
        CatG<T> rC(cs);
        for(int i = 0; i < left.size(); i ++)
          lC.inq(left[i].first);
        for(int i = 0; i < right.size(); i ++)
          rC.inq(right[i].first);
        lC.catg.compute();
        rC.catg.compute();
        result[t] = std::make_pair(std::move(left), std::move(lC.catg));
        result.emplace_back(std::make_pair(std::move(right), std::move(rC.catg)));
      } else
        t ++;
    } else
      t ++;
  }
  return result;
}

template <typename T> std::vector<std::pair<std::vector<std::pair<std::pair<SimpleVector<T>, int>, int> >, Catg<T> > > crushNoContext(const std::vector<SimpleVector<T> >& v, const int& cs) {
  std::vector<std::pair<std::vector<std::pair<std::pair<SimpleVector<T>, int>, int> >, Catg<T> > > result;
  if(! v.size()) return result;
  std::vector<std::vector<std::pair<std::pair<SimpleVector<T>, int>, int> > > vv;
  int t(0);
  T   Mdist(0);
  vv.emplace_back(std::vector<std::pair<std::pair<SimpleVector<T>, int>, int> >());
  vv[0].reserve(v.size());
  for(int i = 0; i < v.size(); i ++)
    vv[0].emplace_back(std::make_pair(std::make_pair(v[i], 0), i));
  std::vector<int> patch;
  patch.resize(v.size(), 0);
  while(t < vv.size()) {
    if(! vv[t].size()) {
      t ++;
      continue;
    }
    CatG<T> cat(cs);
    for(int i = 0; i < vv[t].size(); i ++)
      cat.inqRecur(vv[t][i].first.first);
    std::cerr << "." << std::flush;
    cat.computeRecur();
    if(! t && Mdist == T(0)) Mdist = cat.distance;
    std::cerr << cat.distance << std::flush;
    if(Mdist <= cat.distance) {
      std::vector<std::pair<std::pair<SimpleVector<T>, int>, int> > left;
      std::vector<std::pair<std::pair<SimpleVector<T>, int>, int> > right;
      for(int i = 0; i < vv[t].size(); i ++) {
        const auto lmr(cat.lmrRecur(vv[t][i].first.first));
        patch[vv[t][i].second] = vv[t][i].first.second = lmr.second;
        (lmr.first < 0 ? left : right).emplace_back(vv[t][i]);
      }
      std::vector<std::pair<std::vector<std::pair<SimpleVector<T>, int> >, Catg<T> > > cache;
      std::vector<SimpleVector<T> > ll;
      std::vector<SimpleVector<T> > rr;
      ll.reserve(left.size());
      rr.reserve(right.size());
      for(int i = 0; i < left.size(); i ++)
        ll.emplace_back(left[i].first.first);
      for(int i = 0; i < right.size(); i ++)
        rr.emplace_back(right[i].first.first);
      auto lG(crush<T>(ll, cs));
      auto rG(crush<T>(rr, cs));
      cache.reserve(lG.size() + rG.size());
      for(int i = 0; i < lG.size(); i ++)
        cache.emplace_back(std::move(lG[i]));
      for(int i = 0; i < rG.size(); i ++)
        cache.emplace_back(std::move(rG[i]));
      if(cache.size()) {
        std::vector<std::pair<std::pair<SimpleVector<T>, int>, int> > w0;
        for(int i = 0; i < cache.size(); i ++) {
          std::vector<std::pair<std::pair<SimpleVector<T>, int>, int> > work;
          work.reserve(cache[i].first.size());
          for(int j = 0; j < cache[i].first.size(); j ++) {
            const auto& idx(i < lG.size() ?
              left[cache[i].first[j].second].second :
              right[cache[i].first[j].second].second);
            work.emplace_back(std::make_pair(std::make_pair(std::move(cache[i].first[j].first), patch[idx]), idx));
          }
          if(! i) w0 = std::move(work);
          else vv.emplace_back(std::move(work));
        }
        vv[t] = w0;
      }
      if(cache.size() <= 1)
        t ++;
    } else
      t ++;
  }
  result.reserve(vv.size());
  for(int i = 0; i < vv.size(); i ++) {
    CatG<T> cg(cs);
    for(int j = 0; j < vv[i].size(); j ++)
      cg.inq(vv[i][j].first.first);
    cg.compute();
    result.emplace_back(std::make_pair(std::move(vv[i]), std::move(cg.catg)));
  }
  return result;
}

#define _CATG_
#endif

