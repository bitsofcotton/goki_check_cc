/*
 BSD 3-Clause License

Copyright (c) 2020-2021, bitsofcotton
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
  inline ~Catg();
  inline void inq(const Vec& in);
  Mat compute();
  Mat R;
private:
  std::vector<Vec> cache;
};

template <typename T> inline Catg<T>::Catg() {
  ;
}

template <typename T> inline Catg<T>::~Catg() {
  ;
}

template <typename T> inline void Catg<T>::inq(const Vec& in) {
  if(cache.size()) assert(in.size() == cache[0].size());
  cache.emplace_back(in);
}

template <typename T> inline typename Catg<T>::Mat Catg<T>::compute() {
  Mat At(cache[0].size(), cache.size());
  for(int i = 0; i < At.cols(); i ++)
    At.setCol(i, cache[i]);
  Mat Q(At.rows(), At.cols());
  for(int i = 0; i < Q.rows(); i ++)
    for(int j = 0; j < Q.cols(); j ++)
      Q(i, j) = T(0);
  for(int i = 0; i < At.rows(); i ++) {
    const auto work(At.row(i) - Q.projectionPt(At.row(i)));
    // generally, assert norm > error is needed.
    // in this case, not.
    Q.row(i) = work / sqrt(work.dot(work));
    if(! isfinite(Q(i, 0))) Q.row(i) *= T(0);
  }
  R = (Q * At.transpose()).transpose();
  return Q;
}


template <typename T> class CatG {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleMatrix<T> Mat;
  inline CatG();
  inline CatG(const int& size);
  inline ~CatG();
  inline void inq(const Vec& in, const bool& computer = false);
  inline void inqRecur(const Vec& in, const bool& computer = false);
  inline void compute(const bool& recur = false);
  inline void computeRecur();
  inline T    lmrS(const Vec& in, const bool& computer = false);
  inline int  lmr(const Vec& in, const bool& computer = false);
  inline std::pair<int, int> lmrRecur(const Vec& in, const bool& computer = false);
  Vec cut;
  T   distance;
  T   origin;
  Catg<T> catg;
  std::vector<Vec> cache;
private:
  const vector<Vec>& tayl(const int& in);
  T threshold_p0;
  int size;
};

template <typename T> inline CatG<T>::CatG() {
#if defined(_FLOAT_BITS_)
  const auto epsilon(T(1) >> int64_t(mybits - 2));
#else
  const auto epsilon(std::numeric_limits<T>::epsilon());
#endif
  threshold_p0 = sqrt(epsilon);
  size = 0;
}

template <typename T> inline CatG<T>::CatG(const int& size) {
#if defined(_FLOAT_BITS_)
  const auto epsilon(T(1) >> int64_t(mybits - 2));
#else
  const auto epsilon(std::numeric_limits<T>::epsilon());
#endif
  threshold_p0 = sqrt(epsilon);
  this->size = size;
}

template <typename T> inline CatG<T>::~CatG() {
  ;
}

template <typename T> const vector<typename CatG<T>::Vec>& CatG<T>::tayl(const int& in) {
  static vector<vector<Vec> > t;
  static P0<T> p;
  if(in < t.size()) {
    if(t[in].size())
      return t[in];
  } else
    t.resize(in + 1, vector<Vec>());
  t[in].reserve(size);
  for(int i = 0; i < size; i ++)
    t[in].emplace_back(p.taylor(in, T(i) * T(in) / T(size)));
  return t[in];
}

template <typename T> inline void CatG<T>::inq(const Vec& in, const bool& computer) {
  Vec work(size);
  if(in.size() == size)
    work = in;
  else {
    const auto& t(tayl(in.size()));
    for(int i = 0; i < work.size(); i ++)
      work[i] = t[i].dot(in);
  }
  if(computer) {
    T pd(0);
    for(int i = 0; i < work.size(); i ++)
      pd += log(abs(work[i]));
    pd = exp(pd / T(work.size()));
    for(int i = 0; i < work.size(); i ++)
      work[i] = pd / work[i];
    cache.emplace_back(work);
  } else
    cache.emplace_back(work);
  assert(isfinite(cache[cache.size() - 1][0]));
  catg.inq(cache[cache.size() - 1]);
  return;
}

template <typename T> inline void CatG<T>::inqRecur(const Vec& in, const bool& computer) {
  auto work(in);
  for(int i = 0; i < in.size(); i ++) {
    inq(work, computer);
    if(i == in.size() - 1) break;
    for(int j = 0; j < work.size(); j ++)
      work[j] = in[(j + i * size / in.size()) % in.size()];
  }
  return;
}

template <typename T> inline void CatG<T>::computeRecur() {
  return compute(true);
}

template <typename T> inline void CatG<T>::compute(const bool& recur) {
  Mat Pt(size + 1, cache.size() * 2);
  Vec q(Pt.cols());
  Vec one(q.size());
  SimpleVector<bool> fix(q.size());
  const auto Q(catg.compute());
  for(int i = 0; i < cache.size(); i ++) {
    const auto pp(Q.col(i));
    for(int j = 0; j < Pt.rows() - 1; j ++) {
      Pt(j, 2 * i)     =   pp[j];
      Pt(j, 2 * i + 1) = - pp[j];
    }
    Pt(Pt.rows() - 1, 2 * i) =
      Pt(Pt.rows() - 1, 2 * i + 1) = T(0);
    q[2 * i]       =   T(1);
    q[2 * i + 1]   = - T(1);
    one[2 * i]     = T(1);
    one[2 * i + 1] = T(1);
    fix[2 * i]     = 0;
    fix[2 * i + 1] = 0;
  }
  one /= sqrt(one.dot(one));
  Pt.row(Pt.rows() - 1)  = q - Pt.projectionPt(q);
  Pt.row(Pt.rows() - 1) /= sqrt(Pt.row(Pt.rows() - 1).dot(Pt.row(Pt.rows() - 1)));
  assert(isfinite(Pt(Pt.rows() - 1, 0)));
  distance = origin = T(0);
  cut      = Vec();
  const auto block(recur ? size * 2 : 2);
  // from bitsofcotton/p1/p1.hh
  int  n_fixed;
  Vec  on;
  Mat  Pverb(Pt);
  if(Pt.cols() == Pt.rows()) {
    cut = Pt * one;
    goto pnext;
  }
  for(n_fixed = 0 ; n_fixed < Pverb.rows(); n_fixed ++) {
    on  = Pverb.projectionPt(- one);
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
        const auto score(on[jj]);
        if(score > mm.first || mm.second < 0)
          mm = std::make_pair(score, jj);
      }
      if(0 <= mm.second && (fidx < 0 ||
          (on[fidx] > on[mm.second] &&
           T(0) <= on[mm.second])))
        fidx = mm.second;
    }
    if(fidx < 0)
      break;
    Vec orth(Pverb.col(fidx));
    const auto norm2orth(orth.dot(orth));
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
    for(int j = 0; j < Pverb.cols(); j ++)
      Pverb.setCol(j, Pverb.col(j) - orth * Pverb.col(j).dot(orth) / norm2orth);
    fix[fidx] = true;
  }
  if(n_fixed == Pt.rows()) {
    int j(0);
    Mat F(Pt.rows(), Pt.rows());
    Vec f(F.rows());
    for(int i = 0; i < Pt.cols() && j < f.size(); i ++)
      if(fix[i]) {
        F.row(j) = Pt.col(i) / sqrt(Pt.col(i).dot(Pt.col(i)));
        f[j]     = one[i];;
        j ++;
      }
    assert(j == f.size());
    try {
      cut = F.solve(f);
    } catch (const char* e) {
      std::cerr << e << std::endl;
    }
  } else
    cut = Pt * on;
 pnext:
  {
    Vec rvec(cut.size() - 1);
    for(int i = 0; i < rvec.size(); i ++)
      rvec[i] = cut[i] - Pt.row(i).dot(one) * cut[cut.size() - 1];
    cut = rvec;
  }
  //assert(isfinite(cut[0]));
  cut  = catg.R * cut;
  //assert(isfinite(cut[0]));
  cut /= sqrt(cut.dot(cut));
  //assert(isfinite(cut[0]));
  std::vector<T> s;
  s.reserve(cache.size());
  for(int i = 0; i < cache.size(); i ++)
    s.emplace_back(cache[i].dot(cut));
  std::sort(s.begin(), s.end());
  for(int i = 0; i < s.size() - 1; i ++)
    if(distance <= s[i + 1] - s[i]) {
      distance =  s[i + 1] - s[i];
      origin   = (s[i + 1] + s[i]) / T(2);
    }
  return;
}

template <typename T> inline T CatG<T>::lmrS(const Vec& in, const bool& computer) {
  Vec work(size);
  if(in.size() == size)
    work = in;
  else {
    const auto& t(tayl(in.size()));
    for(int i = 0; i < work.size(); i ++)
      work[i] = t[i].dot(in);
  }
  if(computer) {
    T pd(0);
    for(int i = 0; i < work.size(); i ++)
      pd += log(abs(work[i]));
    pd = exp(pd / T(work.size()));
    for(int i = 0; i < work.size(); i ++)
      work[i] = pd / work[i];
    return work.dot(cut) - origin;
  }
  return work.dot(cut) - origin;
}

template <typename T> inline int CatG<T>::lmr(const Vec& in, const bool& computer) {
  return lmrS(in, computer) < T(0) ? - 1 : 1;
}

template <typename T> inline std::pair<int, int> CatG<T>::lmrRecur(const Vec& in, const bool& computer) {
  std::pair<int, int> res(make_pair(0, 0));
  T    dM(0);
  auto work(in);
  for(int i = 0; i < in.size(); i ++) {
    const auto score(lmrS(work, computer));
    if(abs(dM) < abs(score)) {
      res = make_pair(score < T(0) ? - 1 : 1, i);
      dM  = abs(score);
    }
    auto tmp(work[0]);
    for(int j = 0; j < work.size() - 1; j ++)
      work[j] = work[j + 1];
    work[work.size() - 1] = tmp;
  }
  return res;
}


template <typename T> std::vector<std::pair<std::vector<std::pair<SimpleVector<T>, int> >, Catg<T> > > crush(const std::vector<SimpleVector<T> >& v, const int& cs, T cut = - T(1) / T(2), const int& Mcount = - 1, const bool& computer = false) {
  std::vector<std::pair<std::vector<std::pair<SimpleVector<T>, int> >, Catg<T> > > result;
  if(! v.size()) return result;
  int t(0);
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
      cat.inq(result[t].first[i].first, computer);
    std::cerr << "." << std::flush;
    cat.compute();
    std::cerr << cat.distance << std::flush;
    if(! t && cut <= T(0))
      cut = cat.distance * abs(cut);
    if((cut <= cat.distance && cat.cut.size()) ||
       (0 < Mcount && Mcount < result[t].first.size())) {
      std::vector<std::pair<SimpleVector<T>, int> > left;
      std::vector<std::pair<SimpleVector<T>, int> > right;
      for(int i = 0; i < result[t].first.size(); i ++)
        (cat.lmr(result[t].first[i].first, computer) < 0 ? left : right).emplace_back(result[t].first[i]);
      if(left.size() && right.size()) {
        CatG<T> lC(cs);
        CatG<T> rC(cs);
        for(int i = 0; i < left.size(); i ++)
          lC.inq(left[i].first, computer);
        for(int i = 0; i < right.size(); i ++)
          rC.inq(right[i].first, computer);
        lC.catg.compute();
        rC.catg.compute();
        result[t] = std::make_pair(std::move(left), std::move(lC.catg));
        result.emplace_back(std::make_pair(std::move(right), std::move(rC.catg)));
      } else {
        result[t].second = cat.catg;
        t ++;
      }
    } else
      t ++;
  }
  return result;
}

template <typename T> std::vector<std::pair<std::vector<std::pair<std::pair<SimpleVector<T>, int>, int> >, Catg<T> > > crushNoContext(const std::vector<SimpleVector<T> >& v, const int& cs, T cut = - T(1) / T(2), const int& Mcount = - 1, const bool& computer = false) {
  std::vector<std::pair<std::vector<std::pair<std::pair<SimpleVector<T>, int>, int> >, Catg<T> > > result;
  if(! v.size()) return result;
  std::vector<std::vector<std::pair<std::pair<SimpleVector<T>, int>, int> > > vv;
  int t(0);
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
      cat.inqRecur(vv[t][i].first.first, computer);
    std::cerr << "." << std::flush;
    cat.computeRecur();
    if(! t && cut <= T(0)) cut = cat.distance * abs(cut);
    std::cerr << cat.distance << std::flush;
    if(cut <= cat.distance || (0 < Mcount && Mcount < vv[t].size())) {
      std::vector<std::pair<std::pair<SimpleVector<T>, int>, int> > left;
      std::vector<std::pair<std::pair<SimpleVector<T>, int>, int> > right;
      for(int i = 0; i < vv[t].size(); i ++) {
        const auto lmr(cat.lmrRecur(vv[t][i].first.first, computer));
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
        vv[t] = std::move(w0);
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
      cg.inq(vv[i][j].first.first, computer);
    cg.compute();
    result.emplace_back(std::make_pair(std::move(vv[i]), std::move(cg.catg)));
  }
  return result;
}


template <typename T> class P012L {
public:
  typedef SimpleVector<T> Vec;
  inline P012L();
  inline P012L(const int& d, const int& stat, const int& slide, const T& intensity = - T(1) / T(2));
  inline ~P012L();
  inline T next(const T& in);
  inline T lastAvg() const;
private:
  std::vector<Vec> cache;
  std::vector<Vec> pp;
  Vec work;
  int stat;
  int slide;
  T   inten;
  int t;
};

template <typename T> inline P012L<T>::P012L() {
  t = stat = slide = 0;
}

template <typename T> inline P012L<T>::P012L(const int& d, const int& stat, const int& slide, const T& intensity) {
  work.resize(d);
  cache.reserve(stat);
  this->stat  = stat;
  this->slide = slide;
  this->inten = intensity;
  t = 0;
}

template <typename T> inline P012L<T>::~P012L() {
  ;
}

template <typename T> inline T P012L<T>::next(const T& in) {
  // XXX:
  static Decompose<T> dec(work.size());
  work[work.size() - 1] = in;
  if(t ++ < work.size() - 1) {
    work[(t - 1) % work.size()] = in;
    return in;
  }
  const auto v(dec.next(work));
  cache.emplace_back(v / sqrt(v.dot(v)));
  for(int i = 0; i < work.size() - 1; i ++)
    work[i] = work[i + 1];
  work[work.size() - 2] = in;
  if(stat <= cache.size()) {
    const auto cat(crush<T>(cache, work.size(), inten, - 1, true));
    pp = std::vector<Vec>();
    pp.reserve(cat.size());
    for(int i = 0; i < cat.size(); i ++) {
      pp.emplace_back(cat[i].first[0].first);
      for(int j = 1; j < cat[i].first.size(); j ++)
        pp[i] += cat[i].first[j].first;
      pp[i] /= sqrt(pp[i].dot(pp[i]));
    }
    const auto cache0(cache);
    cache = std::vector<Vec>();
    cache.reserve(stat);
    for(int i = 0; i < slide; i ++)
      cache.emplace_back(cache0[i - slide + cache0.size()]);
  }
  T MM(0);
  T res(0);
  for(int i = 0; i < pp.size(); i ++) {
    const auto& p(pp[i]);
    const auto  vdp(dec.next(work).dot(p));
    const auto  last(p[p.size() - 1] - p[p.size() - 2]);
    if(! isfinite(vdp)) continue;
    if(MM <= abs(vdp) && last != T(0)) {
      MM  = abs(vdp);
      res = last * vdp;
    }
  }
  return in + res;
}

template <typename T> inline T P012L<T>::lastAvg() const{
  T la(0);
  for(int i = 0; i < cache.size(); i ++)
    la += abs(cache[i][cache[i].size() - 1] - cache[i][cache[i].size() - 2]);
  return la /= T(cache.size());
}

#define _CATG_
#endif

