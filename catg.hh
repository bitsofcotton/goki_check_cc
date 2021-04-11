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
  T   MM(0), mm(0);
  for(int i = 0; i < At.cols(); i ++) {
    At.setCol(i, cache[i]);
    if(! i) MM = mm = At(0, 0);
    for(int j = 0; j < At.rows(); j ++) {
      MM = max(MM, At(j, i));
      mm = min(mm, At(j, i));
    }
  }
  for(int i = 0; i < At.rows(); i ++)
    for(int j = 0; j < At.cols(); j ++)
      At(i, j) -= mm + (MM - mm) / T(2);
  Mat Q(At.rows(), At.cols());
  for(int i = 0; i < Q.rows(); i ++)
    for(int j = 0; j < Q.cols(); j ++)
      Q(i, j) = T(0);
  for(int i = 0; i < At.rows(); i ++) {
    const auto work(At.row(i) - Q.projectionPt(At.row(i)));
    // generally, assert norm > error is needed.
    // in this case, not.
    Q.row(i) = work / sqrt(work.dot(work));
    if(! isfinite(Q.row(i).dot(Q.row(i)))) {
      std::cerr << "!" << std::flush;
      Q.row(i) *= T(0);
    }
  }
  R = At * Q.transpose();
  return Q;
}


template <typename T> class CatG {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleMatrix<T> Mat;
  inline CatG();
  inline CatG(const int& size);
  inline ~CatG();
  inline void inq(const Vec& in, const int& computer = - 1);
  inline void inqRecur(const Vec& in, const int& computer = - 1);
         void compute(const bool& recur = false);
  inline void computeRecur();
  inline T    lmrS(const Vec& in, const int& computer = - 1);
  inline int  lmr(const Vec& in, const int& computer = - 1);
  inline std::pair<int, int> lmrRecur(const Vec& in, const int& computer = - 1);
  Vec cut;
  T   distance;
  T   origin;
  std::vector<Vec> cache;
  inline Vec normalizeComputer(const Vec& v, const int& computer);
  Catg<T> catg;
private:
  const Mat& tayl(const int& in);
  int size;
};

template <typename T> inline CatG<T>::CatG() {
  size = 0;
}

template <typename T> inline CatG<T>::CatG(const int& size) {
  this->size = size;
}

template <typename T> inline CatG<T>::~CatG() {
  ;
}

template <typename T> const typename CatG<T>::Mat& CatG<T>::tayl(const int& in) {
  static vector<Mat> t;
  static P0<T> p;
  if(in < t.size()) {
    if(t[in].rows() && t[in].cols())
      return t[in];
  } else
    t.resize(in + 1, Mat());
  t[in].resize(size, in);
  for(int i = 0; i < size; i ++)
    t[in].row(i) = p.taylor(in, T(i) * T(in) / T(size));
  return t[in];
}

template <typename T> inline typename CatG<T>::Vec CatG<T>::normalizeComputer(const Vec& v, const int& computer) {
  if(computer < 0)
    return v;
  T pd(0);
  for(int i = 0; i < v.size(); i ++)
    pd += log(abs(v[i]));
  const auto res(v * exp((computer ? T(computer) * pd : - pd) / T(v.size())));
  assert(isfinite(res.dot(res)) && res.dot(res) != T(0));
  return res;
}

template <typename T> inline void CatG<T>::inq(const Vec& in, const int& computer) {
  cache.emplace_back(
    normalizeComputer(in.size() == size ? in :
      tayl(in.size()) * in, computer));
  catg.inq(cache[cache.size() - 1]);
  return;
}

template <typename T> inline void CatG<T>::inqRecur(const Vec& in, const int& computer) {
  auto work(in);
  for(int i = 1; i <= in.size(); i ++) {
    inq(work, computer);
    if(i == in.size()) break;
    for(int j = 0; j < work.size(); j ++)
      work[j] = in[(j + i * size / in.size()) % in.size()];
  }
  return;
}

template <typename T> inline void CatG<T>::computeRecur() {
  return compute(true);
}

template <typename T> void CatG<T>::compute(const bool& recur) {
  const auto block(recur ? size * 2 : 2);
  const auto Q(catg.compute());
        Mat  Pt(Q.rows() + 1, Q.cols() * 2 - 1);
        Vec  q(Pt.cols());
        Vec  one(q.size());
  SimpleVector<bool> fix(q.size());
  for(int i = 0; i < Pt.cols(); i ++) {
    const auto qq(Q.col(i / 2));
    for(int j = 0; j < Q.rows(); j ++)
      Pt(j, i) = qq[j];
    Pt(Q.rows(), i) = T(0);
    q[i]   = i & 1 ? - T(1) : T(1);
    one[i] = T(1);
    fix[i] = false;
  }
  one /= sqrt(one.dot(one));
  Pt.row(Pt.rows() - 1)  = q - Pt.projectionPt(q);
  Pt.row(Pt.rows() - 1) /= sqrt(Pt.row(Pt.rows() - 1).dot(Pt.row(Pt.rows() - 1)));
  auto Ptt(Pt);
  int  n_fixed;
  int  ntry(0);
  if(! isfinite(Pt.row(Pt.rows() - 1).dot(Pt.row(Pt.rows() - 1)))) {
    ntry ++;
    Pt.resize(Ptt.rows() - 1, Ptt.cols());
    for(int i = 0; i < Pt.rows(); i ++)
      Pt.row(i) = std::move(Ptt.row(i));
  }
 retry:
  for(n_fixed = 0; n_fixed < Pt.rows() - 1; n_fixed ++) {
    const auto on(Pt.projectionPt(- one));
          int  fidx(- 1);
    for(int i = 0; i < Pt.cols() / block; i ++) {
      std::pair<T, int> mm;
      mm = std::make_pair(T(0), - 1);
      for(int j = 0; j < block && i * block + j < Pt.cols(); j ++) {
        const auto jj(i * block + j);
        if(fix[jj]) {
          // no matter other dimensions if one of them is fixed.
          mm.second = - 1;
          break;
        }
        const auto& score(on[jj]);
        if(mm.first < score)
          mm = std::make_pair(score, jj);
      }
      if(fidx < 0 || (0 <= mm.second && on[mm.second] < on[fidx]))
        fidx = mm.second;
    }
    if(fidx < 0)
      break;
    const auto orth(Pt.col(fidx));
    const auto norm2orth(orth.dot(orth));
    if(norm2orth <= T(0)) break;
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
    for(int j = 0; j < Pt.cols(); j ++)
      Pt.setCol(j, Pt.col(j) - orth * Pt.col(j).dot(orth) / norm2orth);
    fix[fidx] = true;
  }
  if(n_fixed) {
    cut = - Pt * one;
    if(catg.R.cols() != cut.size()) {
      Vec rvec(cut.size() - 1);
      for(int i = 0; i < rvec.size(); i ++)
        rvec[i] = cut[i] + Pt.row(i).dot(one) * cut[cut.size() - 1];
      cut = catg.R * rvec;
    } else
      cut = catg.R * cut;
  } else if(ntry) {
    cut = Vec(catg.R.cols());
    for(int i = 0; i < cut.size(); i ++)
      cut[i] = T(0);
    distance = origin = T(0);
    std::cerr << "?" << std::flush;
    return;
  } else {
    ntry ++;
    Pt.resize(Ptt.rows() - 1, Ptt.cols());
    for(int i = 0; i < Pt.rows(); i ++)
      Pt.row(i) = std::move(Ptt.row(i));
    goto retry;
  }
  const auto ncut(cut.dot(cut));
  if(isfinite(ncut) && T(0) < ncut)
    cut /= sqrt(ncut);
  std::vector<T> s;
  s.reserve(cache.size());
  for(int i = 0; i < cache.size(); i ++)
    s.emplace_back(cache[i].dot(cut));
  std::sort(s.begin(), s.end());
  distance = origin = T(0);
  for(int i = 0; i < s.size() - 1; i ++)
    if(distance <= s[i + 1] - s[i]) {
      distance =  s[i + 1] - s[i];
      origin   = (s[i + 1] + s[i]) / T(2);
    }
  return;
}

template <typename T> inline T CatG<T>::lmrS(const Vec& in, const int& computer) {
  return normalizeComputer(in.size() == size ? in :
           tayl(in.size()) * in, computer).dot(cut) - origin;
}

template <typename T> inline int CatG<T>::lmr(const Vec& in, const int& computer) {
  return lmrS(in, computer) < T(0) ? - 1 : 1;
}

template <typename T> inline std::pair<int, int> CatG<T>::lmrRecur(const Vec& in, const int& computer) {
  std::pair<int, int> res(make_pair(0, 0));
  T    dM(0);
  auto work(in);
  for(int i = 1; i <= in.size(); i ++) {
    const auto score(lmrS(work, computer));
    if(abs(dM) < abs(score)) {
      res = make_pair(score < T(0) ? - 1 : 1, i - 1);
      dM  = abs(score);
    }
    if(i == in.size()) break;
    for(int j = 0; j < work.size(); j ++)
      work[j] = in[(j + i * size / in.size()) % in.size()];
  }
  return res;
}


template <typename T> std::vector<std::pair<std::vector<std::pair<SimpleVector<T>, int> >, Catg<T> > > crush(const std::vector<SimpleVector<T> >& v, const int& cs, T cut = - T(1) / T(2), const int& Mcount = - 1, const int& computer = 20) {
  std::vector<std::pair<std::vector<std::pair<SimpleVector<T>, int> >, Catg<T> > > result;
  if(! v.size() || !v[0].size()) return result;
  auto MM(v[0].dot(v[0]));
  for(int i = 1; i < v.size(); i ++)
    MM = max(MM, v[i].dot(v[i]));
  MM = sqrt(MM);
  int t(0);
  result.emplace_back(std::pair<std::vector<std::pair<SimpleVector<T>, int> >, Catg<T> >());
  result[0].first.reserve(v.size());
  for(int i = 0; i < v.size(); i ++)
    result[0].first.emplace_back(std::make_pair(v[i] / MM, i));
  while(t < result.size()) {
    if(! result[t].first.size()) {
      t ++;
      continue;
    }
    CatG<T> cat(cs);
    for(int i = 0; i < result[t].first.size(); i ++)
      cat.inq(result[t].first[i].first, computer);
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
        result[t].second = std::move(cat.catg);
        t ++;
      }
    } else
      t ++;
  }
  for(int i = 0; i < result.size(); i ++)
    for(int j = 0; j < result[i].first.size(); j ++)
      result[i].first[j].first *= MM;
  return result;
}

template <typename T> std::vector<std::pair<std::vector<std::pair<std::pair<SimpleVector<T>, int>, int> >, Catg<T> > > crushNoContext(const std::vector<SimpleVector<T> >& v, const int& cs, T cut = - T(1) / T(2), const int& Mcount = - 1, const int& computer = 20) {
  std::vector<std::pair<std::vector<std::pair<std::pair<SimpleVector<T>, int>, int> >, Catg<T> > > result;
  if(! v.size() || ! v[0].size()) return result;
  auto MM(v[0].dot(v[0]));
  for(int i = 1; i < v.size(); i ++)
    MM = max(MM, v[i].dot(v[i]));
  MM = sqrt(MM);
  std::vector<std::vector<std::pair<std::pair<SimpleVector<T>, int>, int> > > vv;
  int t(0);
  vv.emplace_back(std::vector<std::pair<std::pair<SimpleVector<T>, int>, int> >());
  vv[0].reserve(v.size());
  for(int i = 0; i < v.size(); i ++)
    vv[0].emplace_back(std::make_pair(std::make_pair(v[i] / MM, 0), i));
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
    cat.computeRecur();
    std::cerr << cat.distance << std::flush;
    if(! t && cut <= T(0)) cut = cat.distance * abs(cut);
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
      auto lG(crush<T>(ll, cs, cut, Mcount, computer));
      auto rG(crush<T>(rr, cs, cut, Mcount, computer));
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
              left[i][cache[i].first[j].second].second :
              right[i - lG.size()][cache[i].first[j].second].second);
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
  for(int i = 0; i < result.size(); i ++)
    for(int j = 0; j < result[i].first.size(); j ++)
      result[i].first[j].first.first *= MM;
  return result;
}


template <typename T, bool dec = true> class P012L {
public:
  typedef SimpleVector<T> Vec;
  inline P012L();
  inline P012L(const int& d, const int& stat, const int& slide, const T& intensity = - T(1) / T(2));
  inline ~P012L();
  T next(const T& in, const int& computer = 20);
private:
  std::vector<Vec> cache;
  std::vector<Vec> pp;
  Vec work;
  int stat;
  int slide;
  T   inten;
  int t;
};

template <typename T, bool dec> inline P012L<T,dec>::P012L() {
  t = stat = slide = 0;
}

template <typename T, bool dec> inline P012L<T,dec>::P012L(const int& d, const int& stat, const int& slide, const T& intensity) {
  work.resize(d);
  cache.reserve(stat);
  this->stat  = stat;
  this->slide = slide;
  this->inten = intensity;
  t = 0;
}

template <typename T, bool dec> inline P012L<T,dec>::~P012L() {
  ;
}

template <typename T, bool dec> inline T P012L<T,dec>::next(const T& in, const int& computer) {
  static std::vector<Decompose<T> > decompose;
  static std::vector<bool> isinit;
  if(dec && decompose.size() <= work.size()) {
    decompose.resize(work.size() + 1, Decompose<T>());
    isinit.resize(decompose.size(), false);
  }
  if(dec && ! isinit[work.size()]) {
    decompose[work.size()] = Decompose<T>(work.size());
    isinit[work.size()] = true;
  }
  work[work.size() - 1] = in;
  if(t ++ < work.size() - 1) {
    work[(t - 1) % work.size()] = in;
    return in;
  }
  cache.emplace_back(dec ? decompose[work.size()].mother(work) : work);
  for(int i = 0; i < work.size() - 1; i ++)
    work[i] = work[i + 1];
  if(stat <= cache.size()) {
    const auto cat(crush<T>(cache, work.size(), inten, - 1, computer));
    pp = std::vector<Vec>();
    pp.reserve(cat.size());
    static CatG<T> ncCaller;
    for(int i = 0; i < cat.size(); i ++) {
      pp.emplace_back(ncCaller.normalizeComputer(cat[i].first[0].first, computer));
      for(int j = 1; j < cat[i].first.size(); j ++)
        pp[i] += ncCaller.normalizeComputer(cat[i].first[j].first, computer);
      pp[i] /= sqrt(pp[i].dot(pp[i]));
    }
    auto cache0(cache);
    cache = std::vector<Vec>();
    cache.reserve(stat);
    for(int i = 0; i < slide; i ++)
      cache.emplace_back(std::move(cache0[i - slide + cache0.size()]));
  }
  T MM(0);
  T res(0);
  for(int i = 0; i < pp.size(); i ++) {
    const auto& p(pp[i]);
    const auto  vdp((dec ? decompose[work.size()].mother(work) : work).dot(p));
    const auto  last(p[p.size() - 1] - p[p.size() - 2]);
    if(! isfinite(vdp)) continue;
    if(MM <= abs(vdp) && last != T(0)) {
      MM  = abs(vdp);
      res = last * vdp;
    }
  }
  return in + res;
}

#define _CATG_
#endif

