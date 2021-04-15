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

using std::move;
using std::vector;
using std::pair;
using std::make_pair;
using std::cerr;
using std::flush;

template <typename T> class CatG {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleMatrix<T> Mat;
  inline CatG();
  inline CatG(const int& size, const bool& recur = false, const int& complexity = 8);
  inline ~CatG();
         void compute(const vector<Vec>& in);
  inline pair<T, int> score(const Vec& in);
  Vec cut;
  T   distance;
  T   origin;
  const Mat& tayl(const int& in);
private:
  int  size;
  bool recur;
  int  complexity;
};

template <typename T> inline CatG<T>::CatG() {
  recur = false;
  size  = complexity = 0;
}

template <typename T> inline CatG<T>::CatG(const int& size, const bool& recur, const int& complexity) {
  this->size       = size;
  this->recur      = recur;
  this->complexity = complexity;
}

template <typename T> inline CatG<T>::~CatG() {
  ;
}

template <typename T> const typename CatG<T>::Mat& CatG<T>::tayl(const int& in) {
  static vector<Mat> t;
  if(in < t.size()) {
    if(t[in].rows() && t[in].cols())
      return t[in];
  } else
    t.resize(in + 1, Mat());
  t[in].resize(size, in);
  for(int i = 0; i < size; i ++)
    t[in].row(i) = taylor<T>(in, T(i) * T(in) / T(size));
  return t[in];
}

template <typename T> void CatG<T>::compute(const vector<Vec>& in) {
  const auto block(recur ? size * 2 : 2);
  SimpleMatrix<T> A(in.size() * block - (recur ? size : 1), size + 1);
  for(int i = 0; i < in.size(); i ++) {
    if(recur) {
      Vec inn(in.size());
      for(int k = 0; k < size; k ++) {
        for(int j = 0; j < size; j ++)
          inn[j] = in[i][(j + i * size / in[i].size()) % in[i].size()];
        A.row(i * size * 2 + k) = makeProgramInvariant(inn.size() == size ? inn : tayl(inn.size()) * inn, complexity, - T(1));
      }
    } else
      A.row(i * 2) = makeProgramInvariant(in[i].size() == size ? in[i] : tayl(in[i].size()) * in[i], complexity, - T(1));
    if(in.size() - 1 <= i) break;
    if(recur) {
      for(int k = 0; k < size; k ++)
        A.row(i * size * 2 + size + k) = - A.row(i * size * 2 + k);
    } else
      A.row(i * 2 + 1) = - A.row(i * 2);
  }
        auto Pt(A.QR());
  const auto R(Pt * A);
  Vec  one(Pt.cols());
  SimpleVector<bool> fix(one.size());
  for(int i = 0; i < Pt.cols(); i ++) {
    one[i] = T(1);
    fix[i] = false;
  }
  const auto on(Pt.projectionPt(one));
  vector<pair<T, int> > fidx;
  fidx.reserve(on.size());
  for(int i = 0; i < on.size(); i ++)
    fidx.emplace_back(make_pair(abs(on[i]), i));
  sort(fidx.begin(), fidx.end());
  for(int n_fixed = 0, idx = 0; n_fixed < Pt.rows() - 1 && idx < fidx.size(); n_fixed ++, idx ++) {
    const auto& iidx(fidx[idx].second);
    for(int j = iidx - iidx % block;
            j < min((iidx - iidx % block) + block, fix.size());
            j ++)
      if(fix[j]) {
        fix[iidx] = true;
        break;
      }
    if(fix[iidx]) continue;
    const auto  orth(Pt.col(iidx));
    const auto  n2(orth.dot(orth));
    if(n2 <= Pt.epsilon) continue;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int j = 0; j < Pt.cols(); j ++)
      Pt.setCol(j, Pt.col(j) - orth * Pt.col(j).dot(orth) / n2);
    fix[iidx] = true;
  }
  cut = R.solve(Pt * one);
  std::vector<T> s;
  s.reserve(in.size());
  for(int i = 0; i < in.size(); i ++)
    s.emplace_back(makeProgramInvariant(in[i].size() == size ? in[i] : tayl(in[i].size()) * in[i], complexity, - T(1)).dot(cut));
  std::sort(s.begin(), s.end());
  distance = origin = T(0);
  for(int i = 0; i < s.size() - 1; i ++)
    if(distance <= s[i + 1] - s[i]) {
      distance =  s[i + 1] - s[i];
      origin   = (s[i + 1] + s[i]) / T(2);
    }
  return;
}

template <typename T> inline pair<T, int> CatG<T>::score(const Vec& in) {
  if(! recur)
    return make_pair(makeProgramInvariant<T>(in.size() == size ? in : tayl(in.size()) * in, complexity, - T(1)).dot(cut) - origin, 0);
  pair<T, int> res(make_pair(0, 0));
  for(int i = 0; i < in.size(); i ++) {
    Vec inn(in.size());
    for(int j = 0; j < size; j ++)
      inn[j] = in[(j + i * size / in.size()) % in.size()];
    const auto score(makeProgramInvariant<T>(inn.size() == size ? inn : tayl(inn.size()) * inn, complexity, - T(1)).dot(cut) - origin);
    if(abs(res.first) < abs(score))
      res = make_pair(score, i);
  }
  return res;
}


template <typename T> vector<pair<pair<vector<SimpleVector<T> >, vector<pair<int, int> > >, SimpleMatrix<T> > > crush(const vector<SimpleVector<T> >& v, const int& cs, const bool& recur, T cut = - T(1) / T(2), const int& Mcount = - 1, const int& complexity = 8) {
  vector<pair<pair<vector<SimpleVector<T> >, vector<pair<int, int> > >, SimpleMatrix<T> > > result;
  if(! v.size() || !v[0].size()) return result;
  auto MM(v[0].dot(v[0]));
  for(int i = 1; i < v.size(); i ++)
    MM = max(MM, v[i].dot(v[i]));
  MM = sqrt(MM);
  int t(0);
  result.emplace_back(pair<pair<vector<SimpleVector<T> >, vector<pair<int, int> > >, SimpleMatrix<T>  >());
  result[0].first.first.reserve(v.size());
  result[0].first.second.reserve(v.size());
  for(int i = 0; i < v.size(); i ++) {
    result[0].first.first.emplace_back(v[i] / MM / T(2));
    result[0].first.second.emplace_back(make_pair(0, i));
  }
  while(t < result.size()) {
    if(result[t].first.first.size() < cs + 2) {
      t ++;
      continue;
    }
    CatG<T> catg(cs, recur, complexity);
    catg.compute(result[t].first.first);
    if(! t && cut <= T(0))
      cut = catg.distance * abs(cut);
    if(catg.cut.size() && (cut <= catg.distance ||
       (0 < Mcount && Mcount < result[t].first.first.size())) ) {
      vector<SimpleVector<T> > left;
      vector<SimpleVector<T> > right;
      vector<pair<int, int> >  lidx;
      vector<pair<int, int> >  ridx;
      for(int i = 0; i < result[t].first.first.size(); i ++) {
        const auto score(catg.score(result[t].first.first[i]));
        (score.first < T(0) ? left : right).emplace_back(move(result[t].first.first[i]));
      }
      if(left.size() && right.size()) {
        result[t].first.first  = move(left);
        result[t].first.second = move(lidx);
        result.emplace_back(make_pair(make_pair(move(right), move(ridx)), SimpleMatrix<T>()));
      } else
        result[t ++].first.first = (left.size() ? move(left) : move(right));
    } else
      t ++;
  }
  for(int i = 0; i < result.size(); i ++) {
    if(result[i].first.first.size() < cs + 2) continue;
    SimpleMatrix<T> spec(result[i].first.first.size(), cs);
    CatG<T> catg(cs, recur, complexity);
    for(int j = 0; j < result[i].first.first.size(); j ++) {
      const auto& v(result[i].first.first[j] *= MM);
      spec.row(j) = (v.size() == cs ? v : catg.tayl(v.size()) * v);
    }
    result[i].second = spec.QR() * spec;
  }
  return result;
}


template <typename T, bool dec = true> class P012L {
public:
  typedef SimpleVector<T> Vec;
  inline P012L();
  inline P012L(const int& d, const int& stat, const int& slide, const T& intensity = - T(1) / T(2));
  inline ~P012L();
  T next(const T& in, const int& complexity = 8);
private:
  vector<Vec> cache;
  vector<Vec> pp;
  Vec work;
  int stat;
  int slide;
  T   inten;
  int t;
};

template <typename T, bool dec> inline P012L<T,dec>::P012L() {
  inten = T(t = stat = slide = 0);
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

template <typename T, bool dec> inline T P012L<T,dec>::next(const T& in, const int& complexity) {
  static vector<Decompose<T> > decompose;
  static vector<bool> isinit;
  if(dec && decompose.size() <= work.size()) {
    decompose.resize(work.size() + 1, Decompose<T>());
    isinit.resize(decompose.size(), false);
  }
  if(dec && ! isinit[work.size()]) {
    decompose[work.size()] = Decompose<T>(work.size());
    isinit[work.size()] = true;
  }
  if(t ++ < work.size() - 1) {
    work[(t - 1) % work.size()] = in;
    return in;
  }
  work[work.size() - 1] = in;
  cache.emplace_back(dec ? decompose[work.size()].mother(work) : work);
  for(int i = 0; i < work.size() - 1; i ++)
    work[i] = work[i + 1];
  if(stat <= cache.size()) {
    const auto cat(crush<T>(cache, work.size(), false, inten, work.size() * 4, complexity));
    pp = vector<Vec>();
    pp.reserve(cat.size());
    for(int i = 0; i < cat.size(); i ++) {
      pp.emplace_back(cat[i].first.first[0]);
      for(int j = 1; j < cat[i].first.first.size(); j ++)
        pp[i] += cat[i].first.first[j];
      pp[i] /= sqrt(pp[i].dot(pp[i]));
    }
    auto cache0(cache);
    cache = vector<Vec>();
    cache.reserve(stat);
    for(int i = 0; i < slide; i ++)
      cache.emplace_back(move(cache0[i - slide + cache0.size()]));
  }
  T MM(0);
  T res(0);
  for(int i = 0; i < pp.size(); i ++) {
    const auto& p(pp[i]);
    if(! p.size()) continue;
    const auto  vdp((dec ? decompose[work.size()].mother(work) : work).dot(p));
    const auto  last(p[p.size() - 1] - p[p.size() - 2]);
    if(! isfinite(vdp)) continue;
    if(MM <= abs(vdp)) {
      MM  = abs(vdp);
      res = last * vdp;
    }
  }
  return in + res;
}

#define _CATG_
#endif

