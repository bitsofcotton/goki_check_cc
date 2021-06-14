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
  inline CatG(const int& size0, const vector<Vec>& in, const bool& recur = false);
  inline ~CatG();
  inline pair<T, int> score(const Vec& in);
  Vec  cut;
  T    distance;
  T    origin;
  const Mat& tayl(const int& size, const int& in);
  bool recur;
};

template <typename T> inline CatG<T>::CatG() {
  recur = false;
}

template <typename T> inline CatG<T>::CatG(const int& size0, const vector<Vec>& in, const bool& recur) {
  const auto size(abs(size0));
  this->recur = recur;
  const auto block(recur ? size : 1);
  SimpleMatrix<T> A(in.size() * block, size + 1);
  for(int i = 0; i < in.size(); i ++)
    tayl(size, in[i].size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < in.size(); i ++)
    if(recur) {
      Vec inn(in[i].size());
      for(int k = 0; k < size; k ++) {
        for(int j = 0; j < inn.size(); j ++)
          inn[j] = in[i][(j + k * in[i].size() / size) % in[i].size()];
        A.row(i * size + k) = makeProgramInvariant(inn.size() == size ? inn : tayl(size, inn.size()) * inn);
      }
    } else
      A.row(i) = makeProgramInvariant(in[i].size() == size ? in[i] : tayl(size, in[i].size()) * in[i]);
        auto Pt(A.QR());
  const auto R(Pt * A);
        Vec  one(Pt.cols());
  SimpleVector<bool> fix(one.size());
  one.I(T(1));
  fix.I(false);
  const auto on(Pt.projectionPt(one));
  vector<pair<T, int> > fidx;
  vector<int> pidx;
  fidx.reserve(on.size());
  pidx.resize(on.size(), 0);
  if(0 < size0) {
    for(int i = 0; i < on.size(); i ++)
      fidx.emplace_back(make_pair(abs(on[i]), i));
  } else {
    for(int i = 0; i < on.size(); i ++) {
      T score(0);
      for(int j = 0; j < in.size(); j ++) {
        const auto lscore(abs(on[i] + on[j * block + i % block]));
        if(score == T(0) || lscore < score) {
          score = lscore;
          pidx[i] = j;
        }
      }
      fidx.emplace_back(make_pair(score, i));
    }
  }
  sort(fidx.begin(), fidx.end());
  for(int n_fixed = 0, idx = 0;
          n_fixed < Pt.rows() - 1 && idx < fidx.size();
          n_fixed ++, idx ++) {
    const auto& iidx(fidx[idx].second);
    for(int j = iidx - iidx % block;
            j < min((iidx - iidx % block) + block, fix.size());
            j ++)
      if(fix[j]) {
        for(int k = iidx - iidx % block;
                k < min((iidx - iidx % block) + block, fix.size());
                k ++)
          fix[k] = true;
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
    if(size0 < 0)
      fix[pidx[iidx]] = true;
  }
  cut = R.solve(Pt * one);
  std::vector<T> s;
  s.reserve(in.size());
  for(int i = 0; i < in.size(); i ++)
    s.emplace_back(makeProgramInvariant(in[i].size() == size ? in[i] : tayl(size, in[i].size()) * in[i]).dot(cut));
  std::sort(s.begin(), s.end());
  distance = origin = T(0);
  for(int i = 0; i < s.size() - 1; i ++)
    if(distance >= s[i + 1] - s[i]) {
      distance =  s[i + 1] - s[i];
      origin   = (s[i + 1] + s[i]) / T(2);
    }
  return;
}

template <typename T> inline CatG<T>::~CatG() {
  ;
}

template <typename T> const typename CatG<T>::Mat& CatG<T>::tayl(const int& size, const int& in) {
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

template <typename T> inline pair<T, int> CatG<T>::score(const Vec& in) {
  const auto size(cut.size() - 1);
  assert(0 < size);
  if(! recur)
    return make_pair(makeProgramInvariant<T>(in.size() == size ? in : tayl(size, in.size()) * in).dot(cut) - origin, 0);
  pair<T, int> res(make_pair(0, 0));
  for(int i = 0; i < size; i ++) {
    Vec inn(in.size());
    for(int j = 0; j < inn.size(); j ++)
      inn[j] = in[(j + i * in.size() / size) % in.size()];
    const auto score(makeProgramInvariant<T>(inn.size() == size ? inn : tayl(size, inn.size()) * inn).dot(cut) - origin);
    if(abs(res.first) < abs(score))
      res = make_pair(score, i);
  }
  return res;
}


template <typename T> vector<pair<vector<SimpleVector<T> >, vector<pair<int, int> > > > crush(const vector<SimpleVector<T> >& v, const int& cs, const bool& recur, const int& count) {
  assert(0 <= count);
  vector<pair<vector<SimpleVector<T> >, vector<pair<int, int> > > > result;
  if(! v.size() || !v[0].size()) return result;
  auto MM(v[0].dot(v[0]));
  for(int i = 1; i < v.size(); i ++)
    MM = max(MM, v[i].dot(v[i]));
  MM = sqrt(MM) * T(2);
  int t(0);
  result.emplace_back(pair<vector<SimpleVector<T> >, vector<pair<int, int> > >());
  result[0].first.reserve(v.size());
  result[0].second.reserve(v.size());
  for(int i = 0; i < v.size(); i ++) {
    result[0].first.emplace_back(v[i] / MM);
    result[0].second.emplace_back(make_pair(0, i));
  }
  vector<pair<T, pair<int, bool> > > sidx;
  sidx.emplace_back(make_pair(T(0), make_pair(0, false)));
  while(sidx.size() < (count ? count : v.size())) {
    sort(sidx.begin(), sidx.end());
    int iidx(sidx.size() - 1);
    for( ; - 1 <= iidx; iidx --)
      if(iidx < 0 || (! sidx[iidx].second.second &&
        abs(cs) + 1 < result[sidx[iidx].second.first].first.size()) )
        break;
    if(iidx < 0) break;
    const auto& t(sidx[iidx].second.first);
    CatG<T> catg(cs, result[t].first, recur);
    if(catg.cut.size()) {
      vector<SimpleVector<T> > left;
      vector<SimpleVector<T> > right;
      vector<pair<int, int> >  lidx;
      vector<pair<int, int> >  ridx;
      left.reserve(result[t].first.size());
      right.reserve(result[t].first.size());
      lidx.reserve(result[t].first.size());
      ridx.reserve(result[t].first.size());
      for(int i = 0; i < result[t].first.size(); i ++) {
        const auto score(catg.score(result[t].first[i]));
        (score.first < T(0) ? left : right).emplace_back(move(result[t].first[i]));
        (score.first < T(0) ? lidx : ridx).emplace_back(make_pair(score.second, result[t].second[i].second));
      }
      if((abs(cs) + 1 < left.size() || abs(cs) + 1 < right.size()) && left.size() && right.size()) {
        result[t].first  = move(left);
        result[t].second = move(lidx);
        sidx[iidx].first  = catg.distance;
        sidx[iidx].second = make_pair(t, false);
        result.emplace_back(make_pair(move(right), move(ridx)));
        sidx.emplace_back(make_pair(catg.distance, make_pair(sidx.size(), false)));
      } else {
        result[t].first  = move(left);
        result[t].second = move(lidx);
        result[t].first.reserve(result[t].first.size() + right.size());
        result[t].second.reserve(result[t].second.size() + ridx.size());
        for(int i = 0; i < right.size(); i ++) {
          result[t].first.emplace_back(move(right[i]));
          result[t].second.emplace_back(move(ridx[i]));
        }
        sidx[iidx].first = catg.distance;
        sidx[iidx].second.second = true;
      }
    } else {
      sidx[iidx].first = catg.distance;
      sidx[iidx].second.second = true;
    }
  }
  for(int i = 0; i < result.size(); i ++)
    for(int j = 0; j < result[i].first.size(); j ++)
      result[i].first[j] *= MM;
  return result;
}

template <typename T> static inline vector<pair<vector<SimpleVector<T> >, vector<pair<int, int> > > > crush(const vector<SimpleVector<T> >& v, const int& cs, const bool& recur) {
  return crush<T>(v, cs, recur, max(int(2), int(sqrt(int(v.size())))));
}

template <typename T> vector<pair<vector<SimpleVector<T> >, vector<int> > > crushWithOrder(const vector<T>& v, const int& cs, const int& count) {
  vector<SimpleVector<T> > work;
  vector<int> edge;
  // N.B. it's O(v.size()^3 * cs^2).
  work.reserve((v.size() * v.size() * v.size() - 19 * v.size() + 30) / 6);
  edge.reserve(v.size() - 1);
  edge.emplace_back(0);
  for(int i = 3; i < v.size(); i ++) {
    SimpleVector<T> buf(i);
    for(int j = 0; j <= v.size() - i; j ++) {
      for(int k = j; k < j + i; k ++)
        buf[k - j] = v[k];
      if(count < 0) {
        vector<T> wbuf;
        wbuf.reserve(buf.size());
        for(int k = 0; k < buf.size(); k ++)
          wbuf.emplace_back(buf[k]);
        std::sort(wbuf.begin(), wbuf.end());
        for(int k = 0; k < buf.size(); k ++)
          buf[k] = wbuf[k];
      }
      work.emplace_back(buf);
    }
    edge.emplace_back(work.size());
  }
  auto whole_crush(crush<T>(work, cs, false, count == - 1 ? 0 : abs(count)));
  vector<pair<vector<SimpleVector<T> >, vector<int> > > res;
  res.reserve(whole_crush.size());
  for(int i = 0; i < whole_crush.size(); i ++) {
    vector<int> idx;
    const auto& sec(whole_crush[i].second);
    idx.reserve(sec.size());
    for(int j = 0; j < sec.size(); j ++)
      idx.emplace_back(sec[j].second - *std::lower_bound(edge.begin(), edge.end(), sec[j].second));
    res.emplace_back(make_pair(move(whole_crush[i].first), move(idx)));
  }
  return res;
}

template <typename T> static inline vector<pair<vector<SimpleVector<T> >, vector<int> > > crushWithOrder(const vector<T>& v, const int& cs) {
  return crushWithOrder<T>(v, cs, max(int(2), int(sqrt(int(v.size())))));
}

template <typename T, typename feeder, bool dec = true> class P012L {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleMatrix<T> Mat;
  inline P012L();
  inline P012L(const int& stat, const int& d, const int& comp0 = 0);
  inline ~P012L();
  T next(const T& in);
private:
  vector<Vec> pp;
  feeder f;
  Mat pc;
  T   M;
};

template <typename T, typename feeder, bool dec> inline P012L<T,feeder,dec>::P012L() {
  M = T(0);
}

template <typename T, typename feeder, bool dec> inline P012L<T,feeder,dec>::P012L(const int& stat, const int& d, const int& comp0) {
  M = T(0);
  assert(0 <= comp0);
  auto comp(comp0);
  T    tmp(d + 1);
  while(true) {
    tmp  /= T(2);
    if(tmp < T(1)) break;
    tmp   = ceil(tmp);
    comp += int(tmp);
  }
  pc = Mat(comp, d);
  for(int i = 0; i < pc.rows() - 1; i ++) {
    const auto work(taylor(pc.cols() - 1, T(1) / T(pc.rows() - 2) * T(pc.cols() - 2)));
    for(int j = 0; j < work.size(); j ++)
      pc(i, j) = work[j];
    pc(i, work.size()) = T(0);
  }
  for(int i = 0; i < pc.cols(); i ++)
    pc(pc.rows() - 1, i) = T(i == pc.cols() - 1 ? 1 : 0);
  T MM(0);
  for(int i = 0; i < pc.rows(); i ++)
    for(int j = 0; j < pc.cols(); j ++)
      MM = max(MM, abs(pc(i, j)));
  pc /= MM * T(4);
  f = feeder(stat + comp);
}

template <typename T, typename feeder, bool dec> inline P012L<T,feeder,dec>::~P012L() {
  ;
}

template <typename T, typename feeder, bool dec> inline T P012L<T,feeder,dec>::next(const T& in) {
  if(f.res[f.res.size() - 1] == in) return T(0);
  if(M < abs(in) * T(4) * T(pc.rows()) * atan(T(1)) * T(pc.cols() + 2)) M = abs(in) * T(8 * pc.rows()) * atan(T(1)) * T(pc.cols() + 2);
  const auto d(f.next(in));
  if(! f.full) return T(0);
  vector<SimpleVector<T> > cache;
  cache.reserve(d.size() - pc.cols() + 1);
  Decompose<T> decompose(pc.cols());
  for(int i = 0; i < d.size() - pc.cols() + 1; i ++) {
    auto w(d.subVector(i, pc.cols()) / M);
    cache.emplace_back(pc * (dec ? decompose.mother(move(w)) : move(w)));
  }
  const auto cat(crush<T>(cache, cache[0].size(), false, 0));
  pp = vector<Vec>();
  pp.reserve(cat.size());
  for(int i = 0; i < cat.size(); i ++) {
    if(cat[i].first.size() <= pc.rows()) continue;
    vector<Vec> pw;
    pw.reserve(cat[i].first.size());
    for(int j = 0; j < cat[i].first.size(); j ++)
      pw.emplace_back(makeProgramInvariant<T>(cat[i].first[j]));
    pp.emplace_back(linearInvariant<T>(pw));
  }
  SimpleVector<T> work(pc.cols());
  for(int i = 0; i < work.size() - 1; i ++)
    work[i] = d[i - work.size() + d.size() + 1] / M;
  work[work.size() - 1] = work[work.size() - 2];
  T MM(0);
  T res(0);
  const auto vdp(makeProgramInvariant<T>(
    pc * (dec ? decompose.mother(work) : work)));
  for(int i = 0; i < pp.size(); i ++) {
    const auto& p(pp[i]);
    if(! p.size()) continue;
    const auto vdps((vdp.dot(p) - vdp[pc.rows() - 1] * p[pc.rows() - 1]) / sqrt(vdp.dot(vdp) * p.dot(p)));
    if(! isfinite(vdps)) continue;
    if(MM < abs(vdps) && p[work.size()] != T(0)) {
      const auto v(revertProgramInvariant<T>((p.dot(vdp) - p[pc.rows() - 1] * vdp[pc.rows() - 1]) / p[pc.rows() - 1]));
      if(v != T(0)) {
        MM  = abs(vdps);
        res = v;
      }
    }
  }
  return res * M;
}

#define _CATG_
#endif

