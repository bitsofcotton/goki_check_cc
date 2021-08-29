/* BSD 3-Clause License:
 * Copyright (c) 2018-2021, bitsofcotton.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 *    Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *    Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation or other materials provided with the distribution.
 *    Neither the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * XXX: this could be rewrited by catg.hh .
 */
#if !defined(_MATCH_)

using std::min;
using std::max;
using std::cerr;
using std::endl;
using std::flush;
using std::vector;
using std::pair;
using std::make_pair;
using std::lower_bound;
using std::upper_bound;
using std::distance;
using std::sort;
using std::unique;
using std::istream;
using std::ostream;

template <typename T> class match_t {
public:
  typedef SimpleMatrix<T>   Mat;
  typedef SimpleVector<T>   Vec;
  typedef SimpleVector<int> Veci;
  Mat         rot;
  Vec         offset;
  T           ratio;
  T           rdepth;
  vector<int> dstpoints;
  vector<int> srcpoints;
  T           thresh;
  T           othresh;
  inline match_t() {
    thresh  = T(0);
    othresh = T(0);
    initId();
  }
  inline match_t(const T& thresh, const T& d) {
    this->thresh  = thresh;
    this->othresh = thresh * d;
    initId();
  }
  inline match_t(const match_t<T>& src) {
    *this = src;
  }
  inline void initId() {
    rot       = Mat(3, 3);
    rot(0, 0) = rot(1, 1) = rot(2, 2) = T(1);
    rot(1, 0) = rot(2, 0) = rot(0, 1) = rot(2, 1)
              = rot(0, 2) = rot(1, 2) = T(0);
    offset    = Vec(3);
    offset[0] = offset[1] = offset[2] = T(0);
    ratio     = T(1);
    rdepth    = T(0);
  }
  inline match_t<T>  operator ~ () const {
    match_t<T> result(*this);
    result.rot    = rot.transpose();
    result.ratio  = T(1) / ratio;
    result.offset = - result.rot * offset * result.ratio;
    std::swap(result.srcpoints, result.dstpoints);
    return result;
  }
  inline match_t<T>  operator / (const match_t<T>& src) const {
    match_t<T> result;
    result.rot    = rot   * src.rot.transpose();
    result.ratio  = ratio / src.ratio;
    result.offset = offset - result.rot * src.offset * result.ratio;
    result.rdepth = rdepth + src.rdepth;
    result.dstpoints = vector<int>();
    result.srcpoints = vector<int>();
    for(int i = 0; i < srcpoints.size(); i ++)
      for(int j = 0; j < src.srcpoints.size(); j ++)
        if(srcpoints[i] == src.srcpoints[j]) {
          result.dstpoints.emplace_back(    dstpoints[i]);
          result.srcpoints.emplace_back(src.dstpoints[j]);
        }
    return result;
  }
  inline match_t<T>& operator = (const match_t<T>& other) {
    rot        = other.rot;
    offset     = other.offset;
    ratio      = other.ratio;
    rdepth     = other.rdepth;
    dstpoints  = other.dstpoints;
    srcpoints  = other.srcpoints;
    thresh     = other.thresh;
    othresh    = other.othresh;
    return *this;
  }
  inline T distance(const match_t<T>& other, const Vec& p) {
    const auto d(transform(p) - other.transform(p));
    return sqrt(d.dot(d));
  }
  inline vector<Veci> hullConv(const vector<Veci>& srchull) const {
    assert(srcpoints.size() == dstpoints.size());
    vector<Veci> res;
    res.reserve(srchull.size());
    for(int i = 0; i < srchull.size(); i ++) {
      Veci tmp(3);
      tmp[0] = tmp[1] = tmp[2] = - 1;
      for(int j = 0; j < srchull[i].size(); j ++)
        for(int k = 0; k < srcpoints.size(); k ++)
          if(srcpoints[k] == srchull[i][j]) {
            tmp[j] = dstpoints[k];
            break;
          }
      res.emplace_back(tmp);
    }
    return res;
  }
  inline Vec transform(const Vec& x) const {
    return rot * x * ratio + offset;
  }
  inline vector<Vec> transform(const vector<Vec>& x) const {
    vector<Vec> result(x);
    for(int i = 0; i < result.size(); i ++)
      result[i] = transform(result[i]);
    return result;
  }
  inline bool operator < (const match_t<T>& x1) const {
    const T rratio(max(abs(   ratio), T(1) / abs(   ratio)));
    const T xratio(max(abs(x1.ratio), T(1) / abs(x1.ratio)));
    return rdepth < x1.rdepth || (rdepth == x1.rdepth && rratio < xratio);
  }
  inline bool operator != (const match_t<T>& x) const {
    const auto test(offset - x.offset);
    const auto roterr(rot * x.rot.transpose());
    return !(abs(T(1) - roterr(0, 0)) <= thresh) ||
           !(abs(T(1) - roterr(1, 1)) <= thresh) ||
           !(abs(T(1) - roterr(2, 2)) <= thresh) ||
           !(sqrt(test.dot(test) / 
               (offset.dot(offset) + x.offset.dot(x.offset))) <= othresh) ||
           ratio * x.ratio < T(0) ||
           !(abs(ratio - x.ratio) / sqrt(ratio * x.ratio) <= thresh);
  }
  inline bool operator == (const match_t<T>& x) const {
    return ! (*this != x);
  }
  friend ostream& operator << (ostream& os, const match_t<T>& x) {
    for(int i = 0; i < x.rot.rows(); i ++) {
      for(int j = 0; j < x.rot.cols(); j ++)
        os << x.rot(i, j) << " ";
      os << endl;
    }
    for(int i = 0; i < x.offset.size(); i ++)
      os << x.offset[i] << endl;
    os << x.ratio  << endl;
    os << x.rdepth << endl;
    assert(x.dstpoints.size() == x.srcpoints.size());
    os << x.dstpoints.size() << endl;
    for(int i = 0; i < x.dstpoints.size(); i ++)
      os << x.dstpoints[i] << " ";
    os << endl;
    for(int i = 0; i < x.srcpoints.size(); i ++)
      os << x.srcpoints[i] << " ";
    os << endl;
    os << x.thresh  << endl;
    return os;
  }
  friend istream& operator >> (istream& is, match_t<T>& x) {
    try {
      for(int i = 0; i < x.rot.rows(); i ++)
        for(int j = 0; j < x.rot.cols(); j ++)
          is >> x.rot(i, j);
      for(int i = 0; i < x.offset.size(); i ++)
        is >> x.offset[i];
      is >> x.ratio;
      is >> x.rdepth;
      int size(0);
      is >> size;
      assert(size > 0);
      x.dstpoints.resize(size);
      x.srcpoints.resize(size);
      for(int i = 0; i < size; i ++)
        is >> x.dstpoints[i];
      for(int i = 0; i < size; i ++)
        is >> x.srcpoints[i];
      is >> x.thresh;
    } catch(...) {
      assert(0 && "match_t input failed.");
    }
    return is;
  }
};

template <typename T> static inline pair<SimpleVector<T>, vector<SimpleVector<T> > > makeG(const vector<SimpleVector<T> >& in) {
  pair<SimpleVector<T>, vector<SimpleVector<T> > > res;
  res.second.reserve(in.size());
  assert(in.size() && in[0].size() == 3);
  res.first = in[0];
  for(int i = 1; i < in.size(); i ++)
    res.first += in[i];
  res.first /= T(in.size());
  for(int i = 0; i < in.size(); i ++)
    res.second.emplace_back(in[i] - res.first);
  return move(res);
}

template <typename T> static inline pair<T, vector<SimpleVector<T> > > normalizeG(const vector<SimpleVector<T> >& s) {
  pair<T, vector<SimpleVector<T> > > res;
  res.first = T(0);
  res.second.reserve(s.size());
  for(int i = 0; i < s.size(); i ++)
    res.first = max(max(res.first, abs(s[i][0])), max(abs(s[i][1]), abs(s[i][2])));
  for(int i = 0; i < s.size(); i ++)
    res.second.emplace_back(s[i] / res.first);
  return move(res);
}

template <typename T> pair<vector<pair<T, int> >, vector<int> > makeMatchPrep(const SimpleVector<T>& on, const int& rows, const int& rowcut) {
  vector<int> pidx;
  vector<pair<T, int> > fidx;
  pidx.resize(rows, - 1);
  for(int i = 0; i < rowcut; i ++) {
    T score(0);
    for(int j = rowcut; j < rows; j ++) {
      const auto lscore(abs(on[i] - on[j]));
      if(score == T(int(0)) || lscore < score) {
        score   = lscore;
        pidx[i] = j;
      }
    }
    fidx.emplace_back(make_pair(score, i));
  }
  for(int i = rowcut; i < rows; i ++) {
    T score(0);
    for(int j = 0; j < rowcut; j ++) {
      const auto lscore(abs(on[i] - on[j]));
      if(score == T(int(0)) || lscore < score) {
        score   = lscore;
        pidx[i] = j;
      }
    }
  }
  sort(fidx.begin(), fidx.end());
  return make_pair(move(fidx), move(pidx));
}

template <typename T> vector<match_t<T> > matchPartial(const vector<SimpleVector<T> >& shapebase0, const vector<SimpleVector<T> >& points0, const T& thresh = T(1) / T(100)) {
  const auto  gs(makeG(shapebase0));
  const auto  gp(makeG(points0));
  const auto  ggs(normalizeG(gs.second));
  const auto  ggp(normalizeG(gp.second));
  const auto& shapebase(ggs.second);
  const auto& points(ggs.second);
  std::cerr << "match(" << shapebase.size() << ", " << points.size() << ")" << std::endl;
  SimpleMatrix<T> test(shapebase.size() + points.size(), 3 + 2);
  for(int i = 0; i < shapebase.size(); i ++)
    test.row(i) = makeProgramInvariant<T>(shapebase[i]).first;
  for(int i = 0; i < points.size(); i ++)
    test.row(shapebase.size() + i) = makeProgramInvariant<T>(points[i]).first;
        auto Pt(test.QR());
  const auto R(Pt * test);
  SimpleVector<T>    one(Pt.cols());
  SimpleVector<bool> fix(one.size());
  one.I(T(int(1)));
  fix.I(false);
  const auto on(Pt.projectionPt(one));
  vector<int> pidx;
  vector<pair<T, int> > fidx;
  const auto fpidx(makeMatchPrep(on, test.rows(), shapebase.size()));
  for(int n_fixed = 0, idx = 0;
          n_fixed < Pt.rows() - 1 && idx < fpidx.first.size();
          n_fixed ++, idx ++) {
    const auto& iidx(fpidx.first[idx].second);
    if(fix[iidx]) continue;
    const auto  orth(Pt.col(iidx));
    const auto  n2(orth.dot(orth));
    if(n2 <= Pt.epsilon) continue;
    for(int j = 0; j < Pt.cols(); j ++)
      Pt.setCol(j, Pt.col(j) - orth * Pt.col(j).dot(orth) / n2);
    fix[iidx] = fix[fpidx.second[iidx]] = true;
  }
  const auto cut(R.solve(Pt * one));
  const auto ffpidx(makeMatchPrep(Pt.projectionPt(one), test.rows(), shapebase.size()));
  match_t<T> m(thresh, max(ggs.first, ggp.first));
  vector<SimpleVector<T> > left, right;
  for(int i = 0;
          i < ffpidx.first.size() && ffpidx.first[i].first < thresh;
          i ++) {
    if(ffpidx.second[ffpidx.first[i].second] < 0) continue;
    if(ffpidx.first[i].second < shapebase.size()) {
      m.dstpoints.emplace_back(ffpidx.first[i].second);
      m.srcpoints.emplace_back(ffpidx.second[ffpidx.first[i].second] - shapebase.size());
    } else {
      m.srcpoints.emplace_back(ffpidx.first[i].second - shapebase.size());
      m.dstpoints.emplace_back(ffpidx.second[ffpidx.first[i].second]);
    }
    assert(0 <= m.dstpoints[m.dstpoints.size() - 1] &&
                m.dstpoints[m.dstpoints.size() - 1] < shapebase.size() &&
           0 <= m.srcpoints[m.srcpoints.size() - 1] &&
                m.srcpoints[m.srcpoints.size() - 1] < points.size());
  }
  assert(m.dstpoints.size() == m.srcpoints.size());
  if(m.dstpoints.size() < 4) return vector<match_t<T> >();
  m.rdepth = T(0);
  SimpleMatrix<T> rot(3, 3);
  rot.O();
  for(int k = 0; k < m.dstpoints.size() - 3; k ++) {
    SimpleMatrix<T> rotl(3, 3);
    SimpleMatrix<T> rotr(3, 3);
    rotl.O();
    rotr.O();
    for(int kk = 0; kk < 3; kk ++) {
      rotl.row(kk) = shapebase[m.dstpoints[k + kk]];
      rotr.row(kk) = m.transform(points[m.srcpoints[k + kk]]);
    }
    rot += rotr.QR() * rotl.QR().transpose();
  }
  m.rot = (rot /= T(m.dstpoints.size() - 3));
  T r0(0);
  for(int k = 0; k < m.dstpoints.size(); k ++) {
    const auto& shapek(shapebase[m.dstpoints[k]]);
    const auto  pointk(m.transform(points[m.srcpoints[k]]));
    r0 += shapek.dot(pointk) / sqrt(shapek.dot(shapek) * pointk.dot(pointk));
  }
  m.ratio = r0 / T(m.dstpoints.size());
  for(int k = 0; k < m.dstpoints.size(); k ++) {
    const auto& shapek(shapebase[m.dstpoints[k]]);
    const auto  pointk(m.transform(points[m.srcpoints[k]]));
    const auto  err(shapek - pointk);
    m.rdepth += sqrt(err.dot(err) / sqrt(shapek.dot(shapek) * pointk.dot(pointk)));
  }
  m.rdepth /= m.dstpoints.size() * m.dstpoints.size();
  m.offset += gs.first - m.rot * gp.first * m.ratio;
  vector<match_t<T> > mm;
  mm.emplace_back(m);
  return mm;
}

#define _MATCH_
#endif

