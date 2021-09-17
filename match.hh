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
  vector<int> dst;
  vector<int> src;
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
  inline match_t(const match_t<T>& s) {
    *this = s;
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
    std::swap(result.src, result.dst);
    return result;
  }
  inline match_t<T>  operator / (const match_t<T>& s) const {
    match_t<T> result;
    result.rot    = rot   * s.rot.transpose();
    result.ratio  = ratio / s.ratio;
    result.offset = offset - result.rot * s.offset * result.ratio;
    result.rdepth = rdepth + s.rdepth;
    result.dst = vector<int>();
    result.src = vector<int>();
    for(int i = 0; i < src.size(); i ++)
      for(int j = 0; j < s.src.size(); j ++)
        if(src[i] == s.src[j]) {
          result.dst.emplace_back(  dst[i]);
          result.src.emplace_back(s.dst[j]);
        }
    return result;
  }
  inline match_t<T>& operator = (const match_t<T>& other) {
    rot        = other.rot;
    offset     = other.offset;
    ratio      = other.ratio;
    rdepth     = other.rdepth;
    dst  = other.dst;
    src  = other.src;
    thresh     = other.thresh;
    othresh    = other.othresh;
    return *this;
  }
  inline T distance(const match_t<T>& other, const Vec& p) {
    const auto d(transform(p) - other.transform(p));
    return sqrt(d.dot(d));
  }
  inline vector<Veci> hullConv(const vector<Veci>& srchull) const {
    assert(src.size() == dst.size());
    vector<Veci> res;
    res.reserve(srchull.size());
    for(int i = 0; i < srchull.size(); i ++) {
      Veci tmp(3);
      tmp[0] = tmp[1] = tmp[2] = - 1;
      for(int j = 0; j < srchull[i].size(); j ++)
        for(int k = 0; k < src.size(); k ++)
          if(src[k] == srchull[i][j]) {
            tmp[j] = dst[k];
            break;
          }
      assert(0 <= tmp[0] && 0 <= tmp[1] && 0 <= tmp[2]);
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
    assert(x.dst.size() == x.src.size());
    os << x.dst.size() << endl;
    for(int i = 0; i < x.dst.size(); i ++)
      os << x.dst[i] << " ";
    os << endl;
    for(int i = 0; i < x.src.size(); i ++)
      os << x.src[i] << " ";
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
      x.dst.resize(size);
      x.src.resize(size);
      for(int i = 0; i < size; i ++)
        is >> x.dst[i];
      for(int i = 0; i < size; i ++)
        is >> x.src[i];
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
  assert(res.second.size() == in.size());
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
  assert(res.second.size() == s.size());
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
    fidx.emplace_back(make_pair(score, i));
  }
  for(int i = rows; i < rows + 4; i ++)
    fidx.emplace_back(make_pair(on[i] < T(0) ? - T(on.size() * 2) : T(on.size() * 2), pidx[i] = i));
  sort(fidx.begin(), fidx.end());
  return make_pair(move(fidx), move(pidx));
}

template <typename T> SimpleVector<T> toQuarterNormalize(const SimpleVector<T>& xyz) {
  assert(xyz.size() == 3);
  static const auto twoPi(atan(T(int(1))) * T(int(4)));
  SimpleVector<T> quat(4);
  // XXX: order is broken:
  quat[0] = sqrt(xyz.dot(xyz) / T(int(6)));
  // y-z plane normal vector. <=> xyz.dot(n == [0 y z]) / quat[0] / sqrt(y * y + z * z) == cos theta.
  quat[1] = acos(xyz[1] / sqrt(xyz[1] * xyz[1] + xyz[2] * xyz[2])) / twoPi;
  // same for z-x.
  quat[2] = acos(xyz[2] / sqrt(xyz[2] * xyz[2] + xyz[0] * xyz[0])) / twoPi;
  // same for x-y.
  quat[3] = acos(xyz[0] / sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1])) / twoPi;
  T slog(int(0));
  // XXX : without this, assertion crash.
  /* XXX */ std::cerr << quat;
  for(int i = 0; i < quat.size(); i ++) {
    assert(- T(int(1)) <= quat[i] && quat[i] <= T(int(1)));
    quat[i] += T(int(1));
    if(quat[i] != T(int(0)))
      slog += log(quat[i]);
  }
  return quat /= exp(slog / T(quat.size()));
}

template <typename T> match_t<T> matchPartial(const vector<SimpleVector<T> >& dst0, const vector<SimpleVector<T> >& src0, const T& thresh = T(1) / T(10)) {
  const auto  gs(makeG(dst0));
  const auto  gp(makeG(src0));
  const auto  ggs(normalizeG(gs.second));
  const auto  ggp(normalizeG(gp.second));
  const auto& dst(ggs.second);
  const auto& src(ggp.second);
  assert(dst.size() == dst0.size());
  assert(src.size() == src0.size());
  std::cerr << "match(" << dst.size() << ", " << src.size() << ")" << std::endl;
  SimpleMatrix<T> test(dst.size() + src.size() + 4, 4);
  for(int i = 0; i < dst.size(); i ++)
    test.row(i) = toQuarterNormalize<T>(dst[i]);
  for(int i = 0; i < src.size(); i ++)
    test.row(dst.size() + i) = toQuarterNormalize<T>(src[i]);
  for(int i = 0; i < 4; i ++)
    test.row(i + dst.size() + src.size()).ek(i, thresh);
        auto Pt(test.QR());
  const auto R(Pt * test);
  SimpleVector<T>    one(Pt.cols());
  SimpleVector<bool> fix(one.size());
  one.I(T(int(1)));
  fix.I(false);
  const auto on(Pt.projectionPt(one));
  vector<pair<T, int> > fidx;
  auto fpidx(makeMatchPrep(on, dst.size() + src.size(), dst.size()));
  for(int n_fixed = 0, idx = 0;
          n_fixed < Pt.rows() - 1 && idx < fpidx.first.size();
          n_fixed ++) {
    for(idx = 0; (fix[fpidx.first[idx].second] || fpidx.second[fpidx.first[idx].second] < 0) && idx < fpidx.first.size(); idx ++) ;
    if(idx < 0 || fpidx.first.size() <= idx) break;
    const auto& iidx(fpidx.first[idx].second);
    const auto  orth(Pt.col(iidx));
    const auto  n2(orth.dot(orth));
    if(n2 <= Pt.epsilon) continue;
    for(int j = 0; j < Pt.cols(); j ++)
      Pt.setCol(j, Pt.col(j) - orth * Pt.col(j).dot(orth) / n2);
    fix[iidx] = fix[fpidx.second[iidx]] = true;
    for(int j = 0; j < dst.size(); j ++) {
      if(fix[j]) continue;
      fpidx.second[j] = - 1;
      T score(0);
      for(int jj = dst.size(); jj < dst.size() + src.size(); jj ++) {
        if(fix[jj]) continue;
        const auto lscore(abs(on[j] - on[jj]));
        if(score == T(int(0)) || lscore < score) {
          score   = lscore;
          fpidx.second[j] = jj;
        }
      }
      fpidx.first[j].first = score;
    }
    for(int j = dst.size(); j < dst.size() + src.size(); j ++) {
      if(fix[j]) continue;
      fpidx.second[j] = - 1;
      T score(0);
      for(int jj = 0; jj < dst.size(); jj ++) {
        if(fix[jj]) continue;
        const auto lscore(abs(on[j] - on[jj]));
        if(score == T(int(0)) || lscore < score) {
          score   = lscore;
          fpidx.second[j] = jj;
        }
      }
      fpidx.first[j].first = score;
    }
    for(int j = dst.size() + src.size(); j < Pt.cols(); j ++)
      fpidx.first[j] = make_pair(on[j] < T(0) ? - T(on.size() * 2) : T(on.size() * 2), j);
    sort(fpidx.first.begin(), fpidx.first.end());
  }
  const auto cut(R.solve(Pt * one));
  const auto ffpidx(makeMatchPrep(Pt.projectionPt(one), dst.size() + src.size(), dst.size()));
  // N.B. what we got: <a, [q0 q1 q2 q3] - [q0' q1' q2' q3']> == 0, a != 0.
  //      and 0 <= a condition makes this valid.
  match_t<T> m(thresh, max(ggs.first, ggp.first));
  for(int i = 0;
          i < ffpidx.first.size() && ffpidx.first[i].first < thresh;
          i ++) {
    if(ffpidx.second[ffpidx.first[i].second] < 0) continue;
    if(ffpidx.first[i].second < dst.size()) {
      m.dst.emplace_back(ffpidx.first[i].second);
      m.src.emplace_back(ffpidx.second[ffpidx.first[i].second] - dst.size());
    } else if(ffpidx.first[i].second < dst.size() + src.size()) {
      m.src.emplace_back(ffpidx.first[i].second - dst.size());
      m.dst.emplace_back(ffpidx.second[ffpidx.first[i].second]);
    } else continue;
    assert(0 <= m.dst[m.dst.size() - 1] &&
                m.dst[m.dst.size() - 1] < dst0.size() &&
           0 <= m.src[m.src.size() - 1] &&
                m.src[m.src.size() - 1] < src0.size());
    assert(src0[m.src[m.src.size() - 1]].size() == 3 &&
           dst0[m.dst[m.dst.size() - 1]].size() == 3);
  }
  assert(m.dst.size() == m.src.size());
  SimpleVector<T> off(3);
  off.O();
  for(int k = 0; k < m.dst.size(); k ++)
    off += dst0[m.dst[k]] - m.transform(src0[m.src[k]]);
  m.offset += (off /= T(int(m.dst.size())));
  if(m.dst.size() < 4) return m;
  SimpleMatrix<T> rot(3, 3);
  rot.I();
  for(int k = 0; k < m.dst.size() - 3; k ++) {
    SimpleMatrix<T> rotl(3, 3);
    SimpleMatrix<T> rotr(3, 3);
    rotl.O();
    rotr.O();
    for(int kk = 0; kk < 3; kk ++) {
      rotl.setCol(kk, dst0[m.dst[k + kk]]);
      rotr.setCol(kk, m.transform(src0[m.src[k + kk]]));
    }
    // dst == Q R == P R' == src, dst !~ avg(Q P^t) src'.
    rot *= rotl.QR().transpose() * rotr.QR();
  }
  m.rot = pow(rot, T(1) / T(int(m.dst.size() - 3)));
  off.O();
  for(int k = 0; k < m.dst.size(); k ++)
    off += dst0[m.dst[k]] - m.transform(src0[m.src[k]]);
  m.offset += (off /= T(int(m.dst.size())));
  T r0(int(1));
  for(int k = 0; k < m.dst.size(); k ++) {
    const auto& dstk(dst0[m.dst[k]]);
    const auto  srck(m.transform(src0[m.src[k]]));
    r0 *= dstk.dot(srck) / srck.dot(srck);
  }
  m.ratio *= (r0 < T(int(0)) ? - T(int(1)) : T(int(1))) * pow(abs(r0), T(int(1)) / T(int(m.dst.size())));
  off.O();
  for(int k = 0; k < m.dst.size(); k ++)
    off += dst0[m.dst[k]] - m.transform(src0[m.src[k]]);
  m.offset += (off /= T(int(m.dst.size())));
  for(int k = 0; k < m.dst.size(); k ++) {
    const auto& dstk(dst0[m.dst[k]]);
    const auto  srck(m.transform(src0[m.src[k]]));
    const auto  err(dstk - srck);
    m.rdepth += sqrt(err.dot(err) / sqrt(dstk.dot(dstk) * srck.dot(srck)));
  }
  m.rdepth /= m.dst.size() * m.dst.size();
  return move(m);
}

#define _MATCH_
#endif

