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
    os << x.rot;
    os << x.offset;
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
      is >> x.rot;
      is >> x.offset;
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

template <typename T> static inline SimpleVector<T> toQuarterNormalize(const SimpleVector<T>& xyz) {
  assert(xyz.size() == 3);
  static const auto twoPi(atan(T(int(1))) * T(int(4)));
  SimpleVector<T> quat(4);
  quat[0] = sqrt(xyz.dot(xyz) / T(int(6)));
  // y-z plane
  quat[1] = atan2(xyz[1], xyz[2]) / twoPi;
  // same for z-x.
  quat[2] = atan2(xyz[2], xyz[0]) / twoPi;
  // same for x-y.
  quat[3] = atan2(xyz[0], xyz[1]) / twoPi;
  return quat;
}

template <typename T> match_t<T> reconfigureMatch(match_t<T>& m, const vector<SimpleVector<T> >& dst0, const vector<SimpleVector<T> >& src0) {
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
  T rlog0(int(0));
  for(int k = 0; k < m.dst.size(); k ++) {
    const auto& dstk(dst0[m.dst[k]]);
    const auto  srck(m.transform(src0[m.src[k]]));
    const auto  r(abs(dstk.dot(srck) / srck.dot(srck)));
    if(r != T(0)) rlog0 += log(r);
  }
  m.ratio *= exp(rlog0 / T(int(m.dst.size())));
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
  return m;
}

template <typename T> vector<match_t<T> > matchPartialR(const vector<SimpleVector<T> >& dst0, const vector<SimpleVector<T> >& src0, const T& thresh = T(1) / T(10)) {
  const auto  gs(normalizeG((makeG(dst0).second)));
  const auto  gp(normalizeG((makeG(src0).second)));
  const auto& dst(gs.second);
  const auto& src(gp.second);
  assert(dst.size() == dst0.size());
  assert(src.size() == src0.size());
  std::cerr << "match(" << dst.size() << ", " << src.size() << ")" << std::endl;
  SimpleMatrix<T> qdst(dst.size(), 4);
  SimpleMatrix<T> qsrc(src.size(), 4);
  for(int i = 0; i < qdst.rows(); i ++)
    qdst.row(i) = toQuarterNormalize<T>(dst[i]);
  for(int i = 0; i < qsrc.rows(); i ++)
    qsrc.row(i) = toQuarterNormalize<T>(src[i]);
  vector<SimpleVector<T> > test;
  vector<pair<int, int> > idx;
  test.reserve(qdst.rows() * qsrc.rows());
  idx.reserve(qdst.rows() * qsrc.rows());
  for(int i = 0; i < qdst.rows(); i ++)
    for(int j = 0; j < qsrc.rows(); j ++) {
      idx.emplace_back(make_pair(i, j));
      test.emplace_back(qdst.row(i) - qsrc.row(j));
    }
  const auto cr(crush<T>(test));
  vector<match_t<T> > mm;
  mm.reserve(cr.size());
  for(int i = 0; i < cr.size(); i ++) {
    match_t<T> m(thresh, max(gs.first, gp.first));
    SimpleVector<bool> dfix, sfix;
    dfix.resize(dst.size());
    sfix.resize(src.size());
    dfix.I(false);
    sfix.I(false);
    m.dst.reserve(min(dst.size(), src.size()));
    m.src.reserve(min(dst.size(), src.size()));
    for(int j = 0; j < cr[i].second.size(); j ++) {
      const auto& lidx(idx[cr[i].second[j]]);
      if(dfix[lidx.first] || sfix[lidx.second]) continue;
      dfix[lidx.first] = sfix[lidx.second] = true;
      m.dst.emplace_back(lidx.first);
      m.src.emplace_back(lidx.second);
    }
    if(m.dst.size() < 4) continue;
    m.dst.reserve(m.dst.size());
    m.src.reserve(m.src.size());
    m = reconfigureMatch<T>(m, dst0, src0);
    for(int i = 0; i < m.rot.rows(); i ++)
      for(int j = 0; j < m.rot.cols(); j ++)
        if(! isfinite(m.rot(i, j))) goto nofix;
    for(int i = 0; i < m.offset.size(); i ++)
      if(! isfinite(m.offset[i])) goto nofix;
    if(! isfinite(m.ratio)) goto nofix;
    mm.emplace_back(move(m));
    cerr << mm[mm.size() - 1] << endl;
   nofix:
    ;
  }
  return mm;
}

template <typename T> vector<match_t<T> > matchPartial(const vector<SimpleVector<T> >& dst0, const vector<SimpleVector<T> >& src0, const T& thresh = T(1) / T(10)) {
  auto m(matchPartialR<T>(src0, dst0, thresh));
  for(int i = 0; i < m.size(); i ++)
    m[i] = ~ m[i];
  return m;
}

#define _MATCH_
#endif

