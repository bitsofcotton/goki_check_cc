/* BSD 3-Clause License:
 * Copyright (c) 2018, bitsofcotton.
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

template <typename T> class reDig;

using std::min;
using std::max;
using std::abs;
using std::cerr;
using std::endl;
using std::flush;
using std::vector;
using std::lower_bound;
using std::upper_bound;
using std::distance;
using std::sort;
using std::unique;
using std::isfinite;
using std::istream;
using std::ostream;

template <typename T> class msub_t;
template <typename T> bool lmsublt(const msub_t<T>& x, const msub_t<T>& y) {
  return x.err < y.err;
}

template <typename T> class msub_t {
public:
  T   t;
  int j;
  int k;
  T   err;
  inline msub_t() {
    t = err = T(0);
    j = k = 0;
  }
  inline msub_t(const msub_t<T>& x) {
    *this = x;
  }
  inline ~msub_t() {
    ;
  }
  inline msub_t<T>& operator = (const msub_t<T>& x) {
    t   = T(x.t);
    j   = x.j;
    k   = x.k;
    err = T(x.err);
    return *this;
  }
  inline bool operator < (const msub_t<T>& x) const {
    return  t >  x.t ||
           (t == x.t && err <  x.err) ||
           (t == x.t && err == x.err && j <  x.j) ||
           (t == x.t && err == x.err && j == x.j && k < x.k);
  }
};

template <typename T> class match_t {
public:
  typedef SimpleMatrix<T>   Mat3x3;
  typedef SimpleVector<T>   Vec3;
  typedef SimpleVector<T>   Vec2;
  typedef SimpleVector<int> Vec3i;
  Mat3x3      rot;
  Vec3        offset;
  T           ratio;
  T           rdepth;
  vector<int> dstpoints;
  vector<int> srcpoints;
  T           thresh;
  T           rthresh;
  Vec2        threshsize;
  inline match_t() {
    thresh  = T(0);
    rthresh = T(1) / T(100);
    threshsize    = Vec2(2);
    threshsize[0] = threshsize[1] = T(0);
    initId();
  }
  inline match_t(const T& thresh, const T& h, const T& w) {
    this->thresh        = thresh;
    this->rthresh       = pow(h * w, - T(1) / T(4));
    this->threshsize    = Vec2(2);
    this->threshsize[0] = h;
    this->threshsize[1] = w;
    initId();
  }
  inline void initId() {
    rot       = Mat3x3(3, 3);
    rot(0, 0) = rot(1, 1) = rot(2, 2) = T(1);
    rot(1, 0) = rot(2, 0) = rot(0, 1) = rot(2, 1)
              = rot(0, 2) = rot(1, 2) = T(0);
    offset    = Vec3(3);
    offset[0] = offset[1] = offset[2] = T(0);
    ratio     = T(1);
    rdepth    = T(0);
  }
  inline match_t(const match_t<T>& other) {
    *this = other;
  }
  inline match_t<T>  operator ~ () const {
    match_t<T> result;
    result.rot    = rot.transpose();
    result.ratio  = T(1) / ratio;
    result.offset = - result.rot * offset * result.ratio;
    result.rdepth = rdepth;
    result.dstpoints  = srcpoints;
    result.srcpoints  = dstpoints;
    result.thresh     = thresh;
    result.rthresh    = rthresh;
    result.threshsize = threshsize;
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
    rthresh    = other.rthresh;
    threshsize = other.threshsize;
    return *this;
  }
  inline T distance(const match_t<T>& other, const Vec3& p) {
    const auto d(transform(p) - other.transform(p));
    return sqrt(d.dot(d));
  }
  inline vector<Vec3i> hullConv(const vector<Vec3i>& srchull) const {
    assert(srcpoints.size() == dstpoints.size());
    vector<Vec3i> res;
    res.reserve(srchull.size());
    for(int i = 0; i < srchull.size(); i ++) {
      Vec3i tmp(3);
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
  inline Vec3 transform(const Vec3& x) const {
    return rot * x * ratio + offset;
  }
  inline vector<Vec3> transform(const vector<Vec3>& x) const {
    vector<Vec3> result(x);
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
               (offset.dot(offset) + x.offset.dot(x.offset))) /
                 sqrt(threshsize[0] * threshsize[1]) <= thresh) ||
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
    os << x.rthresh << endl;
    os << x.threshsize[0] << " " << x.threshsize[1] << endl;
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
      is >> x.rthresh;
      is >> x.threshsize[0];
      is >> x.threshsize[1];
    } catch(...) {
      assert(0 && "match_t input failed.");
    }
    return is;
  }
};


template <typename T> class matchPartial {
public:
  typedef SimpleMatrix<T> Mat;
  typedef SimpleMatrix<T> Mat3x3;
  typedef SimpleMatrix<T> Mat2x2;
  typedef SimpleVector<T> Vec3;
  typedef complex<T> U;
  inline matchPartial();
  inline matchPartial(const int& ndiv, const T& threshr, const T& threshp, const T& threshs);
  inline ~matchPartial();
  void init(const int& ndiv, const T& threshr, const T& threshp, const T& threshs);
  
  vector<match_t<T> > match(const vector<Vec3>& shapebase, const vector<Vec3>& points, const bool& norot = false);
  void match(const vector<Vec3>& shapebase0, const vector<Vec3>& points0, vector<match_t<T> >& result, const bool& norot = false);
  
  Vec3 makeG(const vector<Vec3>& in) const;
  vector<match_t<T> > elim(const vector<match_t<T> >& m, const Mat dst[3], const Mat src[3], const Mat& srcbump, const vector<Vec3>& srcpts, const T& thresh = T(4) / T(256));
  
  // theta resolution.
  int ndiv;
  // match points thresh in [0, 1].
  T   threshp;
  // match operator == thresh in [0, 1].
  T   threshs;
  
private:
  T   isElim(const match_t<T>& m, const Mat dst[3], const Mat tsrc[3], const vector<Vec3>& srcpts, const T& thresh);
  U   I;
  T   Pi;
  // match theta  thresh in [0, 1].
  T   thresh;
  // match ratio  thresh in [0, 1].
  T   thresht;
};

template <typename T> inline matchPartial<T>::matchPartial() {
  I  = sqrt(U(- T(1)));
  Pi = atan2(T(1), T(1)) * T(4);
  // rough match.
  init(40, 2, T(1) / T(1000), T(1) / T(10));
}

template <typename T> inline matchPartial<T>::matchPartial(const int& ndiv, const T& threshr, const T& threshp, const T& threshs) {
  I  = sqrt(U(- T(1)));
  Pi = atan2(T(1), T(1)) * T(4);
  init(ndiv, threshr, threshp, threshs);
}

template <typename T> inline matchPartial<T>::~matchPartial() {
  ;
}

template <typename T> void matchPartial<T>::init(const int& ndiv, const T& threshr, const T& threshp, const T& threshs) {
  assert(0 < ndiv && T(1) <= threshr && T(0) <= threshp && threshp <= T(1));
  this->ndiv    = ndiv;
  this->thresh  = sin(T(2) * Pi / T(ndiv)) / T(2) * threshr;
  this->thresht = sin(T(2) * Pi / T(ndiv)) / T(2) * threshr;
  this->threshp = threshp;
  this->threshs = threshs;
  assert(this->thresh < T(1));
  cerr << "e(" << this->thresh << ")" << endl;
  return;
}

template <typename T> vector<match_t<T> > matchPartial<T>::match(const vector<Vec3>& shapebase, const vector<Vec3>& points, const bool& norot) {
  vector<match_t<T> > result;
  match(shapebase, points, result, norot);
  return result;
}

template <typename T> typename matchPartial<T>::Vec3 matchPartial<T>::makeG(const vector<Vec3>& in) const {
  Vec3 result(3);
  result[0] = result[1] = result[2] = T(0);
  for(int i = 0; i < in.size(); i ++)
    result += in[i];
  return result / in.size();
}

template <typename T> void matchPartial<T>::match(const vector<Vec3>& shapebase0, const vector<Vec3>& points0, vector<match_t<T> >& result, const bool& norot) {
  const auto gs(makeG(shapebase0));
  const auto gp(makeG(points0));
  vector<Vec3> shapebase;
  vector<Vec3> points;
  shapebase.reserve(shapebase0.size());
  points.reserve(points0.size());
  Vec3 gd(3);
  gd[0] = gd[1] = gd[2] = T(0);
  for(int i = 0; i < shapebase0.size(); i ++) {
    shapebase.emplace_back(shapebase0[i] - gs);
    const auto& diff(shapebase[i]);
    gd[0] = max(gd[0], abs(diff[0]));
    gd[1] = max(gd[1], abs(diff[1]));
  }
  for(int i = 0; i < points0.size(); i ++) {
    points.emplace_back(points0[i] - gp);
    const auto& diff(points[i]);
    gd[0] = max(gd[0], abs(diff[0]));
    gd[1] = max(gd[1], abs(diff[1]));
  }
  // for each rotation, we can now handle t <= 0:
  // We need large memory for this.
#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic)
#endif
  for(int nd = 0; nd < ndiv; nd ++) {
    vector<T> ddiv;
    ddiv.resize(4, T(0));
    ddiv[0] = cos(T(2) * Pi * T(nd) / T(ndiv));
    ddiv[1] = sin(T(2) * Pi * T(nd) / T(ndiv));
    for(int nd2 = 0; nd2 < (norot ? 1 : ndiv); nd2 ++) {
      ddiv[2] = cos(T(2) * Pi * T(nd2) / T(ndiv));
      ddiv[3] = sin(T(2) * Pi * T(nd2) / T(ndiv));
      match_t<T> work0(threshs, abs(gd[0]), abs(gd[1]));
      for(int k = 0; k < ddiv.size() / 2; k ++) {
        Mat3x3 lrot(3, 3);
        lrot((k    ) % 3, (k    ) % 3) =   ddiv[k * 2 + 0];
        lrot((k + 1) % 3, (k    ) % 3) =   ddiv[k * 2 + 1];
        lrot((k + 2) % 3, (k    ) % 3) = T(0);
        lrot((k    ) % 3, (k + 1) % 3) = - ddiv[k * 2 + 1];
        lrot((k + 1) % 3, (k + 1) % 3) =   ddiv[k * 2 + 0];
        lrot((k + 2) % 3, (k + 1) % 3) = T(0);
        lrot((k    ) % 3, (k + 2) % 3) = T(0);
        lrot((k + 1) % 3, (k + 2) % 3) = T(0);
        lrot((k + 2) % 3, (k + 2) % 3) = T(1);
        work0.rot = lrot * work0.rot;
      }
      // for each near matches:
      vector<msub_t<T> > msub;
      msub.reserve(points.size());
      for(int k = 0; k < points.size(); k ++) {
        const auto bk(work0.transform(points[k]));
        vector<msub_t<T> > lmsub;
        lmsub.reserve(shapebase.size());
        for(int j = 0; j < shapebase.size(); j ++) {
          const auto& aj(shapebase[j]);
          const auto  t(aj.dot(bk) / bk.dot(bk));
          const auto  lerr(aj - bk * t);
          const auto  err(lerr.dot(lerr) / sqrt(aj.dot(aj) * bk.dot(bk) * t * t));
          // if t <= T(0), it's mirrored and this should not match.
          if(T(0) <= t && err <= thresht * thresht &&
             isfinite(t) && isfinite(err)) {
            msub_t<T> work;
            work.t   = t;
            work.j   = j;
            work.k   = k;
            work.err = err;
            lmsub.emplace_back(work);
          }
        }
        if(lmsub.size()) {
          sort(lmsub.begin(), lmsub.end(), lmsublt<T>);
          msub.emplace_back(lmsub[0]);
        }
      }
      sort(msub.begin(), msub.end());
      cerr << msub.size() << ":" << flush;
      if(T(msub.size()) /
           T(min(shapebase.size(), points.size())) < threshp)
        continue;
      // get nearer matches:
      int t0(0);
      for( ; t0 < msub.size(); t0 ++) {
        auto work(work0);
        int tt(t0);
        const auto tr0(msub[t0].t);
        for(int t1 = t0; t1 < msub.size(); t1 ++)
          // N.B. abar:=sum(aj)/n, bbar:=P*sum(bk)/n,
          //   get condition sum||aj-abar||, sum||P*bk-bbar|| -> 0
          //   with this imcomplete set.
          // N.B. t >= 0 and msub is sorted by t0 > t1.
          if((msub[t0].t - msub[t1].t) / msub[t0].t <= thresht) {
            work.dstpoints.emplace_back(msub[t1].j);
            work.srcpoints.emplace_back(msub[t1].k);
            tt = t1;
          }
        // N.B. rough match.
        t0 += max(1, (tt - t0) / 2);
        // if it's good:
        if(threshp <= T(work.dstpoints.size()) /
                        T(min(shapebase.size(), points.size()))) {
          work.ratio  = sqrt(tr0 * msub[tt].t);
          work.rdepth = T(0);
          for(int k = 0; k < work.dstpoints.size(); k ++) {
            const auto& shapek(shapebase[work.dstpoints[k]]);
            const auto  pointk(work.transform(points[work.srcpoints[k]]));
            const auto  err(shapek - pointk);
            work.rdepth += sqrt(err.dot(err) / sqrt(shapek.dot(shapek) * pointk.dot(pointk)));
          }
          work.rdepth /= work.dstpoints.size() * work.dstpoints.size();
          work.offset += gs - work.rot * gp * work.ratio;
#if defined(_OPENMP)
#pragma omp critical
#endif
          {
            int idx(- 1);
            int k0(distance(result.begin(), lower_bound(result.begin(), result.end(), work)));
            for(int k = k0; k < result.size(); k ++)
              if(work < result[k]) {
                if(result[k] == work) {
                  idx = k;
                  break;
                }
              } else if(k - k0)
                break;
            if(idx >= 0)
              result[idx] = work;
            else {
              result.emplace_back(work);
              cerr << "*" << flush;
            }
            sort(result.begin(), result.end());
          }
        }
      }
    }
  }
  return;
}

template <typename T> vector<match_t<T> > matchPartial<T>::elim(const vector<match_t<T> >& m, const Mat dst[3], const Mat src[3], const Mat& srcbump, const vector<Vec3>& srcpts, const T& thresh) {
  vector<match_t<T> > res(m);
  reDig<T> redig;
  for(int i = 0; i < m.size(); i ++) {
    Mat tsrc[3];
    const auto ref(redig.tilt(redig.makeRefMatrix(src[0], 1), srcbump, m[i]));
    for(int j = 0; j < 3; j ++)
      tsrc[j] = redig.pullRefMatrix(ref, 1, src[j]);
    res[i].rdepth *= isElim(m[i], dst, tsrc, srcpts, thresh);
  }
  sort(res.begin(), res.end());
  return res;
}

template <typename T> T matchPartial<T>::isElim(const match_t<T>& m, const Mat dst[3], const Mat tsrc[3], const vector<Vec3>& srcpts, const T& thresh) {
  vector<T> diffs;
  for(int i = 0; i < m.srcpoints.size(); i ++) {
    const auto  g(m.transform(srcpts[m.srcpoints[i]]));
    const auto& y(g[0]);
    const auto& x(g[1]);
    if(T(0) <= y && y < T(tsrc[0].rows()) && T(0) <= x && x < T(tsrc[0].cols())) {
      if(tsrc[0](y, x) == T(0) && tsrc[1](y, x) == T(0) && tsrc[2](y, x) == T(0))
        continue;
      T diff(0);
      for(int j = 0; j < 3; j ++)
        diff += pow(tsrc[j](y, x) - dst[j](y, x), T(2));
      diffs.emplace_back(sqrt(diff));
    }
  }
  sort(diffs.begin(), diffs.end());
  vector<T> ddiffs;
  for(int i = 1; i < diffs.size(); i ++)
    ddiffs.emplace_back(diffs[i] - diffs[i - 1]);
  sort(ddiffs.begin(), ddiffs.end());
  return ddiffs.size() ? T(1) - T(distance(ddiffs.begin(), upper_bound(ddiffs.begin(), ddiffs.end(), thresh))) / T(ddiffs.size()) : T(0);
}


template <typename T> class matchWhole {
public:
  typedef SimpleMatrix<T>   Mat;
  typedef SimpleMatrix<T>   Mat3x3;
  typedef SimpleVector<T>   Vec3;
  typedef SimpleVector<int> Veci4;
  typedef complex<T> U;
  inline matchWhole();
  inline ~matchWhole();
  void init(const int& ndiv, const T& threshp, const T& threshs);

  vector<vector<match_t<T> > > match(const vector<Vec3>& shapebase, const vector<vector<Vec3> >& points, const vector<Vec3>& origins, const vector<vector<Veci4> >& bones, const int& ntry = 6);
private:
  U   I;
  T   Pi;
  int ndiv;
  T   threshp;
  T   threshs;
};

template <typename T> inline matchWhole<T>::matchWhole() {
  I  = sqrt(U(- T(1)));
  Pi = atan2(T(1), T(1)) * T(4);
  init(40, .05, .1);
}

template <typename T> inline matchWhole<T>::~matchWhole() {
  ;
}

template <typename T> void matchWhole<T>::init(const int& ndiv, const T& threshp, const T& threshs) {
  this->ndiv    = ndiv;
  this->threshp = threshp;
  this->threshs = threshs;
  return;
}

template <typename T> vector<vector<match_t<T> > > matchWhole<T>::match(const vector<Vec3>& shapebase, const vector<vector<Vec3> >& points, const vector<Vec3>& origins, const vector<vector<Veci4> >& bones, const int& ntry) {
  assert(points.size() == origins.size());
  vector<vector<match_t<T> > > pmatches;
  matchPartial<T> pmatch(ndiv, threshp, threshs);
  for(int i = 0; i < points.size(); i ++) {
    cerr << "matching partials : " << i << "/" << shapebase.size() << endl;
    pmatches.emplace_back(pmatch.match(shapebase, points[i]));
  }
  int i0idx(0);
  for(int i = 1; i < pmatches.size(); i ++)
    if(pmatches[i0idx].size() < pmatches[i].size())
      i0idx = i;
  vector<vector<match_t<T> > > result;
  for(int i0 = 0; i0 < min(ntry, int(pmatches[i0idx].size())); i0 ++) {
    vector<match_t<T> > lmatch;
    lmatch.resize(pmatches.size(), match_t<T>());
    lmatch[i0idx] = pmatches[i0idx][i0];
    int mcount(1);
    for(int i = 0; i < pmatches.size(); i ++)
      if(i != i0idx && pmatches[i].size()) {
        int idx(0);
        T   m(0);
        T   lratio(1);
        for(int k = 0; k < pmatches[i].size(); k ++) {
          T lm(0);
          T llratio(1);
          vector<int> boneidxs;
          for(int j = 0; j < pmatches[i][0].srcpoints.size(); j ++) {
            const auto& work(bones[i][pmatches[i][0].srcpoints[j]]);
            for(int k = 0; k < work.size(); k ++)
              if(0 <= work[k])
                boneidxs.emplace_back(work[k]);
          }
          sort(boneidxs.begin(), boneidxs.end());
          boneidxs.erase(unique(boneidxs.begin(), boneidxs.end()), boneidxs.end());
          for(int j = 0; j < boneidxs.size(); j ++) {
            lm      += lmatch[boneidxs[j]].distance(pmatches[i][0], origins[i]);
            llratio *= max(abs(pmatches[i][0].ratio) / abs(lmatch[boneidxs[j]].ratio), abs(lmatch[boneidxs[j]].ratio) / abs(pmatches[i][0].ratio));
          }
          if(!k || (m < lm && (llratio <= lratio || llratio < T(1) + threshs))) {
            m      = lm;
            lratio = llratio;
            idx    = k;
          }
        }
        lmatch[i] = pmatches[i][idx];
        mcount ++;
      }
    if(threshp <= mcount / T(pmatches.size()))
      result.emplace_back(lmatch);
  }
  return result;
}

#define _MATCH_
#endif

