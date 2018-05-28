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

#include <Eigen/Core>
#include <Eigen/LU>
#include <cmath>
#include <vector>
#include "tilt.hh"

using std::sqrt;
using std::atan2;
using std::abs;
using std::log;
using std::sin;
using std::cos;
using std::sin;
using std::min;
using std::max;
using std::cerr;
using std::endl;
using std::flush;
using std::vector;
using std::lower_bound;
using std::upper_bound;
using std::distance;
using std::sort;
using std::unique;
using std::complex;
using std::isfinite;

template <typename T> class msub_t {
public:
  T   t;
  int j;
  int k;
  T   err;
  msub_t() {
    t = err = T(0);
    j = k = 0;
  }
  msub_t(const msub_t<T>& x) {
    *this = x;
  }
  ~msub_t() {
    ;
  }
  msub_t<T>& operator = (const msub_t<T>& x) {
    t   = T(x.t);
    j   = x.j;
    k   = x.k;
    err = T(x.err);
    return *this;
  }
  bool operator < (const msub_t<T>& x) const {
    return  t >  x.t ||
           (t == x.t && err <  x.err) ||
           (t == x.t && err == x.err && j <  x.j) ||
           (t == x.t && err == x.err && j == x.j && k < x.k);
  }
};

template <typename T> class match_t {
public:
  Eigen::Matrix<T, 3, 3> rot;
  Eigen::Matrix<T, 3, 1> offset;
  T                      ratio;
  T                      rdepth;
  vector<int>            dstpoints;
  vector<int>            srcpoints;
  T                      thresh;
  Eigen::Matrix<T, 2, 1> threshsize;
  match_t() {
    thresh = T(0);
    threshsize[0] = threshsize[1] = T(0);
  }
  match_t(const T& thresh, const T& h, const T& w) {
    this->thresh        = thresh;
    this->threshsize[0] = h;
    this->threshsize[1] = w;
  }
  match_t(const match_t<T>& other) {
    *this = other;
  }
  match_t<T>  operator ~ () const {
    match_t<T> result;
    result.rot    = rot.transpose();
    result.ratio  = T(1) / ratio;
    result.offset = - result.rot * offset * result.ratio;
    result.rdepth = rdepth;
    result.dstpoints  = srcpoints;
    result.srcpoints  = dstpoints;
    result.thresh     = thresh;
    result.threshsize = threshsize;
    return result;
  }
  match_t<T>  operator / (const match_t<T>& src) const {
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
  match_t<T>& operator = (const match_t<T>& other) {
    rot        = other.rot;
    offset     = other.offset;
    ratio      = other.ratio;
    rdepth     = other.rdepth;
    dstpoints  = other.dstpoints;
    srcpoints  = other.srcpoints;
    thresh     = other.thresh;
    threshsize = other.threshsize;
    return *this;
  }
  vector<int> ptr(const vector<int>& orig, const vector<int>& p) const {
    vector<int> res;
    res.reserve(p.size());
    for(int i = 0; i < p.size(); i ++)
      if(0 <= p[i] && p[i] < orig.size())
        res.emplace_back(orig[p[i]]);
      else
        res.push_back(- 1);
    return res;
  }
  vector<int> reversePtr(const vector<int>& orig, const vector<int>& p) const {
    vector<int> res;
    res.reserve(p.size());
    for(int i = 0; i < p.size(); i ++) {
      int j0(- 1);
      for(int j = 0; j < orig.size(); j ++)
        if(orig[j] == p[i]) {
          j0 = j;
          break;
        }
      res.emplace_back(j0);
    }
    return res;
  }
  vector<Eigen::Matrix<int, 3, 1> > hull(const vector<int>& orig, const vector<Eigen::Matrix<int, 3, 1> >& p) const {
    vector<Eigen::Matrix<int, 3, 1> > res;
    vector<int> work;
    res.reserve(p.size());
    work.reserve(p.size() * 3);
    for(int i = 0; i < p.size(); i ++)
      for(int j = 0; j < 3; j ++)
        work.emplace_back(p[i][j]);
    work = ptr(orig, work);
    for(int i = 0; i < p.size(); i ++) {
      Eigen::Matrix<int, 3, 1> tmp;
      for(int j = 0; j < 3; j ++)
        tmp[j] = work[i * 3 + j];
      res.emplace_back(tmp);
    }
    return res;
  }
  vector<Eigen::Matrix<int, 3, 1> > reverseHull(const vector<int>& orig, const vector<Eigen::Matrix<int, 3, 1> >& p) const {
    vector<Eigen::Matrix<int, 3, 1> > res;
    vector<int> work;
    res.reserve(p.size());
    work.reserve(p.size() * 3);
    for(int i = 0; i < p.size(); i ++)
      for(int j = 0; j < 3; j ++)
        work.emplace_back(p[i][j]);
    work = reversePtr(orig, work);
    for(int i = 0; i < p.size(); i ++) {
      Eigen::Matrix<int, 3, 1> tmp;
      for(int j = 0; j < 3; j ++)
        tmp[j] = work[i * 3 + j];
      res.emplace_back(tmp);
    }
    return res;
  }
  Eigen::Matrix<T, 3, 1> transform(const Eigen::Matrix<T, 3, 1>& x) const {
    return rot * x * ratio + offset;
  }
  bool operator < (const match_t<T>& x1) const {
    const T rratio(max(abs(   ratio), T(1) / abs(   ratio)));
    const T xratio(max(abs(x1.ratio), T(1) / abs(x1.ratio)));
    // N.B. parallel -> max, points ratio -> max.
    return rdepth < x1.rdepth || (rdepth == x1.rdepth && rratio < xratio);
  }
  bool operator != (const match_t<T>& x) const {
    const auto test(offset - x.offset);
    const auto roterr(rot * x.rot.transpose());
    return !(abs(T(1) - roterr(0, 0)) <= thresh) ||
           !(abs(T(1) - roterr(1, 1)) <= thresh) ||
           !(abs(T(1) - roterr(2, 2)) <= thresh) ||
           !(sqrt(test.dot(test) / 
               (offset.dot(offset) + x.offset.dot(x.offset))) /
                 sqrt(threshsize[0] * threshsize[1]) <= thresh) ||
           ratio * x.ratio < 0 ||
           !(abs(ratio - x.ratio) / sqrt(ratio * x.ratio) <= thresh);
  }
  bool operator == (const match_t<T>& x) const {
    return ! (*this != x);
  }
};


template <typename T> class matchPartialPartial {
public:
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T, 3, 3>                           Mat3x3;
  typedef Eigen::Matrix<T, 2, 2>                           Mat2x2;
  typedef Eigen::Matrix<T, 3, 1>                           Vec3;
  typedef complex<T> U;
  matchPartialPartial();
  ~matchPartialPartial();
  void init(const int& ndiv, const T& threshp, const T& threshs);
  
  vector<match_t<T> > match(const vector<Vec3>& shapebase, const vector<Vec3>& points);
  void match(const vector<Vec3>& shapebase, const vector<Vec3>& points, vector<match_t<T> >& result);
  vector<match_t<T> > elim(const vector<match_t<T> >& m, const Mat dst[3], const Mat src[3], const Mat& srcbump, const vector<Vec3>& srcpts, const T& thresh = T(4) / T(256));
  
  // theta resolution.
  int ndiv;
  // match points thresh in [0, 1].
  T   threshp;
  // match operator == thresh in [0, 1].
  T   threshs;
  
private:
  Vec3               makeG(const vector<Vec3>& in) const;
  vector<msub_t<T> > makeMsub(const vector<Vec3>& shapebase, const vector<Vec3>& points, const Vec3& gs, const Vec3& gp, const match_t<T>& m) const;
  void               complementMatch(match_t<T>& work, const vector<Vec3>& shapebase, const vector<Vec3>& points, const Vec3& gs, const Vec3& gp) const;
  T                  isElim(const match_t<T>& m, const Mat dst[3], const Mat tsrc[3], const vector<Vec3>& srcpts, const T& thresh);
  U   I;
  T   Pi;
  // match theta  thresh in [0, 1].
  T   thresh;
  // match ratio  thresh in [0, 1].
  T   thresht;
};

template <typename T> matchPartialPartial<T>::matchPartialPartial() {
  I  = sqrt(U(- T(1)));
  Pi = atan2(T(1), T(1)) * T(4);
  // rough match.
  init(20, .25, .1);
}

template <typename T> matchPartialPartial<T>::~matchPartialPartial() {
  ;
}

template <typename T> void matchPartialPartial<T>::init(const int& ndiv, const T& threshp, const T& threshs) {
  this->ndiv    = ndiv;
  this->thresh  = sin(T(2) * Pi / ndiv) / T(2);
  this->thresht = this->thresh;
  this->threshp = threshp;
  this->threshs = threshs;
  return;
}

template <typename T> vector<match_t<T> > matchPartialPartial<T>::match(const vector<Vec3>& shapebase, const vector<Vec3>& points) {
  vector<match_t<T> > result;
  match(shapebase, points, result);
  return result;
}

template <typename T> Eigen::Matrix<T, 3, 1> matchPartialPartial<T>::makeG(const vector<Vec3>& in) const {
  Vec3 result;
  result[0] = result[1] = result[2] = T(0);
  for(int i = 0; i < in.size(); i ++)
    result += in[i];
  return result / in.size();
}

template <typename T> vector<msub_t<T> > matchPartialPartial<T>::makeMsub(const vector<Vec3>& shapebase, const vector<Vec3>& points, const Vec3& gs, const Vec3& gp, const match_t<T>& m) const {
  vector<msub_t<T> > result;
  vector<Vec3> shapework;
  vector<Vec3> pointswork;
  shapework.reserve(shapebase.size());
  pointswork.reserve(points.size());
  for(int i = 0; i < shapebase.size(); i ++)
    shapework.push_back(shapebase[i] - gs);
  for(int i = 0; i < points.size(); i ++)
    pointswork.push_back(m.transform(points[i] - gp));
  // result.reserve(points.size() * shapebase.size());
  for(int k = 0; k < points.size(); k ++)
    for(int j = 0; j < shapebase.size(); j ++) {
      const Vec3& aj(shapework[j]);
      const Vec3& bk(pointswork[k]);
      const T     t(aj.dot(bk) / bk.dot(bk));
      const Vec3  lerr(aj - bk * t);
      const T     err(lerr.dot(lerr) / sqrt(aj.dot(aj) * bk.dot(bk) * t * t));
      // if t <= T(0), it's mirrored and this should not match.
      if(T(0) <= t && err <= thresh * thresh && isfinite(t) && isfinite(err)) {
        msub_t<T> work;
        work.t   = t;
        work.j   = j;
        work.k   = k;
        work.err = err;
        result.emplace_back(work);
      }
    }
  sort(result.begin(), result.end());
  return result;
}

template <typename T> void matchPartialPartial<T>::complementMatch(match_t<T>& work, const vector<Vec3>& shapebase, const vector<Vec3>& points, const Vec3& gs, const Vec3& gp) const {
  T num(0);
  T denom(0);
  for(int k = 0; k < work.dstpoints.size(); k ++) {
    const auto shapek(shapebase[work.dstpoints[k]] - gs);
    const auto pointk(work.transform(points[work.srcpoints[k]] - gp));
    denom += pointk.dot(pointk);
    num   += shapek.dot(pointk);
  }
  assert(T(0) <= num);
  work.ratio = num / denom;
  auto offset(work.offset);
  for(int k = 0; k < work.dstpoints.size(); k ++)
    offset += shapebase[work.dstpoints[k]] - work.transform(points[work.srcpoints[k]]);
  work.offset = offset / work.dstpoints.size();
  for(int k = 0; k < work.dstpoints.size(); k ++) {
    const auto err(shapebase[work.dstpoints[k]] - work.transform(points[work.srcpoints[k]]));
    work.rdepth += err.dot(err);
  }
  work.rdepth /= work.dstpoints.size();
  return;
}

template <typename T> void matchPartialPartial<T>::match(const vector<Vec3>& shapebase, const vector<Vec3>& points, vector<match_t<T> >& result) {
  // drot0 := I_3.
  Mat3x3 drot0;
  for(int k = 0; k < drot0.rows(); k ++)
    for(int l = 0; l < drot0.cols(); l ++)
      drot0(k, l) = (k == l ? T(1) : T(0));
  const Vec3 gs(makeG(shapebase));
  const Vec3 gp(makeG(points));
  Vec3 gd;
  gd[0] = gd[1] = gd[2] = T(0);
  for(int i = 0; i < shapebase.size(); i ++) {
    const auto diff(shapebase[i] - gs);
    gd[0] = max(gd[0], abs(diff[0]));
    gd[1] = max(gd[1], abs(diff[1]));
  }
  for(int i = 0; i < points.size(); i ++) {
    const auto diff(points[i] - gp);
    gd[0] = max(gd[0], abs(diff[0]));
    gd[1] = max(gd[1], abs(diff[1]));
  }
  // for each rotation, we can now handle t <= 0:
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int nd = 0; nd <= ndiv / 2; nd ++) {
    vector<T> ddiv;
    ddiv.resize(4, T(0));
    ddiv[0] = cos(2 * Pi * nd / ndiv);
    ddiv[1] = sin(2 * Pi * nd / ndiv);
    // with t < 0 match.
    for(int nd2 = 0; nd2 < ndiv; nd2 ++) {
      ddiv[2] = cos(2 * Pi * nd2 / ndiv);
      ddiv[3] = sin(2 * Pi * nd2 / ndiv);
      match_t<T> work0(threshs, abs(gd[0]), abs(gd[1]));
      work0.rot = drot0;
      for(int k = 0; k < ddiv.size() / 2; k ++) {
        Mat3x3 lrot;
        lrot((k    ) % 3, (k    ) % 3) =   ddiv[k * 2 + 0];
        lrot((k + 1) % 3, (k    ) % 3) =   ddiv[k * 2 + 1];
        lrot((k    ) % 3, (k + 1) % 3) = - ddiv[k * 2 + 1];
        lrot((k + 1) % 3, (k + 1) % 3) =   ddiv[k * 2 + 0];
        lrot((k + 2) % 3, (k    ) % 3) = T(0);
        lrot((k + 2) % 3, (k + 1) % 3) = T(0);
        lrot((k + 2) % 3, (k + 2) % 3) = T(1);
        lrot((k    ) % 3, (k + 2) % 3) = T(0);
        lrot((k + 1) % 3, (k + 2) % 3) = T(0);
        work0.rot = lrot * work0.rot;
      }
      work0.rdepth    = T(0);
      work0.offset[0] = work0.offset[1] = work0.offset[2] = T(0);
      work0.ratio     = T(1);
      // for each near matches:
      const auto msub(makeMsub(shapebase, points, gs, gp, work0));
      cerr << msub.size() << ":" << flush;
      if(msub.size() /
           T(min(shapebase.size(), points.size())) < threshp)
        continue;
      // get nearer matches:
      int t0(0);
      for( ; t0 < msub.size(); t0 ++) {
        auto work(work0);
        bool flagj[shapebase.size()];
        for(int kk = 0; kk < shapebase.size(); kk ++)
          flagj[kk] = false;
        bool flagk[points.size()];
        for(int kk = 0; kk < points.size(); kk ++)
          flagk[kk] = false;
        int tt(t0);
        for(int t1 = t0; t1 < msub.size(); t1 ++)
          if(!flagj[msub[t1].j] && !flagk[msub[t1].k] &&
             // N.B. abar:=sum(aj)/n, bbar:=P*sum(bk)/n,
             //   get condition sum||aj-abar||, sum||P*bk-bbar|| -> 0
             //   with this imcomplete set.
             // N.B. t >= 0 and msub is sorted by t0 > t1.
             (msub[t0].t - msub[t1].t) / msub[t0].t <= thresht) {
            work.dstpoints.emplace_back(msub[t1].j);
            work.srcpoints.emplace_back(msub[t1].k);
            flagj[msub[t1].j] = true;
            flagk[msub[t1].k] = true;
            tt = t1;
          }
        t0 = tt;
        // if it's good:
        if(threshp <= work.dstpoints.size() /
                        T(min(shapebase.size(), points.size()))) {
          complementMatch(work, shapebase, points, gs, gp);
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

template <typename T> vector<match_t<T> > matchPartialPartial<T>::elim(const vector<match_t<T> >& m, const Mat dst[3], const Mat src[3], const Mat& srcbump, const vector<Vec3>& srcpts, const T& thresh) {
  vector<match_t<T> > res(m);
  tilter<T> tilt;
  for(int i = 0; i < m.size(); i ++) {
    Mat tsrc[3];
    for(int j = 0; j < 3; j ++)
      tsrc[j] = tilt.tilt(src[j], srcbump, m[i]);
    res[i].rdepth *= isElim(m[i], dst, tsrc, srcpts, thresh);
  }
  sort(res.begin(), res.end());
  return res;
}

template <typename T> T matchPartialPartial<T>::isElim(const match_t<T>& m, const Mat dst[3], const Mat tsrc[3], const vector<Vec3>& srcpts, const T& thresh) {
  assert(dst[0].rows() == tsrc[0].rows() && dst[0].cols() == tsrc[0].cols());
  assert(dst[1].rows() == tsrc[1].rows() && dst[1].cols() == tsrc[2].cols());
  assert(dst[2].rows() == tsrc[2].rows() && dst[2].cols() == tsrc[2].cols());
  assert(dst[0].rows() ==  dst[1].rows() && dst[0].cols() ==  dst[1].cols());
  assert(dst[0].rows() ==  dst[2].rows() && dst[0].cols() ==  dst[2].cols());
  vector<T> diffs;
  for(int i = 0; i < m.srcpoints.size(); i ++) {
    const auto g(m.transform(srcpts[m.srcpoints[i]]));
    const int  y(g[0]);
    const int  x(g[1]);
    if(0 <= y && y < tsrc[0].rows() && 0 <= x && x < tsrc[0].cols()) {
      if(tsrc[0](y, x) == T(0) && tsrc[1](y, x) == T(0) && tsrc[2](y, x) == T(0))
        continue;
      T diff(0);
      for(int j = 0; j < 3; j ++)
        diff += pow(tsrc[j](y, x) - dst[j](y, x), T(2));
      diffs.push_back(sqrt(diff));
    }
  }
  sort(diffs.begin(), diffs.end());
  vector<T> ddiffs;
  for(int i = 1; i < diffs.size(); i ++)
    ddiffs.push_back(diffs[i] - diffs[i - 1]);
  sort(ddiffs.begin(), ddiffs.end());
  return ddiffs.size() ? T(1) - distance(ddiffs.begin(), upper_bound(ddiffs.begin(), ddiffs.end(), thresh)) / T(ddiffs.size()) : T(0);
}


template <typename T> class matchWholePartial {
public:
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T, 3, 3>                           Mat3x3;
  typedef Eigen::Matrix<T, 3, 1>                           Vec3;
  typedef complex<T> U;
  matchWholePartial();
  ~matchWholePartial();
  void init(const int& ndiv, const T& thresh, const T& threshp, const T& threshs);

  vector<vector<match_t<T> > > match(const vector<Vec3>& shapebase, const vector<vector<Vec3> >& points, const vector<Vec3>& origins);
private:
  U I;
  T Pi;
  T thresh;
  T threshp;
  T threshs;
};

template <typename T> matchWholePartial<T>::matchWholePartial() {
  I  = sqrt(U(- T(1)));
  Pi = atan2(T(1), T(1)) * T(4);
}

template <typename T> matchWholePartial<T>::~matchWholePartial() {
  ;
}

template <typename T> void matchWholePartial<T>::init(const int& ndiv, const T& thresh, const T& threshp, const T& threshs) {
  this->thresh  = thresh;
  this->threshp = threshp;
  this->threshs = threshs;
  return;
}

template <typename T> vector<vector<match_t<T> > > matchWholePartial<T>::match(const vector<Vec3>& shapebase, const vector<vector<Vec3> >& points, const vector<Vec3>& origins) {
  assert(points.size() == origins.size());
  vector<vector<match_t<T> > > pmatches;
  matchPartialPartial<T> pmatch(thresh, threshp, threshs);
  for(int i = 0; i < points.size(); i ++) {
    cerr << "matching partials : " << i << "/" << shapebase.size() << endl;
    pmatches.push_back(pmatch.match(shapebase, points[i]));
  }
  cerr << "XXX not implemented: detecting possible whole matches..." << endl;
  vector<vector<match_t<T> > > result;
  return result;
}

#define _MATCH_
#endif

