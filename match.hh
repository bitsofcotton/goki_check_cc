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

#if !defined(_SCAN_CONTEXT_)

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
  // XXX configure us:
  bool operator < (const match_t<T>& x1) const {
    const T rratio(max(abs(   ratio), T(1) / abs(   ratio)));
    const T xratio(max(abs(x1.ratio), T(1) / abs(x1.ratio)));
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
  
  // theta resolution.
  int ndiv;
  // match points thresh in [0, 1].
  T   threshp;
  // match operator == thresh in [0, 1].
  T   threshs;
  
private:
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

template <typename T> void matchPartialPartial<T>::match(const vector<Vec3>& shapebase, const vector<Vec3>& points, vector<match_t<T> >& result) {
  // drot0 := I_3.
  Mat3x3 drot0;
  for(int k = 0; k < drot0.rows(); k ++)
    for(int l = 0; l < drot0.cols(); l ++)
      drot0(k, l) = (k == l ? T(1) : T(0));
  // centering each.
  Vec3 gs, gp, gd;
  gs[0] = gs[1] = gs[2] = T(0);
  gp[0] = gp[1] = gp[2] = T(0);
  gd[0] = gd[1] = gd[2] = T(0);
  for(int i = 0; i < shapebase.size(); i ++)
    gs += shapebase[i];
  gs /= shapebase.size();
  for(int i = 0; i < points.size(); i ++)
    gp += points[i];
  gp /= points.size();
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
  // for each rotation:
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int nd = 0; nd < ndiv; nd ++) {
    Eigen::Matrix<T, 4, 1> ddiv;
    ddiv[0] = cos(2 * Pi * nd / ndiv - Pi);
    ddiv[1] = sin(2 * Pi * nd / ndiv - Pi);
    // with t < 0 match.
    for(int nd2 = 0; nd2 <= ndiv / 2; nd2 ++) {
      ddiv[2] = cos(2 * Pi * nd2 / ndiv - Pi / T(2));
      ddiv[3] = sin(2 * Pi * nd2 / ndiv - Pi / T(2));
      Mat3x3 drot1(drot0);
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
        drot1 = lrot * drot1;
      }
      // for each near matches:
      vector<msub_t<T> > msub;
      for(int k = 0; k < points.size(); k ++)
        for(int j = 0; j < shapebase.size(); j ++) {
          const Vec3 aj(shapebase[j]       - gs );
          const Vec3 bk(drot1 * (points[k] - gp));
          const T    t(aj.dot(bk) / bk.dot(bk));
          const Vec3 lerr(aj - bk * t);
          const T    err(lerr.dot(lerr) / (aj.dot(aj) + bk.dot(bk) * t * t));
          if(err <= thresh * thresh && isfinite(t) && isfinite(err)) {
            msub_t<T> work;
            work.t   = t;
            work.j   = j;
            work.k   = k;
            work.err = err;
            msub.push_back(work);
          }
        } 
      cerr << msub.size() << ":" << flush;
      if(msub.size() /
           T(min(shapebase.size(), points.size())) < threshp)
        continue;
      sort(msub.begin(), msub.end());
      // get nearer matches:
      int t0(0);
      for( ; t0 < msub.size(); t0 ++) {
        match_t<T> work(threshs, abs(gd[0]), abs(gd[1]));
        work.rot = drot1;
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
             abs(msub[t0].t - msub[t1].t) /
               max(abs(msub[t0].t), abs(msub[t1].t)) <= thresht &&
             T(0) <= msub[t0].t * msub[t1].t) {
            work.dstpoints.push_back(msub[t1].j);
            work.srcpoints.push_back(msub[t1].k);
            flagj[msub[t1].j] = true;
            flagk[msub[t1].k] = true;
            tt = t1;
          }
        t0 = tt;
        // if it's good:
        if(threshp <= work.dstpoints.size() /
                        T(min(shapebase.size(), points.size()))) {
          work.rdepth = T(0);
          T num(0), denom(0);
          for(int k = 0; k < work.dstpoints.size(); k ++) {
            const Vec3 shapek(shapebase[work.dstpoints[k]] - gs);
            const Vec3 pointk(work.rot * (points[work.srcpoints[k]] - gp));
            denom       += pointk.dot(pointk);
            num         += shapek.dot(pointk);
            work.rdepth += shapek.dot(shapek);
          }
          work.ratio     = num / denom;
          work.rdepth   /= work.dstpoints.size();
          work.offset[0] = T(0);
          work.offset[1] = T(0);
          work.offset[2] = T(0);
          for(int k = 0; k < work.dstpoints.size(); k ++)
            work.offset += shapebase[work.dstpoints[k]] - (work.rot * points[work.srcpoints[k]] * work.ratio);
          work.offset /= work.dstpoints.size();
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
              result.push_back(work);
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
    cerr << "matching partial polys : " << i << "/" << shapebase.size() << endl;
    pmatches.push_back(pmatch.match(shapebase, points[i]));
  }
  cerr << "detecting possible whole matches (not implemented now.)..." << endl;
  vector<vector<match_t<T> > > result;
  return result;
}


template <typename T> class reDig {
public:
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T, 3, 3>                           Mat3x3;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
  typedef Eigen::Matrix<T,   3, 1>            Vec3;
  typedef Eigen::Matrix<int, 3, 1>            Veci3;
  
  reDig();
  ~reDig();
  void init();
  Mat  emphasis(const Mat& dstimg, const Mat& dstbump, const Mat& srcimg, const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const vector<Veci3>& hull, const T& ratio, tilter<T>& tilt);
  Mat  replace(const Mat& dstimg, const Mat& dstbump, const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const vector<Veci3>& hull, const vector<Vec3>& srcrep, const match_t<T>& match2, const vector<Veci3>& hullrep);
  bool takeShape(vector<Vec3>& points, vector<Veci3>& tris, const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const vector<Veci3>& hull, const T& ratio);
  
  Mat  eliminate(const Mat& dstimg, const vector<Vec3>& dst, const vector<Veci3>& hull);
  Mat  showMatch(const Mat& dstimg, const vector<Vec3>& dst, const vector<Veci3>& hull, const T& emph = T(.2));
  
private:
  void drawMatchLine(Mat& map, const Vec3& lref0, const Vec3& lref1, const T& emph);
};

template <typename T> reDig<T>::reDig() {
  ;
}

template <typename T> reDig<T>::~reDig() {
  ;
}

template <typename T> void reDig<T>::init() {
  return;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> reDig<T>::emphasis(const Mat& dstimg, const Mat& srcimg, const Mat& srcbump, const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const vector<Eigen::Matrix<int, 3, 1> >& hull, const T& ratio, tilter<T>& tilt) {
  cerr << "m" << flush;
  vector<typename tilter<T>::Triangles> triangles;
  for(int i = 0; i < srcimg.rows() - 1; i ++)
    for(int j = 0; j < srcimg.cols() - 1; j ++) {
      triangles.push_back(tilt.makeTriangle(i, j, srcimg, srcbump, false));
      triangles.push_back(tilt.makeTriangle(i, j, srcimg, srcbump, true));
    }
  
  bool *checked;
  checked = new bool[triangles.size() * 3];
  for(int i = 0; i < triangles.size() * 3; i ++)
    checked[i] = false;
  
  cerr << "emphasis(" << hull.size() << ")" << endl;
  const auto rmatch(~ match);
  const auto offset(match.rot.transpose() * match.offset);
  for(int ii = 0; ii < hull.size(); ii ++) {
    const int&  i(hull[ii][0]);
    const int&  j(hull[ii][1]);
    const int&  k(hull[ii][2]);
    const Vec3& p0(dst[match.dstpoints[i]]);
    const Vec3& p1(dst[match.dstpoints[j]]);
    const Vec3& p2(dst[match.dstpoints[k]]);
    const Vec3  src00((src[match.srcpoints[i]] +
                       src[match.srcpoints[j]] +
                       src[match.srcpoints[k]]) / T(3));
    const Vec3  dst00((p0 + p1 + p2) / T(3));
    const Vec3  dst0(rmatch.ratio * rmatch.rot * dst00 + rmatch.offset);
    for(int l = 0; l < triangles.size(); l ++)
      if(!checked[l * 3] || !checked[l * 3 + 1] || !checked[l * 3 + 2])
        for(int ll = 0; ll < 3; ll ++) {
          if(checked[l * 3 + ll])
            continue;
          const Vec3& q(triangles[l].col(ll));
          if(tilt.sameSide2(p0, p1, p2, q) &&
             tilt.sameSide2(p1, p2, p0, q) &&
             tilt.sameSide2(p2, p0, p1, q)) {
            triangles[l].col(ll) += (dst00 * ratio - (src00 * match.ratio + offset) * (T(1) - ratio)) / match.ratio;
            checked[l * 3 + ll] = true;
          }
        }
  }
  delete[] checked;
  
  for(int i = 0; i < triangles.size(); i ++)
    triangles[i](1, 3) = triangles[i].col(4).dot(triangles[i].col(0));
  Mat I3(3, 3);
  Vec zero3(3);
  for(int i = 0; i < 3; i ++) {
    for(int j = 0; j < 3; j ++)
      I3(i, j) = (i == j ? T(1) : T(0));
    zero3[i] = T(0);
  }
  return (dstimg + tilt.tiltsub(dstimg, triangles, match.rot, I3, zero3, match.offset, match.ratio)) / T(2);
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> reDig<T>::eliminate(const Mat& dstimg, const vector<Vec3>& dst, const vector<Veci3>& hull) {
  cerr << "eliminate(" << hull.size() << ")" << endl;
  Mat result(dstimg);
  for(int ii = 0; ii < hull.size(); ii ++) {
    const int& i(hull[ii][0]);
    const int& j(hull[ii][1]);
    const int& k(hull[ii][2]);
    const Vec3 diff(dst[i] - dst[j]);
          Vec3 orth(dst[k] - dst[i]);
    orth -= orth.dot(diff) / diff.dot(diff) * diff;
    const T n2orth(sqrt(orth.dot(orth)));
    for(int l = 0; l < n2orth + 1; l ++)
      drawMatchLine(result, dst[i] + (dst[k] - dst[i]) * l / n2orth, dst[j] + (dst[k] - dst[j]) * l / n2orth);
  }
  return result;
}

template <typename T> void reDig<T>::drawMatchLine(Mat& map, const Vec3& lref0, const Vec3& lref1, const T& emph) {
  int idxm(0);
  int idxM(1);
  if(abs(lref1[idxM] - lref0[idxM]) < abs(lref1[idxm] - lref0[idxm])) {
    idxm = 1;
    idxM = 0;
  }
  Vec3 lm(lref1), lM(lref0);
  if(lM[idxM] < lm[idxM]) {
    lm = lref0;
    lM = lref1;
  }
  for(int i = 0; i <= lM[idxM] - lm[idxM]; i ++) {
    const auto gidx(lm + (lM - lm) * i / (lM[idxM] - lm[idxM]));
    map(max(0, min(int(gidx[0]), int(map.rows() - 1))),
        max(0, min(int(gidx[1]), int(map.cols() - 1)))) = emph;
  }
  return;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> reDig<T>::showMatch(const Mat& dstimg, const vector<Vec3>& dst, const vector<Veci3>& hull, const T& emph) {
  Mat map(dstimg.rows(), dstimg.cols());
  for(int i = 0; i < map.rows(); i ++)
    for(int j = 0; j < map.cols(); j ++)
      map(i, j) = T(0);
#if defined(_OPENMP)
#pragma omp paralell for schedule(static, 1)
#endif
  for(int k = 0; k < hull.size(); k ++) {
    drawMatchLine(map, dst[hull[k][0]], dst[hull[k][1]], emph);
    drawMatchLine(map, dst[hull[k][1]], dst[hull[k][2]], emph);
    drawMatchLine(map, dst[hull[k][2]], dst[hull[k][0]], emph);
  }
  return dstimg + map;
}

#define _SCAN_CONTEXT_
#endif

