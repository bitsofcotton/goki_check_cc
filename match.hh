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
using std::cos;
using std::sin;
using std::sort;
using std::unique;
using std::min;
using std::max;
using std::cerr;
using std::endl;
using std::flush;
using std::vector;
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
    result.offset = - result.rot * offset;
    result.ratio  = T(1) / ratio;
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
    return x1.rdepth < rdepth || (rdepth == x1.rdepth && rratio < xratio);
  }
  bool operator != (const match_t<T>& x) const {
    const auto test(offset - x.offset);
    const auto roterr(rot * x.rot.transpose());
    for(int l = 0; l < rot.rows(); l ++)
      if(!(abs(T(1) - roterr(l, l)) <= thresh))
        return true;
    if(!(sqrt(test.dot(test) / 
              (offset.dot(offset) + x.offset.dot(x.offset))) /
           sqrt(threshsize[0] * threshsize[1])
         <= thresh))
      return true;
    if(ratio * x.ratio < 0 ||
       !(abs((ratio - x.ratio) / sqrt(ratio * x.ratio)) <= thresh))
      return true;
    return false;
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
  void init(const int& ndiv, const T& thresh, const T& threshp, const T& threshs);
  
  vector<match_t<T> > match(const vector<Vec3>& shapebase, const vector<Vec3>& points);
  void match(const vector<Vec3>& shapebase, const vector<Vec3>& points, vector<match_t<T> >& result);
  
  int ndiv;
  T   thresh;
  T   threshp;
  T   thresht;
  T   threshs;
  
private:
  U   I;
  T   Pi;
};

template <typename T> matchPartialPartial<T>::matchPartialPartial() {
  I  = sqrt(U(- T(1)));
  Pi = atan2(T(1), T(1)) * T(4);
  // rough match.
  init(8, .125, .25, .25);
}

template <typename T> matchPartialPartial<T>::~matchPartialPartial() {
  ;
}

template <typename T> void matchPartialPartial<T>::init(const int& ndiv, const T& thresh, const T& threshp, const T& threshs) {
  this->ndiv    = ndiv;
  this->thresh  = thresh;
  this->threshp = threshp;
  this->thresht = thresh;
  this->threshs = threshs;
  return;
}

template <typename T> vector<match_t<T> > matchPartialPartial<T>::match(const vector<Vec3>& shapebase, const vector<Vec3>& points) {
  vector<match_t<T> > result;
  match(shapebase, points, result);
  return result;
}

template <typename T> void matchPartialPartial<T>::match(const vector<Vec3>& shapebase, const vector<Vec3>& points, vector<match_t<T> >& result) {
  Mat3x3 drot0;
  for(int k = 0; k < drot0.rows(); k ++)
    for(int l = 0; l < drot0.cols(); l ++)
      drot0(k, l) = (k == l ? T(1) : T(0));
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
    gd[0] = max(gd[0], abs((shapebase[i] - gs)[0]));
    gd[1] = max(gd[1], abs((shapebase[i] - gs)[1]));
  }
  for(int i = 0; i < points.size(); i ++) {
    gd[0] = max(gd[0], abs((points[i] - gp)[0]));
    gd[1] = max(gd[1], abs((points[i] - gp)[1]));
  }
  for(int i = 0; i < 3; i ++) {
    cerr << " matching" << flush;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int nd = 0; nd < ndiv; nd ++) {
      cerr << "." << flush;
      Eigen::Matrix<T, 4, 1> ddiv;
      ddiv[0] = cos(2 * Pi * nd / ndiv);
      ddiv[1] = sin(2 * Pi * nd / ndiv);
      for(int nd2 = 0; nd2 < ndiv; nd2 ++) {
        ddiv[2] = cos(2 * Pi * nd2 / ndiv);
        ddiv[3] = sin(2 * Pi * nd2 / ndiv);
        Mat3x3 drot1(drot0);
        for(int k = 0; k < ddiv.size() / 2; k ++) {
          Mat3x3 lrot;
          lrot((k + i    ) % 3, (k + i    ) % 3) =   ddiv[k * 2 + 0];
          lrot((k + i + 1) % 3, (k + i    ) % 3) =   ddiv[k * 2 + 1];
          lrot((k + i    ) % 3, (k + i + 1) % 3) = - ddiv[k * 2 + 1];
          lrot((k + i + 1) % 3, (k + i + 1) % 3) =   ddiv[k * 2 + 0];
          lrot((k + i + 2) % 3, (k + i    ) % 3) = T(0);
          lrot((k + i + 2) % 3, (k + i + 1) % 3) = T(0);
          lrot((k + i + 2) % 3, (k + i + 2) % 3) = T(1);
          lrot((k + i    ) % 3, (k + i + 2) % 3) = T(0);
          lrot((k + i + 1) % 3, (k + i + 2) % 3) = T(0);
          drot1 = lrot * drot1;
        }
        vector<msub_t<T> > msub;
        for(int k = 0; k < points.size(); k ++)
          for(int j = 0; j < shapebase.size(); j ++) {
            const Vec3 aj(shapebase[j]       - gs);
            const Vec3 bk(drot1 * (points[k] - gp));
            const T    t(aj.dot(bk) / bk.dot(bk));
            const Vec3 lerr(aj - bk * t);
            const T    err(lerr.dot(lerr) / (aj.dot(aj) + bk.dot(bk) * t * t));
            if(err <= thresh * thresh && isfinite(t) &&
               isfinite(err) && T(0) <= t) {
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
        for(int t0 = 0; t0 < msub.size(); t0 ++) {
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
               //   get condition sum||aj-abar||, sum||bk-bbar|| -> 0
               //   with this imcomplete set.
               abs(msub[t0].t - msub[t1].t) / abs(msub[t0].t) <= thresht) {
              work.dstpoints.push_back(msub[t1].j);
              work.srcpoints.push_back(msub[t1].k);
              flagj[msub[t1].j] = true;
              flagk[msub[t1].k] = true;
              tt = t1;
            }
          t0 = tt;
          if(threshp <= work.dstpoints.size() /
                          T(min(shapebase.size(), points.size()))) {
            Vec3 sbar, pbar;
            sbar[0] = sbar[1] = sbar[2] = T(0);
            pbar[0] = pbar[1] = pbar[2] = T(0);
            for(int k = 0; k < work.dstpoints.size(); k ++) {
              sbar += shapebase[work.dstpoints[k]];
              pbar += work.rot * points[work.srcpoints[k]];
            }
            sbar /= work.dstpoints.size();
            pbar /= work.srcpoints.size();
            T num(0), denom(0);
            for(int k = 0; k < work.dstpoints.size(); k ++) {
              const Vec3 shapek(shapebase[work.dstpoints[k]] - sbar);
              const Vec3 pointk(work.rot * points[work.srcpoints[k]] - pbar);
              num   += shapek.dot(pointk);
              denom += pointk.dot(pointk);
            }
            work.ratio = num / denom;
            work.offset[0] = T(0);
            work.offset[1] = T(0);
            work.offset[2] = T(0);
            for(int k = 0; k < work.dstpoints.size(); k ++)
              work.offset += shapebase[work.dstpoints[k]] - (work.rot * points[work.srcpoints[k]] * work.ratio);
            work.offset /= work.dstpoints.size();
            work.rdepth  = T(0);
            for(int k = 0; k < work.dstpoints.size(); k ++) {
              const Vec3 aj(shapebase[work.dstpoints[k]] - sbar);
              const Vec3 bk(work.rot * points[work.srcpoints[k]] * work.ratio + work.offset - sbar);
              work.rdepth += aj.dot(aj);
            }
            work.rdepth /= work.dstpoints.size();
#if defined(_OPENMP)
#pragma omp critical
#endif
            {
              int  idx(- 1);
              // we can't use STL because of operator < implementation.
              for(int k = 0; k < result.size(); k ++)
                if(result[k] == work && work < result[k]) {
                  idx = k;
                  break;
                }
              if(idx >= 0)
                result[idx] = work;
              else {
                result.push_back(work);
                cerr << "*" << flush;
              }
            }
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
  Vec3 emphasis0(const Vec3& dst, const Vec3& refdst, const Vec3& src, const match_t<T>& match, const T& ratio);
  Mat  emphasis(const Mat& dstimg, const Mat& dstbump, const Mat& srcimg, const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const vector<Veci3>& hull, const T& ratio);
  Mat  replace(const Mat& dstimg, const Mat& dstbump, const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const vector<Veci3>& hull, const vector<Vec3>& srcrep, const match_t<T>& match2, const vector<Veci3>& hullrep);
  bool takeShape(vector<Vec3>& points, vector<Veci3>& tris, const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const vector<Veci3>& hull, const T& ratio);
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

template <typename T> Eigen::Matrix<T, 3, 1> reDig<T>::emphasis0(const Vec3& dst, const Vec3& refdst, const Vec3& src, const match_t<T>& match, const T& ratio) {
  const Vec3& a(refdst);
  const Vec3  b(match.rot * src * match.ratio + match.offset);
  return dst + (b - a) * ratio;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> reDig<T>::emphasis(const Mat& dstimg, const Mat& dstbump, const Mat& srcimg, const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const vector<Eigen::Matrix<int, 3, 1> >& hull, const T& ratio) {
  cerr << " making triangles" << flush;
  tilter<T> tilt;
  vector<typename tilter<T>::Triangles> triangles;
  for(int i = 0; i < dstimg.rows() - 1; i ++)
    for(int j = 0; j < dstimg.cols() - 1; j ++) {
      triangles.push_back(tilt.makeTriangle(i, j, dstimg, dstbump, false));
      triangles.push_back(tilt.makeTriangle(i, j, dstimg, dstbump, true));
    }
  bool *checked;
  checked = new bool[triangles.size() * 3];
  for(int i = 0; i < triangles.size() * 3; i ++)
    checked[i] = false;
  
  cerr << "emphasis(" << hull.size() << ")" << flush;
  for(int ii = 0; ii < hull.size(); ii ++) {
    cerr << "." << flush;
    const int  i(hull[ii][0]);
    const int  j(hull[ii][1]);
    const int  k(hull[ii][2]);
    const Vec3 dst0((dst[match.dstpoints[i]] + dst[match.dstpoints[j]] + dst[match.dstpoints[k]]) / T(3));
    const Vec3 src0((src[match.srcpoints[i]] + src[match.srcpoints[j]] + src[match.srcpoints[k]]) / T(3));
    const Vec3& p0(src[match.srcpoints[i]]);
    const Vec3& p1(src[match.srcpoints[j]]);
    const Vec3& p2(src[match.srcpoints[k]]);
    for(int l = 0; l < triangles.size(); l ++) {
      if(!checked[l * 3] || !checked[l * 3 + 1] || !checked[l * 3 + 2]) {
        for(int ll = 0; ll < 3; ll ++) {
          if(checked[l * 3 + ll])
            continue;
          const Vec3& q(triangles[l].col(ll));
          if(tilt.sameSide2(p0, p1, p2, q) &&
             tilt.sameSide2(p1, p2, p0, q) &&
             tilt.sameSide2(p2, p0, p1, q)) {
            triangles[l].col(ll) = emphasis0(triangles[l].col(ll), dst0, src0, match, ratio);
            checked[l * 3 + ll]  = true;
          }
        }
      }
    }
  }
  const T zr(sqrt(T(dstimg.rows() * dstimg.cols())));
  for(int i = 0; i < triangles.size(); i ++) {
    cerr << "." << flush;
    for(int k = 0; k < 3; k ++)
      triangles[i](2, k) *= zr;
    triangles[i].col(4) = tilt.solveN(triangles[i].col(0), triangles[i].col(1), triangles[i].col(2));
    triangles[i](1, 3)  = triangles[i].col(4).dot(triangles[i].col(0));
  }
  
  Mat I3(3, 3);
  Vec zero3(3);
  for(int i = 0; i < 3; i ++) {
    for(int j = 0; j < 3; j ++)
      I3(i, j) = (i == j ? T(1) : T(0));
    zero3[i] = T(0);
  }
  delete[] checked;
  auto rmatch(~ match);
  // XXX inverse ratio, but better works.
  return (srcimg * ratio + tilt.tiltsub(dstimg, triangles, rmatch.rot, I3, zero3, rmatch.offset, rmatch.ratio) * (T(1) - ratio)) / T(2);
}

template <typename T> void drawMatchLine(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& map, const Eigen::Matrix<T, 3, 1>& lref0, const Eigen::Matrix<T, 3, 1>& lref1, const T& emph, const T& epsilon = T(1e-4)) {
  if(abs(lref1[0] - lref0[0]) <= epsilon) {
    int sgndelta(1);
    if(lref1[1] < lref0[1])
      sgndelta = - 1;
    for(int i = lref0[1]; (lref1[1] - i) * sgndelta > 0; i += sgndelta)
      map(max(0, min(int(lref0[0]), int(map.rows() - 1))),
          max(0, min(i,             int(map.cols() - 1)))) = emph;
  } else if(abs(lref1[1] - lref0[1]) <= epsilon) {
    int sgndelta(1);
    if(lref1[0] < lref0[0])
      sgndelta = - 1;
    for(int j = lref0[0]; (lref1[0] - j) * sgndelta > 0; j += sgndelta) {
      map(max(0, min(j,             int(map.rows() - 1))),
          max(0, min(int(lref0[1]), int(map.cols() - 1)))) = emph;
    }
  } else if(abs(lref1[1] - lref0[1]) > abs(lref1[0] - lref0[0])) {
    const T tilt((lref1[1] - lref0[1]) / (lref1[0] - lref0[0]));
    int sgndelta(1);
    if(lref1[0] < lref0[0])
      sgndelta = - 1;
    int i = lref0[0], ii = 0;
    for(; i != int(lref1[0]) &&
          (i < 0 || map.rows() <= i ||
           lref0[1] + tilt * ii < 0 ||
           map.rows() <= lref0[1] + tilt * ii) &&
          abs(tilt * ii) <= abs(lref0[1]) + T(map.rows());
          i += sgndelta, ii += sgndelta);
    for(; i != int(lref1[0]) &&
          0 <= lref0[1] + tilt * ii &&
               lref0[1] + tilt * ii <= map.cols();
          i += sgndelta, ii += sgndelta)
      for(int jj = 0; jj < abs(tilt); jj ++) {
        int j = tilt * (ii - sgndelta) +
                jj * (tilt * sgndelta < T(0) ? - 1 : 1);
        map(max(0, min(i, int(map.rows() - 1))),
            max(0, min(int(lref0[1] + j), int(map.cols() - 1)))) = emph;
      }
  } else {
    const T tilt((lref1[0] - lref0[0]) / (lref1[1] - lref0[1]));
    int sgndelta(1);
    if(lref1[1] < lref0[1])
      sgndelta = - 1;
    int j = lref0[1], jj = 0;
    for(; j != int(lref1[1]) &&
          (j < 0 || map.cols() <= j ||
           lref0[0] + tilt * jj < 0 ||
           map.cols() <= lref0[0] + tilt * jj) &&
          abs(tilt * jj) <= abs(lref0[0]) + T(map.rows());
          j += sgndelta, jj += sgndelta);
    for(; j != int(lref1[1]) &&
          0 <= lref0[0] + tilt * jj &&
               lref0[0] + tilt * jj <= map.rows();
          j += sgndelta, jj += sgndelta)
      for(int ii = 0; ii < abs(tilt); ii ++) {
        int i = tilt * (jj - sgndelta) +
                ii * (tilt * sgndelta < T(0) ? - 1 : 1);
        map(max(0, min(int(lref0[0] + i), int(map.rows() - 1))),
            max(0, min(j, int(map.cols() - 1)))) = emph;
      }
  }
  return;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> showMatch(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& input, const vector<Eigen::Matrix<T, 3, 1> >& refpoints, const vector<Eigen::Matrix<int, 3, 1> >& prefpoints, const T& emph = T(.8)) {
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> map(input.rows(), input.cols());
  for(int i = 0; i < map.rows(); i ++)
    for(int j = 0; j < map.cols(); j ++)
      map(i, j) = T(0);
#if defined(_OPENMP)
#pragma omp paralell for schedule(static, 1)
#endif
  for(int k = 0; k < prefpoints.size(); k ++) {
    drawMatchLine<T>(map,
                     refpoints[prefpoints[k][0]],
                     refpoints[prefpoints[k][1]],
                     emph);
    drawMatchLine<T>(map,
                     refpoints[prefpoints[k][1]],
                     refpoints[prefpoints[k][2]],
                     emph);
    drawMatchLine<T>(map,
                     refpoints[prefpoints[k][2]],
                     refpoints[prefpoints[k][0]],
                     emph);
  }
  return input + map;
}

#define _SCAN_CONTEXT_
#endif

