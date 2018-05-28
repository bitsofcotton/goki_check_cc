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

#if !defined(_REDIG_)

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

template <typename T> bool less0(const T& x, const T& y) {
  return x.first[0] < y.first[0] || (x.first[0] == y.first[0] && x.first[1] < y.first[1]);
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
  Mat  emphasis(const Mat& dstimg, const Mat& dstbump, const Mat& srcimg, const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const vector<Veci3>& hulldst, const vector<Veci3>& hullsrc, const T& ratio, tilter<T>& tilt);
  Mat  replace(const Mat& dstimg, const vector<Vec3>& src, const match_t<T>& match, const vector<Veci3>& hullsrc);
  Mat  replace(const Mat& dstimg, const Mat& srcimg, const Mat& dstbump, const Mat& srcbump, const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const vector<Veci3>& hullsrc);
  vector<Vec3> takeShape(const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const vector<Veci3>& hulldst, const vector<Veci3>& hullsrc, const T& ratio);
  Mat  showMatch(const Mat& dstimg, const vector<Vec3>& dst, const vector<Veci3>& hull, const T& emph = T(.2));
  Mat  makeRefMatrix(const Mat& orig, const int& start) const;
  Mat  pullRefMatrix(const Mat& ref, const int& start, const Mat& orig) const;
  vector<Veci3> delaunay2(const vector<Vec3>& p, const vector<int>& pp, const T& epsilon = T(1e-5), const int& mdiv = 300) const;
  void maskVectors(vector<Vec3>& points, vector<Veci3>& polys, const Mat& mask);
  vector<vector<int> > getEdges(const Mat& mask, const vector<Vec3>& points, const int& vbox);
  Mat  rgb2l(const Mat rgb[3]);
  Mat  rgb2xz(const Mat rgb[3]);
  Mat  tilt45(const Mat& in, const bool& invert, const Mat& orig = Mat());
  void normalize(Mat data[3], const T& upper);

private:
  void drawMatchLine(Mat& map, const Vec3& lref0, const Vec3& lref1, const T& emph);
  void drawMatchTriangle(Mat& map, const Vec3& lref0, const Vec3& lref1, const Vec3& lref2);
  bool isDelaunay2(T& cw, const Vec3 p[4], const T& epsilon) const;
  bool isCrossing(const Vec3& p0, const Vec3& p1, const Vec3& q0, const Vec3& q1, const T& err = T(1e-4)) const;
  void floodfill(Mat& checked, vector<pair<int, int> >& store, const Mat& mask, const int& y, const int& x);
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

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> reDig<T>::emphasis(const Mat& dstimg, const Mat& srcimg, const Mat& srcbump, const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const vector<Veci3>& hulldst, const vector<Veci3>& hullsrc, const T& ratio, tilter<T>& tilt) {
  cerr << "m" << flush;
  assert(hulldst.size() == hullsrc.size());
  vector<typename tilter<T>::Triangles> triangles;
  triangles.reserve((srcimg.rows() - 1) * (srcimg.cols() - 1) * 2);
  for(int i = 0; i < srcimg.rows() - 1; i ++)
    for(int j = 0; j < srcimg.cols() - 1; j ++) {
      triangles.push_back(tilt.makeTriangle(i, j, srcimg, srcbump, false));
      triangles[triangles.size() - 1].c = srcimg(i, j);
      triangles.push_back(tilt.makeTriangle(i, j, srcimg, srcbump, true));
      triangles[triangles.size() - 1].c = srcimg(i, j);
    }
  
  bool *checked = new bool[triangles.size() * 3];
  assert(checked);
  for(int i = 0; i < triangles.size() * 3; i ++)
    checked[i] = false;
  
  cerr << "e(" << hulldst.size() << ")" << endl;
  const auto rmatch(~ match);
  for(int i = 0; i < hulldst.size(); i ++) {
    // XXX matching places inversion needed?
    const Vec3 p0(rmatch.transform(dst[hulldst[i][0]]));
    const Vec3 p1(rmatch.transform(dst[hulldst[i][1]]));
    const Vec3 p2(rmatch.transform(dst[hulldst[i][2]]));
    const Vec3 src0((src[hullsrc[i][0]] +
                     src[hullsrc[i][1]] +
                     src[hullsrc[i][2]]) / T(3));
    const Vec3 dst0((p0 + p1 + p2) / T(3));
    for(int l = 0; l < triangles.size(); l ++)
      for(int ll = 0; ll < 3; ll ++) {
        const Vec3& q(triangles[l].p.col(ll));
        if(!checked[l * 3 + ll] &&
           tilt.sameSide2(p0, p1, p2, q) &&
           tilt.sameSide2(p1, p2, p0, q) &&
           tilt.sameSide2(p2, p0, p1, q)) {
          triangles[l].p.col(ll) += (dst0 - src0) * ratio / match.ratio;
          checked[l * 3 + ll] = true;
        }
      }
  }
  vector<typename tilter<T>::Triangles> wt;
  wt.reserve(triangles.size());
  for(int i = 0; i < triangles.size(); i ++)
    if(checked[i * 3 + 0] || checked[i * 3 + 1] || checked[i * 3 + 2])
      wt.push_back(triangles[i]);
  delete[] checked;
  return tilt.tilt(dstimg, wt, match);
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> reDig<T>::replace(const Mat& dstimg, const vector<Vec3>& src, const match_t<T>& match, const vector<Veci3>& hullsrc) {
  Mat result(dstimg);
  vector<Vec3> tsrc;
  T            M(0);
  T            m(0);
  for(int i = 0; i < src.size(); i ++) {
    tsrc.push_back(match.transform(src[i]));
    if(i) {
      M = max(M, tsrc[i][2]);
      m = min(m, tsrc[i][2]);
    } else
      M = m = tsrc[i][2];
  }
  if(M - m != T(0))
    for(int i = 0; i < tsrc.size(); i ++)
      tsrc[i][2] = (tsrc[i][2] - m) / (M - m);
  for(int ii = 0; ii < hullsrc.size(); ii ++)
    drawMatchTriangle(result, tsrc[hullsrc[ii][0]],
                              tsrc[hullsrc[ii][1]],
                              tsrc[hullsrc[ii][2]]);
  return result;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> reDig<T>::replace(const Mat& dstimg, const Mat& srcimg, const Mat& dstbump, const Mat& srcbump, const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const vector<Veci3>& hullsrc) {
  Mat result(dstimg);
  // stub.
  return result;
}

template <typename T> vector<Eigen::Matrix<T, 3, 1> > reDig<T>::takeShape(const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const vector<Veci3>& hulldst, const vector<Veci3>& hullsrc, const T& ratio) {
  assert(hulldst.size() == hullsrc.size());
  bool* checked = new bool[hullsrc.size() * 3];
  assert(checked);
  for(int i = 0; i < hullsrc.size() * 3; i ++)
    checked[i] = false;
  vector<Vec3> result(dst);
  tilter<T> tilt;
  const auto rmatch(~ match);
  for(int i = 0; i < hullsrc.size(); i ++) {
    const auto& p0(src[hullsrc[i][0]]);
    const auto& p1(src[hullsrc[i][1]]);
    const auto& p2(src[hullsrc[i][2]]);
    for(int j = 0; j < 3; j ++) {
      const auto q(rmatch.transform(dst[hulldst[i][j]]));
      if(!checked[i * 3 + j] &&
         tilt.sameSide2(p0, p1, p2, q) &&
         tilt.sameSide2(p1, p2, p0, q) &&
         tilt.sameSide2(p2, p0, p1, q)) {
        result[hulldst[i][j]] += (match.transform(src[hullsrc[i][j]]) - dst[hulldst[i][j]]) * ratio;
        checked[i * 3 + j] = true;
      }
    }
  }
  delete[] checked;
  return result;
}

template <typename T> void reDig<T>::drawMatchLine(Mat& map, const Vec3& lref0, const Vec3& lref1, const T& emph) {
  int idxm(0);
  int idxM(1);
  if(abs(lref1[idxM] - lref0[idxM]) < abs(lref1[idxm] - lref0[idxm])) {
    idxm = 1;
    idxM = 0;
  }
  for(int i = 0; i <= abs(lref0[idxM] - lref1[idxM]); i ++) {
    const auto gidx(lref0 + (lref1 - lref0) * i / abs(lref0[idxM] - lref1[idxM]));
    map(max(0, min(int(gidx[0]), int(map.rows() - 1))),
        max(0, min(int(gidx[1]), int(map.cols() - 1)))) = emph;
  }
  return;
}

template <typename T> void reDig<T>::drawMatchTriangle(Mat& map, const Vec3& lref0, const Vec3& lref1, const Vec3& lref2) {
  int idxm(0);
  int idxM(1);
  if(abs(lref1[idxM] - lref0[idxM]) < abs(lref1[idxm] - lref0[idxm])) {
    idxm = 1;
    idxM = 0;
  }
  const Vec3 ldiff0(lref0 - lref1);
        Vec3 ldiff(lref2 - lref0);
  ldiff -= ldiff.dot(ldiff0) * ldiff0 / ldiff0.dot(ldiff0);
  const T    lnum(sqrt(ldiff.dot(ldiff)) + 1);
  for(int k = 0; k < lnum; k ++) {
    const Vec3 l0(lref0 + (lref2 - lref0) * k / int(lnum));
    const Vec3 l1(lref1 + (lref2 - lref1) * k / int(lnum));
    for(int i = 0; i <= int(abs(l0[idxM] - l1[idxM]) + 1); i ++) {
      const auto gidx(l0 + (l1 - l0) * i / int(abs(l0[idxM] - l1[idxM]) + 1));
      map(max(0, min(int(gidx[0]), int(map.rows() - 1))),
          max(0, min(int(gidx[1]), int(map.cols() - 1)))) = gidx[2];
    }
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

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> reDig<T>::makeRefMatrix(const Mat& orig, const int& start) const {
  Mat result(orig.rows(), orig.cols());
  for(int i = 0; i < orig.rows() * orig.cols(); i ++)
    result(i % orig.rows(), i / orig.rows()) = i + start;
  return result;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> reDig<T>::pullRefMatrix(const Mat& ref, const int& start, const Mat& orig) const {
  assert(orig.rows() == ref.rows() && orig.cols() == ref.cols());
  Mat result(orig.rows(), orig.cols());
  for(int i = 0; i < ref.rows() * ref.cols(); i ++) {
    const int ly(i % ref.rows());
    const int lx(i / ref.rows());
    const int v(ref(ly, lx) - start);
    if(0 <= v && v < orig.rows() * orig.cols())
      result(ly, lx) = orig(v % orig.rows(), v / orig.rows());
    else
      result(ly, lx) = T(0);
  }
  return result;
}

template <typename T> vector<Eigen::Matrix<int, 3, 1> > reDig<T>::delaunay2(const vector<Vec3>& p, const vector<int>& pp, const T& epsilon, const int& mdiv) const {
  vector<Veci3> res;
  tilter<T> tilt;
  cerr << pp.size() << ":" << flush;
  if(pp.size() > mdiv) {
    vector<pair<Vec3, int> > div;
    div.reserve(pp.size());
    for(int i = 0; i < pp.size(); i ++)
      div.emplace_back(make_pair(p[pp[i]], pp[i]));
    sort(div.begin(), div.end(), less0<pair<Vec3, int> >);
    vector<int> lo, mid, hi;
    lo.reserve( div.size() / 2);
    mid.reserve(div.size() / 2);
    hi.reserve( div.size() / 2);
    for(int i = 0; i < div.size() / 2; i ++)
      lo.emplace_back(div[i].second);
    for(int i = div.size() / 4; i < div.size() * 3 / 4; i ++)
      mid.emplace_back(div[i].second);
    for(int i = div.size() / 2; i < div.size(); i ++)
      hi.emplace_back(div[i].second);
    const auto left(  delaunay2(p, lo,  epsilon));
    const auto middle(delaunay2(p, mid, epsilon));
    const auto right( delaunay2(p, hi,  epsilon));
    vector<Veci3> work;
    work.reserve(left.size() + right.size());
    res.reserve(middle.size());
    for(int i = 0; i < left.size(); i ++) {
      for(int ii = 0; ii < 3; ii ++) {
        const auto itr(upper_bound(div.begin(), div.end(), make_pair(p[left[i][ii]], left[i][ii]), less0<pair<Vec3, int> >));
        if(itr != div.end() && itr->second == left[i][ii] &&
           div.size() / 2 < distance(div.begin(), itr) < div.size())
          goto nextl;
      }
      work.emplace_back(left[i]);
     nextl:
      ;
    }
    for(int i = 0; i < right.size(); i ++) {
      for(int ii = 0; ii < 3; ii ++) {
        const auto itr(upper_bound(div.begin(), div.end(), make_pair(p[right[i][ii]], right[i][ii]), less0<pair<Vec3, int> >));
        if(itr != div.end() && itr->second == right[i][ii] &&
           distance(div.begin(), itr) < div.size() / 2)
          goto nextr;
      }
      work.emplace_back(right[i]);
     nextr:
      ;
    }
    for(int i = 0; i < middle.size(); i ++) {
      for(int ii = 0; ii < 3; ii ++) {
        const auto itr(upper_bound(div.begin(), div.end(), make_pair(p[middle[i][ii]], middle[i][ii]), less0<pair<Vec3, int> >));
        if(itr != div.end() && itr->second == middle[i][ii] &&
           (distance(div.begin(), itr) < div.size() / 4 ||
            div.size() * 3 / 4 < distance(div.begin(), itr) ) )
          goto next;
      }
      res.emplace_back(middle[i]);
     next:
      ;
    }
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int ii = 0; ii < work.size(); ii ++) {
      const int& i(work[ii][0]);
      const int& j(work[ii][1]);
      const int& k(work[ii][2]);
      T cw;
      Eigen::Matrix<int, 3, 1> idx;
      Vec3 q[4];
      for(int jj = 0; jj < middle.size(); jj ++)
        for(int i0 = 0; i0 < 3; i0 ++)
          for(int j0 = 0; j0 < 3; j0 ++)
            if(isCrossing(p[work[ii][ i0      % 3]],
                          p[work[ii][(i0 + 1) % 3]],
                          p[middle[jj][ j0      % 3]],
                          p[middle[jj][(j0 + 1) % 3]]))
              goto fixnext0;
      q[0] = p[i]; q[1] = p[j]; q[2] = p[k];
      q[3] = p[middle[0][0]];
      isDelaunay2(cw, q, epsilon);
      idx[0] = i;
      if(cw < 0) {
        idx[1] = k;
        idx[2] = j;
      } else {
        idx[1] = j;
        idx[2] = k;
      }
#if defined(_OPENMP)
#pragma omp atomic
#endif
      res.emplace_back(idx);
     fixnext0:
      ;
    }
  } else {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < pp.size(); i ++)
      for(int j = i + 1; j < pp.size(); j ++)
        for(int k = j + 1; k < pp.size(); k ++) {
          T cw;
          Eigen::Matrix<int, 3, 1> idx;
          Vec3 q[4];
          q[0] = p[pp[i]]; q[1] = p[pp[j]]; q[2] = p[pp[k]];
          for(int l = 0; l < pp.size(); l ++) {
            q[3] = p[pp[l]];
            if(!isDelaunay2(cw, q, epsilon))
              goto fixnext;
          }
          idx[0] = pp[i];
          if(cw < 0) {
            idx[1] = pp[k];
            idx[2] = pp[j];
          } else {
            idx[1] = pp[j];
            idx[2] = pp[k];
          }
#if defined(_OPENMP)
#pragma omp critical
#endif
          {
            for(int i = 0; i < res.size(); i ++)
              for(int j = 0; j < res[i].size(); j ++)
                for(int k = 0; k < idx.size(); k ++)
                  if(isCrossing(p[res[i][(j + 0) % 3]],
                                p[res[i][(j + 1) % 3]],
                                p[idx[(k + 0) % 3]],
                                p[idx[(k + 1) % 3]]))
                    goto fixnext;
            res.emplace_back(idx);
          }
         fixnext:
          ;
        }
  }
  return res;
}

template <typename T> bool reDig<T>::isDelaunay2(T& cw, const Vec3 p[4], const T& epsilon) const {
  // sameline?
  Vec3 bcn(p[1] - p[2]);
  bcn[2] = T(0);
  Vec3 err(p[0]);
  err[2] = T(0);
  err   -= err.dot(bcn) * bcn / bcn.dot(bcn);
  if(err.dot(err) <= epsilon)
    return false;
  Eigen::Matrix<T, 4, 4> dc;
  const auto g((p[0] + p[1] + p[2]) / T(3));
  for(int i = 0; i < 4; i ++) {
    dc(i, 0) = T(1);
    dc(i, 1) = p[i][0] - g[0];
    dc(i, 2) = p[i][1] - g[1];
    dc(i, 3) = dc(i, 1) * dc(i, 1) + dc(i, 2) * dc(i, 2);
  }
  Eigen::Matrix<T, 3, 3> dc0;
  for(int i = 0; i < 3; i ++) {
    dc0(i, 0) = T(1);
    dc0(i, 1) = p[i][0] - g[0];
    dc0(i, 2) = p[i][1] - g[1];
  }
  cw = dc0.determinant();
  if(abs(cw) <= pow(epsilon, T(3)))
    cw = T(0);
  else if(cw < T(0))
    cw = - T(1);
  else if(T(0) < cw)
    cw =   T(1);
  tilter<T> tilt;
  if(cw * dc.determinant() < - epsilon ||
     (tilt.sameSide2(p[0], p[1], p[2], p[3], false) &&
      tilt.sameSide2(p[1], p[2], p[0], p[3], false) &&
      tilt.sameSide2(p[2], p[0], p[1], p[3], false)) )
    return false;
  return true;
}

template <typename T> bool reDig<T>::isCrossing(const Vec3& p0, const Vec3& p1, const Vec3& q0, const Vec3& q1, const T& err) const {
  // t * p0 + (1 - t) * p1 == s * q0 + (1 - s) * q1
  // <=> p1 + (p0 - p1) t == q1 + (q0 - q1) s
  // <=> [(p0 - p1), (q1 - q0)][t, s] == q1 - p1.
  // <=> Ax==b.
  Eigen::Matrix<T, 2, 2> A;
  Eigen::Matrix<T, 2, 1> b;
  A(0, 0) = p0[0] - p1[0];
  A(1, 0) = p0[1] - p1[1];
  A(0, 1) = q1[0] - q0[0];
  A(1, 1) = q1[1] - q0[1];
  b[0]    = q1[0] - p1[0];
  b[1]    = q1[1] - p1[1];
  if(abs(A.determinant()) <= pow(err, T(2)))
    return false;
  auto x(A.inverse() * b);
  return (err <= x[0] && x[0] <= T(1) - err) &&
         (err <= x[1] && x[1] <= T(1) - err);
}

template <typename T> void reDig<T>::maskVectors(vector<Vec3>& points, vector<Veci3>& polys, const Mat& mask) {
  std::vector<int> elim, elimp, after;
  for(int i = 0, ii = 0; i < points.size(); i ++) {
    const int y(std::max(std::min(int(points[i][0]), int(mask.rows() - 1)), 0));
    const int x(std::max(std::min(int(points[i][1]), int(mask.cols() - 1)), 0));
    if(mask(y, x) > .5) {
      elim.emplace_back(i);
      after.push_back(- 1);
    } else
      after.push_back(ii ++);
  }
  for(int i = 0; i < polys.size(); i ++)
    if(std::binary_search(elim.begin(), elim.end(), polys[i][0]) ||
       std::binary_search(elim.begin(), elim.end(), polys[i][1]) ||
       std::binary_search(elim.begin(), elim.end(), polys[i][2]))
      elimp.emplace_back(i);
  for(int i = 0, j = 0; i < elim.size(); i ++)
    points.erase(points.begin() + (elim[i] - i));
  for(int i = 0; i < polys.size(); i ++)
    for(int j = 0; j < polys[i].size(); j ++)
      polys[i][j] = after[polys[i][j]];
  for(int i = 0; i < elimp.size(); i ++)
    polys.erase(polys.begin() + (elimp[i] - i));
  return;
}

template <typename T> void reDig<T>::floodfill(Mat& checked, vector<pair<int, int> >& store, const Mat& mask, const int& y, const int& x) {
  assert(mask.rows() == checked.rows() && mask.cols() == checked.cols());
  vector<pair<int, int> > tries;
  tries.emplace_back(make_pair(+ 1,   0));
  tries.emplace_back(make_pair(  0, + 1));
  tries.emplace_back(make_pair(- 1,   0));
  tries.emplace_back(make_pair(  0, - 1));
  vector<pair<int, int> > stack;
  stack.emplace_back(make_pair(y, x));
  while(stack.size()) {
    const auto pop(stack[stack.size() - 1]);
    stack.pop_back();
    const int& yy(pop.first);
    const int& xx(pop.second);
    if(! (0 <= yy && yy < checked.rows() && 0 <= xx && xx < checked.cols()) )
      store.emplace_back(pop);
    else if(!checked(yy, xx)) {
      checked(yy, xx) = true;
      if(T(.5) < mask(yy, xx))
        store.emplace_back(pop);
      else
        for(int i = 0; i < tries.size(); i ++)
          stack.emplace_back(make_pair(yy + tries[i].first, xx + tries[i].second));
    }
  }
  return;
}

template <typename T> vector<vector<int> > reDig<T>::getEdges(const Mat& mask, const vector<Vec3>& points, const int& vbox) {
  cerr << "getEdges" << flush;
  Mat checked(mask.rows(), mask.cols());
  for(int i = 0; i < checked.rows(); i ++)
    for(int j = 0; j < checked.cols(); j ++)
      checked(i, j) = false;
  vector<pair<int, int> > store;
  for(int i = 0; i < checked.rows(); i ++)
    for(int j = 0; j < checked.cols(); j ++)
      if(mask(i, j) < T(.5))
        floodfill(checked, store, mask, i, j);
  sort(store.begin(), store.end());
  store.erase(unique(store.begin(), store.end()), store.end());
  cerr << " with " << store.size() << " edge points " << flush;
  
  vector<vector<int> > result;
  // stored indices.
  vector<int>          se;
  for( ; se.size() < store.size(); ) {
    // tree index.
    vector<int> e;
    int         i(0);
    for( ; i < store.size(); i ++)
      if(!binary_search(se.begin(), se.end(), i))
        break;
    if(store.size() <= i)
      break;
    // get edge point tree.
    for( ; i < store.size(); ) {
      // get nearest points.
      vector<int> si;
      for(int j = 0; j < store.size(); j ++)
        if(!binary_search(se.begin(), se.end(), j) &&
           abs(store[i].first  - store[j].first)  <= 1 &&
           abs(store[i].second - store[j].second) <= 1)
          si.emplace_back(j);
      if(!si.size())
        break;
      // normal vector direction.
      sort(si.begin(), si.end());
      int j(0);
      for( ; j < si.size(); j ++) {
        const auto& sti(store[i]);
        const auto& stj(store[si[j]]);
        const T y(stj.first  - sti.first);
        const T x(stj.second - stj.second);
        const T x2(x - y + sti.second);
        const T y2(x + y + sti.first);
        if(0 <= x2 && x2 < mask.cols() && 0 <= y2 && y2 < mask.rows() &&
           mask(y2, x2) < T(.5))
          break;
      }
      if(si.size() <= j)
        j = 0;
      // store.
      i = si[j];
      e.emplace_back(i);
      se.emplace_back(i);
      sort(se.begin(), se.end());
    }
    // N.B. almost bruteforce...
    if(1 < e.size()) {
      result.push_back(vector<int>());
      // apply to points index.
      vector<int> pj;
      for(int i = 0; i < e.size(); i ++) {
        const auto& s(store[e[i]]);
        vector<pair<T, int> > distances;
        for(int j = 0; j < points.size(); j ++)
          distances.emplace_back(make_pair(
                                  pow(T(s.first  - points[j][0]), T(2)) +
                                    pow(T(s.second - points[j][1]), T(2)),
                                  j));
        sort(distances.begin(), distances.end());
        if(distances.size() && !binary_search(pj.begin(), pj.end(), distances[0].second)) {
          result[result.size() - 1].emplace_back(distances[0].second);
          pj.emplace_back(distances[0].second);
        }
      }
      result[result.size() - 1].emplace_back(result[result.size() - 1][0]);
    } else
      se.emplace_back(i);
    cerr << "." << flush;
  }
  return result;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> reDig<T>::rgb2l(const Mat rgb[3]) {
  // CIE 1931 XYZ from wikipedia.org
  auto X(.49000/.17697 * rgb[0] + .31000/.17697  * rgb[1] + .30000/.17697  * rgb[2]);
  auto Y(.17697/.17697 * rgb[0] + .81240/.17697  * rgb[1] + .010630/.17697 * rgb[2]);
  auto Z(                         .010000/.17697 * rgb[1] + .99000/.17697  * rgb[2]);
  return Y;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> reDig<T>::rgb2xz(const Mat rgb[3]) {
  // CIE 1931 XYZ from wikipedia.org
  auto X(.49000/.17697 * rgb[0] + .31000/.17697  * rgb[1] + .30000/.17697  * rgb[2]);
  auto Y(.17697/.17697 * rgb[0] + .81240/.17697  * rgb[1] + .010630/.17697 * rgb[2]);
  auto Z(                         .010000/.17697 * rgb[1] + .99000/.17697  * rgb[2]);
  Mat result(rgb[0].rows(), rgb[0].cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int j = 0; j < rgb[0].rows(); j ++)
    for(int k = 0; k < rgb[0].cols(); k ++)
      result(j, k) = sqrt(X(j, k) * X(j, k) + Z(j, k) * Z(j, k));
  return result;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> reDig<T>::tilt45(const Mat& in, const bool& invert, const Mat& orig) {
  Mat res;
  if(invert) {
    assert(orig.rows() + orig.cols() == in.rows() &&
           orig.rows() + orig.cols() == in.cols());
    res = Mat(orig.rows(), orig.cols());
    for(int i = 0; i < orig.rows(); i ++)
      for(int j = 0; j < orig.cols(); j ++)
        res(i, j) = in(i + j, orig.rows() - i + j);
  } else {
    res = Mat(in.rows() + in.cols(), in.cols() + in.rows());
    for(int i = - in.rows() - in.cols(); i < 2 * in.rows() + in.cols(); i ++)
      for(int j = - in.rows() - in.cols(); j < 2 * in.cols() + in.rows(); j ++) {
        const int y(i + j);
        const int x(in.rows() - i + j);
        if(0 <= y && y < res.rows() && 0 <= x && x < res.cols())
          res(y, x) = in(abs(abs((i + in.rows()) % (2 * in.rows())) - in.rows()) % in.rows(),
                         abs(abs((j + in.cols()) % (2 * in.cols())) - in.cols()) % in.cols());
      }
    for(int i = - in.rows() - in.cols(); i < 2 * in.rows() + in.cols(); i ++)
      for(int j = - in.rows() - in.cols(); j < 2 * in.cols() + in.rows(); j ++) {
        const int y(i + j + 1);
        const int x(in.rows() - i + j);
        if(!(0 <= y && y < res.rows() && 0 <= x && x < res.cols()))
          continue;
        int cnt(0);
        T   sum(0);
        const int y0(abs(abs((i + in.rows()) % (2 * in.rows())) - in.rows()) % in.rows());
        const int x0(abs(abs((j + in.cols()) % (2 * in.cols())) - in.cols()) % in.cols());
        if(0 <= y0 - 1) {
          sum += in(y0 - 1, x0);
          cnt ++;
        }
        if(0 <= x0 - 1) {
          sum += in(y0, x0 - 1);
          cnt ++;
        }
        if(y0 + 1 < in.rows()) {
          sum += in(y0 + 1, x0);
          cnt ++;
        }
        if(x0 + 1 < in.cols()) {
          sum += in(y0, x0 + 1);
          cnt ++;
        }
        res(y, x) = sum / cnt;
      }
  }
  return res;
}

template <typename T> void reDig<T>::normalize(Mat data[3], const T& upper) {
  T MM(data[0](0, 0)), mm(data[0](0, 0));
  for(int k = 0; k < 3; k ++)
    for(int i = 0; i < data[k].rows(); i ++)
      for(int j = 0; j < data[k].cols(); j ++) {
        MM = max(MM, data[k](i, j));
        mm = min(mm, data[k](i, j));
      }
  if(MM == mm)
    return;
  for(int k = 0; k < 3; k ++) {
    for(int i = 0; i < data[k].rows(); i ++)
      for(int j = 0; j < data[k].cols(); j ++)
        data[k](i, j) -= mm;
    data[k] *= upper / (MM - mm);
  }
  return;
}

#define _REDIG_
#endif

