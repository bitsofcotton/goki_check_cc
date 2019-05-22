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

template <typename T> class match_t;

using std::sqrt;
using std::atan2;
using std::abs;
using std::log;
using std::sin;
using std::cos;
using std::sin;
using std::min;
using std::max;
using std::ceil;
using std::floor;
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

template <typename T> class triangles_t {
public:
#if defined(_WITHOUT_EIGEN_)
  typedef SimpleMatrix<T> Mat3x3;
  typedef SimpleVector<T> Vec3;
#else
  typedef Eigen::Matrix<T, 3, 3> Mat3x3;
  typedef Eigen::Matrix<T, 3, 1> Vec3;
#endif
  Mat3x3 p;
  Vec3   n;
  T      c;
  T      z;
  triangles_t() {
    p = Mat3x3(3, 3);
    n = Vec3(3);
  }
  triangles_t<T>& rotate(const Mat3x3& R, const Vec3& origin) {
    for(int i = 0; i < 3; i ++)
#if defined(_WITHOUT_EIGEN_)
      p.setCol(i, R * (p.col(i) - origin) + origin);
#else
      p.col(i) = R * (p.col(i) - origin) + origin;
#endif
    return *this;
  }
  triangles_t<T>& solveN() {
    const auto pq(p.col(1) - p.col(0));
    const auto pr(p.col(2) - p.col(0));
    n[0] =   (pq[1] * pr[2] - pq[2] * pr[1]);
    n[1] = - (pq[0] * pr[2] - pq[2] * pr[0]);
    n[2] =   (pq[0] * pr[1] - pq[1] * pr[0]);
    if(n.dot(n) > 0)
      n /= sqrt(n.dot(n));
    z = n.dot(p.col(0));
    return *this;
  }
};

template <typename T> class reDig {
public:
#if defined(_WITHOUT_EIGEN_)
  typedef SimpleMatrix<T> Mat;
  typedef SimpleMatrix<complex<T> > MatU;
  typedef SimpleMatrix<T> Mat4x4;
  typedef SimpleMatrix<T> Mat3x3;
  typedef SimpleMatrix<T> Mat2x2;
  typedef SimpleVector<T> Vec;
  typedef SimpleVector<T> Vec3;
  typedef SimpleVector<T> Vec2;
  typedef SimpleVector<int> Veci3;
#else
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<complex<T>, Eigen::Dynamic, Eigen::Dynamic> MatU;
  typedef Eigen::Matrix<T, 4, 4>                           Mat4x4;
  typedef Eigen::Matrix<T, 3, 3>                           Mat3x3;
  typedef Eigen::Matrix<T, 2, 2>                           Mat2x2;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
  typedef Eigen::Matrix<T,   3, 1> Vec3;
  typedef Eigen::Matrix<int, 3, 1> Veci3;
  typedef Eigen::Matrix<T,   2, 1> Vec2;
#endif
  typedef triangles_t<T>           Triangles;
  
  reDig();
  ~reDig();
  void initialize(const int& vbox, const T& rz = - T(1));
  Mat  emphasis(const Mat& dstimg, const Mat& srcimg, const Mat& srcbump, const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const vector<Veci3>& hulldst, const vector<Veci3>& hullsrc, const T& ratio);
  Mat  replace(const Mat& dstimg, const vector<Vec3>& src, const vector<Veci3>& hullsrc, const bool& elim = false);
  Mat  replace(const Mat& dstimg, const Mat& srcimg, const Mat& srcbump, const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const vector<Veci3>& hulldst, const vector<Veci3>& hullsrc);
  vector<Vec3> takeShape(const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const vector<Veci3>& hulldst, const vector<Veci3>& hullsrc, const T& ratio);
  Mat  showMatch(const Mat& dstimg, const vector<Vec3>& dst, const vector<Veci3>& hull, const T& emph = T(.2));
  Mat  makeRefMatrix(const Mat& orig, const int& start) const;
  Mat  pullRefMatrix(const Mat& ref, const int& start, const Mat& orig) const;
  vector<Veci3> delaunay2(const vector<Vec3>& p, const vector<int>& pp, const T& epsilon = T(1e-5), const int& mdiv = 400) const;
  void maskVectors(vector<Vec3>& points, const vector<Veci3>& polys, const Mat& mask);
  void maskVectors(vector<Vec3>& points, vector<Veci3>& polys, const Mat& mask);
  Mat  reShape(const Mat& cbase, const Mat& vbase, const int& count = 20);
  vector<vector<int> > getEdges(const Mat& mask, const vector<Vec3>& points);
  Mat  rgb2l(const Mat rgb[3]);
  Mat  rgb2d(const Mat rgb[3]);
  Mat  rgb2xz(const Mat rgb[3]);
  void rgb2xyz(Mat xyz[3], const Mat rgb[3]);
  Mat  contrast(const Mat& in, const T& intensity, const T& thresh = T(.5));
  Mat  applytilt(const Mat& in, const int& dx, const int& dy);
  Mat  normalize(const Mat& data, const T& upper);
  void normalize(Mat data[3], const T& upper);
  Mat  autoLevel(const Mat& data, const int& count = 0);
  void getTileVec(const Mat& in, vector<Vec3>& geoms, vector<Veci3>& delaunay);
  match_t<T> tiltprep(const Mat& in, const int& idx, const int& samples, const T& psi);
  vector<Triangles> tiltprep(const vector<Vec3>& points, const vector<Veci3>& polys, const Mat& in, const match_t<T>& m);
  Mat  tilt(const Mat& in, const vector<Triangles>& triangles, const T& z0 = - T(1e8));
  Mat  tilt(const Mat& in, const Mat& bump, const match_t<T>& m, const T& z0 = - T(1e8));

private:
  void drawMatchLine(Mat& map, const Vec3& lref0, const Vec3& lref1, const T& emph);
  void drawMatchTriangle(Mat& map, const Vec3& lref0, const Vec3& lref1, const Vec3& lref2);
  bool isDelaunay2(T& cw, const Vec3 p[4], const T& epsilon) const;
  bool isCrossing(const Vec3& p0, const Vec3& p1, const Vec3& q0, const Vec3& q1, const T& err = T(1e-2)) const;
  void floodfill(Mat& checked, vector<pair<int, int> >& store, const Mat& mask, const int& y, const int& x);
  bool onTriangle(T& z, const Triangles& tri, const Vec2& geom);
  Triangles makeTriangle(const int& u, const int& v, const Mat& in, const Mat& bump, const int& flg);
  bool sameSide2(const Vec2& p0, const Vec2& p1, const Vec2& p, const Vec2& q, const bool& extend = true, const T& err = T(1e-5)) const;
  bool sameSide3(const Vec3& p0, const Vec3& p1, const Vec3& p, const Vec3& q, const bool& extend = true, const T& err = T(1e-5)) const;
  Mat  tilt(const Mat& in, const vector<Triangles>& triangles0, const match_t<T>& m, const T& z0 = - T(1e8));
  int  getImgPtRecursive(const T& h, const T& y) const;
  
  T   Pi;
  int vbox;
  T   rz;
};

template <typename T> reDig<T>::reDig() {
  initialize(3, .15);
}

template <typename T> reDig<T>::~reDig() {
  ;
}

template <typename T> void reDig<T>::initialize(const int& vbox, const T& rz) {
  assert(0 < vbox);
  Pi         = T(4) * atan2(T(1), T(1));
  this->vbox = vbox;
  if(T(0) < rz)
    this->rz = rz;
  else
    this->rz = T(.3);
  return;
}

template <typename T> typename reDig<T>::Mat reDig<T>::emphasis(const Mat& dstimg, const Mat& srcimg, const Mat& srcbump, const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const vector<Veci3>& hulldst, const vector<Veci3>& hullsrc, const T& ratio) {
  cerr << "m" << flush;
  assert(hulldst.size() == hullsrc.size());
  vector<Triangles> triangles;
  triangles.reserve((srcimg.rows() - 1) * (srcimg.cols() - 1) * 2);
  for(int i = 0; i < srcimg.rows() - 1; i ++)
    for(int j = 0; j < srcimg.cols() - 1; j ++) {
      triangles.push_back(makeTriangle(i, j, srcimg, srcbump, false));
      triangles[triangles.size() - 1].c = srcimg(i, j);
      triangles.push_back(makeTriangle(i, j, srcimg, srcbump, true));
      triangles[triangles.size() - 1].c = srcimg(i, j);
    }
  
  cerr << "e(" << hulldst.size() << ")" << endl;
  const auto rmatch(~ match);
  vector<vector<int> > emphs;
  emphs.resize(triangles.size() * 3);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < hulldst.size(); i ++) {
    const auto& p0(dst[hulldst[i][0]]);
    const auto& p1(dst[hulldst[i][1]]);
    const auto& p2(dst[hulldst[i][2]]);
    const auto  dp0(rmatch.transform(p0));
    const auto  dp1(rmatch.transform(p1));
    const auto  dp2(rmatch.transform(p2));
    for(int l = 0; l < triangles.size(); l ++)
      for(int ll = 0; ll < 3; ll ++) {
        const auto dq(triangles[l].p.col(ll));
        const auto q(match.transform(dq));
        if((sameSide3(p0, p1, p2, q) &&
            sameSide3(p1, p2, p0, q) &&
            sameSide3(p2, p0, p1, q)) ||
           (sameSide3(dp0, dp1, dp2, dq) &&
            sameSide3(dp1, dp2, dp0, dq) &&
            sameSide3(dp2, dp0, dp1, dq)) ) {
#if defined(_OPENMP)
#pragma omp critical
#endif
            {
              emphs[l * 3 + ll].push_back(hulldst[i][ll]);
            }
        }
      }
  }
  for(int i = 0; i < emphs.size(); i ++) if(emphs[i].size()) {
    Vec3 diff(3);
    diff[0] = diff[1] = diff[2] = T(0);
    for(int j = 0; j < emphs[i].size(); j ++)
      diff += rmatch.transform(dst[emphs[i][j]]) - triangles[i / 3].p.col(i % 3);
#if defined(_WITHOUT_EIGEN_)
    triangles[i / 3].p.setCol(i % 3, triangles[i / 3].p.col(i % 3) + diff * ratio / emphs[i].size());
#else
    triangles[i / 3].p.col(i % 3) += diff * ratio / emphs[i].size();
#endif
  }
  return tilt(dstimg, triangles, match);
}

template <typename T> typename reDig<T>::Mat reDig<T>::replace(const Mat& dstimg, const vector<Vec3>& src, const vector<Veci3>& hullsrc, const bool& elim) {
  Mat result(dstimg);
  T   M(0);
  T   m(0);
  for(int i = 0; i < src.size(); i ++) {
    if(i) {
      M = max(M, src[i][2]);
      m = min(m, src[i][2]);
    } else
      M = m = src[i][2];
  }
  auto tsrc(src);
  if(M - m != T(0))
    for(int i = 0; i < tsrc.size(); i ++)
      tsrc[i][2] = elim ? T(0) : (tsrc[i][2] - m) / (M - m);
  for(int ii = 0; ii < hullsrc.size(); ii ++)
    drawMatchTriangle(result, tsrc[hullsrc[ii][0]],
                              tsrc[hullsrc[ii][1]],
                              tsrc[hullsrc[ii][2]]);
  return result;
}

template <typename T> typename reDig<T>::Mat reDig<T>::replace(const Mat& dstimg, const Mat& srcimg, const Mat& srcbump, const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const vector<Veci3>& hulldst, const vector<Veci3>& hullsrc) {
  const Mat repl(replace(dstimg, src, match, hullsrc, true));
  return repl + emphasis(repl, srcimg, srcbump, dst, src, match,
                         hulldst, hullsrc, T(0));
}

template <typename T> vector<typename reDig<T>::Vec3> reDig<T>::takeShape(const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const vector<Veci3>& hulldst, const vector<Veci3>& hullsrc, const T& ratio) {
  assert(hulldst.size() == hullsrc.size());
  vector<Vec3> result(dst);
  const auto rmatch(~ match);
  vector<vector<int> > emphs;
  emphs.resize(result.size() * 3);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < hullsrc.size(); i ++) {
    const auto& p0(src[hullsrc[i][0]]);
    const auto& p1(src[hullsrc[i][1]]);
    const auto& p2(src[hullsrc[i][2]]);
    const auto  dp0(match.transform(p0));
    const auto  dp1(match.transform(p1));
    const auto  dp2(match.transform(p2));
    for(int j = 0; j < 3; j ++) {
      const auto& dq(dst[hulldst[i][j]]);
      const auto  q(rmatch.transform(dq));
      if((sameSide3(p0, p1, p2, q) &&
          sameSide3(p1, p2, p0, q) &&
          sameSide3(p2, p0, p1, q)) ||
         (sameSide3(dp0, dp1, dp2, dq) &&
          sameSide3(dp1, dp2, dp0, dq) &&
          sameSide3(dp2, dp0, dp1, dq)) ) {
#if defined(_OPENMP)
#pragma omp critical
#endif
          {
            emphs[hulldst[i][j]].push_back(hullsrc[i][j]);
          }
        }
    }
  }
  for(int i = 0; i < emphs.size(); i ++) if(emphs[i].size()) {
    Vec3 diff(3);
    diff[0] = diff[1] = diff[2] = T(0);
    for(int j = 0; j < emphs[i].size(); j ++)
      diff += match.transform(src[emphs[i][j]]) - dst[i];
    result[i] += diff * ratio / emphs[i].size();
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
  ldiff -= ldiff0 * ldiff.dot(ldiff0) / ldiff0.dot(ldiff0);
  const T    lnum(sqrt(ldiff.dot(ldiff)) + 1);
  // XXX : tan theta depend loop num, this have glitches.
  for(int k = 0; k < lnum * 4; k ++) {
    const Vec3 l0(lref0 + (lref2 - lref0) * k / (lnum * 4));
    const Vec3 l1(lref1 + (lref2 - lref1) * k / (lnum * 4));
    for(int i = 0; i <= int(abs(l0[idxM] - l1[idxM]) + 1); i ++) {
      const auto gidx(l0 + (l1 - l0) * i / int(abs(l0[idxM] - l1[idxM]) + 1));
      map(max(0, min(int(gidx[0]), int(map.rows() - 1))),
          max(0, min(int(gidx[1]), int(map.cols() - 1)))) = gidx[2];
    }
  }
  return;
}

template <typename T> typename reDig<T>::Mat reDig<T>::showMatch(const Mat& dstimg, const vector<Vec3>& dst, const vector<Veci3>& hull, const T& emph) {
  Mat map(dstimg.rows(), dstimg.cols());
  for(int i = 0; i < map.rows(); i ++)
    for(int j = 0; j < map.cols(); j ++)
      map(i, j) = T(0);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int k = 0; k < hull.size(); k ++) {
    drawMatchLine(map, dst[hull[k][0]], dst[hull[k][1]], emph);
    drawMatchLine(map, dst[hull[k][1]], dst[hull[k][2]], emph);
    drawMatchLine(map, dst[hull[k][2]], dst[hull[k][0]], emph);
  }
  return dstimg + map;
}

template <typename T> typename reDig<T>::Mat reDig<T>::makeRefMatrix(const Mat& orig, const int& start) const {
  Mat result(orig.rows(), orig.cols());
  for(int i = 0; i < orig.rows() * orig.cols(); i ++)
    result(i % orig.rows(), i / orig.rows()) = i + start;
  return result;
}

template <typename T> typename reDig<T>::Mat reDig<T>::pullRefMatrix(const Mat& ref, const int& start, const Mat& orig) const {
  assert(orig.rows() == ref.rows() && orig.cols() == ref.cols());
  Mat result(ref.rows(), ref.cols());
  for(int i = 0; i < ref.rows() * ref.cols(); i ++) {
    const int ly(i % ref.rows());
    const int lx(i / ref.rows());
    const int v(ref(ly, lx) - start);
    if(0 <= v && v < orig.rows() * orig.cols())
      result(ly, lx) = orig(v % orig.rows(), v / orig.rows());
    else
      result(ly, lx) = T(0);
  }
  for(int i = 0; i < result.rows(); i ++) {
    T c(0);
    for(int j = 0; j < result.cols(); j ++)
      if(result(i, j) != T(0)) {
        c = result(i, j);
        break;
      }
    for(int j = 0; j < result.cols(); j ++)
      if(result(i, j) == T(0))
        result(i, j) = c;
      else
        break;
    c = T(0);
    for(int j = 0; j < result.cols(); j ++)
      if(result(i, result.cols() - 1 - j) != T(0)) {
        c = result(i, result.cols() - 1 - j);
        break;
      }
    for(int j = 0; j < result.cols(); j ++)
      if(result(i, result.cols() - 1 - j) == T(0))
        result(i, result.cols() - 1 - j) = c;
      else
        break;
  }
  for(int i = 0; i < result.cols(); i ++) {
    T c(0);
    for(int j = 0; j < result.rows(); j ++)
      if(result(j, i) != T(0)) {
        c = result(j, i);
        break;
      }
    for(int j = 0; j < result.rows(); j ++)
      if(result(j, i) == T(0)) {
        result(j, i) = c;
      } else
        break;
    c = T(0);
    for(int j = 0; j < result.rows(); j ++)
      if(result(result.rows() - 1 - j, i) != T(0)) {
        c = result(result.rows() - 1 - j, i);
        break;
      }
    for(int j = 0; j < result.rows(); j ++)
      if(result(result.rows() - 1 - j, i) == T(0))
        result(result.rows() - 1 - j, i) = c;
      else
        break;
  }
  return result;
}

// XXX this have glitches.
template <typename T> vector<typename reDig<T>::Veci3> reDig<T>::delaunay2(const vector<Vec3>& p, const vector<int>& pp, const T& epsilon, const int& mdiv) const {
  vector<Veci3> res;
  cerr << pp.size() << ":" << flush;
  if(pp.size() > mdiv) {
    vector<pair<Vec3, int> > div;
    div.reserve(pp.size());
    for(int i = 0; i < pp.size(); i ++)
      div.push_back(make_pair(p[pp[i]], pp[i]));
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
        const auto itr(upper_bound(div.begin(), div.end(),
                         make_pair(p[left[i][ii]], left[i][ii]),
                         less0<pair<Vec3, int> >));
        if(div.size() * 7 / 16 < distance(div.begin(), itr))
          goto nextl;
      }
      work.emplace_back(left[i]);
     nextl:
      ;
    }
    for(int i = 0; i < right.size(); i ++) {
      for(int ii = 0; ii < 3; ii ++) {
        const auto itr(upper_bound(div.begin(), div.end(),
                         make_pair(p[right[i][ii]], right[i][ii]),
                         less0<pair<Vec3, int> >));
        if(distance(div.begin(), itr) < div.size() * 9 / 16)
          goto nextr;
      }
      work.emplace_back(right[i]);
     nextr:
      ;
    }
    for(int i = 0; i < middle.size(); i ++) {
      for(int ii = 0; ii < 3; ii ++) {
        const auto itr(upper_bound(div.begin(), div.end(),
                         make_pair(p[middle[i][ii]], middle[i][ii]),
                         less0<pair<Vec3, int> >));
        if(distance(div.begin(), itr) < div.size() * 5 / 16 ||
                 div.size() * 11 / 16 < distance(div.begin(), itr) )
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
      Veci3 idx;
      Vec3 q[4];
      assert(middle.size());
      for(int jj = 0; jj < middle.size(); jj ++)
        for(int i0 = 0; i0 < 3; i0 ++)
          for(int j0 = 0; j0 < 3; j0 ++)
            if(isCrossing(p[work[ii][   i0      % 3]],
                          p[work[ii][  (i0 + 1) % 3]],
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
#pragma omp critical
#endif
      {
        res.emplace_back(idx);
      }
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
          Veci3 idx;
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
            for(int ii = 0; ii < res.size(); ii ++)
              for(int jj = 0; jj < res[ii].size(); jj ++)
                for(int kk = 0; kk < idx.size(); kk ++)
                  if(isCrossing(p[res[ii][(jj + 0) % 3]],
                                p[res[ii][(jj + 1) % 3]],
                                p[idx[(kk + 0) % 3]],
                                p[idx[(kk + 1) % 3]]))
                    goto fixnextcr;
            res.emplace_back(idx);
           fixnextcr:
            ;
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
  Vec3 err(p[0] - p[2]);
  err[2] = T(0);
  err   -= bcn * err.dot(bcn) / bcn.dot(bcn);
  if(err.dot(err) <= epsilon)
    return false;
  Mat4x4 dc(4, 4);
  const auto g((p[0] + p[1] + p[2]) / T(3));
  for(int i = 0; i < 4; i ++) {
    dc(i, 0) = T(1);
    dc(i, 1) = p[i][0] - g[0];
    dc(i, 2) = p[i][1] - g[1];
    dc(i, 3) = dc(i, 1) * dc(i, 1) + dc(i, 2) * dc(i, 2);
  }
  Mat3x3 dc0(3, 3);
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
  if(cw * dc.determinant() < - epsilon ||
     (sameSide3(p[0], p[1], p[2], p[3], false) &&
      sameSide3(p[1], p[2], p[0], p[3], false) &&
      sameSide3(p[2], p[0], p[1], p[3], false)) )
    return false;
  return true;
}

template <typename T> bool reDig<T>::isCrossing(const Vec3& p0, const Vec3& p1, const Vec3& q0, const Vec3& q1, const T& err) const {
  // t * p0 + (1 - t) * p1 == s * q0 + (1 - s) * q1
  // <=> p1 + (p0 - p1) t == q1 + (q0 - q1) s
  // <=> [(p0 - p1), (q1 - q0)][t, s] == q1 - p1.
  // <=> Ax==b.
  Mat2x2 A(2, 2);
  Vec2   b(2);
  A(0, 0) = p0[0] - p1[0];
  A(1, 0) = p0[1] - p1[1];
  A(0, 1) = q1[0] - q0[0];
  A(1, 1) = q1[1] - q0[1];
  b[0]    = q1[0] - p1[0];
  b[1]    = q1[1] - p1[1];
  if(abs(A.determinant()) <= pow(err, T(2)))
    return false;
#if defined(_WITHOUT_EIGEN_)
  auto x(A.solve(b));
#else
  auto x(A.inverse() * b);
#endif
  return (err <= x[0] && x[0] <= T(1) - err) &&
         (err <= x[1] && x[1] <= T(1) - err);
}

template <typename T> void reDig<T>::maskVectors(vector<Vec3>& points, const vector<Veci3>& polys, const Mat& mask) {
  vector<Veci3> tpoly(polys);
  return maskVectors(points, tpoly, mask);
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
  for(int i = 0; i < elimp.size(); i ++)
    polys.erase(polys.begin() + (elimp[i] - i));
  for(int i = 0; i < polys.size(); i ++)
    for(int j = 0; j < polys[i].size(); j ++)
      polys[i][j] = after[polys[i][j]];
  return;
}

template <typename T> typename reDig<T>::Mat reDig<T>::reShape(const Mat& cbase, const Mat& vbase, const int& count) {
  assert(cbase.rows() && cbase.cols() && vbase.rows() && cbase.cols());
  assert(vbase.rows() % cbase.rows() == 0 && vbase.cols() % cbase.cols() == 0);
  assert(vbase.rows() / cbase.rows() == vbase.cols() / cbase.cols());
  Mat res(vbase.rows(), vbase.cols());
  vector<pair<T, pair<int, int> > > vpoints;
  vpoints.reserve(vbase.rows() * vbase.cols());
  for(int i = 0; i < vbase.rows(); i ++)
    for(int j = 0; j < vbase.cols(); j ++)
      vpoints.push_back(make_pair(vbase(i, j), make_pair(i, j)));
  sort(vpoints.begin(), vpoints.end());
  vector<T> vvpoints;
  vvpoints.reserve(vpoints.size());
  vvpoints.push_back(T(0));
  for(int i = 1; i < vpoints.size(); i ++)
    vvpoints.push_back(vvpoints[i - 1] + vpoints[i].first - vpoints[i - 1].first);
  T avg(0);
  const int ratio(vbase.rows() / cbase.rows());
  for(int i = 0, ii = 0; i < vpoints.size(); i ++) {
    if(abs(vvpoints[i] - vvpoints[ii]) < abs(vvpoints[vvpoints.size() - 1]) / count && i < vpoints.size() - 1)
      avg += cbase(vpoints[i].second.first / ratio, vpoints[i].second.second / ratio);
    else {
      if(i == vpoints.size() - 1) {
        avg += cbase(vpoints[i].second.first / ratio, vpoints[i].second.second / ratio);
        i ++;
      }
      avg /= (i - ii + 1);
      for(int j = ii; j < i; j ++)
        res(vpoints[j].second.first, vpoints[j].second.second) = avg;
      avg  = T(0);
      ii   = i;
    }
  }
  return res;
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

template <typename T> vector<vector<int> > reDig<T>::getEdges(const Mat& mask, const vector<Vec3>& points) {
  cerr << "getEdges" << flush;
  vector<vector<int> > result;
  if(mask.rows() <= 0 || mask.cols() <= 0)
    return result;
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
      // XXX: sorting is needed.
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

template <typename T> typename reDig<T>::Mat reDig<T>::rgb2l(const Mat rgb[3]) {
  Mat work[3];
  rgb2xyz(work, rgb);
  return work[0];
}

template <typename T> typename reDig<T>::Mat reDig<T>::rgb2d(const Mat rgb[3]) {
  Mat xyz[3];
  rgb2xyz(xyz, rgb);
  Mat result(rgb[0].rows(), rgb[0].cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int j = 0; j < rgb[0].rows(); j ++)
    for(int k = 0; k < rgb[0].cols(); k ++)
      result(j, k) = sqrt(xyz[0](j, k) * xyz[0](j, k) + xyz[1](j, k) * xyz[1](j, k) + xyz[2](j, k) * xyz[2](j, k));
  return result;
}

template <typename T> typename reDig<T>::Mat reDig<T>::rgb2xz(const Mat rgb[3]) {
  Mat xyz[3];
  rgb2xyz(xyz, rgb);
  Mat result(rgb[0].rows(), rgb[0].cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int j = 0; j < rgb[0].rows(); j ++)
    for(int k = 0; k < rgb[0].cols(); k ++)
      result(j, k) = sqrt(xyz[1](j, k) * xyz[1](j, k) + xyz[2](j, k) * xyz[2](j, k));
  return result;
}

template <typename T> void reDig<T>::rgb2xyz(Mat xyz[3], const Mat rgb[3]) {
  // CIE 1931 XYZ from wikipedia.org
  xyz[0] = rgb[0] * .49000/.17697 + rgb[1] * .31000/.17697  + rgb[2] * .30000/.17697;
  xyz[1] = rgb[0] * .17697/.17697 + rgb[1] * .81240/.17697  + rgb[2] * .010630/.17697;
  xyz[2] =                          rgb[1] * .010000/.17697 + rgb[2] * .99000/.17697;
  return;
}

template <typename T> typename reDig<T>::Mat reDig<T>::contrast(const Mat& in, const T& intensity, const T& thresh) {
  Mat result(in);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < result.rows(); i ++)
    for(int j = 0; j < result.cols(); j ++)
      result(i, j) = min(abs(thresh) + T(.5), max(- abs(thresh) + T(.5), intensity * (result(i, j) - T(.5)) + T(.5)));
  return result;
}

template <typename T> typename reDig<T>::Mat reDig<T>::applytilt(const Mat& in, const int& dx, const int& dy) {
  assert(dx == 1 || dy == 1);
  if(abs(dx) < abs(dy))
    return applytilt(in.transpose(), dy, dx).transpose();
  Mat res(in.rows(), in.cols());
  for(int i = 0; i < in.rows(); i ++)
    for(int j = 0; j < in.cols(); j ++)
      res(i, j) = in(getImgPtRecursive(in.rows(), i + j / dx), j);
  return res;
}

template <typename T> typename reDig<T>::Mat reDig<T>::normalize(const Mat& data, const T& upper) {
  T MM(data(0, 0)), mm(data(0, 0));
  for(int i = 0; i < data.rows(); i ++)
    for(int j = 0; j < data.cols(); j ++)
      if(isfinite(data(i, j))) {
        MM = max(MM, data(i, j));
        mm = min(mm, data(i, j));
      }
  if(MM == mm)
    return data;
  Mat result(data);
  for(int i = 0; i < data.rows(); i ++)
    for(int j = 0; j < data.cols(); j ++) {
      if(isfinite(result(i, j)))
        result(i, j) -= mm;
      else
        result(i, j)  = T(0);
      assert(T(0) <= result(i, j) && result(i, j) <= MM - mm);
    }
  return result * upper / (MM - mm);
}

template <typename T> void reDig<T>::normalize(Mat data[3], const T& upper) {
  T MM(data[0](0, 0)), mm(data[0](0, 0));
  for(int k = 0; k < 3; k ++)
    for(int i = 0; i < data[k].rows(); i ++)
      for(int j = 0; j < data[k].cols(); j ++)
        if(isfinite(data[k](i, j))) {
          MM = max(MM, data[k](i, j));
          mm = min(mm, data[k](i, j));
        }
  if(MM == mm)
    return;
  for(int k = 0; k < 3; k ++) {
    for(int i = 0; i < data[k].rows(); i ++)
      for(int j = 0; j < data[k].cols(); j ++) {
        if(isfinite(data[k](i, j)))
          data[k](i, j) -= mm;
        else
          data[k](i, j)  = T(0);
        // assert(T(0) <= data[k](i, j) && data[k](i, j) <= MM - mm);
      }
    data[k] *= upper / (MM - mm);
  }
  return;
}

template <typename T> typename reDig<T>::Mat reDig<T>::autoLevel(const Mat& data, const int& count) {
  vector<T> res;
  res.reserve(data.rows() * data.cols());
  for(int i = 0; i < data.rows(); i ++)
    for(int j = 0; j < data.cols(); j ++)
      res.push_back(data(i, j));
  sort(res.begin(), res.end());
  Mat result(data.rows(), data.cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < data.rows(); i ++)
    for(int j = 0; j < data.cols(); j ++)
      result(i, j) = max(min(data(i, j), res[res.size() - count - 1]), res[count]);
  return result;
}

// get bump with multiple scale and vectorized result.
template <typename T> void reDig<T>::getTileVec(const Mat& in, vector<Vec3>& geoms, vector<Veci3>& delaunay) {
  // get vectorize.
  geoms = vector<Vec3>();
  T aavg(0);
  for(int i = 0; i < in.rows(); i ++)
    for(int j = 0; j < in.cols(); j ++)
      aavg += in(i, j);
  aavg /= in.rows() * in.cols();
  for(int i = 0; i < in.rows() / vbox + 1; i ++)
    for(int j = 0; j < in.cols() / vbox + 1; j ++) {
      if(in.rows() < (i + 1) * vbox ||
         in.cols() < (j + 1) * vbox) {
        if(geoms.size() >= 1) {
          Vec3 gbuf(3);
          gbuf[0] = i * vbox;
          gbuf[1] = j * vbox;
          gbuf[2] = geoms[geoms.size() - 1][2];
          geoms.push_back(gbuf);
        }
        continue;
      }
      T avg(0);
      for(int ii = i * vbox; ii < (i + 1) * vbox; ii ++)
        for(int jj = j * vbox; jj < (j + 1) * vbox; jj ++)
          avg += in(ii, jj);
      Vec3 work(3);
      work[0] = i * vbox;
      work[1] = j * vbox;
      work[2] = sqrt(T(in.rows() * in.cols())) * rz * (avg / vbox / vbox - aavg);
      geoms.push_back(work);
    }
  Vec3 avg(3);
  avg[0] = avg[1] = avg[2] = T(0);
  for(int i = 0; i < geoms.size(); i ++)
    avg += geoms[i];
  avg /= geoms.size();
  for(int i = 0; i < geoms.size(); i ++)
    geoms[i][2] -= avg[2];
  delaunay = vector<Veci3>();
  for(int i = 1; i < in.rows() / vbox + 1; i ++)
    for(int j = 0; j < in.cols() / vbox; j ++) {
      Veci3 work(3), work2(3);
      work[0]  = (i - 1) * (in.cols() / vbox + 1) + j;
      work[1]  =  i      * (in.cols() / vbox + 1) + j;
      work[2]  =  i      * (in.cols() / vbox + 1) + j + 1;
      work2[0] = (i - 1) * (in.cols() / vbox + 1) + j;
      work2[2] = (i - 1) * (in.cols() / vbox + 1) + j + 1;
      work2[1] =  i      * (in.cols() / vbox + 1) + j + 1;
      delaunay.push_back(work);
      delaunay.push_back(work2);
    }
  return;
}

template <typename T> triangles_t<T> reDig<T>::makeTriangle(const int& u, const int& v, const Mat& in, const Mat& bump, const int& flg) {
  Triangles work;
  if(flg) {
    work.p(0, 0) = u;
    work.p(1, 0) = v;
    work.p(0, 1) = u + 1;
    work.p(1, 1) = v;
    work.p(0, 2) = u + 1;
    work.p(1, 2) = v + 1;
  } else {
    work.p(0, 0) = u;
    work.p(1, 0) = v;
    work.p(0, 1) = u;
    work.p(1, 1) = v + 1;
    work.p(0, 2) = u + 1;
    work.p(1, 2) = v + 1;
  }
  work.c = in(u, v);
  for(int i = 0; i < 3;  i ++)
    work.p(2, i) = bump(int(work.p(0, i)), int(work.p(1, i))) * sqrt(T(bump.rows() * bump.cols())) * rz;
  return work.solveN();
}

template <typename T> bool reDig<T>::sameSide2(const Vec2& p0, const Vec2& p1, const Vec2& p2, const Vec2& q, const bool& extend, const T& err) const {
  const Vec2 dlt(p1 - p0);
  Vec2 dp(p2 - p0);
  dp -= dlt * dlt.dot(dp) / dlt.dot(dlt);
  // N.B. dp.dot(p1 - p0) >> 0.
  return dp.dot(q - p0) >= (extend ? - T(1) : T(1)) * (abs(dp[0]) + abs(dp[1])) * err;
}

template <typename T> bool reDig<T>::sameSide3(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& q, const bool& extend, const T& err) const {
  Vec2 q0(2), q1(2), q2(2), qq(2);
  q0[0] = p0[0]; q0[1] = p0[1]; q1[0] = p1[0]; q1[1] = p1[1];
  q2[0] = p2[0]; q2[1] = p2[1]; qq[0] = q[0];  qq[1] = q[1];
  return sameSide2(q0, q1, q2, qq, extend, err);
}

// <[x, y, t], triangle.n> == triangle.z
template <typename T> bool reDig<T>::onTriangle(T& z, const Triangles& tri, const Vec2& geom) {
  Vec3 v0(3);
  Vec3 camera(3);
  v0[0] = 0;
  v0[1] = 0;
  v0[2] = 1;
  camera[0] = geom[0];
  camera[1] = geom[1];
  camera[2] = 0;
  // <v0 t + camera, tri.n> = tri.z
  const T t((tri.z - tri.n.dot(camera)) / (tri.n.dot(v0)));
  z = camera[2] + v0[2] * t;
  return (sameSide3(tri.p.col(0), tri.p.col(1), tri.p.col(2), camera, true, T(.125)) &&
          sameSide3(tri.p.col(1), tri.p.col(2), tri.p.col(0), camera, true, T(.125)) &&
          sameSide3(tri.p.col(2), tri.p.col(0), tri.p.col(1), camera, true, T(.125)));
}

template <typename T> match_t<T> reDig<T>::tiltprep(const Mat& in, const int& idx, const int& samples, const T& psi) {
  const T theta(2. * Pi * idx / samples);
  const T lpsi(Pi * psi);
  Mat3x3 R0(3, 3);
  Mat3x3 R1(3, 3);
  R0(0, 0) =   cos(theta);
  R0(0, 1) =   sin(theta);
  R0(0, 2) = 0.;
  R0(1, 0) = - sin(theta);
  R0(1, 1) =   cos(theta);
  R0(1, 2) = 0.;
  R0(2, 0) = 0.;
  R0(2, 1) = 0.;
  R0(2, 2) = 1.;
  R1(0, 0) = 1.;
  R1(0, 1) = 0.;
  R1(0, 2) = 0.;
  R1(1, 0) = 0.;
  R1(1, 1) = 1.;
  R1(1, 2) = 0.;
  R1(2, 0) = 0.;
  R1(2, 1) = 0.;
  R1(2, 2) = 1.;
  R1(0, 0) =   cos(lpsi);
  R1(0, 2) = - sin(lpsi);
  R1(2, 0) =   sin(lpsi);
  R1(2, 2) =   cos(lpsi);
  Vec3 pcenter(3);
  pcenter[0] = T(in.rows() - 1.) / 2;
  pcenter[1] = T(in.cols() - 1.) / 2;
  pcenter[2] = T(0);
  match_t<T> m;
  m.rot    = R0.transpose() * R1 * R0;
  // x -> m.rot * x, same center
  // x - pcenter -> m.rot * (x - pcenter)
  // x -> m.rot * x - m.rot * pcenter + pcenter.
  m.offset = pcenter - m.rot * pcenter;
  m.ratio  = T(1);
  return m;
}

template <typename T> vector<typename reDig<T>::Triangles> reDig<T>::tiltprep(const vector<Vec3>& points, const vector<Veci3>& polys, const Mat& in, const match_t<T>& m) {
  vector<Triangles> result;
  for(int i = 0; i < polys.size(); i ++) {
    Triangles work;
    for(int j = 0; j < 3; j ++) {
#if defined(_WITHOUT_EIGEN_)
      work.p.setCol(j, m.transform(points[polys[i][j]]));
#else
      work.p.col(j) = m.transform(points[polys[i][j]]);
#endif
    }
    if(T(0) <= points[polys[i][0]][0] && points[polys[i][0]][0] < in.rows() &&
       T(0) <= points[polys[i][0]][1] && points[polys[i][0]][1] < in.cols())
      work.c = in(int(points[polys[i][0]][0]),
                  int(points[polys[i][0]][1]));
    else
      work.c = T(0);
    result.push_back(work.solveN());
  }
  return result;
}

template <typename T> typename reDig<T>::Mat reDig<T>::tilt(const Mat& in, const Mat& bump, const match_t<T>& m, const T& z0) {
  assert(in.rows() == bump.rows() && in.cols() == bump.cols());
  vector<Triangles> triangles;
  triangles.reserve((in.rows() - 1) * (in.cols() - 1) * 2);
  for(int i = 0; i < in.rows() - 1; i ++)
    for(int j = 0; j < in.cols() - 1; j ++) {
      triangles.push_back(makeTriangle(i, j, in, bump, false));
      triangles.push_back(makeTriangle(i, j, in, bump, true));
    }
  return tilt(in, triangles, m, z0);
}

template <typename T> typename reDig<T>::Mat reDig<T>::tilt(const Mat& in, const vector<Triangles>& triangles0, const match_t<T>& m, const T& z0) {
  vector<Triangles> triangles(triangles0);
  for(int j = 0; j < triangles.size(); j ++) {
    for(int k = 0; k < 3; k ++)
#if defined(_WITHOUT_EIGEN_)
      triangles[j].p.setCol(k, m.transform(triangles[j].p.col(k)));
#else
      triangles[j].p.col(k) = m.transform(triangles[j].p.col(k));
#endif
    triangles[j].solveN();
  }
  return tilt(in, triangles, z0);
}

template <typename T> typename reDig<T>::Mat reDig<T>::tilt(const Mat& in, const vector<Triangles>& triangles, const T& z0) {
  Mat result(in.rows(), in.cols());
  for(int i = 0; i < in.rows(); i ++)
    for(int j = 0; j < in.cols(); j ++)
      result(i, j) = 0.;
  cerr << "t" << flush;
  Mat zb(in.rows(), in.cols());
  for(int j = 0; j < zb.rows(); j ++)
    for(int k = 0; k < zb.cols(); k ++)
      zb(j, k) = z0;
  // able to boost with divide and conquer.
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int j = 0; j < triangles.size(); j ++) {
    const Triangles& tri(triangles[j]);
    int ll = int( min(min(tri.p(0, 0), tri.p(0, 1)), tri.p(0, 2)));
    int rr = ceil(max(max(tri.p(0, 0), tri.p(0, 1)), tri.p(0, 2))) + 1;
    int bb = int( min(min(tri.p(1, 0), tri.p(1, 1)), tri.p(1, 2)));
    int tt = ceil(max(max(tri.p(1, 0), tri.p(1, 1)), tri.p(1, 2))) + 1;
    for(int y = max(0, ll); y < min(rr, int(in.rows())); y ++)
      for(int x = max(0, bb); x < min(tt, int(in.cols())); x ++) {
        T z;
        Vec2 midgeom(2);
        midgeom[0] = y;
        midgeom[1] = x;
        if(onTriangle(z, tri, midgeom) && isfinite(z)) {
#if defined(_OPENMP)
#pragma omp critical
#endif
          {
            if(zb(y, x) < z) {
              result(y, x) = tri.c;
              zb(y, x)     = z;
            }
          }
        }
      }
  }
  return result;
}

template <typename T> int reDig<T>::getImgPtRecursive(const T& h, const T& y) const {
  const int lgyh(log(y) / log(h));
  return int(y - (lgyh ? pow(h, lgyh) : 0) + h * h) % int(h);
}

#define _REDIG_
#endif

