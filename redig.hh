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
template <typename T> class Filter;

using std::abs;
using std::min;
using std::max;
using std::isinf;
using std::cerr;
using std::endl;
using std::flush;
using std::vector;
using std::lower_bound;
using std::upper_bound;
using std::distance;
using std::sort;
using std::unique;

template <typename T> static inline bool less0(const T& x, const T& y) {
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
  inline triangles_t() {
    p = Mat3x3(3, 3);
    n = Vec3(3);
  }
  inline triangles_t<T>& rotate(const Mat3x3& R, const Vec3& origin) {
    for(int i = 0; i < 3; i ++)
#if defined(_WITHOUT_EIGEN_)
      p.setCol(i, R * (p.col(i) - origin) + origin);
#else
      p.col(i) = R * (p.col(i) - origin) + origin;
#endif
    return *this;
  }
  inline triangles_t<T>& solveN() {
    const auto pq(p.col(1) - p.col(0));
    const auto pr(p.col(2) - p.col(0));
    n[0] =   (pq[1] * pr[2] - pq[2] * pr[1]);
    n[1] = - (pq[0] * pr[2] - pq[2] * pr[0]);
    n[2] =   (pq[0] * pr[1] - pq[1] * pr[0]);
    if(n.dot(n) > T(0))
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
  
  inline reDig();
  inline ~reDig();
  inline void initialize(const int& vbox, const T& rz = - T(1));
  Mat  draw(const Mat& img, const vector<Vec3>& shape, const vector<Vec3>& emph, const vector<Veci3>& hull);
  Mat  draw(const Mat& img, const vector<Vec3>& shape, const vector<Veci3>& hull, const bool& elim = false);
  Mat  drawBone(const vector<Vec3>& center, const vector<T>& r, const int& rows, const int& cols);
  vector<Vec3> takeShape(const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const T& ratio);
  vector<Vec3> takeShape(const vector<Vec3>& shape, const vector<Vec3>& center, const vector<Vec3>& outcenter, const vector<vector<int> >& attend, const T& ratio);
  void complement(vector<Mat>& dstimg, vector<Vec3>& dstcenter, const vector<vector<Mat> >& srcimg, const vector<vector<Vec3> >& srccenter, const vector<vector<int> >& attend, const vector<vector<pair<int, int> > >& a2xy, const T& iemph);
  Mat  showMatch(const Mat& dstimg, const vector<Vec3>& dst, const vector<Veci3>& hull, const T& emph = T(1));
  Mat  makeRefMatrix(const Mat& orig, const int& start) const;
  Mat  pullRefMatrix(const Mat& ref, const int& start, const Mat& orig) const;
  vector<Veci3> mesh2(const vector<Vec3>& p, const vector<int>& pp) const;
  void maskVectors(vector<Vec3>& points, const vector<Veci3>& polys, const Mat& mask);
  void maskVectors(vector<Vec3>& points, vector<Veci3>& polys, const Mat& mask);
  Mat  reShape(const Mat& cbase, const Mat& vbase, const int& count = 20);
  Mat  reColor(const Mat& cbase, const Mat& vbase, const int& count = 20);
  Mat  reColor3(const Mat& cbase, const Mat& vbase, const int& count = 20);
  Mat  reColor(const Mat& cbase, const int& count = 20, const T& intensity = T(1));
  Mat  reTrace(const Mat& dst, const Mat& src, const T& intensity, const int& count = 20);
  Mat  reTrace(const Mat& dst, const T& intensity, const int& count = 20);
  Mat  reImage(const Mat& dst, const Mat& src, const T& intensity, const int& count = 20);
  Mat  reImage(const Mat& dst, const T& intensity, const int& count = 20);
  Mat  catImage(const vector<Mat>& imgs);
  vector<vector<int> > getEdges(const Mat& mask, const vector<Vec3>& points);
  Mat  rgb2d(const Mat rgb[3]);
  void rgb2xyz(Mat xyz[3], const Mat rgb[3]);
  void xyz2rgb(Mat rgb[3], const Mat xyz[3]);
  Mat  contrast(const Mat& in, const T& intensity, const T& thresh = T(1) / T(2));
  Mat  normalize(const Mat& data, const T& upper);
  void normalize(Mat data[3], const T& upper);
  Mat  autoLevel(const Mat& data, const int& count = 0);
  void autoLevel(Mat data[3], const int& count = 0);
  void getTileVec(const Mat& in, vector<Vec3>& geoms, vector<Veci3>& delaunay);
  void getBone(const Mat& in, const vector<Vec3>& geoms, vector<Vec3>& center, vector<T>& r, vector<vector<int> >& attend, const T& thresh = T(1) / T(20));
  vector<vector<std::pair<int, int> > > getReverseLookup(const vector<vector<int> >& attend, const Mat& refimg);
  vector<Vec3> copyBone(const vector<Vec3>& centerdst, const vector<T>& rdst, const vector<Vec3>& centersrc, const vector<T>& rsrc);
  match_t<T> tiltprep(const Mat& in, const int& idx, const int& samples, const T& psi);
  vector<Triangles> tiltprep(const vector<Vec3>& points, const vector<Veci3>& polys, const Mat& in, const match_t<T>& m);
  Mat  tilt(const Mat& in, const vector<Triangles>& triangles, const T& z0 = - T(1000000));
  Mat  tilt(const Mat& in, const Mat& bump, const match_t<T>& m, const T& z0 = - T(1000000));
  void floodfill(Mat& checked, vector<pair<int, int> >& store, const Mat& mask, const int& y, const int& x);

private:
  void drawMatchLine(Mat& map, const Vec3& lref0, const Vec3& lref1, const T& emph);
  void drawMatchTriangle(Mat& map, const Vec3& lref0, const Vec3& lref1, const Vec3& lref2);
  inline bool isClockwise(const Vec3 p[3]) const;
  inline bool onTriangle(T& z, const Triangles& tri, const Vec2& geom);
  inline Triangles makeTriangle(const int& u, const int& v, const Mat& in, const Mat& bump, const int& flg);
  inline bool sameSide2(const Vec2& p0, const Vec2& p1, const Vec2& p, const Vec2& q, const bool& extend = true, const T& err = T(1) / T(100000)) const;
  inline bool sameSide3(const Vec3& p0, const Vec3& p1, const Vec3& p, const Vec3& q, const bool& extend = true, const T& err = T(1) / T(100000)) const;
  Mat  tilt(const Mat& in, const vector<Triangles>& triangles0, const match_t<T>& m, const T& z0 = - T(100000));
  inline int  getImgPtRecursive(const int& h, const int& y) const;
  void prepTrace(pair<Vec, Vec>& v, pair<pair<int, int>, pair<int, int> >& hw, const Mat& mask);
  Mat  applyTrace(const pair<Vec, Vec>& v, const pair<pair<pair<int, int>, pair<int, int> >, pair<pair<int, int>, pair<int, int> > >& hw);
  
  T   Pi;
  int vbox;
  T   rz;
};

template <typename T> inline reDig<T>::reDig() {
  initialize(3, T(15) / T(100));
}

template <typename T> inline reDig<T>::~reDig() {
  ;
}

template <typename T> inline void reDig<T>::initialize(const int& vbox, const T& rz) {
  assert(0 < vbox);
  Pi         = T(4) * atan2(T(1), T(1));
  this->vbox = vbox;
  if(T(0) < rz)
    this->rz = rz;
  else
    this->rz = T(03) / T(10);
  return;
}

template <typename T> typename reDig<T>::Mat reDig<T>::draw(const Mat& img, const vector<Vec3>& shape, const vector<Vec3>& emph, const vector<Veci3>& hull) {
  assert(shape.size() == emph.size());
  vector<Triangles> tris;
  tris.reserve(hull.size());
  for(int i = 0; i < hull.size(); i ++) {
    assert(hull[i].size() == 3);
    assert(0 <= hull[i][0] && hull[i][0] < shape.size());
    assert(0 <= hull[i][1] && hull[i][1] < shape.size());
    assert(0 <= hull[i][2] && hull[i][2] < shape.size());
    Triangles work;
    for(int j = 0; j < 3; j ++) {
#if defined(_WITHOUT_EIGEN_)
      work.p.setCol(j, emph[hull[i][j]]);
#else
      work.p.col(j) = emph[hull[i][j]];
#endif
    }
    work.c = img(max(0, min(int(img.rows() - 1),
                   int(shape[hull[i][0]][0]))),
                 max(0, min(int(img.cols() - 1),
                   int(shape[hull[i][0]][1]))));
    tris.emplace_back(work.solveN());
  }
  return tilt(img * T(0), tris);
}

template <typename T> typename reDig<T>::Mat reDig<T>::draw(const Mat& img, const vector<Vec3>& shape, const vector<Veci3>& hull, const bool& elim) {
  Mat result(img);
  T   M(0);
  T   m(0);
  for(int i = 0; i < shape.size(); i ++) {
    if(i) {
      M = max(M, shape[i][2]);
      m = min(m, shape[i][2]);
    } else
      M = m = shape[i][2];
  }
  auto tsrc(shape);
  if(M - m != T(0))
    for(int i = 0; i < tsrc.size(); i ++)
      tsrc[i][2] = elim ? T(0) : (tsrc[i][2] - m) / (M - m);
  for(int ii = 0; ii < hull.size(); ii ++)
    drawMatchTriangle(result, tsrc[hull[ii][0]],
                              tsrc[hull[ii][1]],
                              tsrc[hull[ii][2]]);
  return result;
}

template <typename T> typename reDig<T>::Mat reDig<T>::drawBone(const vector<Vec3>& center, const vector<T>& r, const int& rows, const int& cols) {
  Mat result(rows, cols);
  assert(center.size() == r.size());
  for(int i = 0; i < result.rows(); i ++)
    for(int j = 0; j < result.cols(); j ++)
      result(i, j) = T(0);
  Mat count(result);
  for(int i = 0; i < result.rows(); i ++)
    for(int j = 0; j < result.cols(); j ++)
      for(int l = 0; l < center.size(); l ++) {
        const auto score(pow(num_t(i) - center[l][0], num_t(2)) +
                         pow(num_t(j) - center[l][1], num_t(2)) +
                         pow(num_t(0) - center[l][2], num_t(2)) -
                         pow(r[l], num_t(2)));
        if(score <= num_t(0)) {
          result(i, j) += sqrt(abs(score));
          count(i, j)  += num_t(1);
        }
      }
  for(int i = 0; i < result.rows(); i ++)
    for(int j = 0; j < result.cols(); j ++)
      if(count(i, j) != num_t(0))
        result(i, j) /= count(i, j);
  return normalize(result, 1.);
}

template <typename T> vector<typename reDig<T>::Vec3> reDig<T>::takeShape(const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const T& ratio) {
  vector<Vec3> result(dst);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < match.srcpoints.size(); i ++)
    result[match.dstpoints[i]] += (match.transform(src[match.srcpoints[i]]) - dst[match.dstpoints[i]]) * ratio;
  return result;
}

template <typename T> vector<typename reDig<T>::Vec3> reDig<T>::takeShape(const vector<Vec3>& shape, const vector<Vec3>& center, const vector<Vec3>& outcenter, const vector<vector<int> >& attend, const T& ratio) {
  assert(center.size() == outcenter.size());
  assert(center.size() == attend.size());
  vector<Vec3> result(shape);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < center.size(); i ++) {
    const auto delta((outcenter[i] - center[i]) * ratio);
    for(int j = 0; j < attend[i].size(); j ++)
      result[attend[i][j]] += delta;
  }
  return result;
}

template <typename T> void reDig<T>::complement(vector<Mat>& dstimg, vector<Vec3>& dstcenter, const vector<vector<Mat> >& srcimg, const vector<vector<Vec3> >& srccenter, const vector<vector<int> >& attend, const vector<vector<pair<int, int> > >& a2xy, const T& iemph) {
  cerr << "p" << flush;
  static P0<T,false> p;
  static Filter<T> filter;
  const int  idx(max(0, min(int(srccenter.size()) - 1, int(floor(iemph)))));
  const auto comp(p.taylor(srccenter.size(), iemph));
  assert(srcimg.size() == srccenter.size());
  assert(srccenter[0].size() == attend.size());
  dstcenter = vector<Vec3>();
  dstimg    = vector<Mat>();
  dstcenter.resize(srccenter[0].size(), Vec3(3));
  dstimg.reserve(3);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < srccenter[idx].size(); i ++) {
    dstcenter[i][0] = dstcenter[i][1] = dstcenter[i][2] = T(0);
    for(int j = 0; j < srccenter.size(); j ++)
      for(int k = 0; k < 3; k ++)
        dstcenter[i][k] += srccenter[j][i][k] * comp[j];
  }
  for(int i = 0; i < 3; i ++)
    dstimg.emplace_back(Mat(srcimg[0][i].rows(), srcimg[0][i].cols()));
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < attend.size(); i ++)
    for(int j = 0; j < attend[i].size(); j ++) {
      const auto yy(filter.getImgPt(a2xy[i][j].first,  srcimg[idx][0].rows()));
      const auto xx(filter.getImgPt(a2xy[i][j].second, srcimg[idx][0].cols()));
      dstimg[0](yy, xx) = dstimg[1](yy, xx) = dstimg[2](yy, xx) = T(0);
      for(int k = 0; k < srccenter.size(); k ++) {
        const auto yf(filter.getImgPt(a2xy[i][j].first  + int(srccenter[k][i][0] - srccenter[idx][i][0]), srcimg[idx][0].rows()));
        const auto xf(filter.getImgPt(a2xy[i][j].second + int(srccenter[k][i][1] - srccenter[idx][i][1]), srcimg[idx][0].cols()));
        for(int kk = 0; kk < 3; kk ++)
          dstimg[kk](yy, xx) += srcimg[k][kk](yf, xf) * comp[k];
      }
      const auto n0(dstimg[0](yy, xx));
      const auto n1(dstimg[1](yy, xx));
      const auto n2(dstimg[2](yy, xx));
#if defined(_OPENMP)
#pragma omp critical
#endif
      {
        for(int y0 = yy; y0 < min(yy + vbox, int(srcimg[idx][0].rows()) - 1); y0 ++)
          for(int x0 = xx; x0 < min(xx + vbox, int(srcimg[idx][0].cols()) - 1); x0 ++) {
            dstimg[0](y0, x0) = n0;
            dstimg[1](y0, x0) = n1;
            dstimg[2](y0, x0) = n2;
          }
      }
    }
  return;
}

template <typename T> void reDig<T>::drawMatchLine(Mat& map, const Vec3& lref0, const Vec3& lref1, const T& emph) {
  int idxm(0);
  int idxM(1);
  if(abs(lref1[idxM] - lref0[idxM]) < abs(lref1[idxm] - lref0[idxm])) {
    idxm = 1;
    idxM = 0;
  }
  for(int i = 0; i <= int(abs(lref0[idxM] - lref1[idxM])); i ++) {
    const auto gidx(lref0 + (lref1 - lref0) * T(i) / abs(lref0[idxM] - lref1[idxM]));
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
  const T    lnum(sqrt(ldiff.dot(ldiff)) + T(1));
  // XXX : tan theta depend loop num, this have glitches.
  for(int k = 0; k < int(lnum * T(4)); k ++) {
    const Vec3 l0(lref0 + (lref2 - lref0) * T(k) / (lnum * T(4)));
    const Vec3 l1(lref1 + (lref2 - lref1) * T(k) / (lnum * T(4)));
    for(int i = 0; i <= int(abs(l0[idxM] - l1[idxM])) + 1; i ++) {
      const auto gidx(l0 + (l1 - l0) * T(i) / T(int(abs(l0[idxM] - l1[idxM])) + 1));
      map(max(0, min(int(gidx[0]), int(map.rows() - 1))),
          max(0, min(int(gidx[1]), int(map.cols() - 1)))) = gidx[2];
    }
  }
  return;
}

template <typename T> typename reDig<T>::Mat reDig<T>::showMatch(const Mat& dstimg, const vector<Vec3>& dst, const vector<Veci3>& hull, const T& emph) {
  Mat map(dstimg);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int k = 0; k < hull.size(); k ++) {
    drawMatchLine(map, dst[hull[k][0]], dst[hull[k][1]], emph);
    drawMatchLine(map, dst[hull[k][1]], dst[hull[k][2]], emph);
    drawMatchLine(map, dst[hull[k][2]], dst[hull[k][0]], emph);
  }
  return map;
}

template <typename T> typename reDig<T>::Mat reDig<T>::makeRefMatrix(const Mat& orig, const int& start) const {
  Mat result(orig.rows(), orig.cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < orig.rows() * orig.cols(); i ++)
    result(i % orig.rows(), i / orig.rows()) = i + start;
  return result;
}

template <typename T> typename reDig<T>::Mat reDig<T>::pullRefMatrix(const Mat& ref, const int& start, const Mat& orig) const {
  assert(orig.rows() == ref.rows() && orig.cols() == ref.cols());
  Mat result(ref.rows(), ref.cols());
  for(int i = 0; i < result.rows(); i ++)
    for(int j = 0; j < result.cols(); j ++)
      result(i, j) = T(0);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < ref.rows() * ref.cols(); i ++) {
    const int ly(i % ref.rows());
    const int lx(i / ref.rows());
    const int v(int(ref(ly, lx)) - start);
    if(0 <= v && v < orig.rows() * orig.cols())
      result(ly, lx) = orig(v % orig.rows(), v / orig.rows());
    else
      result(ly, lx) = T(0);
  }
  return result;
}

template <typename T> vector<typename reDig<T>::Veci3> reDig<T>::mesh2(const vector<Vec3>& p, const vector<int>& pp) const {
  vector<pair<Vec3, int> > sp;
  vector<pair<Vec3, int> > sp2;
  T m0(0);
  T m1(0);
  T M0(0);
  T M1(0);
  for(int i = 0; i < pp.size(); i ++) {
    sp.emplace_back(make_pair(p[pp[i]], pp[i]));
    m0 = min(m0, sp[i].first[0]);
    m1 = min(m1, sp[i].first[1]);
    M0 = max(M0, sp[i].first[0]);
    M1 = max(M1, sp[i].first[1]);
  }
  m0 -= T(1);
  m1 -= T(1);
  M0 += T(1);
  M1 += T(1);
  sp.emplace_back(make_pair(Vec3(3), p.size()));
  sp.emplace_back(make_pair(Vec3(3), p.size() + 1));
  sp.emplace_back(make_pair(Vec3(3), p.size() + 2));
  sp.emplace_back(make_pair(Vec3(3), p.size() + 3));
  sp[sp.size() - 4].first[0] = m0;
  sp[sp.size() - 4].first[1] = m1;
  sp[sp.size() - 4].first[2] = T(0);
  sp[sp.size() - 3].first    = sp[sp.size() - 4].first;
  sp[sp.size() - 3].first[0] = M0;
  sp[sp.size() - 2].first    = sp[sp.size() - 4].first;
  sp[sp.size() - 2].first[1] = M1;
  sp[sp.size() - 1].first    = sp[sp.size() - 3].first;
  sp[sp.size() - 1].first[1] = M1;
  sp[sp.size() - 4].first[0] -= T(1);
  sp[sp.size() - 3].first[0] += T(1);
  sort(sp.begin(), sp.end(), less0<pair<Vec3, int> >);
  sp2.reserve(sp.size());
  for(int i = 0; i < sp.size(); i ++)
    sp2.emplace_back(sp[sp.size() - i - 1]);
  vector<Veci3> res0;
  for(int i = 2; i < sp.size(); i ++) {
    Veci3 lres(3);
    lres[0] = lres[1] = lres[2] = i;
    for(int j = i - 1; 0 <= j; j --)
      if(sp[j].first[1] <= sp[i].first[1] &&
         (lres[1] == i ||
          abs(sp[j].first[1] - sp[i].first[1]) <
            abs(sp[lres[1]].first[1] - sp[i].first[1]) ) )
        lres[1] = j;
    for(int j = i - 1; 0 <= j; j --)
      if(sp[j].first[1] >= sp[i].first[1] &&
         ! (sp[j].first[0] == sp[lres[0]].first[0] &&
            sp[j].first[0] == sp[lres[1]].first[0]) &&
         ! (sp[j].first[1] == sp[lres[0]].first[1] &&
            sp[j].first[1] == sp[lres[1]].first[1]) &&
         (lres[2] == i ||
          abs(sp[j].first[1] - sp[i].first[1]) <
            abs(sp[lres[2]].first[1] - sp[i].first[1]) ) )
        lres[2] = j;
    res0.emplace_back(lres);
    lres[0] = lres[1] = lres[2] = i;
    for(int j = i - 1; 0 <= j; j --)
      if(sp2[j].first[1] >= sp2[i].first[1] &&
         (lres[1] == i ||
          abs(sp2[j].first[1] - sp2[i].first[1]) <
            abs(sp2[lres[1]].first[1] - sp2[i].first[1]) ) )
        lres[1] = j;
    for(int j = i - 1; 0 <= j; j --)
      if(sp2[j].first[1] <= sp2[i].first[1] &&
         ! (sp2[j].first[0] == sp2[lres[0]].first[0] &&
            sp2[j].first[0] == sp2[lres[1]].first[0]) &&
         ! (sp2[j].first[1] == sp2[lres[0]].first[1] &&
            sp2[j].first[1] == sp2[lres[1]].first[1]) &&
         (lres[2] == i ||
          abs(sp2[j].first[1] - sp2[i].first[1]) <
            abs(sp2[lres[2]].first[1] - sp2[i].first[1]) ) )
        lres[2] = j;
    for(int j = 0; j < lres.size(); j ++)
      lres[j] = sp.size() - 1 - lres[j];
    std::swap(lres[0], lres[2]);
    res0.emplace_back(lres);
  }
  vector<Veci3> res;
  for(int i = 0; i < res0.size(); i ++)
    if(sp[res0[i][0]].second < p.size() &&
       sp[res0[i][1]].second < p.size() &&
       sp[res0[i][2]].second < p.size()) {
      for(int j = 0; j < res0[i].size(); j ++)
        res0[i][j] = sp[res0[i][j]].second;
      for(int j = 0; j < res.size(); j ++)
        if(res[j] == res0[i])
          goto nofix;
      res.emplace_back(res0[i]);
     nofix:
      ;
    }
  return res;
}

template <typename T> inline bool reDig<T>::isClockwise(const Vec3 p[3]) const {
  Mat3x3 dc(3, 3);
  for(int i = 0; i < 3; i ++) {
    dc(i, 0) = T(1);
    dc(i, 1) = p[i][0];
    dc(i, 2) = p[i][1];
  }
  return dc.determinant() <= T(0);
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
    if(mask(y, x) > T(1) / T(2)) {
      elim.emplace_back(i);
      after.emplace_back(- 1);
    } else
      after.emplace_back(ii ++);
  }
  for(int i = 0; i < polys.size(); i ++)
    if(std::binary_search(elim.begin(), elim.end(), polys[i][0]) ||
       std::binary_search(elim.begin(), elim.end(), polys[i][1]) ||
       std::binary_search(elim.begin(), elim.end(), polys[i][2]))
      elimp.emplace_back(i);
  for(int i = 0; i < elim.size(); i ++)
    points.erase(points.begin() + (elim[i] - i));
  for(int i = 0; i < elimp.size(); i ++)
    polys.erase(polys.begin() + (elimp[i] - i));
  for(int i = 0; i < polys.size(); i ++)
    for(int j = 0; j < polys[i].size(); j ++)
      polys[i][j] = after[polys[i][j]];
  return;
}

template <typename T> typename reDig<T>::Mat reDig<T>::reShape(const Mat& cbase, const Mat& vbase, const int& count) {
  assert(cbase.rows() && cbase.cols() && vbase.rows() && vbase.cols());
  assert(cbase.rows() == vbase.rows() && cbase.cols() == vbase.cols());
  vector<pair<T, pair<int, int> > > vpoints;
  vpoints.reserve(vbase.rows() * vbase.cols());
  for(int i = 0; i < vbase.rows(); i ++)
    for(int j = 0; j < vbase.cols(); j ++)
      vpoints.emplace_back(make_pair(vbase(i, j), make_pair(i, j)));
  sort(vpoints.begin(), vpoints.end());
  Mat res(vbase.rows(), vbase.cols());
  T   avg(0);
  for(int i = 0, ii = 0; i < vpoints.size(); i ++) {
    if(abs(vpoints[i].first - vpoints[ii].first) < abs(vpoints[vpoints.size() - 1].first - vpoints[0].first) / T(count) && i < vpoints.size() - 1)
      avg += cbase(vpoints[i].second.first, vpoints[i].second.second);
    else if(i != ii) {
      if(i == vpoints.size() - 1) {
        avg += cbase(vpoints[i].second.first, vpoints[i].second.second);
        i ++;
      }
      avg /= i - ii;
      for(int j = ii; j < i; j ++)
        res(vpoints[j].second.first, vpoints[j].second.second) = avg;
      avg  = T(0);
      ii   = i;
    }
  }
  return res;
}

template <typename T> typename reDig<T>::Mat reDig<T>::reColor(const Mat& cbase, const Mat& vbase, const int& count) {
  assert(cbase.rows() && cbase.cols() && vbase.rows() && vbase.cols());
  vector<pair<T, pair<int, int> > > vpoints;
  vector<pair<T, pair<int, int> > > cpoints;
  vpoints.reserve(vbase.rows() * vbase.cols());
  cpoints.reserve(cbase.rows() * cbase.cols());
  for(int i = 0; i < vbase.rows(); i ++)
    for(int j = 0; j < vbase.cols(); j ++)
      vpoints.emplace_back(make_pair(vbase(i, j), make_pair(i, j)));
  for(int i = 0; i < cbase.rows(); i ++)
    for(int j = 0; j < cbase.cols(); j ++)
      cpoints.emplace_back(make_pair(cbase(i, j), make_pair(i, j)));
  sort(vpoints.begin(), vpoints.end());
  sort(cpoints.begin(), cpoints.end());
  Vec vv(vpoints.size());
  Vec cc(cpoints.size());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < vpoints.size(); i ++)
    vv[i] = vpoints[i].first;
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < cpoints.size(); i ++)
    cc[i] = cpoints[i].first;
  const auto ccc(Decompose<T>(count).mimic(cc, vv));
  Mat res(cbase);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < cpoints.size(); i ++)
    res(cpoints[i].second.first, cpoints[i].second.second) = ccc[i];
  return res;
}

template <typename T> typename reDig<T>::Mat reDig<T>::reColor3(const Mat& cbase, const Mat& vbase, const int& count) {
  assert(cbase.rows() && cbase.cols() && vbase.rows() && vbase.cols());
  vector<pair<T, pair<int, int> > > vpoints;
  vector<pair<T, pair<int, int> > > cpoints;
  vpoints.reserve(vbase.rows() * vbase.cols());
  cpoints.reserve(cbase.rows() * cbase.cols());
  for(int i = 0; i < vbase.rows(); i ++)
    for(int j = 0; j < vbase.cols(); j ++)
      vpoints.emplace_back(make_pair(vbase(i, j), make_pair(i, j)));
  for(int i = 0; i < cbase.rows(); i ++)
    for(int j = 0; j < cbase.cols(); j ++)
      cpoints.emplace_back(make_pair(cbase(i, j), make_pair(i, j)));
  sort(vpoints.begin(), vpoints.end());
  sort(cpoints.begin(), cpoints.end());
  Mat res(cbase.rows(), cbase.cols());
  for(int i = 0; i < count; i ++) {
    T scorev(0);
    T scorec(0);
    const auto vstart(i * int(vpoints.size() / count));
    const auto vend(min((i + 1) * int(vpoints.size() / count), int(vpoints.size())));
    const auto cstart(i * int(cpoints.size() / count));
    const auto cend(min((i + 1) * int(cpoints.size() / count), int(cpoints.size())));
    for(int j = vstart; j < vend; j ++)
      scorev += vpoints[j].first;
    for(int j = cstart; j < cend; j ++)
      scorec += cpoints[j].first;
    scorev /= T(vend - vstart);
    scorec /= T(cend - cstart);
    for(int j = cstart; j < cend; j ++)
      res(cpoints[j].second.first, cpoints[j].second.second) =
        cpoints[j].first * scorev / scorec;
  }
  return res;
}

template <typename T> typename reDig<T>::Mat reDig<T>::reColor(const Mat& cbase, const int& count, const T& intensity) {
  assert(cbase.rows() && cbase.cols());
  vector<pair<T, pair<int, int> > > cpoints;
  cpoints.reserve(cbase.rows() * cbase.cols());
  for(int i = 0; i < cbase.rows(); i ++)
    for(int j = 0; j < cbase.cols(); j ++)
      cpoints.emplace_back(make_pair(cbase(i, j), make_pair(i, j)));
  sort(cpoints.begin(), cpoints.end());
  Vec cc(cpoints.size());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < cpoints.size(); i ++)
    cc[i] = cpoints[i].first;
  const auto ccc(Decompose<T>(count).emphasis(cc, intensity));
  Mat res(cbase);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < cpoints.size(); i ++)
    res(cpoints[i].second.first, cpoints[i].second.second) = ccc[i];
  return res;
}

template <typename T> typename reDig<T>::Mat reDig<T>::reTrace(const Mat& dst, const Mat& src, const T& intensity, const int& count) {
  pair<Vec, Vec> pdst, psrc;
  pair<pair<int, int>, pair<int, int> > dsthw, srchw;
  prepTrace(pdst, dsthw, dst);
  prepTrace(psrc, srchw, src);
  Decompose<T> decom(count);
  return applyTrace(make_pair(
    decom.mimic(pdst.first,  psrc.first,  intensity),
    decom.mimic(pdst.second, psrc.second, intensity)), make_pair(dsthw, srchw));
}

template <typename T> typename reDig<T>::Mat reDig<T>::reTrace(const Mat& dst, const T& intensity, const int& count) {
  pair<Vec, Vec> pdst;
  pair<pair<int, int>, pair<int, int> > hw;
  prepTrace(pdst, hw, dst);
  Decompose<T> decom(count);
  return applyTrace(make_pair(
    decom.emphasis(pdst.first,  intensity),
    decom.emphasis(pdst.second, intensity)), make_pair(hw, hw));
}

template <typename T> void reDig<T>::prepTrace(pair<Vec, Vec>& v, pair<pair<int, int>, pair<int, int> >& hw, const Mat& mask) {
  vector<Vec3>  pdst;
  vector<Veci3> facets;
  getTileVec(mask, pdst, facets);
  const auto idsts(getEdges(mask, pdst));
  assert(idsts.size());
  int iidst(0);
  for(int i = 1; i < idsts.size(); i ++)
    if(idsts[iidst].size() < idsts[i].size()) iidst = i;
  const auto& idst(idsts[iidst]);
  cerr << iidst << "/" << idsts.size() << ":" << idst.size() << flush;
  assert(idst.size());
  auto& yMd(hw.first.first   = pdst[idst[0]][0]);
  auto& ymd(hw.second.first  = pdst[idst[0]][0]);
  auto& xMd(hw.first.second  = pdst[idst[0]][1]);
  auto& xmd(hw.second.second = pdst[idst[0]][1]);
  for(int i = 1; i < idst.size(); i ++) {
    yMd = max(yMd, int(pdst[idst[i]][0]));
    ymd = min(ymd, int(pdst[idst[i]][0]));
    xMd = max(xMd, int(pdst[idst[i]][1]));
    xmd = min(xMd, int(pdst[idst[i]][1]));
  }
  v.first.resize(idst.size());
  v.second.resize(idst.size());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < idst.size(); i ++) {
    v.first[i]  = pdst[idst[i]][0];
    v.second[i] = pdst[idst[i]][1];
  }
  return;
}

template <typename T> typename reDig<T>::Mat reDig<T>::applyTrace(const pair<Vec, Vec>& v, const pair<pair<pair<int, int>, pair<int, int> >, pair<pair<int, int>, pair<int, int> > >& hw) {
  const auto& vy(v.first);
  const auto& vx(v.second);
  pair<int, int> MM(make_pair(int(vy[0]), int(vx[0])));
  auto mm(MM);
  for(int i = 1; i < vy.size(); i ++) {
    MM.first  = max(MM.first,  int(vy[i]));
    MM.second = max(MM.second, int(vx[i]));
    mm.first  = min(mm.first,  int(vy[i]));
    mm.second = min(mm.second, int(vx[i]));
  }
  const auto& dsthw(hw.first);
  const auto& srchw(hw.second);
  const auto yy(max(dsthw.first.first - dsthw.second.first,
                    srchw.first.first - srchw.second.first) + 1);
  const auto xx(max(dsthw.first.second - dsthw.second.second,
                    srchw.first.second - srchw.second.second) + 1);
  const auto yyy(MM.first - mm.first + 1);
  const auto xxx(MM.second - mm.second + 1);
  cerr << ":" << yy << ":" << xx << ":" << yyy << ":" << xxx << flush;
  Mat res(yy, xx);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < res.rows(); i ++)
    for(int j = 0; j < res.cols(); j ++)
      res(i, j) = T(0);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < vy.size(); i ++) {
    Vec3 v0(3);
    Vec3 v1(3);
    v0[2] = v1[2] = T(0);
    v0[0] = (vy[i] - T(mm.first))  / T(yyy) * T(yy);
    v0[1] = (vx[i] - T(mm.second)) / T(xxx) * T(xx);
    if(i == vy.size() - 1) {
      v1[0] = (vy[0] - T(mm.first))  / T(yyy) * T(yy);
      v1[1] = (vx[0] - T(mm.second)) / T(xxx) * T(xx);
    } else {
      v1[0] = (vy[i + 1] - T(mm.first))  / T(yyy) * T(yy);
      v1[1] = (vx[i + 1] - T(mm.second)) / T(xxx) * T(xx);
    }
    drawMatchLine(res, v0, v1, T(1));
  }
  return res;
}

template <typename T> typename reDig<T>::Mat reDig<T>::reImage(const Mat& dst, const Mat& src, const T& intensity, const int& count) {
  Decompose<T> decom(count);
  assert(dst.rows() == src.rows() && dst.cols() == src.cols());
  Mat res(dst.rows(), dst.cols());
  for(int i = 0; i < res.rows(); i ++)
    res.row(i) = decom.mimic(dst.row(i), src.row(i), intensity);
  return res;
}

template <typename T> typename reDig<T>::Mat reDig<T>::reImage(const Mat& dst, const T& intensity, const int& count) {
  Decompose<T> decom(count);
  Mat res(dst.rows(), dst.cols());
  for(int i = 0; i < res.rows(); i ++)
    res.row(i) = decom.emphasis(dst.row(i), intensity);
  return res;
}

template <typename T> typename reDig<T>::Mat reDig<T>::catImage(const vector<Mat>& imgs) {
  assert(imgs.size());
  for(int i = 1; i < imgs.size(); i ++)
    assert(imgs[i].rows() == imgs[0].rows() && imgs[i].cols() == imgs[0].cols());
  Decompose<T> dec(imgs[0].cols());
  Catg<T> cat(imgs[0].cols());
  for(int i = 0; i < imgs.size(); i ++) {
    for(int j = 0; j < imgs[i].rows(); j ++) {
      cat.inq(dec.next(imgs[i].row(j)));
      std::cerr << "." << std::flush;
    }
  }
  cat.compute();
  std::vector<std::pair<T, int> > scat;
  scat.resize(cat.lambda.size());
  for(int i = 0; i < cat.lambda.size(); i ++) {
    scat[i].first  = abs(cat.lambda[i]);
    scat[i].second = i;
    assert(isfinite(scat[i].first));
  }
  std::sort(scat.begin(), scat.end());
  Mat res(imgs[0].cols(), imgs[0].cols());
  for(int i = 0; i < res.rows(); i ++) {
    const auto& ii(scat[scat.size() - 1 - i].second);
    res.row(i) = cat.Left.col(ii);
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
      if(T(1) / T(2) < mask(yy, xx))
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
      if(mask(i, j) < T(1) / T(2))
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
        const T x2(x - y + T(sti.second));
        const T y2(x + y + T(sti.first));
        if(T(0) <= x2 && x2 < T(mask.cols()) &&
           T(0) <= y2 && y2 < T(mask.rows()) &&
           mask(y2, x2) < T(1) / T(2))
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
      result.emplace_back(vector<int>());
      // apply to points index.
      vector<int> pj;
      for(int i = 0; i < e.size(); i ++) {
        const auto& s(store[e[i]]);
        vector<pair<T, int> > distances;
        for(int j = 0; j < points.size(); j ++)
          distances.emplace_back(make_pair(
                                  pow(T(s.first)  - points[j][0], T(2)) +
                                    pow(T(s.second) - points[j][1], T(2)),
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

template <typename T> typename reDig<T>::Mat reDig<T>::rgb2d(const Mat rgb[3]) {
  Mat xyz[3];
  rgb2xyz(xyz, rgb);
  Mat result(rgb[0].rows(), rgb[0].cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int j = 0; j < rgb[0].rows(); j ++) {
    for(int k = 0; k < rgb[0].cols(); k ++)
      result(j, k) = sqrt(xyz[0](j, k) * xyz[0](j, k) + xyz[1](j, k) * xyz[1](j, k) + xyz[2](j, k) * xyz[2](j, k));
  }
  return result;
}

template <typename T> void reDig<T>::rgb2xyz(Mat xyz[3], const Mat rgb[3]) {
  // CIE 1931 XYZ from wikipedia.org
  Mat3x3 mRGB2XYZ(3, 3);
  mRGB2XYZ(0, 0) = T(49000);
  mRGB2XYZ(0, 1) = T(31000);
  mRGB2XYZ(0, 2) = T(20000);
  mRGB2XYZ(1, 0) = T(17697);
  mRGB2XYZ(1, 1) = T(81240);
  mRGB2XYZ(1, 2) = T( 1063);
  mRGB2XYZ(2, 0) = T(0);
  mRGB2XYZ(2, 1) = T( 1000);
  mRGB2XYZ(2, 2) = T(99000);
  mRGB2XYZ /= T(17697);
  xyz[0] = rgb[0] * mRGB2XYZ(0, 0) + rgb[1] * mRGB2XYZ(0, 1) + rgb[2] * mRGB2XYZ(0, 2);
  xyz[1] = rgb[0] * mRGB2XYZ(1, 0) + rgb[1] * mRGB2XYZ(1, 1) + rgb[2] * mRGB2XYZ(1, 2);
  xyz[2] = rgb[0] * mRGB2XYZ(2, 0) + rgb[1] * mRGB2XYZ(2, 1) + rgb[2] * mRGB2XYZ(2, 2);
  return;
}

template <typename T> void reDig<T>::xyz2rgb(Mat rgb[3], const Mat xyz[3]) {
  // CIE 1931 XYZ from wikipedia.org
  Mat3x3 mRGB2XYZ(3, 3);
  mRGB2XYZ(0, 0) = T(49000);
  mRGB2XYZ(0, 1) = T(31000);
  mRGB2XYZ(0, 2) = T(20000);
  mRGB2XYZ(1, 0) = T(17697);
  mRGB2XYZ(1, 1) = T(81240);
  mRGB2XYZ(1, 2) = T( 1063);
  mRGB2XYZ(2, 0) = T(0);
  mRGB2XYZ(2, 1) = T( 1000);
  mRGB2XYZ(2, 2) = T(99000);
  mRGB2XYZ /= T(17697);
  Mat3x3 mXYZ2RGB(3, 3);
  // XXX: very slow with simplelin.hh
  mXYZ2RGB = mRGB2XYZ.inverse();
  rgb[0] = xyz[0] * mXYZ2RGB(0, 0) + xyz[1] * mXYZ2RGB(0, 1) + xyz[2] * mXYZ2RGB(0, 2);
  rgb[1] = xyz[0] * mXYZ2RGB(1, 0) + xyz[1] * mXYZ2RGB(1, 1) + xyz[2] * mXYZ2RGB(1, 2);
  rgb[2] = xyz[0] * mXYZ2RGB(2, 0) + xyz[1] * mXYZ2RGB(2, 1) + xyz[2] * mXYZ2RGB(2, 2);
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

template <typename T> typename reDig<T>::Mat reDig<T>::normalize(const Mat& data, const T& upper) {
  T MM(0), mm(0);
  bool fixed(false);
  for(int i = 0; i < data.rows(); i ++)
    for(int j = 0; j < data.cols(); j ++)
      if(! fixed || (isfinite(data(i, j)) && ! isinf(data(i, j)) && ! isnan(data(i, j)))) {
        if(! fixed)
          MM = mm = data(i, j);
        else {
          MM = max(MM, data(i, j));
          mm = min(mm, data(i, j));
        }
        fixed = true;
      }
  if(MM == mm || ! fixed)
    return data;
  Mat result(data);
  for(int i = 0; i < data.rows(); i ++)
    for(int j = 0; j < data.cols(); j ++) {
      if(isfinite(result(i, j)) && ! isinf(data(i, j)) && ! isnan(result(i, j)))
        result(i, j) -= mm;
      else
        result(i, j)  = T(0);
      assert(T(0) <= result(i, j) && result(i, j) <= MM - mm);
    }
  return result * upper / (MM - mm);
}

template <typename T> void reDig<T>::normalize(Mat data[3], const T& upper) {
  T MM(0), mm(0);
  bool fixed(false);
  for(int k = 0; k < 3; k ++)
    for(int i = 0; i < data[k].rows(); i ++)
      for(int j = 0; j < data[k].cols(); j ++)
        if(! fixed || (isfinite(data[k](i, j)) && ! isinf(data[k](i, j)) && ! isnan(data[k](i, j)))) {
          if(! fixed)
            MM = mm = data[k](i, j);
          else {
            MM = max(MM, data[k](i, j));
            mm = min(mm, data[k](i, j));
          }
          fixed = true;
        }
  if(MM == mm || ! fixed)
    return;
  for(int k = 0; k < 3; k ++) {
    for(int i = 0; i < data[k].rows(); i ++)
      for(int j = 0; j < data[k].cols(); j ++) {
        if(isfinite(data[k](i, j)) && ! isinf(data[k](i, j)) && ! isnan(data[k](i, j)))
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
      res.emplace_back(data(i, j));
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

template <typename T> void reDig<T>::autoLevel(Mat data[3], const int& count) {
  vector<T> res;
  res.reserve(data[0].rows() * data[0].cols() * 3);
  for(int k = 0; k < 3; k ++)
    for(int i = 0; i < data[k].rows(); i ++)
      for(int j = 0; j < data[k].cols(); j ++)
        res.emplace_back(data[k](i, j));
  sort(res.begin(), res.end());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int k = 0; k < 3; k ++)
    for(int i = 0; i < data[k].rows(); i ++)
      for(int j = 0; j < data[k].cols(); j ++)
        data[k](i, j) = max(min(data[k](i, j), res[res.size() - count - 1]), res[count]);
  return;
}

// get bump with multiple scale and vectorized result.
template <typename T> void reDig<T>::getTileVec(const Mat& in, vector<Vec3>& geoms, vector<Veci3>& delaunay) {
  // get vectorize.
  geoms = vector<Vec3>();
  T aavg(0);
  for(int i = 0; i < in.rows(); i ++)
    for(int j = 0; j < in.cols(); j ++)
      aavg += in(i, j);
  aavg /= T(in.rows() * in.cols());
  for(int i = 0; i < in.rows() / vbox + 1; i ++)
    for(int j = 0; j < in.cols() / vbox + 1; j ++) {
      if(in.rows() < (i + 1) * vbox ||
         in.cols() < (j + 1) * vbox) {
        Vec3 gbuf(3);
        gbuf[0] = T(i * vbox);
        gbuf[1] = T(j * vbox);
        gbuf[2] = geoms[geoms.size() - 1][2];
        geoms.emplace_back(gbuf);
      } else {
        T avg(0);
        for(int ii = i * vbox; ii < (i + 1) * vbox; ii ++)
          for(int jj = j * vbox; jj < (j + 1) * vbox; jj ++)
            avg += in(ii, jj);
        Vec3 work(3);
        work[0] = T(i * vbox);
        work[1] = T(j * vbox);
        work[2] = sqrt(T(in.rows() * in.cols())) * rz * abs(avg / T(vbox) / T(vbox) - aavg);
        geoms.emplace_back(work);
      }
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
      work[0]  = T((i - 1) * (in.cols() / vbox + 1) + j);
      work[1]  = T( i      * (in.cols() / vbox + 1) + j);
      work[2]  = T( i      * (in.cols() / vbox + 1) + j + 1);
      work2[0] = T((i - 1) * (in.cols() / vbox + 1) + j);
      work2[2] = T((i - 1) * (in.cols() / vbox + 1) + j + 1);
      work2[1] = T( i      * (in.cols() / vbox + 1) + j + 1);
      delaunay.emplace_back(work);
      delaunay.emplace_back(work2);
    }
  return;
}

template <typename T> void reDig<T>::getBone(const Mat& in, const vector<Vec3>& geoms, vector<Vec3>& center, vector<T>& r, vector<vector<int> >& attend, const T& thresh) {
  int idx(0);
  r      = vector<T>();
  center = vector<Vec3>();
  attend = vector<vector<int> >();
  for(int i = 0; i < in.rows() / vbox + 1; i ++) {
    vector<T> workx;
    vector<T> workz;
    for(int j = 0; j < in.cols() / vbox + 1; j ++) {
      assert(0 <= idx && idx < geoms.size());
      workx.emplace_back(geoms[idx][1]);
      workz.emplace_back(geoms[idx][2]);
      // N.B. for any k, on the line (workx[k] - x)^2 + (workz[k] - z)^2 = r^2.
      T x(0);
      T z(0);
      for(int k = 0; k < workx.size(); k ++) {
        x += workx[k];
        z += workz[k];
      }
      x /= T(workx.size());
      z /= T(workz.size());
      T r(0);
      for(int k = 0; k < workx.size(); k ++)
        r += (workx[k] - x) * (workx[k] - x) + (workz[k] - z) * (workz[k] - z);
      r  = sqrt(r / T(workx.size()));
      T err(0);
      for(int k = 0; k < workx.size(); k ++)
        err += abs((workx[k] - x) * (workx[k] - x) +
                   (workz[k] - z) * (workz[k] - z) - r * r);
      err = sqrt(err / T(workx.size()));
      if(j && err < thresh * thresh * sqrt(T(in.rows() * in.cols())) && center.size()) {
        center[center.size() - 1][1] = x;
        center[center.size() - 1][2] = z;
      } else {
        auto gc(geoms[idx]);
        gc[1] = x;
        gc[2] = z;
        center.emplace_back(gc);
        attend.emplace_back(vector<int>());
        workx = vector<T>();
        workz = vector<T>();
      }
      attend[attend.size() - 1].emplace_back(idx);
      idx ++;
    }
  }
  if(attend.size() && !attend[attend.size() - 1].size()) {
    center.resize(center.size() - 1);
    attend.resize(attend.size() - 1);
  }
  vector<Vec3>         newcenter;
  vector<vector<int> > newattend;
  assert(attend.size() == center.size());
  int i;
  for(i = 0; i < attend.size(); i ++) {
    auto z(center[i] * attend[i].size());
    int  zc(attend[i].size());
    T    rr(0);
    int  lastj(i);
    for(int j = i + 1; j < attend.size(); j ++) {
      // N.B. for any k, on the line ||center[k] - z||^2 = r^2.
      // <=> sum_k (ck_0 - z_0)^2 + (ck_1 - z_1)^2 + (ck_2 - z_2)^2 == r^2.
      const auto newz(z + center[j] * T(attend[j].size()));
      const auto newz0(newz / T(zc + attend[j].size()));
      T   newr(0);
      int cnt(0);
      for(int jj = i; jj <= j; jj ++)
        for(int k = 0; k < attend[jj].size(); k ++) {
          const auto diff(geoms[attend[jj][k]] - newz0);
          newr += diff.dot(diff);
          cnt ++;
        }
      newr = sqrt(newr / T(cnt));
      T err(0);
      for(int jj = i; jj <= j; jj ++)
        for(int k = 0; k < attend[jj].size(); k ++) {
          const auto diff(geoms[attend[jj][k]] - newz0);
          err += abs(diff.dot(diff) - newr * newr);
        }
      err = sqrt(err / T(cnt));
      if(err < thresh * sqrt(T(in.rows() * in.cols()))) {
        z   = newz;
        zc += attend[j].size();
        rr  = newr;
        lastj = j;
      } else
        break;
    }
    r.emplace_back(rr);
    newcenter.emplace_back(z / T(zc));
    newattend.emplace_back(attend[i]);
    for(int k = i + 1; k <= lastj; k ++)
      newattend[newattend.size() - 1].insert(
        newattend[newattend.size() - 1].end(),
        attend[k].begin(), attend[k].end());
    i = lastj;
  }
  center = newcenter;
  attend = newattend;
  return;
}

template <typename T> vector<vector<std::pair<int, int> > > reDig<T>::getReverseLookup(const vector<vector<int> >& attend, const Mat& refimg) {
  vector<vector<std::pair<int, int> > > res;
  res.resize(attend.size());
  const auto h(refimg.rows() / vbox + 1);
  const auto w(refimg.cols() / vbox + 1);
  for(int i = 0; i < attend.size(); i ++) {
    res[i].resize(attend[i].size());
    for(int j = 0; j < attend[i].size(); j ++) {
      const auto& a0ij(attend[i][j]);
            auto& a2xyij(res[i][j]);
      a2xyij.first  = (a0ij / w) * vbox;
      a2xyij.second = (a0ij % w) * vbox;
    }
  }
  return res;
}

template <typename T> vector<typename reDig<T>::Vec3> reDig<T>::copyBone(const vector<Vec3>& centerdst, const vector<T>& rdst, const vector<Vec3>& centersrc, const vector<T>& rsrc) {
  assert(centerdst.size() == rdst.size());
  assert(centersrc.size() == rsrc.size());
  assert(centerdst.size());
  vector<Vec3> result(centerdst);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < centerdst.size(); i ++) {
    int midx(i ? 0 : 1);
    for(int j = 0; j < centersrc.size(); j ++) if(i != j) {
      const auto diff0(centerdst[i] - centersrc[j]);
      const auto diff1(centerdst[i] - centersrc[midx]);
      if(diff0.dot(diff0) + (rdst[i] - rsrc[j]) * (rdst[i] - rsrc[j]) <=
         diff1.dot(diff1) + (rdst[i] - rsrc[midx]) * (rdst[i] - rsrc[midx]))
        midx = j;
    }
    result[i] = centersrc[midx];
  }
  return result;
}

template <typename T> inline typename reDig<T>::Triangles reDig<T>::makeTriangle(const int& u, const int& v, const Mat& in, const Mat& bump, const int& flg) {
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

template <typename T> inline bool reDig<T>::sameSide2(const Vec2& p0, const Vec2& p1, const Vec2& p2, const Vec2& q, const bool& extend, const T& err) const {
  const Vec2 dlt(p1 - p0);
        Vec2 dp(p2 - p0);
  if(T(0) < dlt.dot(dlt))
    dp -= dlt * dlt.dot(dp) / dlt.dot(dlt);
  // N.B. dp.dot(p1 - p0) >> 0.
  return dp.dot(q - p0) >= (extend ? - T(1) : T(1)) * (abs(dp[0]) + abs(dp[1])) * err;
}

template <typename T> inline bool reDig<T>::sameSide3(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& q, const bool& extend, const T& err) const {
  Vec2 q0(2), q1(2), q2(2), qq(2);
  q0[0] = p0[0]; q0[1] = p0[1]; q1[0] = p1[0]; q1[1] = p1[1];
  q2[0] = p2[0]; q2[1] = p2[1]; qq[0] = q[0];  qq[1] = q[1];
  return sameSide2(q0, q1, q2, qq, extend, err);
}

// <[x, y, t], triangle.n> == triangle.z
template <typename T> inline bool reDig<T>::onTriangle(T& z, const Triangles& tri, const Vec2& geom) {
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
  return (sameSide3(tri.p.col(0), tri.p.col(1), tri.p.col(2), camera, true, T(0125) / T(1000)) &&
          sameSide3(tri.p.col(1), tri.p.col(2), tri.p.col(0), camera, true, T(0125) / T(1000)) &&
          sameSide3(tri.p.col(2), tri.p.col(0), tri.p.col(1), camera, true, T(0125) / T(1000)));
}

template <typename T> match_t<T> reDig<T>::tiltprep(const Mat& in, const int& idx, const int& samples, const T& psi) {
  const T theta(T(2) * Pi * T(idx) / T(samples));
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
  pcenter[0] = T(in.rows() - 1) / T(2);
  pcenter[1] = T(in.cols() - 1) / T(2);
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
    if(T(0) <= points[polys[i][0]][0] && points[polys[i][0]][0] < T(in.rows()) &&
       T(0) <= points[polys[i][0]][1] && points[polys[i][0]][1] < T(in.cols()))
      work.c = in(int(points[polys[i][0]][0]),
                  int(points[polys[i][0]][1]));
    else
      work.c = T(0);
    result.emplace_back(work.solveN());
  }
  return result;
}

template <typename T> typename reDig<T>::Mat reDig<T>::tilt(const Mat& in, const Mat& bump, const match_t<T>& m, const T& z0) {
  assert(in.rows() == bump.rows() && in.cols() == bump.cols());
  vector<Triangles> triangles;
  triangles.reserve((in.rows() - 1) * (in.cols() - 1) * 2);
  for(int i = 0; i < in.rows() - 1; i ++)
    for(int j = 0; j < in.cols() - 1; j ++) {
      triangles.emplace_back(makeTriangle(i, j, in, bump, false));
      triangles.emplace_back(makeTriangle(i, j, in, bump, true));
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
#if defined(_OPENMP)
  SimpleMatrix<omp_lock_t> lock(result.rows(), result.cols());
  for(int i = 0; i < lock.rows(); i ++)
    for(int j = 0; j < lock.cols(); j ++)
      omp_init_lock(&lock(i, j));
#endif
  // able to boost with divide and conquer.
#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic)
#endif
  for(int j = 0; j < triangles.size(); j ++) {
    const Triangles& tri(triangles[j]);
    int ll = int( min(min(tri.p(0, 0), tri.p(0, 1)), tri.p(0, 2)));
    int rr = int(ceil(max(max(tri.p(0, 0), tri.p(0, 1)), tri.p(0, 2)))) + 1;
    int bb = int( min(min(tri.p(1, 0), tri.p(1, 1)), tri.p(1, 2)));
    int tt = int(ceil(max(max(tri.p(1, 0), tri.p(1, 1)), tri.p(1, 2)))) + 1;
    for(int y = max(0, ll); y < min(rr, int(in.rows())); y ++)
      for(int x = max(0, bb); x < min(tt, int(in.cols())); x ++) {
        T z;
        Vec2 midgeom(2);
        midgeom[0] = y;
        midgeom[1] = x;
        if(onTriangle(z, tri, midgeom) && isfinite(z)) {
          {
#if defined(_OPENMP)
            omp_set_lock(&lock(y, x));
#endif
            if(zb(y, x) < z) {
              result(y, x) = tri.c;
              zb(y, x)     = z;
            }
#if defined(_OPENMP)
            omp_unset_lock(&lock(y, x));
#endif
          }
        }
      }
  }
#if defined(_OPENMP)
  for(int i = 0; i < lock.rows(); i ++)
    for(int j = 0; j < lock.cols(); j ++)
      omp_destroy_lock(&lock(i, j));
#endif
  return result;
}

template <typename T> inline int reDig<T>::getImgPtRecursive(const int& h, const int& y) const {
  return ((y % h) + h) % h;
}

#define _REDIG_
#endif

