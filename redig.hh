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

#if !defined(_REDIG_)

template <typename T> class match_t;
template <typename T> static inline T getImgPt(const T& y, const T& h);
template <typename T> SimpleMatrix<T> rotate(const SimpleMatrix<T>& d, const T& theta);
template <typename T> static inline SimpleMatrix<T> center(const SimpleMatrix<T>& dr, const SimpleMatrix<T>& d);

using std::swap;
using std::abs;
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
using std::binary_search;
using std::unique;

template <typename T> static inline bool less0(const T& x, const T& y) {
  return x.first[0] < y.first[0] || (x.first[0] == y.first[0] && x.first[1] < y.first[1]);
}

template <typename T> static inline bool lessf(const T& x, const T& y) {
  return x.first < y.first;
}

template <typename T> class triangles_t {
public:
  typedef SimpleMatrix<T> Mat;
  typedef SimpleVector<T> Vec;
  Mat p;
  Vec n;
  T   c;
  T   z;
  inline triangles_t() {
    p = Mat(3, 3);
    n = Vec(3);
  }
  inline triangles_t<T>& rotate(const Mat& R, const Vec& origin) {
    for(int i = 0; i < 3; i ++)
      p.setCol(i, R * (p.col(i) - origin) + origin);
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
  typedef SimpleMatrix<T>   Mat;
  typedef SimpleVector<T>   Vec;
  typedef SimpleVector<int> Veci;
  typedef triangles_t<T>    Triangles;
  
  inline reDig();
  inline ~reDig();
  inline void initialize(const int& vbox, const T& rz = - T(1));
  Mat  draw(const Mat& img, const vector<Vec>& shape, const vector<Vec>& emph, const vector<Veci>& hull);
  Mat  draw(const Mat& img, const vector<Vec>& shape, const vector<Veci>& hull, const bool& elim = false);
  vector<Vec> takeShape(const vector<Vec>& dst, const vector<Vec>& src, const match_t<T>& match, const T& ratio);
  Mat  showMatch(const Mat& dstimg, const vector<Vec>& dst, const vector<Veci>& hull, const T& emph = T(1));
  Mat  makeRefMatrix(const Mat& orig, const int& start) const;
  Mat  pullRefMatrix(const Mat& ref, const int& start, const Mat& orig) const;
  vector<Veci> mesh2(const vector<Vec>& p, const vector<int>& pp) const;
  vector<int>  edge(const vector<Vec>& p, const vector<int>& pp) const;
  void maskVectors(vector<Vec>& points, const vector<Veci>& polys, const Mat& mask);
  void maskVectors(vector<Vec>& points, vector<Veci>& polys, const Mat& mask);
  Mat  reShape(const Mat& cbase, const Mat& vbase, const int& count = 20, const T& thresh = T(1) / T(128));
  Mat  reColor(const Mat& cbase, const Mat& vbase, const int& count = 20, const T& intensity = T(1));
  Mat  reColor3(const Mat& cbase, const Mat& vbase, const int& count = 20);
  Mat  reColor(const Mat& cbase, const int& count = 20, const T& intensity = T(1));
  Mat  reTrace(const Mat& dst, const Mat& src, const T& intensity, const int& count = 20);
  Mat  reTrace(const Mat& dst, const T& intensity, const int& count = 20);
  Mat  reImage(const Mat& dst, const Mat& src, const T& intensity, const int& count = 20);
  Mat  reImage(const Mat& dst, const T& intensity, const int& count = 20);
  Mat  optImage(const vector<pair<Mat, Mat> >& img, const int& comp) const;
  Mat  compImage(const Mat& in, const Mat& opt, const int& comp) const;
  vector<Mat> catImage(const vector<Mat>& rep, const vector<Mat>& imgs, const int& cs = 40);
  vector<Mat> compositeImage(const vector<Mat>& imgs);
  Mat  bump(const Mat& color, const Mat& bumpm, const T& psi, const int& n = 0, const int& origin = - 1) const;
  vector<vector<int> > floodfill(const Mat& mask, const vector<Vec>& points);
  Mat  rgb2d(const Mat rgb[3]);
  void rgb2xyz(Mat xyz[3], const Mat rgb[3]);
  void xyz2rgb(Mat rgb[3], const Mat xyz[3]);
  Mat  contrast(const Mat& in, const T& intensity, const T& thresh = T(1) / T(2));
  Mat  normalize(const Mat& data, const T& upper) const;
  void normalize(Mat data[3], const T& upper) const;
  Mat  autoLevel(const Mat& data, const int& count = 0);
  void autoLevel(Mat data[3], const int& count = 0);
  void getTileVec(const Mat& in, vector<Vec>& geoms, vector<Veci>& delaunay) const;
  match_t<T> tiltprep(const Mat& in, const int& idx, const int& samples, const T& psi, const int& origin = 0) const;
  vector<Triangles> tiltprep(const vector<Vec>& points, const vector<Veci>& polys, const Mat& in, const match_t<T>& m) const;
  Mat  tilt(const Mat& in, const Mat& bump, const match_t<T>& m, const T& depth = - T(1000)) const;
  Mat  tilt(const Mat& in, vector<Triangles>& triangles, const T& depth = - T(1000)) const;
  Mat  tilt(const Mat& in, const vector<Triangles>& triangles, const T& depth = - T(1000)) const;
  Mat  applyTrace(const pair<Vec, Vec>& v, const pair<pair<pair<int, int>, pair<int, int> >, pair<pair<int, int>, pair<int, int> > >& hw);

private:
  void drawMatchLine(Mat& map, const Vec& lref0, const Vec& lref1, const T& c) const;
  void drawMatchTriangle(Mat& map, Vec lref0, Vec lref1, Vec lref2, const T& c) const;
  inline Triangles makeTriangle(const int& u, const int& v, const Mat& in, const Mat& bump, const int& flg) const;
  void prepTrace(pair<Vec, Vec>& v, pair<pair<int, int>, pair<int, int> >& hw, const Mat& mask);
  
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

template <typename T> typename reDig<T>::Mat reDig<T>::draw(const Mat& img, const vector<Vec>& shape, const vector<Vec>& emph, const vector<Veci>& hull) {
  assert(shape.size() == emph.size());
  vector<Triangles> tris;
  tris.reserve(hull.size());
  for(int i = 0; i < hull.size(); i ++) {
    assert(hull[i].size() == 3);
    assert(0 <= hull[i][0] && hull[i][0] < shape.size());
    assert(0 <= hull[i][1] && hull[i][1] < shape.size());
    assert(0 <= hull[i][2] && hull[i][2] < shape.size());
    Triangles work;
    for(int j = 0; j < 3; j ++)
      work.p.setCol(j, emph[hull[i][j]]);
    work.c = img(max(int(0), min(int(img.rows() - 1),
                   int(shape[hull[i][0]][0]))),
                 max(int(0), min(int(img.cols() - 1),
                   int(shape[hull[i][0]][1]))));
    tris.emplace_back(work.solveN());
  }
  return tilt(img * T(0), tris);
}

template <typename T> typename reDig<T>::Mat reDig<T>::draw(const Mat& img, const vector<Vec>& shape, const vector<Veci>& hull, const bool& elim) {
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
                              tsrc[hull[ii][2]],
                              (tsrc[hull[ii][0]][2] + tsrc[hull[ii][1]][2] +
                               tsrc[hull[ii][2]][2]) / T(3));
  return result;
}

template <typename T> vector<typename reDig<T>::Vec> reDig<T>::takeShape(const vector<Vec>& dst, const vector<Vec>& src, const match_t<T>& match, const T& ratio) {
  vector<Vec> result(dst);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < match.srcpoints.size(); i ++)
    result[match.dstpoints[i]] += (match.transform(src[match.srcpoints[i]]) - dst[match.dstpoints[i]]) * ratio;
  return result;
}

template <typename T> void reDig<T>::drawMatchLine(Mat& map, const Vec& lref0, const Vec& lref1, const T& c) const {
  int idxm(0);
  int idxM(1);
  if(abs(lref1[idxM] - lref0[idxM]) < abs(lref1[idxm] - lref0[idxm])) {
    idxm = 1;
    idxM = 0;
  }
  const auto d10(lref1 - lref0);
  const auto dlt(abs(lref0[idxM] - lref1[idxM]));
  const auto denom(T(1) / dlt);
  for(int i = 0; i <= int(dlt); i ++) {
    const auto gidx(lref0 + d10 * T(i) * denom);
    map(max(int(0), min(int(gidx[0]), int(map.rows() - 1))),
        max(int(0), min(int(gidx[1]), int(map.cols() - 1)))) = c;
  }
  return;
}

template <typename T> void reDig<T>::drawMatchTriangle(Mat& map, Vec lref0, Vec lref1, Vec lref2, const T& c) const {
  int idxm(0);
  int idxM(1);
  if(abs(lref1[idxM] - lref0[idxM]) < abs(lref1[idxm] - lref0[idxm])) {
    idxm = 1;
    idxM = 0;
  }
  lref0[2] = lref1[2] = lref2[2] = T(0);
  const auto ldiff0(lref1 - lref0);
        auto ldiff(lref2 - lref0);
  ldiff -= ldiff0 * ldiff.dot(ldiff0) / ldiff0.dot(ldiff0);
  const auto lnum((sqrt(ldiff.dot(ldiff)) + T(1)) * T(2));
  // XXX : tan theta depend loop num, this have glitches.
  const auto d20(lref2 - lref0);
  const auto d21(lref2 - lref1);
  for(int k = 0; k < int(lnum); k ++)
    drawMatchLine(map, lref0 + d20 * T(k) / lnum,
                       lref1 + d21 * T(k) / lnum, c);
  return;
}

template <typename T> typename reDig<T>::Mat reDig<T>::showMatch(const Mat& dstimg, const vector<Vec>& dst, const vector<Veci>& hull, const T& emph) {
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
  result.O();
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

template <typename T> vector<typename reDig<T>::Veci> reDig<T>::mesh2(const vector<Vec>& p, const vector<int>& pp) const {
  vector<pair<Vec, int> > sp;
  vector<pair<Vec, int> > sp2;
  sp.reserve(pp.size());
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
  sp.emplace_back(make_pair(Vec(3), p.size()));
  sp.emplace_back(make_pair(Vec(3), p.size() + 1));
  sp.emplace_back(make_pair(Vec(3), p.size() + 2));
  sp.emplace_back(make_pair(Vec(3), p.size() + 3));
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
  sort(sp.begin(), sp.end(), less0<pair<Vec, int> >);
  sp2.reserve(sp.size());
  for(int i = 0; i < sp.size(); i ++)
    sp2.emplace_back(sp[sp.size() - i - 1]);
  vector<Veci> res0;
  for(int i = 2; i < sp.size(); i ++) {
    Veci lres(3);
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
    swap(lres[0], lres[2]);
    res0.emplace_back(lres);
  }
  vector<Veci> res;
  for(int i = 0; i < res0.size(); i ++)
    if(sp[res0[i][0]].second < p.size() &&
       sp[res0[i][1]].second < p.size() &&
       sp[res0[i][2]].second < p.size()) {
      for(int j = 0; j < res0[i].size(); j ++)
        res0[i][j] = sp[res0[i][j]].second;
      for(int j = 0; j < res0[i].size(); j ++)
        for(int k = j + 1; k < res0[i].size(); k ++)
          if(res0[i][j] == res0[i][k])
            goto nofix;
      for(int j = 0; j < res.size(); j ++)
        if(res[j] == res0[i])
          goto nofix;
      res.emplace_back(res0[i]);
     nofix:
      ;
    }
  for(int i = 0; i < res.size(); i ++) {
    Mat d(3, 3);
    for(int j = 0; j < 3; j ++) {
      d(j, 0) = T(1);
      d(j, 1) = p[res[i][j]][0];
      d(j, 2) = p[res[i][j]][1];
    }
    const auto det(d.determinant());
    assert(det != T(0));
    if(det < T(0))
      swap(res[i][0], res[i][1]);
  }
  return res;
}

template <typename T> vector<int> reDig<T>::edge(const vector<Vec>& p, const vector<int>& pp) const {
  vector<pair<Vec, int> > sp;
  sp.reserve(pp.size());
  for(int i = 0; i < pp.size(); i ++)
    sp.emplace_back(make_pair(p[pp[i]], pp[i]));
  sort(sp.begin(), sp.end(), less0<pair<Vec, int> >);
  vector<int> resl, resr;
  int i;
  for(i = 0; i < sp.size(); i ++) {
    const auto y(sp[i].first[0]);
    int jj = i;
    for(int j = i; j < sp.size() && y == sp[j].first[0]; j ++)
      if(sp[j].first[1] <= sp[jj].first[1])
        jj = j;
    resl.emplace_back(sp[jj].second);
    jj = i;
    for(int j = i; j < sp.size() && y == sp[j].first[0]; j ++)
      if(sp[j].first[1] >= sp[jj].first[1])
        jj = j;
    resr.emplace_back(sp[jj].second);
    i = jj + 1;
  }
  vector<int> res;
  res.reserve(resl.size() + resr.size());
  for(int i = 0; i < resl.size(); i ++)
    res.emplace_back(std::move(resl[i]));
  for(int i = 0; i < resr.size(); i ++)
    res.emplace_back(std::move(resr[resr.size() - i - 1]));
  return res;
}

template <typename T> void reDig<T>::maskVectors(vector<Vec>& points, const vector<Veci>& polys, const Mat& mask) {
  vector<Veci> tpoly(polys);
  return maskVectors(points, tpoly, mask);
}

template <typename T> void reDig<T>::maskVectors(vector<Vec>& points, vector<Veci>& polys, const Mat& mask) {
  vector<int> elim, elimp, after;
  for(int i = 0, ii = 0; i < points.size(); i ++) {
    const int y(max(min(int(points[i][0]), int(mask.rows() - 1)), int(0)));
    const int x(max(min(int(points[i][1]), int(mask.cols() - 1)), int(0)));
    if(mask(y, x) > T(1) / T(2)) {
      elim.emplace_back(i);
      after.emplace_back(- 1);
    } else
      after.emplace_back(ii ++);
  }
  for(int i = 0; i < polys.size(); i ++)
    if(binary_search(elim.begin(), elim.end(), polys[i][0]) ||
       binary_search(elim.begin(), elim.end(), polys[i][1]) ||
       binary_search(elim.begin(), elim.end(), polys[i][2]))
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

template <typename T> typename reDig<T>::Mat reDig<T>::reShape(const Mat& cbase, const Mat& vbase, const int& count, const T& thresh) {
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
  Mat masker(res);
  Mat mask(res.rows(), res.cols());
  for(int i = 0; i < mask.rows(); i ++)
    for(int j = 0; j < mask.cols(); j ++)
      mask(i, j) = false;
  for(int i = 0; i < mask.rows(); i ++)
    for(int j = 0; j < mask.cols(); j ++)
      if(! mask(i, j)) {
        vector<pair<int, int> > store;
        vector<pair<int, int> > tries;
        tries.emplace_back(make_pair(+ 1,   0));
        tries.emplace_back(make_pair(  0, + 1));
        tries.emplace_back(make_pair(- 1,   0));
        tries.emplace_back(make_pair(  0, - 1));
        vector<pair<int, int> > stack;
        stack.emplace_back(make_pair(i, j));
        T   avg(0);
        int cnt(0);
        while(stack.size()) {
          const auto pop(stack[stack.size() - 1]);
          stack.pop_back();
          const int& yy(pop.first);
          const int& xx(pop.second);
          if(0 <= yy && yy < mask.rows() && 0 <= xx && xx < mask.cols() &&
             masker(i, j) == masker(yy, xx) &&
             abs(cbase(i, j) - cbase(yy, xx)) <= thresh &&
             ! mask(yy, xx)) {
            mask(yy, xx) = true;
            store.emplace_back(make_pair(yy, xx));
            for(int ii = 0; ii < tries.size(); ii ++)
              stack.emplace_back(make_pair(yy + tries[ii].first, xx + tries[ii].second));
            avg += cbase(yy, xx);
            cnt ++;
          }
        }
        avg /= T(cnt);
        for(int i = 0; i < store.size(); i ++)
          res(store[i].first, store[i].second) = avg;
      }
  return res;
}

template <typename T> typename reDig<T>::Mat reDig<T>::reColor(const Mat& cbase, const Mat& vbase, const int& count, const T& intensity) {
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
  SimpleVector<T> vv(vpoints.size());
  SimpleVector<T> cc(cpoints.size());
  for(int i = 0; i < vv.size(); i ++)
    vv[i] = vpoints[i].first;
  for(int i = 0; i < cc.size(); i ++)
    cc[i] = cpoints[i].first;
  const auto ccc(Decompose<T>(count).mimic(cc, vv, intensity));
  Mat res(cbase);
  for(int i = 0; i < ccc.size(); i ++)
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
  SimpleVector<T> cc(cpoints.size());
  for(int i = 0; i < cc.size(); i ++)
    cc[i] = cpoints[i].first;
  const auto ccc(Decompose<T>(count).emphasis(cc, intensity));
  Mat res(cbase);
  for(int i = 0; i < ccc.size(); i ++)
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
  vector<Vec>  pdst;
  vector<Veci> facets;
  getTileVec(mask, pdst, facets);
  const auto idsts(floodfill(mask, pdst));
  assert(idsts.size());
  int iidst(0);
  for(int i = 1; i < idsts.size(); i ++)
    if(idsts[iidst].size() < idsts[i].size()) iidst = i;
  const auto idst(edge(pdst, idsts[iidst]));
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
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < idst.size(); i ++) {
    v.first[ i] = pdst[idst[i]][0];
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
  Mat res(yy, xx);
  res.O();
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < vy.size(); i ++) {
    Vec v0(3);
    Vec v1(3);
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

template <typename T> typename reDig<T>::Mat reDig<T>::optImage(const vector<pair<Mat, Mat> >& img, const int& comp) const {
  assert(0 <= comp && 0 < img.size());
  const auto ss0(img[0].first.rows() * img[0].first.cols());
  const auto ss1(img[0].second.rows() * img[0].second.cols());
        Mat  tayl0(ss0 + comp, ss0);
        Mat  tayl1(ss1 + comp, ss1);
  for(int i = 0; i < tayl0.rows(); i ++)
    tayl0.row(i) = taylor<T>(ss0, T(i) / T(tayl0.rows() - 1) * T(ss0 - 1));
  for(int i = 0; i < tayl1.rows(); i ++)
    tayl1.row(i) = taylor<T>(ss1, T(i) / T(tayl1.rows() - 1) * T(ss1 - 1));
  vector<Vec> serialize;
  serialize.reserve(img.size());
  for(int i = 0; i < img.size(); i ++) {
    assert(img[i].first.rows() == img[0].first.rows());
    assert(img[i].first.cols() == img[0].first.cols());
    assert(img[i].second.rows() == img[0].second.rows());
    assert(img[i].second.cols() == img[0].second.cols());
    Vec lseri0(tayl0.cols());
    Vec lseri1(tayl1.cols());
    for(int j = 0; j < img[i].first.rows(); j ++)
      for(int k = 0; k < img[i].first.cols(); k ++)
        lseri0[j * img[i].first.cols() + k] = img[i].first(j, k);
    for(int j = 0; j < img[i].second.rows(); j ++)
      for(int k = 0; k < img[i].second.cols(); k ++)
        lseri1[j * img[i].second.cols() + k] = img[i].second(j, k);
    serialize.emplace_back(makeProgramInvariant<T>(Vec(tayl0.rows() + tayl1.rows()).O().setVector(0, tayl0 * lseri0).setVector(tayl0.rows(), tayl1 * lseri1)));
  }
  Mat res(serialize[0].size(), serialize[0].size());
  for(int i = 0; i < res.rows(); i ++) {
    res.row(i) = linearInvariant(serialize);
    if(res.row(i).dot(res.row(i)) <= T(1) / T(10000)) {
      vector<int> idx;
      idx.reserve(res.rows() - i + 1);
      for(int j = i; j < res.rows(); j ++) {
        idx.emplace_back(j);
        for(int k = 0; k < res.cols(); k ++)
          res(j, k) = T(0);
      }
      return res.fillP(idx);
    }
    res.row(i) /= sqrt(res.row(i).dot(res.row(i)));
    for(int j = 0; j < serialize.size(); j ++)
      serialize[j] -= res.row(i) * serialize[j].dot(res.row(i));
  }
  return res;
}

template <typename T> typename reDig<T>::Mat reDig<T>::compImage(const Mat& in, const Mat& opt, const int& comp) const {
  const auto ss(in.rows() * in.cols());
        Mat  tayl(ss + comp, ss);
        Mat  taylr((in.rows() + 1) * (in.cols() + 1), opt.cols() - tayl.rows() - 1);
  for(int i = 0; i < tayl.rows(); i ++)
    tayl.row(i) = taylor<T>(ss, T(i) / T(tayl.rows() - 1) * T(ss - 1));
  for(int i = 0; i < taylr.rows(); i ++)
    taylr.row(i) = taylor<T>(taylr.cols(), T(i) / T(taylr.rows() - 1) * T(taylr.cols() - 1));
  Vec work(ss);
  for(int j = 0; j < in.rows(); j ++)
    for(int k = 0; k < in.cols(); k ++)
      work[j * in.cols() + k] = in(j, k);
  work = makeProgramInvariant<T>(tayl * work);
  work = opt.subMatrix(0, 0, opt.cols(), opt.cols()).solve(opt.subMatrix(0, opt.cols()- work.size(), opt.cols(), work.size()) * work).subVector(0, taylr.cols());
  for(int i = 0; i < work.size(); i ++)
    work[i] = atan(work[i]);
  work = taylr * work;
  Mat res(sqrt(work.size()), sqrt(work.size()));
  for(int j = 0; j < res.rows(); j ++)
    for(int k = 0; k < res.cols(); k ++)
      res(j, k) = work[(j * res.cols() + k) % work.size()];
  return res;
}

template <typename T> vector<typename reDig<T>::Mat> reDig<T>::catImage(const vector<Mat>& rep, const vector<Mat>& imgs, const int& cs) {
  assert(imgs.size() && rep.size() == imgs.size());
  for(int i = 1; i < rep.size(); i ++) {
    assert(rep[i].cols() == rep[0].cols());
    assert(imgs[i].rows() == imgs[0].rows());
    assert(imgs[i].cols() == imgs[0].cols());
  }
  vector<Vec> work;
  vector<int> workidx;
  vector<int> workidx2;
  work.reserve(imgs.size());
  for(int i = 0; i < rep.size(); i ++)
    for(int j = 0; j < min(rep[i].rows(), int(1 + 5 + 1)); j ++) {
      work.emplace_back(rep[i].row(j));
      workidx.emplace_back(i);
      workidx2.emplace_back(j);
    }
  const auto cg(crush<T>(work, cs, false, 0));
  vector<Mat> res;
  res.reserve(cg.size());
  for(int i = 0; i < cg.size(); i ++) {
    if(! cg[i].first.size()) continue;
    res.emplace_back(Mat(cg[i].first.size() * imgs[0].rows(), imgs[0].cols()));
    for(int j = 0; j < cg[i].first.size(); j ++) {
      for(int k = 0; k < imgs[0].rows(); k ++)
        res[i].row(j * imgs[0].rows() + k) = imgs[workidx[cg[i].second[j].second]].row(k);
      const auto widx2(workidx2[cg[i].second[j].second]);
      if(!widx2)
        for(int k = 0; k < imgs[0].rows(); k ++)
          res[i].row(j * imgs[0].rows() + k) = - res[i].row(j * imgs[0].rows() + k);
      else if(widx2 < 1 + 5 + 1)
        for(int k = 0; k < imgs[0].rows() / 2; k ++)
          for(int kk = 0; kk < imgs[0].cols() / 2; kk ++) {
            auto& rr(res[i](j * imgs[0].rows() + k +
              (widx2 == 2 || widx2 == 4 ? imgs[0].rows() / 2 :
                (widx2 == 5 ? imgs[0].rows() / 4 : 0)), kk +
              (widx2 == 3 || widx2 == 4 ? imgs[0].cols() / 2 :
                (widx2 == 5 ? imgs[0].cols() / 4 : 0)) ) );
            rr = - rr;
          }
      else
        assert(0 && "Should not be reached.");
    }
  }
  return res;
}

template <typename T> vector<typename reDig<T>::Mat> reDig<T>::compositeImage(const vector<Mat>& imgs) {
  assert(imgs.size());
  vector<SimpleVector<T> > work;
  for(int i = 0; i < imgs.size(); i ++) {
    assert(imgs[i].rows() == imgs[0].rows());
    assert(imgs[i].cols() == imgs[0].cols());
    work.emplace_back(SimpleVector<T>(imgs[i].rows() * imgs[i].cols()));
    for(int ii = 0; ii < imgs[i].rows(); ii ++)
      for(int jj = 0; jj < imgs[i].cols(); jj ++)
        work[i][ii * imgs[i].cols() + jj] =
          imgs[i](ii, ii & 1 ? imgs[i].cols() - 1 - jj : jj);
    assert(imgs[i].rows() * imgs[i].cols() == work[i].size());
  }
  const auto cg(crush<T>(work, work[0].size(), false, 0));
  vector<Mat> res;
  res.reserve(cg.size());
  for(int i = 0; i < cg.size(); i ++) {
    res.emplace_back(SimpleMatrix<T>(imgs[0].rows(), imgs[0].cols()));
    auto work(cg[i].first[0]);
    for(int j = 1; j < cg[i].first.size(); j ++)
      work += cg[i].first[j];
    assert(res[i].rows() * res[i].cols() == work.size());
    for(int ii = 0; ii < res[i].rows(); ii ++)
      for(int jj = 0; jj < res[i].cols(); jj ++)
        res[i](ii, ii & 1 ? res[i].cols() - 1 - jj : jj) =
          atan(work[ii * res[i].cols() + jj]);
  }
  return res;
}

template <typename T> typename reDig<T>::Mat reDig<T>::bump(const Mat& color, const Mat& bumpm, const T& psi, const int& n, const int& origin) const {
  assert(color.rows() == bumpm.rows() && color.cols() == bumpm.cols());
  if(n == 0) {
    if(origin < 0) {
      auto res(bump(color, bumpm, psi, n, 0));
      for(int i = 1; i < 5; i ++)
        res += bump(color, bumpm, psi, n, i);
      return res /= T(5);
    }
    assert(0 <= origin && origin < 5);
    const auto color0(tilt(color, bumpm, tiltprep(bumpm, 1, 2, - abs(psi), origin)));
    const auto color1(tilt(color, bumpm, tiltprep(bumpm, 1, 2,   abs(psi), origin)));
          Mat  result(color.rows(), color.cols());
          auto zscore(result);
    result.O();
    zscore.I(- T(1));
    const T   rxy(T(color0.rows() - 1) / T(2));
    const int dratio(sqrt(sqrt(rxy)));
          Vec camera(2);
    camera[0] = T(0);
    camera[1] = T(1);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int j = 0; j < result.rows(); j ++) {
      for(int zi = 0; zi < dratio; zi ++) {
        Vec cpoint(2);
        cpoint[0] = T(1) / T(2 * dratio);
        cpoint[1] = T(zi) / T(dratio);
        const auto t(- camera[1] / (cpoint[1] - camera[1]));
        const auto x0((camera + (cpoint - camera) * t)[0] * rxy);
        if(int(x0) < 3 * vbox || rxy < abs(x0) * T(2)) continue;
        Vec work(result.cols());
        work.O();
        SimpleVector<T> Dop0;
#if defined(_OPENMP)
#pragma omp critical
#endif
        {
          const auto Dop(diff<T>(abs(int(x0 * T(2))) & ~ int(1)));
          Dop0 = Dop.row(Dop.rows() / 2) + Dop.row(Dop.rows() / 2 + 1);
        }
        for(int k = 0; k < Dop0.size(); k ++) {
          // N.B. projection scale is linear.
          T xorigin(0);
          switch(origin) {
          case 0:
            cpoint[0] = (T(j + k) - T(Dop0.size() + result.rows() - 2) / T(2)) / T(2 * dratio) / rxy;
            xorigin   = T(result.rows() - 1) / T(2);
            break;
          case 1:
          case 3:
            cpoint[0] = (T(j + k) - T(Dop0.size() - 1) / T(2)) / T(2 * dratio) / rxy;
            break;
          case 2:
          case 4:
            cpoint[0] = (T(j + k) - T(Dop0.size() - 1) / T(2) - T(result.rows() - 1)) / T(2 * dratio) / rxy;
            xorigin   = T(result.rows() - 1);
            break;
          default:
            assert(0 && "Should not be reached in redig::bump");
          }
          // x-z plane projection of point p with camera geometry c to z=0.
          // c := camera, p := cpoint.
          // <c + (p - c) * t, [0, 1]> = 0
          const auto t(- camera[1] / (cpoint[1] - camera[1]));
          const auto x0(getImgPt<int>(int((camera + (cpoint - camera) * t)[0] * rxy + xorigin), result.rows()));
          work += (k < Dop0.size() / 2 ? color0.row(x0) : color1.row(x0)) * Dop0[k];
        }
        for(int i = 0; i < work.size(); i ++)
          if(zscore(j, i) < abs(work[i])) {
            result(j, i) = T(zi + 1);
            zscore(j, i) = abs(work[i]);
          }
      }
    }
    return result;
  }
  static const auto Pi(atan2(T(1), T(1)) * T(4));
  Mat res(color * T(0));
  const auto ct(color.transpose());
  const auto bt(bumpm.transpose());
  for(int i = 0; i < n; i ++) {
    const auto theta((T(i) - T(n - 1) / T(2)) * atan(T(1)) / (T(n) / T(2)));
    res += center<T>(rotate<T>(bump(rotate<T>(color, theta),
             rotate<T>(bumpm, theta), psi), - theta), color);
    res += center<T>(rotate<T>(bump(rotate<T>(ct, theta),
             rotate<T>(bt, theta), psi), - theta).transpose(), color);
  }
  return res /= T(n * 2);
}

template <typename T> vector<vector<int> > reDig<T>::floodfill(const Mat& mask, const vector<Vec>& points) {
  vector<vector<int> > result;
  SimpleMatrix<bool> checked(mask.rows(), mask.cols());
  for(int i = 0; i < checked.rows(); i ++)
    for(int j = 0; j < checked.cols(); j ++)
      checked(i, j) = false;
  for(int i = 0; i < points.size(); i ++) {
    const auto& pi(points[i]);
    const int   y0(pi[0]);
    const int   x0(pi[1]);
    if(0 <= y0 && y0 < mask.rows() &&
       0 <= x0 && x0 < mask.cols() &&
       mask(y0, x0) < T(1) / T(2)) {
      vector<pair<T, T> > tries;
      tries.emplace_back(make_pair(  T(vbox),         0));
      tries.emplace_back(make_pair(        0,   T(vbox)));
      tries.emplace_back(make_pair(- T(vbox),         0));
      tries.emplace_back(make_pair(        0, - T(vbox)));
      vector<pair<int, int> > stack;
      vector<int> store;
      stack.emplace_back(make_pair(i, i));
      while(stack.size()) {
        const auto pop(stack[stack.size() - 1]);
        stack.pop_back();
        if(pop.first < 0 || points.size() <= pop.first)
          store.emplace_back(pop.second);
        const int& yy(points[pop.first][0]);
        const int& xx(points[pop.first][1]);
        if(! (0 <= yy && yy < checked.rows() &&
              0 <= xx && xx < checked.cols() &&
              mask(yy, xx) < T(1) / T(2) &&
              !checked(yy, xx) ) )
          store.emplace_back(pop.first);
        else if(!checked(yy, xx)) {
          checked(yy, xx) = true;
          for(int ii = 0; ii < tries.size(); ii ++) {
            const auto ny(T(yy) + tries[ii].first);
            const auto nx(T(xx) + tries[ii].second);
            for(int k = 0; k < points.size(); k ++)
              if(abs(points[k][0] - ny) < T(1) / T(2) &&
                 abs(points[k][1] - nx) < T(1) / T(2)) {
                stack.emplace_back(make_pair(k, pop.first));
                break;
              }
          }
        }
      }
      if(! store.size()) continue;
      result.emplace_back(std::move(store));
    }
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
  Mat mRGB2XYZ(3, 3);
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
  Mat mRGB2XYZ(3, 3);
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
  Mat mXYZ2RGB(3, 3);
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
      result(i, j) = min(abs(thresh) + T(.5), max(- abs(thresh) + T(.5), intensity * (result(i, j) - T(.5)) )) + T(.5);
  return result;
}

template <typename T> typename reDig<T>::Mat reDig<T>::normalize(const Mat& data, const T& upper) const {
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

template <typename T> void reDig<T>::normalize(Mat data[3], const T& upper) const {
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
template <typename T> void reDig<T>::getTileVec(const Mat& in, vector<Vec>& geoms, vector<Veci>& delaunay) const {
  // get vectorize.
  geoms = vector<Vec>();
  T aavg(0);
  for(int i = 0; i < in.rows(); i ++)
    for(int j = 0; j < in.cols(); j ++)
      aavg += in(i, j);
  aavg /= T(in.rows() * in.cols());
  for(int i = 0; i < in.rows() / vbox + 1; i ++)
    for(int j = 0; j < in.cols() / vbox + 1; j ++) {
      if(in.rows() < (i + 1) * vbox ||
         in.cols() < (j + 1) * vbox) {
        Vec gbuf(3);
        gbuf[0] = T(i * vbox);
        gbuf[1] = T(j * vbox);
        gbuf[2] = geoms[geoms.size() - 1][2];
        geoms.emplace_back(gbuf);
      } else {
        T avg(0);
        for(int ii = i * vbox; ii < (i + 1) * vbox; ii ++)
          for(int jj = j * vbox; jj < (j + 1) * vbox; jj ++)
            avg += in(ii, jj);
        Vec work(3);
        work[0] = T(i * vbox);
        work[1] = T(j * vbox);
        work[2] = sqrt(T(in.rows() * in.cols())) * rz * abs(avg / T(vbox) / T(vbox) - aavg);
        geoms.emplace_back(work);
      }
    }
  Vec avg(3);
  avg[0] = avg[1] = avg[2] = T(0);
  for(int i = 0; i < geoms.size(); i ++)
    avg += geoms[i];
  avg /= geoms.size();
  for(int i = 0; i < geoms.size(); i ++)
    geoms[i][2] -= avg[2];
  delaunay = vector<Veci>();
  for(int i = 1; i < in.rows() / vbox + 1; i ++)
    for(int j = 0; j < in.cols() / vbox; j ++) {
      Veci work(3), work2(3);
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

template <typename T> inline typename reDig<T>::Triangles reDig<T>::makeTriangle(const int& u, const int& v, const Mat& in, const Mat& bump, const int& flg) const {
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

template <typename T> match_t<T> reDig<T>::tiltprep(const Mat& in, const int& idx, const int& samples, const T& psi, const int& origin) const {
  const T theta(T(2) * Pi * T(idx) / T(samples));
  const T lpsi(Pi * psi);
  Mat R0(3, 3);
  Mat R1(3, 3);
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
  Vec pcenter(3);
  switch(origin) {
  case 0:
    pcenter[0] = T(in.rows() - 1) / T(2);
    pcenter[1] = T(in.cols() - 1) / T(2);
    break;
  case 1:
    pcenter[0] = pcenter[1] = T(0);
    break;
  case 2:
    pcenter[0] = T(in.rows() - 1);
    pcenter[1] = T(0);
    break;
  case 3:
    pcenter[0] = T(0);
    pcenter[1] = T(in.cols() - 1);
    break;
  case 4:
    pcenter[0] = T(in.rows() - 1);
    pcenter[1] = T(in.cols() - 1);
    break;
  default:
    assert(0 && "Should not be reached in tiltprep.");
  }
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

template <typename T> vector<typename reDig<T>::Triangles> reDig<T>::tiltprep(const vector<Vec>& points, const vector<Veci>& polys, const Mat& in, const match_t<T>& m) const {
  vector<Triangles> result;
  for(int i = 0; i < polys.size(); i ++) {
    Triangles work;
    for(int j = 0; j < 3; j ++)
      work.p.setCol(j, m.transform(points[polys[i][j]]));
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

template <typename T> typename reDig<T>::Mat reDig<T>::tilt(const Mat& in, const Mat& bump, const match_t<T>& m, const T& depth) const {
  assert(in.rows() == bump.rows() && in.cols() == bump.cols());
  vector<Vec>  points;
  vector<Veci> facets;
  getTileVec(bump, points, facets);
  return tilt(in, tiltprep(points, facets, in, m), depth);
}

template <typename T> typename reDig<T>::Mat reDig<T>::tilt(const Mat& in, const vector<Triangles>& triangles, const T& depth) const {
  auto tris(triangles);
  return tilt(in, tris, depth);
}

template <typename T> typename reDig<T>::Mat reDig<T>::tilt(const Mat& in, vector<Triangles>& triangles, const T& depth) const {
  cerr << "t" << flush;
  Mat result(in.rows(), in.cols());
  Vec vz(3);
  result.O();
  vz.ek(2);
  // XXX: patent???
  vector<pair<T, Triangles> > zbuf;
  zbuf.resize(triangles.size(), make_pair(T(0), Triangles()));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int j = 0; j < triangles.size(); j ++) {
          auto& tri(triangles[j]);
    const auto  p0(tri.p.col(0));
    const auto  p1(tri.p.col(1));
    const auto  p2(tri.p.col(2));
          auto  camera((p0 + p1 + p2) / T(3));
    camera[2] = T(0);
    const auto t((tri.z - tri.n.dot(camera)) / (tri.n.dot(vz)));
    zbuf[j].first  = camera[2] + vz[2] * t;
    zbuf[j].second = move(tri);
  }
  sort(zbuf.begin(), zbuf.end(), lessf<pair<T, Triangles> >);
  int i;
  for(i = 0; i < zbuf.size() && zbuf[i].first < depth; i ++) ;
  for( ; i < zbuf.size(); i ++) {
    const auto& zbi(zbuf[i].second);
    drawMatchTriangle(result, zbi.p.col(0), zbi.p.col(1), zbi.p.col(2), zbi.c);
  }
  return result;
}

#define _REDIG_
#endif

