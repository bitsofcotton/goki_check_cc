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

template <typename T> void drawMatchLine(SimpleMatrix<T>& map, const SimpleVector<T>& lref0, const SimpleVector<T>& lref1, const T& c) {
  int idxm(0);
  int idxM(1);
  if(abs(lref1[idxM] - lref0[idxM]) < abs(lref1[idxm] - lref0[idxm])) {
    idxm = 1;
    idxM = 0;
  }
  const auto d10(lref1 - lref0);
  const auto dlt(abs(lref0[idxM] - lref1[idxM]));
  if(dlt == T(0)) return;
  const auto denom(T(1) / dlt);
  for(int i = 0; i <= int(dlt); i ++) {
    const auto gidx(lref0 + d10 * T(i) * denom);
    map(max(int(0), min(int(gidx[0]), int(map.rows() - 1))),
        max(int(0), min(int(gidx[1]), int(map.cols() - 1)))) = c;
  }
  return;
}

template <typename T> void drawMatchTriangle(SimpleMatrix<T>& map, SimpleVector<T> lref0, SimpleVector<T> lref1, SimpleVector<T> lref2, const T& c) {
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
    drawMatchLine<T>(map, lref0 + d20 * T(k) / lnum,
                          lref1 + d21 * T(k) / lnum, c);
  return;
}

template <typename T> void addMeshTri(vector<SimpleVector<int> >& res, vector<pair<SimpleVector<T>, int> >& scan, const vector<SimpleVector<T> >& p, const int& idx) {
  assert(0 <= idx && idx < scan.size());
  vector<int> elim;
  if(0 <= idx - 1 &&
     scan[idx].first[0] < scan[idx - 1].first[0] &&
     scan[idx].first[0] < scan[idx + 1].first[0]) {
    elim.emplace_back(idx);
    SimpleVector<int> lres(3);
    lres[0] = scan[idx - 1].second;
    lres[1] = scan[idx].second;
    lres[2] = scan[idx + 1].second;
    bool psize(false);
    for(int k = 0; k < 3; k ++)
      psize = psize || p.size() <= lres[k];
    if(! psize)
      res.emplace_back(move(lres));
  }
  if(idx + 3 < scan.size() &&
     scan[idx + 2].first[0] < scan[idx + 1].first[0] &&
     scan[idx + 2].first[0] < scan[idx + 3].first[0]) {
    elim.emplace_back(idx + 2);
    SimpleVector<int> lres(3);
    lres[0] = scan[idx + 1].second;
    lres[1] = scan[idx + 2].second;
    lres[2] = scan[idx + 3].second;
    bool psize(false);
    for(int k = 0; k < 3; k ++)
      psize = psize || p.size() <= lres[k];
    if(! psize)
      res.emplace_back(move(lres));
  }
  {
    SimpleVector<int> lres(3);
    lres[0] = scan[idx].second;
    lres[1] = scan[idx + 1].second;
    lres[2] = scan[idx + 2].second;
    bool psize(false);
    for(int k = 0; k < 3; k ++)
      psize = psize || p.size() <= lres[k];
    if(! psize)
      res.emplace_back(move(lres));
  }
  sort(elim.begin(), elim.end());
  for(int j = 0; j < elim.size(); j ++)
    scan.erase(scan.begin() + elim[j] - j);
  return;
}

template <typename T> vector<SimpleVector<int> > mesh2(const vector<SimpleVector<T> >& p, const vector<int>& pp) {
  vector<pair<SimpleVector<T>, int> > sp;
  sp.reserve(pp.size());
  SimpleMatrix<T> lrot(3, 3);
  lrot.I();
  lrot(0, 0) =    lrot(1, 1) = cos(T(int(1)) / T(p.size()));
  lrot(0, 1) = - (lrot(1, 0) = sin(T(int(1)) / T(p.size())));
  T    m1((lrot * p[pp[0]])[1]);
  auto M1(m1);
  for(int i = 0; i < pp.size(); i ++) {
    sp.emplace_back(make_pair(lrot * p[pp[i]], pp[i]));
    sp[i].first[2] = T(0);
    m1 = min(m1, sp[i].first[1]);
    M1 = max(M1, sp[i].first[1]);
  }
  sort(sp.begin(), sp.end(), less0<pair<SimpleVector<T>, int> >);
  vector<pair<SimpleVector<T>, int> > scan;
  scan.reserve(sp.size() + 2);
  scan.emplace_back(sp[0]);
  scan[scan.size() - 1].first[0] -= T(2);
  scan[scan.size() - 1].first[1]  = m1 - T(1);
  scan.emplace_back(sp[0]);
  scan[scan.size() - 1].first[0] -= T(1);
  scan[scan.size() - 1].first[1]  = M1 + T(1);
  vector<SimpleVector<int> > res;
  res.reserve(sp.size());
  for(int i = 0; i < sp.size(); i ++) {
    // N.B. lrot support this on lattice.
    assert(! i || (sp[i].first[0] != sp[i - 1].first[0] &&
                   sp[i].first[1] != sp[i - 1].first[1]) );
    // scanline update
    int idx;
    for(idx = 0; idx < scan.size(); idx ++)
      if(sp[i].first[1] < scan[idx].first[1]) break;
    idx = max(0, min(int(scan.size()) - 2, idx - 1));
    assert(scan[idx].first[1] < sp[i].first[1]);
    assert(sp[i].first[1] < scan[idx + 1].first[1]);
    assert(scan[idx].first[0] < sp[i].first[0]);
    assert(scan[idx + 1].first[0] < sp[i].first[0]);
    scan.insert(scan.begin() + idx + 1, pair<SimpleVector<T>, int>(sp[i]));
    assert(scan[idx].first[1] < scan[idx + 1].first[1] &&
           scan[idx + 1].first[1] < scan[idx + 2].first[1]);
    addMeshTri<T>(res, scan, p, idx);
  }
  while(6 < scan.size()) {
    const auto before(scan.size());
    for(int i = 2; i < scan.size() - 3; i ++)
      addMeshTri<T>(res, scan, p, i);
    if(before == scan.size()) break;
  }
  res.reserve(res.size());
  for(int i = 0; i < res.size(); i ++)
    if(p[res[i][0]][0] * p[res[i][1]][1]
     + p[res[i][1]][0] * p[res[i][2]][1]
     + p[res[i][2]][0] * p[res[i][0]][1]
     - p[res[i][0]][1] * p[res[i][1]][0]
     - p[res[i][1]][1] * p[res[i][2]][0]
     - p[res[i][2]][1] * p[res[i][0]][0] < T(0))
      swap(res[i][0], res[i][1]);
  return res;
}

template <typename T> vector<SimpleVector<int> > mesh2(const vector<SimpleVector<T> >& p) {
  vector<int> pp;
  pp.reserve(p.size());
  for(int i = 0; i < p.size(); i ++) pp.emplace_back(i);
  return mesh2<T>(p, pp);
}

// get bump with multiple scale and vectorized result.
template <typename T> vector<SimpleVector<T> > getTileVec(const SimpleMatrix<T>& in, const int& vbox = 3) {
  vector<SimpleVector<T> > geoms;
  T aavg(0);
  for(int i = 0; i < in.rows(); i ++)
    for(int j = 0; j < in.cols(); j ++)
      aavg += in(i, j);
  aavg /= T(in.rows() * in.cols());
  for(int i = 0; i < in.rows() / vbox + 1; i ++)
    for(int j = 0; j < in.cols() / vbox + 1; j ++) {
      if(in.rows() < (i + 1) * vbox ||
         in.cols() < (j + 1) * vbox) {
        SimpleVector<T> gbuf(3);
        gbuf[0] = T(i * vbox);
        gbuf[1] = T(j * vbox);
        gbuf[2] = geoms[geoms.size() - 1][2];
        geoms.emplace_back(gbuf);
      } else {
        T avg(0);
        for(int ii = i * vbox; ii < (i + 1) * vbox; ii ++)
          for(int jj = j * vbox; jj < (j + 1) * vbox; jj ++)
            avg += in(ii, jj);
        SimpleVector<T> work(3);
        work[0] = T(i * vbox);
        work[1] = T(j * vbox);
        work[2] = T(min(in.rows(), in.cols())) * abs(avg / T(vbox) / T(vbox) - aavg);
        geoms.emplace_back(work);
      }
    }
  SimpleVector<T> avg(3);
  avg.O();
  for(int i = 0; i < geoms.size(); i ++)
    avg += geoms[i];
  avg /= geoms.size();
  for(int i = 0; i < geoms.size(); i ++)
    geoms[i][2] -= avg[2];
  return geoms;
}

template <typename T> vector<SimpleVector<T> > getHesseVec(const SimpleMatrix<T>& in, const int& vbox = 300) {
  const auto guard(max(1, int(sqrt(T(in.rows() * in.cols() / vbox)))));
  vector<SimpleVector<T> > geoms;
  geoms.reserve(vbox + 4);
  const auto xx(in * diff<T>(in.cols()).transpose() * diff<T>(in.cols()).transpose());
  const auto xy(diff<T>(in.rows()) * in * diff<T>(in.cols()).transpose());
  const auto yy(diff<T>(in.rows()) * diff<T>(in.rows()) * in);
  vector<pair<T, pair<int, int> > > score;
  score.reserve(in.rows() * in.cols());
  for(int i = 0; i < in.rows(); i ++)
    for(int j = 0; j < in.cols(); j ++)
      score.emplace_back(make_pair(abs(xx(i, j) * yy(i, j) - xy(i, j) * xy(i, j)), make_pair(i, j)));
  sort(score.begin(), score.end());
  vector<pair<int, int> > cache;
  cache.reserve(score.size());
  for(int i = score.size() - 1; 0 <= i && geoms.size() < abs(vbox); i --)
    if(! binary_search(cache.begin(), cache.end(),
           make_pair(score[i].second.first / guard,
                     score[i].second.second / guard)) ) {
      SimpleVector<T> g(3);
      g[0] = T(int(score[i].second.first));
      g[1] = T(int(score[i].second.second));
      g[2] = T(min(in.rows(), in.cols())) *
        in(score[i].second.first, score[i].second.second);
      geoms.emplace_back(move(g));
      cache.emplace_back(make_pair(score[i].second.first / guard,
                                   score[i].second.second / guard));
      sort(cache.begin(), cache.end());
    }
  SimpleVector<T> g(3);
  g[0] = T(int(0));
  g[1] = T(int(0));
  g[2] = T(min(in.rows(), in.cols())) * in(0, 0);
  geoms.emplace_back(g);
  g[0] = T(int(in.rows() - 1));
  g[1] = T(int(0));
  g[2] = T(min(in.rows(), in.cols())) * in(in.rows() - 1, 0);
  geoms.emplace_back(g);
  g[0] = T(int(0));
  g[1] = T(int(in.cols() - 1));
  g[2] = T(min(in.rows(), in.cols())) * in(0, in.cols() - 1);
  geoms.emplace_back(g);
  g[0] = T(int(in.rows() - 1));
  g[1] = T(int(in.cols() - 1));
  g[2] = T(min(in.rows(), in.cols())) * in(in.rows() - 1, in.cols() - 1);
  geoms.emplace_back(g);
  return geoms;
}

template <typename T> SimpleMatrix<T> tilt(const SimpleMatrix<T>& in, vector<triangles_t<T> >& triangles, const T& depth = - T(10000)) {
  cerr << "t" << flush;
  SimpleMatrix<T> result(in.rows(), in.cols());
  SimpleVector<T> vz(3);
  result.O();
  vz.ek(2);
  // XXX: patent???
  vector<pair<T, triangles_t<T>> > zbuf;
  zbuf.resize(triangles.size(), make_pair(T(0), triangles_t<T>()));
  assert(zbuf.size() == triangles.size());
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
    zbuf[j].first  = - (camera[2] + vz[2] * t);
    if(! isfinite(zbuf[j].first)) zbuf[j].first = depth;
    zbuf[j].second = move(tri);
  }
  sort(zbuf.begin(), zbuf.end(), lessf<pair<T, triangles_t<T>> >);
  int i;
  for(i = 0; i < zbuf.size() && zbuf[i].first < depth; i ++) ;
  for( ; i < zbuf.size(); i ++) {
    const auto& zbi(zbuf[i].second);
    drawMatchTriangle<T>(result, zbi.p.col(0), zbi.p.col(1), zbi.p.col(2), zbi.c);
  }
  return result;
}

template <typename T> SimpleMatrix<T> tilt(const SimpleMatrix<T>& in, const vector<triangles_t<T> >& triangles, const T& depth = - T(10000)) {
  auto tris(triangles);
  return tilt<T>(in, tris, depth);
}

template <typename T> SimpleMatrix<T> tilt(const SimpleMatrix<T>& in, const SimpleMatrix<T>& bump, const match_t<T>& m, const T& depth = - T(10000)) {
  assert(in.rows() == bump.rows() && in.cols() == bump.cols());
  auto points(getTileVec<T>(bump));
  auto facets(mesh2<T>(points));
  vector<triangles_t<T> > triangles;
  triangles.reserve(facets.size());
  for(int i = 0; i < facets.size(); i ++) {
    triangles_t<T> work;
    for(int j = 0; j < 3; j ++) {
      assert(0 <= facets[i][j] && facets[i][j] < points.size());
      work.p.setCol(j, m.transform(points[facets[i][j]]));
    }
    if(T(0) <= points[facets[i][0]][0] && points[facets[i][0]][0] < T(in.rows()) &&
       T(0) <= points[facets[i][0]][1] && points[facets[i][0]][1] < T(in.cols()))
      work.c = in(int(points[facets[i][0]][0]),
                  int(points[facets[i][0]][1]));
    else
      work.c = T(0);
    triangles.emplace_back(work.solveN());
  }
  return tilt<T>(in, triangles, depth);
}

template <typename T> SimpleMatrix<T> draw(const SimpleMatrix<T>& img, const vector<SimpleVector<T> >& shape, const vector<SimpleVector<T> >& emph, const vector<SimpleVector<int> >& hull) {
  assert(shape.size() == emph.size());
  vector<triangles_t<T> > tris;
  tris.reserve(hull.size());
  for(int i = 0; i < hull.size(); i ++) {
    assert(hull[i].size() == 3);
    assert(0 <= hull[i][0] && hull[i][0] < shape.size());
    assert(0 <= hull[i][1] && hull[i][1] < shape.size());
    assert(0 <= hull[i][2] && hull[i][2] < shape.size());
    triangles_t<T> work;
    for(int j = 0; j < 3; j ++)
      work.p.setCol(j, emph[hull[i][j]]);
    work.c = img(max(int(0), min(int(img.rows() - 1),
                   int(shape[hull[i][0]][0]))),
                 max(int(0), min(int(img.cols() - 1),
                   int(shape[hull[i][0]][1]))));
    tris.emplace_back(work.solveN());
  }
  return tilt<T>(img * T(0), tris);
}

template <typename T> SimpleMatrix<T> draw(const SimpleMatrix<T>& img, const vector<SimpleVector<T> >& shape, const vector<SimpleVector<int> >& hull, const bool& elim = false) {
  auto result(img);
  T    M(0);
  T    m(0);
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
    drawMatchTriangle<T>(result, tsrc[hull[ii][0]],
                                 tsrc[hull[ii][1]],
                                 tsrc[hull[ii][2]],
                                 (tsrc[hull[ii][0]][2] + tsrc[hull[ii][1]][2] +
                                  tsrc[hull[ii][2]][2]) / T(3));
  return result;
}

template <typename T> vector<SimpleVector<T> > takeShape(const vector<SimpleVector<T> >& dst, const vector<SimpleVector<T> >& src, const match_t<T>& match, const T& ratio) {
  auto result(dst);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < match.src.size(); i ++)
    result[match.dst[i]] += (match.transform(src[match.src[i]]) - dst[match.dst[i]]) * ratio;
  return result;
}

template <typename T> SimpleMatrix<T> showMatch(const SimpleMatrix<T>& dstimg, const vector<SimpleVector<T> >& dst, const vector<SimpleVector<int> >& hull, const T& emph = T(1)) {
  auto map(dstimg);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int k = 0; k < hull.size(); k ++) {
    drawMatchLine<T>(map, dst[hull[k][0]], dst[hull[k][1]], emph);
    drawMatchLine<T>(map, dst[hull[k][1]], dst[hull[k][2]], emph);
    drawMatchLine<T>(map, dst[hull[k][2]], dst[hull[k][0]], emph);
  }
  return map;
}

template <typename T> SimpleMatrix<T> makeRefMatrix(const SimpleMatrix<T>& orig, const int& start) {
  SimpleMatrix<T> result(orig.rows(), orig.cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < orig.rows() * orig.cols(); i ++)
    result(i % orig.rows(), i / orig.rows()) = i + start;
  return result;
}

template <typename T> SimpleMatrix<T> pullRefMatrix(const SimpleMatrix<T>& ref, const int& start, const SimpleMatrix<T>& orig) {
  assert(orig.rows() == ref.rows() && orig.cols() == ref.cols());
  SimpleMatrix<T> result(ref.rows(), ref.cols());
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

template <typename T> vector<int> edge(const vector<SimpleVector<T> >& p, const vector<int>& pp) {
  vector<pair<SimpleVector<T>, int> > sp;
  sp.reserve(pp.size());
  for(int i = 0; i < pp.size(); i ++)
    sp.emplace_back(make_pair(p[pp[i]], pp[i]));
  sort(sp.begin(), sp.end(), less0<pair<SimpleVector<T>, int> >);
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

template <typename T> SimpleMatrix<T> reShape(const SimpleMatrix<T>& cbase, const SimpleMatrix<T>& vbase, const int& count, const T& thresh) {
  assert(cbase.rows() && cbase.cols() && vbase.rows() && vbase.cols());
  assert(cbase.rows() == vbase.rows() && cbase.cols() == vbase.cols());
  vector<pair<T, pair<int, int> > > vpoints;
  vpoints.reserve(vbase.rows() * vbase.cols());
  for(int i = 0; i < vbase.rows(); i ++)
    for(int j = 0; j < vbase.cols(); j ++)
      vpoints.emplace_back(make_pair(vbase(i, j), make_pair(i, j)));
  sort(vpoints.begin(), vpoints.end());
  SimpleMatrix<T> res(vbase.rows(), vbase.cols());
  T avg(0);
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
  SimpleMatrix<T> masker(res);
  SimpleMatrix<T> mask(res.rows(), res.cols());
  mask.O(false);
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

template <typename T> SimpleMatrix<T> reColor(const SimpleMatrix<T>& cbase, const SimpleMatrix<T>& vbase, const int& count, const T& intensity) {
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
        auto res(cbase);
  for(int i = 0; i < ccc.size(); i ++)
    res(cpoints[i].second.first, cpoints[i].second.second) = ccc[i];
  return res;
}

template <typename T> SimpleMatrix<T> reColor3(const SimpleMatrix<T>& cbase, const SimpleMatrix<T>& vbase, const int& count) {
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
  SimpleMatrix<T> res(cbase.rows(), cbase.cols());
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

template <typename T> SimpleMatrix<T> reColor(const SimpleMatrix<T>& cbase, const int& count, const T& intensity) {
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
        auto res(cbase);
  for(int i = 0; i < ccc.size(); i ++)
    res(cpoints[i].second.first, cpoints[i].second.second) = ccc[i];
  return res;
}

template <typename T> vector<vector<int> > floodfill(const SimpleMatrix<T>& mask, const vector<SimpleVector<T> >& points, const int& vbox = 3) {
  vector<vector<int> > result;
  SimpleMatrix<bool> checked(mask.rows(), mask.cols());
  checked.O(false);
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

template <typename T> void prepTrace(pair<SimpleVector<T>, SimpleVector<T> >& v, pair<pair<int, int>, pair<int, int> >& hw, const SimpleMatrix<T>& mask, const int& vbox = 3) {
        auto pdst(getTileVec<T>(mask, vbox));
  const auto idsts(floodfill<T>(mask, pdst, vbox));
  assert(idsts.size());
  int iidst(0);
  for(int i = 1; i < idsts.size(); i ++)
    if(idsts[iidst].size() < idsts[i].size()) iidst = i;
  const auto idst(edge<T>(pdst, idsts[iidst]));
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

template <typename T> SimpleMatrix<T> applyTrace(const pair<SimpleVector<T>, SimpleVector<T> >& v, const pair<pair<pair<int, int>, pair<int, int> >, pair<pair<int, int>, pair<int, int> > >& hw) {
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
  SimpleMatrix<T> res(yy, xx);
  res.O();
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < vy.size(); i ++) {
    SimpleVector<T> v0(3);
    SimpleVector<T> v1(3);
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
    drawMatchLine<T>(res, v0, v1, T(1));
  }
  return res;
}

template <typename T> SimpleMatrix<T> reTrace(const SimpleMatrix<T>& dst, const SimpleMatrix<T>& src, const T& intensity, const int& count) {
  pair<SimpleVector<T>, SimpleVector<T> > pdst, psrc;
  pair<pair<int, int>, pair<int, int> > dsthw, srchw;
  prepTrace<T>(pdst, dsthw, dst);
  prepTrace<T>(psrc, srchw, src);
  Decompose<T> decom(count);
  return applyTrace<T>(make_pair(
    decom.mimic(pdst.first,  psrc.first,  intensity),
    decom.mimic(pdst.second, psrc.second, intensity)), make_pair(dsthw, srchw));
}

template <typename T> SimpleMatrix<T> reTrace(const SimpleMatrix<T>& dst, const T& intensity, const int& count) {
  pair<SimpleVector<T>, SimpleVector<T> > pdst;
  pair<pair<int, int>, pair<int, int> > hw;
  prepTrace<T>(pdst, hw, dst);
  Decompose<T> decom(count);
  return applyTrace<T>(make_pair(
    decom.emphasis(pdst.first,  intensity),
    decom.emphasis(pdst.second, intensity)), make_pair(hw, hw));
}

template <typename T> SimpleMatrix<T> reImage(const SimpleMatrix<T>& dst, const SimpleMatrix<T>& src, const T& intensity, const int& count) {
  Decompose<T> decom(count);
  assert(dst.rows() == src.rows() && dst.cols() == src.cols());
  SimpleMatrix<T> res(dst.rows(), dst.cols());
  for(int i = 0; i < res.rows(); i ++)
    res.row(i) = decom.mimic(dst.row(i), src.row(i), intensity);
  return res;
}

template <typename T> SimpleMatrix<T> reImage(const SimpleMatrix<T>& dst, const T& intensity, const int& count) {
  Decompose<T> decom(count);
  SimpleMatrix<T> res(dst.rows(), dst.cols());
  for(int i = 0; i < res.rows(); i ++)
    res.row(i) = decom.emphasis(dst.row(i), intensity);
  return res;
}

template <typename T> SimpleMatrix<T> optImage(const vector<pair<SimpleMatrix<T>, SimpleMatrix<T> > >& img) {
  assert(0 < img.size());
  const auto ss0(img[0].first.rows() * img[0].first.cols());
  const auto ss1(img[0].second.rows() * img[0].second.cols());
  SimpleMatrix<T> serialize(img.size(), ss0 + ss1 + 2);
  for(int i = 0; i < img.size(); i ++) {
    assert(img[i].first.rows() == img[0].first.rows());
    assert(img[i].first.cols() == img[0].first.cols());
    assert(img[i].second.rows() == img[0].second.rows());
    assert(img[i].second.cols() == img[0].second.cols());
    SimpleVector<T> buf(ss0 + ss1);
    for(int j = 0; j < img[i].first.rows(); j ++)
      buf.setVector(j * img[i].first.cols(), img[i].first.row(j));
    for(int j = 0; j < img[i].second.rows(); j ++)
      buf.setVector(ss0 + j * img[i].second.cols(), img[i].second.row(j));
    serialize.row(i) = makeProgramInvariant<T>(buf).first;
  }
  SimpleMatrix<T> res(serialize.rows() - 1, serialize.cols());
  for(int i = 0; i < res.rows(); i ++) {
    res.row(i) = linearInvariant(serialize);
    if(res.row(i).dot(res.row(i)) <= T(1) / T(10000)) {
      SimpleMatrix<T> ares(i, res.cols());
      for(int j = 0; j < ares.rows(); j ++)
        ares.row(j) = move(res.row(j));
      return ares;
    }
    const auto orth(res.row(i) /= sqrt(res.row(i).dot(res.row(i))));
    for(int j = 0; j < serialize.rows(); j ++)
      serialize.row(j) -= orth * orth.dot(serialize.row(j));
  }
  return res;
}

template <typename T> SimpleMatrix<T> compImage(const SimpleMatrix<T>& in, const SimpleMatrix<T>& opt) {
  const auto ss(opt.cols() - 2);
  SimpleVector<T> work(ss);
  work.O();
  for(int j = 0; j < in.rows(); j ++)
    work.setVector(j * in.cols(), in.row(j));
  auto work2(makeProgramInvariant<T>(work));
  work  = move(work2.first);
  work -= opt.projectionPt(work);
  SimpleMatrix<T> res(sqrt(work.size() - in.rows() * in.cols()),
                      sqrt(work.size() - in.rows() * in.cols()));
  for(int j = 0; j < res.rows(); j ++)
    for(int k = 0; k < res.cols(); k ++)
      res(j, k) = revertProgramInvariant<T>(make_pair(work[(j * res.cols() + k) % work.size()], work2.second));
  return res;
}

template <typename T> vector<vector<int> > catImage(const vector<SimpleMatrix<T> >& imgs, const int& cs = 40) {
  for(int i = 1; i < imgs.size(); i ++) {
    assert(imgs[i].rows() == imgs[0].rows());
    assert(imgs[i].cols() == imgs[0].cols());
  }
  vector<SimpleVector<T> > work;
  work.reserve(imgs.size());
  for(int i = 0; i < imgs.size(); i ++) {
    work.emplace_back(SimpleVector<T>(imgs[i].rows() * imgs[i].cols()).O());
    for(int j = 0; j < imgs[i].rows(); j ++)
      work[i].setVector(imgs[i].cols() * j, imgs[i].row(j));
  }
  auto cg(crush<T>(work, cs, 0));
  vector<vector<int> > res;
  res.reserve(cg.size());
  for(int i = 0; i < cg.size(); i ++)
    res.emplace_back(std::move(cg[i].second));
  return res;
}

template <typename T> vector<SimpleMatrix<T> > compositeImage(const vector<SimpleMatrix<T> >& imgs, const int& cs = 40) {
  assert(imgs.size());
  vector<SimpleVector<T> > work;
  work.reserve(imgs.size());
  for(int i = 0; i < imgs.size(); i ++) {
    work.emplace_back(SimpleVector<T>(imgs[i].rows() * imgs[i].cols()).O());
    for(int j = 0; j < imgs[i].rows(); j ++)
      work[i].setVector(imgs[i].cols() * j, imgs[i].row(j));
  }
  auto cg(crush<T>(work, cs, 0));
  vector<SimpleMatrix<T> > res;
  res.reserve(cg.size());
  for(int i = 0; i < cg.size(); i ++) {
    res.emplace_back(SimpleMatrix<T>(imgs[0].rows(), imgs[0].cols()).O());
    auto work(cg[i].first[0]);
    for(int j = 1; j < cg[i].first.size(); j ++)
      work += cg[i].first[j];
    assert(res[i].rows() * res[i].cols() == work.size());
    for(int j = 0; j < res[i].rows(); j ++)
      res[i].row(j) = work.subVector(res[i].cols() * j, res[i].cols());
  }
  return res;
}

template <typename T> vector<SimpleMatrix<T> > rgb2xyz(const vector<SimpleMatrix<T> >& rgb) {
  // CIE 1931 XYZ from wikipedia.org
  SimpleMatrix<T> mRGB2XYZ(3, 3);
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
  assert(rgb.size() == 3);
  assert(rgb[0].rows() == rgb[1].rows() && rgb[1].rows() == rgb[2].rows());
  assert(rgb[0].cols() == rgb[1].cols() && rgb[1].cols() == rgb[2].cols());
  auto xyz(rgb);
  xyz[0] = rgb[0] * mRGB2XYZ(0, 0) + rgb[1] * mRGB2XYZ(0, 1) + rgb[2] * mRGB2XYZ(0, 2);
  xyz[1] = rgb[0] * mRGB2XYZ(1, 0) + rgb[1] * mRGB2XYZ(1, 1) + rgb[2] * mRGB2XYZ(1, 2);
  xyz[2] = rgb[0] * mRGB2XYZ(2, 0) + rgb[1] * mRGB2XYZ(2, 1) + rgb[2] * mRGB2XYZ(2, 2);
  assert(xyz.size() == 3);
  assert(xyz[0].rows() == xyz[1].rows() && xyz[1].rows() == xyz[2].rows());
  assert(xyz[0].cols() == xyz[1].cols() && xyz[1].cols() == xyz[2].cols());
  return xyz;
}

template <typename T> vector<SimpleMatrix<T> > xyz2rgb(const vector<SimpleMatrix<T> >& xyz) {
  // CIE 1931 XYZ from wikipedia.org
  SimpleMatrix<T> mRGB2XYZ(3, 3);
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
  const auto mXYZ2RGB(mRGB2XYZ.inverse());
  assert(xyz.size() == 3);
  assert(xyz[0].rows() == xyz[1].rows() && xyz[1].rows() == xyz[2].rows());
  assert(xyz[0].cols() == xyz[1].cols() && xyz[1].cols() == xyz[2].cols());
  auto rgb(xyz);
  rgb[0] = xyz[0] * mXYZ2RGB(0, 0) + xyz[1] * mXYZ2RGB(0, 1) + xyz[2] * mXYZ2RGB(0, 2);
  rgb[1] = xyz[0] * mXYZ2RGB(1, 0) + xyz[1] * mXYZ2RGB(1, 1) + xyz[2] * mXYZ2RGB(1, 2);
  rgb[2] = xyz[0] * mXYZ2RGB(2, 0) + xyz[1] * mXYZ2RGB(2, 1) + xyz[2] * mXYZ2RGB(2, 2);
  assert(rgb.size() == 3);
  assert(rgb[0].rows() == rgb[1].rows() && rgb[1].rows() == rgb[2].rows());
  assert(rgb[0].cols() == rgb[1].cols() && rgb[1].cols() == rgb[2].cols());
  return rgb;
}

template <typename T> SimpleMatrix<T> rgb2d(const vector<SimpleMatrix<T> > rgb) {
  auto xyz(rgb2xyz<T>(rgb));
  SimpleMatrix<T> result(rgb[0].rows(), rgb[0].cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int j = 0; j < rgb[0].rows(); j ++) {
    for(int k = 0; k < rgb[0].cols(); k ++)
      result(j, k) = sqrt(xyz[0](j, k) * xyz[0](j, k) + xyz[1](j, k) * xyz[1](j, k) + xyz[2](j, k) * xyz[2](j, k));
  }
  return result;
}

template <typename T> vector<SimpleMatrix<T> > contrast(const vector<SimpleMatrix<T> >& in, const T& intensity, const T& thresh) {
  auto result(in);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int k = 0; k < result.size(); k ++)
    for(int i = 0; i < result[k].rows(); i ++)
      for(int j = 0; j < result[k].cols(); j ++)
        result[k](i, j) = min(abs(thresh) + T(.5), max(- abs(thresh) + T(.5), intensity * (result[k](i, j) - T(.5)) )) + T(.5);
  return result;
}

template <typename T> SimpleMatrix<T> contrast(const SimpleMatrix<T>& in, const T& intensity, const T& thresh) {
  vector<SimpleMatrix<T> > work;
  work.emplace_back(in);
  return contrast<T>(work, intensity, thresh)[0];
}

template <typename T> vector<SimpleMatrix<T> > normalize(const vector<SimpleMatrix<T> >& data, const T& upper = T(1)) {
  T MM(0), mm(0);
  bool fixed(false);
  for(int k = 0; k < data.size(); k ++)
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
    return data;
  auto result(data);
  for(int k = 0; k < data.size(); k ++)
    for(int i = 0; i < data[k].rows(); i ++)
      for(int j = 0; j < data[k].cols(); j ++) {
        if(isfinite(result[k](i, j)) && ! isinf(data[k](i, j)) && ! isnan(result[k](i, j)))
          result[k](i, j) -= mm;
        else
          result[k](i, j)  = T(0);
        assert(T(0) <= result[k](i, j) && result[k](i, j) <= MM - mm);
        result[k](i, j) *= upper / (MM - mm);
      }
  return result;
}

template <typename T> SimpleMatrix<T> normalize(const SimpleMatrix<T>& data, const T& upper = T(1)) {
  vector<SimpleMatrix<T> > work;
  work.emplace_back(data);
  return normalize<T>(work, upper)[0];
}

template <typename T> vector<SimpleMatrix<T> > autoLevel(const vector<SimpleMatrix<T> >& data, const int& count) {
  vector<T> res;
  res.reserve(data[0].rows() * data[0].cols() * data.size());
  for(int k = 0; k < data.size(); k ++)
    for(int i = 0; i < data[k].rows(); i ++)
      for(int j = 0; j < data[k].cols(); j ++)
        res.emplace_back(data[k](i, j));
  sort(res.begin(), res.end());
  auto result(data);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int k = 0; k < data.size(); k ++)
    for(int i = 0; i < data[k].rows(); i ++)
      for(int j = 0; j < data[k].cols(); j ++)
        result[k](i, j) = max(min(data[k](i, j), res[res.size() - count - 1]), res[count]);
  return result;
}

template <typename T> SimpleMatrix<T> autoLevel(const SimpleMatrix<T>& data, const int& count) {
  vector<SimpleMatrix<T> > work;
  work.emplace_back(data);
  return autoLevel(work, count)[0];
}

template <typename T> match_t<T> tiltprep(const SimpleMatrix<T>& in, const int& idx, const int& samples, const T& psi) {
  const auto Pi(atan2(T(int(1)), T(int(1))) * T(int(4)));
  const auto theta(T(2) * Pi * T(idx) / T(samples));
  const auto lpsi(Pi * psi);
  SimpleMatrix<T> R0(3, 3);
  SimpleMatrix<T> R1(3, 3);
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
  match_t<T> m;
  m.rot    = R0.transpose() * R1 * R0;
  SimpleVector<T> pcenter(3);
  pcenter[0] = T(in.rows() - 1) / T(2);
  pcenter[1] = T(in.cols() - 1) / T(2);
  pcenter[2] = T(0);
  // x -> m.rot * x, same center
  // x - origin -> m.rot * (x - origin)
  // x -> m.rot * x - m.rot * origin + origin.
  m.offset = pcenter - m.rot * pcenter;
  m.ratio  = T(1);
  return m;
}

#define _REDIG_
#endif

