#if !defined(_SCAN_CONTEXT_)

#include <Eigen/Core>
#include <Eigen/LU>
#include <cmath>
#include <vector>
#include "tilt.hh"

using std::sqrt;
using std::atan2;
using std::abs;
using std::cos;
using std::acos;
using std::sin;
using std::pow;
using std::sort;
using std::unique;
using std::min;
using std::max;
using std::floor;
using std::cerr;
using std::endl;
using std::fflush;
using std::vector;

template <typename T> class lfmatch_t {
public:
  Eigen::Matrix<T, 3, 1> pt;
  T                      score;
  lfmatch_t<T>& operator = (const lfmatch_t<T>& other) {
    pt    = other.pt;
    score = other.score;
    return *this;
  }
  bool operator < (const lfmatch_t<T>& x1) const {
    return score > x1.score;
  }
};

template <typename T> class lowFreq {
public:
  typedef complex<T> U;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T, 3, 1>                           Vec3;
  lowFreq();
  ~lowFreq();
  void init(const T& zr);
  
  vector<Vec3> getLowFreq(const Mat& data, const int& npoints = 60);
  Mat getLowFreqImage(const Mat& data, const int& npoints = 60);
private:
  vector<lfmatch_t<T> > prepareCost(const Mat& data, const int& npoints);
  T zr;
  T Pi;
  U I;
};

template <typename T> lowFreq<T>::lowFreq() {
  Pi          = atan2(T(1), T(1)) * T(4);
  I           = sqrt(U(- 1));
  init(60);
}

template <typename T> lowFreq<T>::~lowFreq() {
  ;
}

template <typename T> void lowFreq<T>::init(const T& zr) {
  this->zr = zr;
  return;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> lowFreq<T>::getLowFreqImage(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& data, const int& npoints) {
  Mat result(data);
  vector<lfmatch_t<T> > match(prepareCost(data, npoints));
  result *= T(0);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int k = 0; k < match.size(); k ++) {
    const Eigen::Matrix<T, 3, 1>& pt(match[k].pt);
    result(int(pt[0]), int(pt[1])) = pt[2];
  }
  return result;
}

template <typename T> vector<Eigen::Matrix<T, 3, 1> > lowFreq<T>::getLowFreq(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& data, const int& npoints) {
  vector<Eigen::Matrix<T, 3, 1> > result;
  vector<lfmatch_t<T> > match(prepareCost(data, npoints));
  for(int k = 0; k < match.size(); k ++) {
    const Eigen::Matrix<T, 3, 1>& pt(match[k].pt);
    result.push_back(pt);
  }
  return result;
}

template <typename T> vector<lfmatch_t<T> > lowFreq<T>::prepareCost(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& data, const int& npoints) {
  vector<lfmatch_t<T> > match;
  const T guard(T(npoints) / T(data.rows() * data.cols()));
  for(int i = 0; i < data.rows() * sqrt(guard); i ++)
    for(int j = 0; j < data.cols() * sqrt(guard); j ++) {
      lfmatch_t<T> m;
      Eigen::Matrix<T, 3, 1>& pt(m.pt);
      T csum(0);
      pt[0] = pt[1] = pt[2] = T(0);
      int count(0);
      for(int ii = i / sqrt(guard); ii < min(T(data.rows()), (i + 1) / sqrt(guard)); ii ++)
        for(int jj = j / sqrt(guard); jj < min(T(data.cols()), (j + 1) / sqrt(guard)); jj ++) {
          csum  += data(ii, jj);
          pt[0] += ii;
          pt[1] += jj;
          pt[2] += data(ii, jj);
          count ++;
        }
      m.score = csum / count;
      pt    /= count;
      pt[2] *= zr;
      match.push_back(m);
    }
  sort(match.begin(), match.end());
  return match;
  // XXX add me lowpoly:
}

template <typename T> class match_t {
public:
  Eigen::Matrix<T, 3, 3> rot;
  Eigen::Matrix<T, 3, 1> offset;
  T                      ratio;
  T                      rdepth;
  T                      rpoints;
  vector<int>            dstpoints;
  vector<int>            srcpoints;
  T                      threshs;
  T                      threshc;
  match_t() {
    threshs = T(0);
    threshc = T(0);
  }
  match_t(const T& threshs, const T& threshc) {
    this->threshs = threshs;
    this->threshc = threshc;
  }
  match_t(const match_t<T>& other) {
    *this = other;
    return;
  }
  match_t<T>& operator = (const match_t<T>& other) {
    rot       = other.rot;
    offset    = other.offset;
    ratio     = other.ratio;
    rdepth    = other.rdepth;
    rpoints   = other.rpoints;
    dstpoints = other.dstpoints;
    srcpoints = other.srcpoints;
    threshs   = other.threshs;
    threshc   = other.threshc;
    return *this;
  }
  // XXX configure us:
  bool operator < (const match_t<T>& x1) const {
    return rpoints > x1.rpoints || (rpoints == x1.rpoints && rdepth < x1.rdepth) || (rpoints == x1.rpoints && rdepth == x1.rdepth && ratio > x1.ratio);
    // return rdepth < x1.rdepth || (rdepth == x1.rdepth && rpoints > x1.rpoints) || (rdepth == x1.rdepth && rpoints == x1.rpoints && ratio > x1.ratio);
  }
  bool operator == (const match_t<T>& x) const {
    const auto test(offset - x.offset);
    const auto roterr(rot * x.rot.transpose());
    const T    Pi(T(4) * atan2(T(1), T(1)));
          T    rotdnorm(0);
    for(int l = 0; l < rot.rows(); l ++)
      if(!(abs(acos(roterr(l, l))) / Pi <= threshs))
        return false;
    if(!(sqrt(test.dot(test) / 
              (offset.dot(offset) * x.offset.dot(x.offset)))
         <= threshc))
      return false;
    if(!(max(T(1) - ratio   / x.ratio,
             T(1) - x.ratio / ratio)
         <= threshc))
      return false;
    return true;
  }
};

template <typename T> class msub_t {
public:
  int mbufj;
  int mbufk;
  T   mbufN;
  msub_t() {
    ;
  }
  msub_t(const msub_t<T>& other) {
    *this = other;
    return;
  }
  msub_t<T>& operator = (const msub_t<T>& other) {
    mbufj = other.mbufj;
    mbufk = other.mbufk;
    mbufN = other.mbufN;
    return *this;
  }
  bool operator < (const msub_t<T>& x1) const {
    return mbufN > x1.mbufN || (mbufN == x1.mbufN && mbufj < x1.mbufj) || (mbufN == x1.mbufN && mbufj == x1.mbufj && mbufk < x1.mbufk);
  }
};

template <typename T> class matchPartialPartial {
public:
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T, 3, 3>                           Mat3x3;
  typedef Eigen::Matrix<T, 2, 2>                           Mat2x2;
  typedef Eigen::Matrix<T, 3, 1>              Vec3;
  typedef complex<T> U;
  matchPartialPartial();
  ~matchPartialPartial();
  void init(const int& ndiv, const T& thresh, const T& threshp, const T& threshr, const T& threshs);
  
  vector<match_t<T> > match(const vector<Vec3>& shapebase, const vector<Vec3>& points);
  void match(const vector<Vec3>& shapebase, const vector<Vec3>& points, vector<match_t<T> >& result);
private:
  U   I;
  T   Pi;
  int ndiv;
  T   thresh;
  T   threshp;
  T   threshr;
  T   threshs;
};

template <typename T> matchPartialPartial<T>::matchPartialPartial() {
  I  = sqrt(U(- T(1)));
  Pi = atan2(T(1), T(1)) * T(4);
  // rough match.
  init(16, .25, .125, .125, .0625);
}

template <typename T> matchPartialPartial<T>::~matchPartialPartial() {
  ;
}

template <typename T> void matchPartialPartial<T>::init(const int& ndiv, const T& thresh, const T& threshp, const T& threshr, const T& threshs) {
  this->ndiv        = ndiv;
  this->thresh      = thresh;
  this->threshp     = threshp;
  this->threshr     = threshr;
  this->threshs     = threshs;
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
  Vec3 gs, gp;
  gs[0] = gs[1] = gs[2] = T(0);
  gp[0] = gp[1] = gp[2] = T(0);
  for(int i = 0; i < shapebase.size(); i ++)
    gs += shapebase[i];
  gs /= shapebase.size();
  for(int i = 0; i < points.size(); i ++)
    gp += points[i];
  gp /= points.size();
  for(int i = 0; i < 3; i ++) {
    Eigen::Matrix<Eigen::Matrix<T, 2, 1>, Eigen::Dynamic, Eigen::Dynamic> table(shapebase.size(), points.size());
    // init table.
    cerr << "making table (" << i << "/3)" << endl;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int j = 0; j < shapebase.size(); j ++) {
      fflush(stderr);
      for(int k = 0; k < points.size(); k ++) {
        Vec3 aj(shapebase[j] - gs);
        Vec3 bk(points[k]    - gp);
        aj /= sqrt(aj.dot(aj));
        bk /= sqrt(bk.dot(bk));
        for(int l = 0; l < table(j, k).size(); l ++) {
          const T a(bk[(l + i    ) % 3]);
          const T b(bk[(l + i + 1) % 3]);
          const T c(aj[(l + i    ) % 3]);
          table(j, k)[l] = sqrt(- T(1));
          if(a * a + b * b < c * c)
            continue;
          T theta(   T(2) * atan2(sqrt(a * a + b * b - c * c) - b, a + c));
          T theta1(- T(2) * atan2(sqrt(a * a + b * b - c * c) + b, a + c));
          if(!isfinite(theta) || (isfinite(theta1) && 
               abs(cos(theta1) * a - sin(theta1) * b - c) <
               abs(cos(theta)  * a - sin(theta)  * b - c)))
            theta = theta1;
          theta -= floor(theta / (T(2) * Pi) + T(.5)) * T(2) * Pi;
          table(j, k)[l] = theta;
          bk[(l + i    ) % 3] = cos(theta) * a - sin(theta) * b;
          bk[(l + i + 1) % 3] = sin(theta) * a + cos(theta) * b;
        }
      }
    }
    // matches.
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int nd = 0; nd < ndiv; nd ++) {
      cerr << "matching table (" << nd << "/" << ndiv << ")";
      Eigen::Matrix<T, 2, 1> ddiv;
      ddiv[0] = cos(2. * Pi * nd / T(ndiv));
      ddiv[1] = sin(2. * Pi * nd / T(ndiv));
      vector<msub_t<T> > msub;
      for(int k = 0; k < points.size(); k ++)
        for(int j = 0; j < shapebase.size(); j ++) {
          const T lnorm(sqrt(table(j, k).dot(table(j, k))));
          const T ldepth(min(cos(table(j, k)[0] / lnorm) * cos(ddiv[0]) +
                             sin(table(j, k)[0] / lnorm) * sin(ddiv[0]),
                             cos(table(j, k)[1] / lnorm) * cos(ddiv[1]) +
                             sin(table(j, k)[1] / lnorm) * sin(ddiv[1])) - T(1));
          if(isfinite(lnorm) && (lnorm <= Pi / ndiv ||
              (isfinite(ldepth) && abs(ldepth) <= T(1) / ndiv) ) ) {
            msub_t<T> workm;
            workm.mbufj = j;
            workm.mbufk = k;
            workm.mbufN = lnorm;
            msub.push_back(workm);
          }
        }
      cerr << " : " << msub.size() << endl;
      if(!msub.size())
        continue;
      sort(msub.begin(), msub.end());
      for(int k0 = 0; k0 < msub.size(); k0 ++) {
        match_t<T> work(threshs, threshs);
        work.rot = drot0;
        for(int k = 0; k < ddiv.size(); k ++) {
          Mat3x3 lrot;
          T theta(ddiv[k] * msub[k0].mbufN);
          lrot((k + i    ) % 3, (k + i    ) % 3) =   cos(theta);
          lrot((k + i + 1) % 3, (k + i    ) % 3) =   sin(theta);
          lrot((k + i    ) % 3, (k + i + 1) % 3) = - sin(theta);
          lrot((k + i + 1) % 3, (k + i + 1) % 3) =   cos(theta);
          lrot((k + i + 2) % 3, (k + i    ) % 3) = T(0);
          lrot((k + i + 2) % 3, (k + i + 1) % 3) = T(0);
          lrot((k + i + 2) % 3, (k + i + 2) % 3) = T(1);
          lrot((k + i    ) % 3, (k + i + 2) % 3) = T(0);
          lrot((k + i + 1) % 3, (k + i + 2) % 3) = T(0);
          work.rot = lrot * work.rot;
        }
        bool flagk[points.size()];
        bool flagj[shapebase.size()];
        for(int kk = 0; kk < points.size(); kk ++)
          flagk[kk] = false;
        for(int kk = 0; kk < shapebase.size(); kk ++)
          flagj[kk] = false;
        const Vec3 ajk0(shapebase[msub[k0].mbufj]          - gs);
        const Vec3 bkk0(work.rot * (points[msub[k0].mbufk] - gp));
        const T    t0(ajk0.dot(bkk0) / bkk0.dot(bkk0));
        const Vec3 lerr0(ajk0 - bkk0 * t0);
        if(!isfinite(t0) || thresh * thresh < lerr0.dot(lerr0) / (ajk0.dot(ajk0) + bkk0.dot(bkk0) * t0 * t0))
          continue;
        for(int kk = k0 + 1; kk < msub.size(); kk ++)
          if(!flagj[msub[kk].mbufj] &&
             !flagk[msub[kk].mbufk]) {
            const Vec3 aj(shapebase[msub[kk].mbufj]          - gs);
            const Vec3 bk(work.rot * (points[msub[kk].mbufk] - gp));
            const Vec3 lerr(aj - bk * t0);
            if(thresh * thresh < lerr.dot(lerr) / (aj.dot(aj) + bk.dot(bk) * t0 * t0))
              break;
            work.dstpoints.push_back(msub[kk].mbufj);
            work.srcpoints.push_back(msub[kk].mbufk);
            flagj[msub[kk].mbufj] = true;
            flagk[msub[kk].mbufk] = true;
          }
        work.rpoints = work.dstpoints.size() /
                         T(min(shapebase.size(), points.size()));
        if(threshp <= work.rpoints) {
          Vec3 sbar, pbar;
          sbar[0] = sbar[1] = sbar[2] = T(0);
          pbar[0] = pbar[1] = pbar[2] = T(0);
          for(int k = 0; k < work.dstpoints.size(); k ++) {
            sbar += shapebase[work.dstpoints[k]]          - gs;
            pbar += work.rot * (points[work.srcpoints[k]] - gp);
          }
          sbar /= work.dstpoints.size();
          pbar /= work.srcpoints.size();
          T num(0), denom(0);
          for(int k = 0; k < work.dstpoints.size(); k ++) {
            num   += (shapebase[work.dstpoints[k]] - gs - sbar).dot(work.rot * (points[work.srcpoints[k]] - gp) - pbar);
            denom += (work.rot * (points[work.srcpoints[k]] - gp) - pbar).dot(work.rot * (points[work.srcpoints[k]] - gp) - pbar);
          }
          work.ratio = num / denom;
          work.offset[0] = T(0);
          work.offset[1] = T(0);
          work.offset[2] = T(0);
          for(int k = 0; k < work.dstpoints.size(); k ++)
            work.offset += shapebase[work.dstpoints[k]] - work.rot * points[work.srcpoints[k]] * work.ratio;
          work.offset /= work.dstpoints.size();
          if(abs(work.ratio) < threshr || T(1) / threshr < abs(work.ratio))
            continue;
          work.rdepth = T(0);
          for(int k = 0; k < work.dstpoints.size(); k ++) {
            const Vec3& aj(shapebase[msub[k].mbufj]);
            const Vec3  bk(work.rot * points[msub[k].mbufk] * work.ratio + work.offset);
            work.rdepth += (aj - bk).dot(aj - bk) / (aj.dot(aj) + bk.dot(bk));
          }
          work.rdepth /= work.dstpoints.size();
          if(work.rdepth < thresh) {
            cerr << ".";
            fflush(stderr);
#if defined(_OPENMP)
#pragma omp critical
#endif
            {
              int  idx(- 1);
              bool flag(false);
              for(int k = 0; k < result.size(); k ++)
                if(result[k] == work) {
                  flag = true;
                  if(work < result[k]) {
                    idx = k;
                    break;
                  }
               }
              if(flag) {
                if(idx >= 0)
                  result[idx] = work;
              } else {
                result.push_back(work);
                cerr << "*";
                fflush(stderr);
              }
              sort(result.begin(), result.end());
            }
          }
        }
      }
    }
  }
  sort(result.begin(), result.end());
  return;
}


template <typename T> class matchWholePartial {
public:
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T, 3, 3>                           Mat3x3;
  typedef Eigen::Matrix<T, 3, 1>              Vec3;
  typedef complex<T> U;
  matchWholePartial();
  ~matchWholePartial();
  void init(const vector<vector<Vec3> >& shapebase, const T& thresh, const T& threshp);

  vector<match_t<T> > match(const vector<Vec3>& points);
private:
  U I;
  T Pi;
  T thresh;
  vector<matchPartialPartial<T> > matches;
  vector<vector<Vec3> >           shapebase;
};

template <typename T> matchWholePartial<T>::matchWholePartial() {
  I  = sqrt(U(- T(1)));
  Pi = atan2(T(1), T(1)) * T(4);
}

template <typename T> matchWholePartial<T>::~matchWholePartial() {
  ;
}

template <typename T> void matchWholePartial<T>::init(const vector<vector<Vec3> >& shapebase, const T& thresh, const T& threshp) {
  for(int i = 0; i < shapebase.size(); i ++) {
    matchPartialPartial<T> work;
    work.init(shapebase[i]);
    matches.push_back(work);
  }
  this->thresh  = thresh;
  this->threshp = threshp;
  return;
}

template <typename T> vector<match_t<T> > matchWholePartial<T>::match(const vector<Vec3>& points) {
  cerr << "matching partial polys..." << endl;
  vector<vector<match_t<T> > > pmatches;
  for(int i = 0; i < matches.size(); i ++)
    pmatches.push_back(matches[i].match(points));
  cerr << "detecting possible whole matches..." << endl;
  vector<match_t<T> > result;
  return result;
}


template <typename T> class tilter;
template <typename T> class reDig {
public:
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T, 3, 3>                           Mat3x3;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
  typedef Eigen::Matrix<T, 3, 1>              Vec3;
  typedef Eigen::Matrix<int, 3, 1>            Veci3;
  
  reDig();
  ~reDig();
  void init();
  Vec3 emphasis0(const Vec3& dst, const Vec3& refdst, const Vec3& src, const match_t<T>& match, const T& ratio);
  Mat  emphasis(const Mat& dstimg, const Mat& dstbump, const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const vector<Veci3>& hull, const T& ratio);
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
  const Vec3 a(refdst);
  const Vec3 b(match.rot * src * match.ratio + match.offset);
  return dst + (b - a) * (exp(ratio) - exp(T(1))) / exp(T(1));
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> reDig<T>::emphasis(const Mat& dstimg, const Mat& dstbump, const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const vector<Eigen::Matrix<int, 3, 1> >& hull, const T& ratio) {
  cerr << " making triangles";
  fflush(stderr);
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
  
  for(int ii = 0; ii < hull.size(); ii ++) {
    cerr << "emphasis: " << ii << "/" << hull.size() << endl;
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
        triangles[l].col(4) = tilt.solveN(triangles[l].col(0), triangles[l].col(1), triangles[l].col(2));
      }
      triangles[l](1, 3)  = triangles[l].col(4).dot(triangles[l].col(0));
    }
  }
  
  Mat I3(3, 3);
  Vec zero3(3);
  for(int i = 0; i < 3; i ++) {
    for(int j = 0; j < 3; j ++)
      I3(i, j) = (i == j ? T(1) : T(0));
    zero3[i] = T(0);
  }
  delete[] checked;
  return tilt.tiltsub(dstimg, triangles, I3, I3, zero3, zero3, T(1));
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
           map.rows() <= lref0[1] + tilt * ii);
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
           map.cols() <= lref0[0] + tilt * jj);
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
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result(input), map(input.rows(), input.cols());
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
  return result + map;
}

#define _SCAN_CONTEXT_
#endif

