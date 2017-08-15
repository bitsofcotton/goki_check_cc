#if !defined(_SCAN_CONTEXT_)

#include <Eigen/Core>
#include <Eigen/LU>
#include <cmath>
#include <vector>
#include "edgedetect.hh"
#include "tilt.hh"

using std::sqrt;
using std::atan2;
using std::abs;
using std::cos;
using std::sin;
using std::pow;
using std::sort;
using std::min;
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
};

template <typename T> int cmplfwrap(const lfmatch_t<T>& x0, const lfmatch_t<T>& x1) {
  return x0.score < x1.score;
}

template <typename T> class lowFreq {
public:
  typedef complex<T> U;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1>              Vec;
  typedef Eigen::Matrix<T, 3, 1>                           Vec3;
  lowFreq(const T& guard = T(.125));
  ~lowFreq();
  
  vector<Eigen::Matrix<T, 3, 1> > getLowFreq(const Mat& data, const int& npoints = 600);
  Mat getLowFreqImage(const Mat& data, const int& npoints = 120);
private:
  vector<lfmatch_t<T> > prepareCost(const Mat& data, const int& npoints);
  T Pi;
  U I;
  T guard;
};

template <typename T> lowFreq<T>::lowFreq(const T& guard) {
  Pi          = atan2(T(1), T(1)) * T(4);
  I           = sqrt(U(- 1));
  this->guard = guard;
}

template <typename T> lowFreq<T>::~lowFreq() {
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> lowFreq<T>::getLowFreqImage(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& data, const int& npoints) {
  Mat result(data);
  vector<lfmatch_t<T> > match(prepareCost(data, npoints));
  result *= T(0);
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
  edgedetect<T> differ;
  Mat costs(differ.detect(data, edgedetect<T>::COLLECT_BOTH));
  for(int i = 0; i < data.rows(); i ++)
    for(int j = 0; j < data.cols(); j ++) {
      lfmatch_t<T> m;
      m.score = T(1) / costs(i, j);
      Eigen::Matrix<T, 3, 1>& pt(m.pt);
      pt[0] = i;
      pt[1] = j;
      // XXX check me:
      pt[2] = data(i, j) * sqrt(T(data.rows() * data.cols()));
      match.push_back(m);
    }
  sort(match.begin(), match.end(), cmplfwrap<T>);
  match.resize(match.size() * guard);
  vector<lfmatch_t<T> > result;
  while(true) {
    result = vector<lfmatch_t<T> >();
    for(int i = 0; i < match.size(); i ++) {
      cerr << "lpoly : " << match.size() << " : " << i << endl;
      vector<T> lmbuf;
      for(int j = 0; j < match.size(); j ++) {
        Eigen::Matrix<T, 3, 1> diff(match[i].pt - match[j].pt);
        lmbuf.push_back(diff.dot(diff));
      }
      vector<T> slmbuf(lmbuf);
      sort(slmbuf.begin(), slmbuf.end());
      if(!slmbuf.size()) break;
      Vec3 p[3];
      int  pidx[3];
      int  k0(1);
      for(int k = 0; k < p[0].size(); k ++)
        for(int j = 0; j < lmbuf.size(); j ++)
          if(slmbuf[k + k0] == lmbuf[j]) {
            bool flag = false;
            for(int kk = 0; kk < k; kk ++)
              if(j == pidx[kk])
                flag = true;
            if(flag) continue;
            pidx[k] = j;
            p[k]    = match[j].pt;
            break;
          }
      Vec3 q(p[1] - p[0]), r(p[2] - p[0]);
      r -= r.dot(q) * q / q.dot(q);
      Vec3 ek;
      for(int k = 0; k < p[0].size(); k ++) {
        for(int kk = 0; kk < p[0].size(); kk ++)
          ek[kk] = k == kk ? T(1) : T(0);
        ek -= ek.dot(q) * q / q.dot(q) + ek.dot(r) * r / r.dot(r);
        if(ek.dot(ek) > T(1) / T(9))
          break;
      }
      result.push_back(match[i]);
      result[i].score = abs(ek.dot(result[i].pt - p[0]) / sqrt(ek.dot(ek)));
    }
    sort(result.begin(), result.end(), cmplfwrap<T>);
    result.resize(result.size() * guard);
    if(result.size() == match.size() || result.size() <= npoints)
      break;
    match = result;
  }
  result.resize(min(int(result.size()), npoints));
  return result;
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
  match_t() {
    dstpoints = vector<int>();
    srcpoints = vector<int>();
  }
  ~match_t() {
    ;
  }
  match_t& operator = (const match_t& other) {
    rot       = other.rot;
    offset    = other.offset;
    ratio     = other.ratio;
    rdepth    = other.rdepth;
    rpoints   = other.rpoints;
    dstpoints = other.dstpoints;
    srcpoints = other.srcpoints;
    return *this;
  }
};

template <typename T> class msub_t {
public:
  int mbufj;
  int mbufk;
  T   mbufN;
  msub_t& operator = (const msub_t& other) {
    mbufj   = other.mbufj;
    mbufk   = other.mbufk;
    mbufN   = other.mbufN;
    return *this;
  }
};

template <typename T> int cmpwrap(const match_t<T>& x0, const match_t<T>& x1) {
  // XXX configure me:
  // return abs(x0.ratio) > abs(x1.ratio) || (x0.ratio == x1.ratio && x0.rpoints > x1.rpoints);
  // return abs(x0.ratio) > abs(x1.ratio) || (x0.ratio == x1.ratio && x0.rdepth > x1.rdepth) || (x0.ratio == x1.ratio && x0.rdepth == x1.rdepth && x0.rpoints > x1.rpoints);
  return x0.rpoints > x1.rpoints || (x0.rpoints == x1.rpoints && x0.rdepth > x1.rdepth);
}

template <typename T> int cmpsubwrap(const msub_t<T>& x0, const msub_t<T>& x1) {
  return x0.mbufN > x1.mbufN || (x0.mbufN == x1.mbufN && x0.mbufj < x1.mbufj) || (x0.mbufN == x1.mbufN && x0.mbufj == x1.mbufj && x0.mbufk < x1.mbufk);
}


template <typename T> class matchPartialPartial {
public:
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T, 3, 3>                           Mat3x3;
  typedef Eigen::Matrix<T, 2, 2>                           Mat2x2;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
  typedef Eigen::Matrix<T, 3, 1>              Vec3;
  typedef Eigen::Matrix<T, 2, 1>              Vec2;
  typedef complex<T> U;
  matchPartialPartial();
  ~matchPartialPartial();
  void init(const vector<Vec3>& shapebase, const T& thresh, const T& threshp, const T& threshr);
  
  vector<match_t<T> > match(const vector<Vec3>& points, const int& ndiv = 120, const T& r_max_theta = T(0));
private:
  U I;
  T Pi;
  // thresh, threshp in [0, 1]
  T thresh;
  T threshp;
  T threshr;
  vector<Vec3> shapebase;
};

template <typename T> matchPartialPartial<T>::matchPartialPartial() {
  I  = sqrt(U(- T(1)));
  Pi = atan2(T(1), T(1)) * T(4);
}

template <typename T> matchPartialPartial<T>::~matchPartialPartial() {
  ;
}

template <typename T> void matchPartialPartial<T>::init(const vector<Vec3>& shapebase, const T& thresh, const T& threshp, const T& threshr) {
  this->shapebase = shapebase;
  this->thresh    = T(1) - thresh;
  this->threshp   = threshp;
  this->threshr   = threshr;
  return;
}

template <typename T> vector<match_t<T> > matchPartialPartial<T>::match(const vector<Vec3>& points, const int& ndiv, const T& r_max_theta) {
  vector<match_t<T> > result;
  for(int i = 0; i < 3; i ++) {
    Eigen::Matrix<Eigen::Matrix<T, 2, 1>, Eigen::Dynamic, Eigen::Dynamic> table(shapebase.size(), points.size());
    // init table.
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for
#endif
    for(int j = 0; j < shapebase.size(); j ++) {
      cerr << "making table (" << i << "/3)" << " : " << j << "/" << shapebase.size() << endl;
      for(int k = 0; k < points.size(); k ++) {
        Vec3 aj(shapebase[j]), bk(points[k]);
        aj /= sqrt(aj.dot(aj));
        bk /= sqrt(bk.dot(bk));
        for(int l = 0; l < 2; l ++) {
          const T a(bk[(l + i) % 3]);
          const T b(bk[(l + i + 1) % 3]);
          const T c(aj[(l + i) % 3]);
          if(a * a + b * b < c * c) {
            table(j, k)[l] = sqrt(- T(1));
            continue;
          }
          const T theta0(T(2) * atan2(sqrt(a * a + b * b - c * c) - b, a + c));
          const T theta1(T(2) * atan2(- sqrt(a * a + b * b - c * c) - b, a + c));
          table(j, k)[l] = sqrt(- T(1));
          if(isfinite(theta0) &&
             abs((cos(theta0) * a - sin(theta0) * b) - c) < thresh)
            table(j, k)[l] = theta0;
          else if(isfinite(theta1) &&
             abs((cos(theta1) * a - sin(theta1) * b) - c) < thresh)
            table(j, k)[l] = theta1;
          const T& theta(table(j, k)[l]);
          bk[(l + i    ) % 3] = cos(theta) * a - sin(theta) * b;
          bk[(l + i + 1) % 3] = sin(theta) * a + cos(theta) * b;
        }
        const Vec3 err(aj - aj.dot(bk) * bk / bk.dot(bk));
        if(thresh * thresh < err.dot(err))
          table(j, k)[0] = table(j, k)[1] = sqrt(- T(1));
      }
    }
    // matches.
#if defined(_OPENMP)
#pragma omp for
#endif
    for(int nd = 0; nd < ndiv + 1; nd ++) {
      cerr << "matching table (" << i << "/3)" << " : " << nd << "/" << ndiv + 1 << endl;
      Eigen::Matrix<T, 2, 1> ddiv;
      ddiv[0] = cos(2. * Pi * nd / ndiv);
      ddiv[1] = sin(2. * Pi * nd / ndiv);
      if(nd == ndiv)
        ddiv *= T(0);
      Mat3x3 drot0;
      for(int k = 0; k < drot0.rows(); k ++)
        for(int l = 0; l < drot0.cols(); l ++)
          drot0(k, l) = (k == l ? T(1) : T(0));
      vector<msub_t<T> > msub;
      for(int j = 0; j < shapebase.size(); j ++) {
        for(int k = 0; k < points.size(); k ++) {
          const T ldepth(table(j, k).dot(ddiv) / sqrt(table(j, k).dot(table(j, k))) - T(1));
          if(isfinite(ldepth) && abs(ldepth) <= thresh) {
            msub_t<T> work;
            work.mbufj = j;
            work.mbufk = k;
            work.mbufN = sqrt(table(j, k).dot(table(j, k)));
            msub.push_back(work);
          }
        }
      }
      sort(msub.begin(), msub.end(), cmpsubwrap<T>);
      for(int k0 = 0; k0 < msub.size(); k0 ++) {
        match_t<T> work;
        work.rot = drot0;
        for(int k = 0; k < ddiv.size(); k ++) {
          Mat3x3 lrot;
          const T theta(ddiv[k] * msub[k0].mbufN);
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
        if(abs(work.rot(2, 2)) < r_max_theta)
          continue;
        bool flagk[points.size()];
        bool flagj[shapebase.size()];
        for(int kk = 0; kk < points.size(); kk ++)
          flagk[kk] = false;
        for(int kk = 0; kk < shapebase.size(); kk ++)
          flagj[kk] = false;
        for(int k = k0; k < msub.size(); k ++) {
          int kfix = k;
          T   err(thresh * T(2));
          for(int kk = k;
                  kk < msub.size();
                  kk ++)
            if(!flagj[msub[kk].mbufj] &&
               !flagk[msub[kk].mbufk]) {
              Vec3 aj(shapebase[msub[kk].mbufj]);
              Vec3 bk(work.rot * points[msub[kk].mbufk]);
              aj /= sqrt(aj.dot(aj));
              aj -= aj.dot(bk) * bk / bk.dot(bk);
              if(!isfinite(aj.dot(aj)) || thresh < aj.dot(aj))
                continue;
              if(aj.dot(aj) < err) {
                err = aj.dot(aj);
                kfix = kk;
              }
              // XXX confirm me. ;?
              break;
            }
          if(thresh < err)
            break;
          work.dstpoints.push_back(msub[kfix].mbufj);
          work.srcpoints.push_back(msub[kfix].mbufk);
          flagj[msub[kfix].mbufj] = true;
          flagk[msub[kfix].mbufk] = true;
        }
        // if there's a match.
        work.rpoints = work.dstpoints.size() / T(min(shapebase.size(), points.size()));
        if(threshp <= work.rpoints) {
          cerr << "*";
          fflush(stderr);
          work.offset[0] = T(0);
          work.offset[1] = T(0);
          work.offset[2] = T(0);
          for(int k = 0; k < work.dstpoints.size(); k ++)
            work.offset += shapebase[work.dstpoints[k]];
          // maximize parallel parts.
          Mat2x2 A;
          Vec2   b;
          b[0]       = b[1] = T(0);
          A(0, 0)    = T(0);
          A(0, 1)    = T(0);
          A(1, 1)    = T(0);
          for(int jj = 0; jj < work.dstpoints.size(); jj ++) {
            A(0, 0) += work.offset.dot(work.offset);
            A(0, 1) += work.offset.dot(work.rot * points[work.srcpoints[jj]]);
            A(1, 1) += points[work.srcpoints[jj]].dot(points[work.srcpoints[jj]]);
            b[0]    += shapebase[work.dstpoints[jj]].dot(work.offset);
            b[1]    += shapebase[work.dstpoints[jj]].dot(work.rot * points[work.srcpoints[jj]]);
          }
          A(1, 0) = A(0, 1);
          Eigen::RealSchur<Mat2x2> schur;
          schur.compute(A, true);
          const Mat2x2 U(schur.matrixU());
          const Mat2x2 L(schur.matrixT());
          b     = U.transpose() * b;
          b[0] /= L(0, 0);
          b[1] /= L(1, 1);
          b     = U * b;
          work.offset *= b[0];
          work.ratio   = b[1];
          if(abs(work.ratio) < threshr)
            continue;
          T err(0), n2(0);
          for(int k = 0; k < work.dstpoints.size(); k ++) {
            const Vec3 a(shapebase[work.dstpoints[k]]);
            const Vec3 b(work.rot * points[work.srcpoints[k]] * work.ratio + work.offset);
            n2  += a.dot(a) * b.dot(b);
            err += (a - b).dot(a - b);
          }
          work.rdepth = T(1) / sqrt(err / n2);
          if(sqrt(err / n2) < thresh) {
#if defined(_OPENMP)
#pragma omp critical
#endif
            {
              bool flag = false;
              // eliminate similar matches.
              for(int kk = 0; kk < result.size(); kk ++) {
                Vec3   test(work.offset - result[kk].offset);
                Mat3x3 roterr(work.rot * result[kk].rot.transpose());
                T      rotdnorm(0);
                for(int l = 0; l < work.rot.rows(); l ++)
                  rotdnorm += pow(roterr(l, l) - T(1), T(2));
                // XXX magic number:
                if(sqrt(test.dot(test) / work.offset.dot(work.offset) +
                        rotdnorm +
                        pow(max(T(1) - work.ratio / result[kk].ratio,
                                T(1) - result[kk].ratio / work.ratio), T(2)))
                    <= T(7) * T(.2)) {
                  flag = true;
                  break;
                }
              }
              if(!flag)
                result.push_back(work);
            }
          }
        }
      }
    }
  }
  sort(result.begin(), result.end(), cmpwrap<T>);
  return result;
}


template <typename T> class matchWholePartial {
public:
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T, 3, 3>                           Mat3x3;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
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
  typedef Eigen::Matrix<T, 2, 1>              Vec2;
  
  reDig();
  ~reDig();
  void init();
  Vec3 emphasis0(const Vec3& dst, const Vec3& src, const match_t<T>& match, const T& ratio);
  vector<Vec3> emphasis(const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const T& ratio);
  Mat          emphasis(const Mat& dstimg, const Mat& dstbump, const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const T& ratio);
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

template <typename T> Eigen::Matrix<T, 3, 1> reDig<T>::emphasis0(const Vec3& dst, const Vec3& src, const match_t<T>& match, const T& ratio) {
  const Vec3 a(dst);
  const Vec3 b(match.rot * src * match.ratio + match.offset);
  const T    r0(sqrt((a - b).dot(a - b)));
  return a + (b - a) * (ratio - T(1)) / r0;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> reDig<T>::emphasis(const Mat& dstimg, const Mat& dstbump, const vector<Vec3>& dst, const vector<Vec3>& src, const match_t<T>& match, const T& ratio) {
  cerr << " making triangles";
  fflush(stderr);
  tilter<T> tilt;
  vector<typename tilter<T>::Triangles> triangles;
  for(int i = 0; i < dstimg.rows() - 1; i ++)
    for(int j = 0; j < dstimg.cols() - 1; j ++) {
      triangles.push_back(tilt.makeTriangle(i, j, dstimg, dstbump, false));
      triangles.push_back(tilt.makeTriangle(i, j, dstimg, dstbump, true));
    }
  bool checked[triangles.size()];
  for(int i = 0; i < triangles.size(); i ++)
    checked[i] = false;
  
  for(int i = 0; i < match.dstpoints.size(); i ++)
    for(int j = i + 1; j < match.dstpoints.size(); j ++)
      for(int k = max(i, j) + 1; k < match.dstpoints.size(); k ++) {
        const Vec3 n(tilt.solveN(dst[match.dstpoints[i]],
                                 dst[match.dstpoints[j]],
                                 dst[match.dstpoints[k]]));
        if(n[2] < T(0))
          continue;
        bool flag = true;
        for(int l = 1; l < match.dstpoints.size(); l ++)
          if(n.dot(dst[match.dstpoints[0]]) *
             n.dot(dst[match.dstpoints[l]]) < T(0)) {
            flag = false;
            break;
          }
        if(!flag)
          continue;
        Vec2 p0, p1, p2;
        p0[0] = dst[match.dstpoints[i]][0];
        p0[1] = dst[match.dstpoints[i]][1];
        p1[0] = dst[match.dstpoints[j]][0];
        p1[1] = dst[match.dstpoints[j]][1];
        p2[0] = dst[match.dstpoints[k]][0];
        p2[1] = dst[match.dstpoints[k]][1];
        for(int l = 0; l < triangles.size(); l ++)
          if(!checked[l]) {
            Vec2 q;
            for(int ll = 0; ll < 3; ll ++) {
              q[0] = triangles[l](0, ll);
              q[1] = triangles[l](1, ll);
              if(tilt.sameSide2(p0, p1, p2, q) &&
                 tilt.sameSide2(p1, p2, p0, q) &&
                 tilt.sameSide2(p2, p0, p1, q)) {
                triangles[l].col(ll) = emphasis0(triangles[l].col(ll), src[match.srcpoints[i]], match, ratio);
              }
            }
            triangles[l].col(4) = tilt.solveN(triangles[l].col(0), triangles[l].col(1), triangles[l].col(2));
            triangles[l](1, 3)  = triangles[l].col(4).dot(triangles[l].col(0));
            checked[l] = true;
          }
      }
  
  // for each pixels extends.
  Mat I3(3, 3);
  Vec zero3(3);
  for(int i = 0; i < 3; i ++) {
    for(int j = 0; j < 3; j ++)
      I3(i, j) = (i == j ? T(1) : T(0));
    zero3[i] = T(0);
  }
  return tilt.tiltsub(dstimg, triangles, I3, I3, zero3, zero3, T(1));
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> showMatch(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& input, const vector<Eigen::Matrix<T, 3, 1> >& refpoints, const vector<int>& prefpoints, const T& emph = T(.8)) {
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result(input), map(input.rows(), input.cols());
  for(int i = 0; i < map.rows(); i ++)
    for(int j = 0; j < map.cols(); j ++)
      map(i, j) = T(0);
  for(int k = 0; k < prefpoints.size(); k ++) {
    const Eigen::Matrix<T, 3, 1>& lref0(refpoints[prefpoints[k]]);
    for(int l = k + 1; l < prefpoints.size(); l ++) {
      const Eigen::Matrix<T, 3, 1>& lref1(refpoints[prefpoints[l]]);
      if(lref0 == lref1) continue;
      if(abs(lref1[1] - lref0[1]) > abs(lref1[0] - lref0[0])) {
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
                   lref0[1] + tilt * ii < map.cols();
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
                   lref0[0] + tilt * jj < map.rows();
              j += sgndelta, jj += sgndelta)
          for(int ii = 0; ii < abs(tilt); ii ++) {
            int i = tilt * (jj - sgndelta) +
                    ii * (tilt * sgndelta < T(0) ? - 1 : 1);
            map(max(0, min(int(lref0[0] + i), int(map.rows() - 1))),
                max(0, min(j, int(map.cols() - 1)))) = emph;
          }
      }
    }
  }
  return result + map;
}

#define _SCAN_CONTEXT_
#endif

