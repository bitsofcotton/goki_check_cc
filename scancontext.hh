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
using std::max;
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
  return x0.score > x1.score;
}

template <typename T> int cmplfrwrap(const lfmatch_t<T>& x0, const lfmatch_t<T>& x1) {
  return x0.pt[0] < x1.pt[0] || (x0.pt[0] == x1.pt[0] && x0.pt[1] < x1.pt[1]);
}

template <typename T> class lowFreq {
public:
  typedef complex<T> U;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1>              Vec;
  typedef Eigen::Matrix<T, 3, 1>                           Vec3;
  lowFreq();
  ~lowFreq();
  
  vector<Eigen::Matrix<T, 3, 1> > getLowFreq(const Mat& data, const int& npoints = 60);
  Mat getLowFreqImage(const Mat& data, const int& npoints = 60);
private:
  vector<lfmatch_t<T> > prepareCost(const Mat& data, const int& npoints);
  T Pi;
  U I;
};

template <typename T> lowFreq<T>::lowFreq() {
  Pi          = atan2(T(1), T(1)) * T(4);
  I           = sqrt(U(- 1));
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
      pt /= count;
      // XXX magic number:
      pt[2] *= T(60);
      match.push_back(m);
    }
  sort(match.begin(), match.end(), cmplfwrap<T>);
  // XXX add me lowpoly:
  vector<lfmatch_t<T> > result(match);
  sort(result.begin(), result.end(), cmplfrwrap<T>);
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
    ;
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
    return *this;
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
  void init(const int& ndiv, const T& thresh, const T& threshl, const T& threshp, const T& threshr, const T& threshN, const T& threshc, const T& r_max_theta);
  
  vector<match_t<T> > match(const vector<Vec3>& shapebase, const vector<Vec3>& points);
private:
  U   I;
  T   Pi;
  int ndiv;
  T   thresh;
  T   threshl;
  T   threshN;
  T   threshp;
  T   threshr;
  T   threshc;
  T   r_max_theta;
};

template <typename T> matchPartialPartial<T>::matchPartialPartial() {
  I  = sqrt(U(- T(1)));
  Pi = atan2(T(1), T(1)) * T(4);
  init(12, .95, .8, .0625, .125, 2., .75, .1);
}

template <typename T> matchPartialPartial<T>::~matchPartialPartial() {
  ;
}

template <typename T> void matchPartialPartial<T>::init(const int& ndiv, const T& thresh, const T& threshl, const T& threshp, const T& threshr, const T& threshN, const T& threshc, const T& r_max_theta) {
  this->ndiv        = ndiv;
  this->thresh      = abs(T(1) - thresh);
  this->threshl     = abs(T(1) - threshl);
  this->threshN     = threshN;
  this->threshp     = threshp;
  this->threshr     = threshr;
  this->threshc     = threshc;
  this->r_max_theta = r_max_theta;
  return;
}

template <typename T> vector<match_t<T> > matchPartialPartial<T>::match(const vector<Vec3>& shapebase, const vector<Vec3>& points) {
  vector<match_t<T> > result;
  Mat3x3 drot0;
  for(int k = 0; k < drot0.rows(); k ++)
    for(int l = 0; l < drot0.cols(); l ++)
      drot0(k, l) = (k == l ? T(1) : T(0));
  for(int i = 0; i < 3; i ++) {
    Eigen::Matrix<Eigen::Matrix<T, 2, 1>, Eigen::Dynamic, Eigen::Dynamic> table(shapebase.size(), points.size());
    // init table.
    cerr << "making table (" << i << "/3)" << endl;
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for
#endif
    for(int j = 0; j < shapebase.size(); j ++) {
      cerr << ".";
      fflush(stderr);
      for(int k = 0; k < points.size(); k ++) {
        Vec3 aj(shapebase[j]), bk(points[k]);
        aj /= sqrt(aj.dot(aj));
        bk /= sqrt(bk.dot(bk));
        for(int l = 0; l < table(j, k).size(); l ++) {
          const T a(bk[(l + i) % 3]);
          const T b(bk[(l + i + 1) % 3]);
          const T c(aj[(l + i) % 3]);
          table(j, k)[l] = sqrt(- T(1));
          if(a * a + b * b < c * c)
            continue;
          const T theta0(T(2) * atan2(  sqrt(a * a + b * b - c * c) - b, a + c));
          const T theta1(T(2) * atan2(- sqrt(a * a + b * b - c * c) - b, a + c));
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
      cerr << "matching table (" << i << "/3)" << " : " << nd << "/" << ndiv + 1;
      Eigen::Matrix<T, 2, 1> ddiv;
      ddiv[0] = cos(2. * Pi * nd / ndiv);
      ddiv[1] = sin(2. * Pi * nd / ndiv);
      if(nd == ndiv)
        ddiv *= T(0);
      vector<msub_t<T> > msub;
      for(int k = 0; k < points.size(); k ++) {
        for(int j = 0; j < shapebase.size(); j ++) {
          T ldepth(1);
          if(nd < ndiv)
            ldepth = table(j, k).dot(ddiv) / sqrt(table(j, k).dot(table(j, k))) - T(1);
          else
            ldepth = sqrt(table(j, k).dot(table(j, k)));
          if(isfinite(ldepth) && abs(ldepth) <= thresh) {
            msub_t<T> workm;
            workm.mbufj = j;
            workm.mbufk = k;
            workm.mbufN = sqrt(table(j, k).dot(table(j, k)));
            msub.push_back(workm);
          }
        }
      }
      cerr << " : " << msub.size() << endl;
      if(!msub.size())
        continue;
      sort(msub.begin(), msub.end(), cmpsubwrap<T>);
      int k0(0);
      for(; k0 < msub.size(); k0 ++) {
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
        if(abs(work.rot((i + 2) % 3, (i + 2) % 3)) < r_max_theta)
          continue;
        bool flagk[points.size()];
        bool flagj[shapebase.size()];
        for(int kk = 0; kk < points.size(); kk ++)
          flagk[kk] = false;
        for(int kk = 0; kk < shapebase.size(); kk ++)
          flagj[kk] = false;
        const Vec3 ajk0(shapebase[msub[k0].mbufj]);
        const Vec3 bkk0(work.rot * points[msub[k0].mbufk]);
        const T    t0(ajk0.dot(bkk0) / bkk0.dot(bkk0));
        int  kk(k0 + 1);
        for(; kk < msub.size(); kk ++) {
          bool flagt(false);
          for(; kk < msub.size(); kk ++)
            if(!flagj[msub[kk].mbufj] &&
               !flagk[msub[kk].mbufk] &&
               abs(msub[k0].mbufN - msub[kk].mbufN) / msub[kk].mbufN < threshN) {
              const Vec3 aj0(shapebase[msub[kk].mbufj]);
              const Vec3 bk(work.rot * points[msub[kk].mbufk]);
              Vec3 aj(aj0);
              aj /= sqrt(aj.dot(aj));
              aj -= aj.dot(bk) * bk / bk.dot(bk);
              if(!isfinite(aj.dot(aj)) || thresh * thresh < aj.dot(aj))
                continue;
              const T t(aj0.dot(bk) / bk.dot(bk));
              if(!isfinite(t) || threshl < max(abs(t0 / t - T(1)),
                                               abs(t / t0 - T(1))))
                continue;
              flagt = true;
              break;
            }
          if(!flagt || msub.size() <= kk)
            break;
          work.dstpoints.push_back(int(msub[kk].mbufj));
          work.srcpoints.push_back(int(msub[kk].mbufk));
          flagj[msub[kk].mbufj] = true;
          flagk[msub[kk].mbufk] = true;
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
          work.offset /= work.dstpoints.size();
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
          T err(0);
          for(int k = 0; k < work.dstpoints.size(); k ++) {
            const Vec3 a(shapebase[work.dstpoints[k]]);
            const Vec3 b(work.rot * points[work.srcpoints[k]] * work.ratio + work.offset);
            err += (a - b).dot(a - b) / a.dot(a) / b.dot(b);
          }
          err /= work.dstpoints.size();
          work.rdepth = T(1) / sqrt(err);
          if(sqrt(err) < thresh) {
#if defined(_OPENMP)
#pragma omp critical
#endif
            {
              int idxf(- 1);
              // eliminate similar matches.
              for(int kk = 0; kk < result.size(); kk ++) {
                Vec3   test(work.offset - result[kk].offset);
                Mat3x3 roterr(work.rot * result[kk].rot.transpose());
                T      rotdnorm(0);
                for(int l = 0; l < work.rot.rows(); l ++)
                  rotdnorm += pow(roterr(l, l) - T(1), T(2));
                // XXX checkme:
                rotdnorm  = sqrt(sqrt(rotdnorm)) / T(3);
                rotdnorm += sqrt(test.dot(test) / 
                                 (work.offset.dot(work.offset) *
                                  result[kk].offset.dot(result[kk].offset)))
                            / T(3);
                rotdnorm +=
                   min(min(pow(T(1) - work.ratio / result[kk].ratio, T(2)),
                           pow(T(1) - result[kk].ratio / work.ratio, T(2))),
                       T(1));
                if(rotdnorm <= threshc) {
                  vector<match_t<T>> workbuf;
                  workbuf.push_back(result[kk]);
                  workbuf.push_back(work);
                  sort(workbuf.begin(), workbuf.end(), cmpwrap<T>);
                  if(workbuf[0].rpoints == work.rpoints &&
                     workbuf[0].rdepth  == work.rdepth) {
                    idxf = kk;
                    break;
                  }
                }
              }
              if(idxf >= 0)
                result[idxf] = work;
              else
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

