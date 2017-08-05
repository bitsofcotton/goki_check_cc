#if !defined(_SCAN_CONTEXT_)

#include <Eigen/Core>
#include <Eigen/LU>
#include <cmath>
#include <vector>

using std::sqrt;
using std::atan2;
using std::abs;
using std::cos;
using std::sin;
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
  return x0.score > x1.score;
  // return x0.score < x1.score;
}

template <typename T> class lowFreq {
public:
  typedef complex<T> U;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1>              Vec;
  lowFreq();
  ~lowFreq();
  
  vector<Eigen::Matrix<T, 3, 1> > getLowFreq(const Mat& data, const int& npoints = 600);
  Mat getLowFreqImage(const Mat& data, const int& npoints = 600);
private:
  vector<lfmatch_t<T> > prepareCost(const Mat& data, const int& npoints);
  T Pi;
  U I;
  T guard;
};

template <typename T> lowFreq<T>::lowFreq() {
  Pi     = atan2(T(1), T(1)) * T(4);
  I      = sqrt(U(- 1));
  guard  = T(20);
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
  for(int i = 1; i < data.rows() - 1; i ++)
    for(int j = 1; j < data.cols() - 1; j ++) {
      lfmatch_t<T> m;
      T& cost(m.score);
      cost  = abs(data(i - 1, j - 1) - data(i, j));
      cost += abs(data(i - 1, j)     - data(i, j));
      cost += abs(data(i - 1, j + 1) - data(i, j));
      cost += abs(data(i, j - 1) - data(i, j));
      cost += abs(data(i, j + 1) - data(i, j));
      cost += abs(data(i + 1, j - 1) - data(i, j));
      cost += abs(data(i + 1, j)     - data(i, j));
      cost += abs(data(i + 1, j + 1) - data(i, j));
      cost /= T(8);
      Eigen::Matrix<T, 3, 1>& pt(m.pt);
      pt[0] = i;
      pt[1] = j;
      pt[2] = data(i, j);
      match.push_back(m);
    }
  sort(match.begin(), match.end(), cmplfwrap<T>);
  vector<lfmatch_t<T> > result;
  for(int i = 0; i < match.size() && result.size() < npoints; i ++) {
    bool flag = false;
    for(int j = 0; j < result.size(); j ++) {
      const Eigen::Matrix<T, 3, 1> diff(match[i].pt - result[j].pt);
      const T    r2(diff.dot(diff));
      if(r2 < guard) {
        flag = true;
        break;
      }
    }
    if(!flag)
      result.push_back(match[i]);
  }
  return result;
}

template <typename T> class match_t {
public:
  Eigen::Matrix<T, 3, 3> rot;
  Eigen::Matrix<T, 3, 1> offset;
  T                      ratio;
  T                      rdepth;
  T                      rpoints;
  vector<int>       dstpoints;
  vector<int>       srcpoints;
  match_t() {
    ;
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
  return x0.rdepth > x1.rdepth || (x0.rdepth == x1.rdepth && x0.rpoints > x1.rpoints);
//  return x0.rpoints > x1.rpoints || (x0.rpoints == x1.rpoints && x0.rdepth > x1.rdepth);
}

template <typename T> int cmpsubwrap(const msub_t<T>& x0, const msub_t<T>& x1) {
  return x0.mbufN < x1.mbufN || (x0.mbufN == x1.mbufN && x0.mbufj < x1.mbufj) || (x0.mbufN == x1.mbufN && x0.mbufj == x1.mbufj && x0.mbufk < x1.mbufk);
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
  void init(const vector<Vec3>& shapebase, const T& thresh, const T& threshp);
  
  vector<match_t<T> > match(const vector<Vec3>& points, const int& ndiv = 100);
private:
  U I;
  T Pi;
  // thresh, threshp in [0, 1]
  T thresh;
  T threshp;
  vector<Vec3> shapebase;
};

template <typename T> matchPartialPartial<T>::matchPartialPartial() {
  I  = sqrt(U(- T(1)));
  Pi = atan2(T(1), T(1)) * T(4);
}

template <typename T> matchPartialPartial<T>::~matchPartialPartial() {
  ;
}

template <typename T> void matchPartialPartial<T>::init(const vector<Vec3>& shapebase, const T& thresh, const T& threshp) {
  this->shapebase = shapebase;
  this->thresh    = T(1) - thresh;
  this->threshp   = threshp;
  return;
}

template <typename T> vector<match_t<T> > matchPartialPartial<T>::match(const vector<Vec3>& points, const int& ndiv) {
  vector<match_t<T> > result;
  // init tables.
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
          T theta(theta0);
          T b0(cos(theta) * a - sin(theta) * b);
          if(!(isfinite(theta) && abs(b0 - c) < thresh) && isfinite(theta1))
            theta = theta1;
          table(j, k)[l] = theta;
          b0 = cos(theta) * a - sin(theta) * b;
          const T b1(sin(theta) * a + cos(theta) * b);
          bk[(l + i    ) % 3] = b0;
          bk[(l + i + 1) % 3] = b1;
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
    for(int nd = 0; nd < ndiv; nd ++) {
      cerr << "matching table (" << i << "/3)" << " : " << nd << "/" << ndiv << endl;
      Eigen::Matrix<T, 2, 1> ddiv;
      ddiv[0] = cos(2. * Pi * nd / ndiv);
      ddiv[1] = sin(2. * Pi * nd / ndiv);
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
        int kl = k0;
        Eigen::Matrix<T, 2, 1> theta0;
        theta0[0] = theta0[1] = T(0);
        T                 mmdepth(0);
        vector<int>  jmatch;
        vector<int>  kmatch;
        bool              flag[points.size()];
        for(int kk = 0; kk < points.size(); kk ++)
          flag[kk] = false;
        for(int kk = k0, k = k0; kk < msub.size(); kk ++) {
          const T diff(abs(msub[k0].mbufN - msub[kk].mbufN) / msub[msub.size() - 1].mbufN);
          if(diff < thresh &&
             (kk == k0 || msub[k].mbufj != msub[kk].mbufj)) {
            for(; k + 1 < msub.size() && kk < msub.size() &&
                  flag[msub[kk].mbufk] &&
                  msub[k + 1].mbufj == msub[kk].mbufj; kk ++);
            if(msub.size() <= k + 1 || msub.size() <= kk ||
               flag[msub[kk].mbufk] ||
               msub[k + 1].mbufj != msub[kk].mbufj) {
              k = kk;
              continue;
            }
            jmatch.push_back(msub[kk].mbufj);
            kmatch.push_back(msub[kk].mbufk);
            theta0  += ddiv * msub[kk].mbufN;
            mmdepth = max(diff, mmdepth);
            flag[msub[kk].mbufk] = true;
            // skips through.
            kl = kk;
          }
          k  = kk;
        }
        k0 = kl;
        if(threshp <= jmatch.size() / T(min(shapebase.size(), points.size()))) {
#if defined(_OPENMP)
#pragma omp critical
#endif
         {
          cerr << "*";
          fflush(stderr);
          match_t<T> work;
          work.dstpoints = jmatch;
          work.srcpoints = kmatch;
          work.rpoints   = jmatch.size() / min(T(shapebase.size()), T(points.size()));
          work.rdepth    = T(1) / mmdepth;
          theta0        /= jmatch.size();
          work.offset[0] = T(0);
          work.offset[1] = T(0);
          work.offset[2] = T(0);
          for(int k = 0; k < jmatch.size(); k ++)
            work.offset += shapebase[jmatch[k]];
          for(int k = 0; k < work.rot.rows(); k ++)
            for(int l = 0; l < work.rot.cols(); l ++)
              work.rot(k, l) = (k == l ? T(1) : T(0));
          for(int k = 0; k < theta0.size(); k ++) {
            Mat3x3 lrot;
            lrot((k + i    ) % 3, (k + i    ) % 3) =   cos(theta0[k]);
            lrot((k + i + 1) % 3, (k + i    ) % 3) =   sin(theta0[k]);
            lrot((k + i    ) % 3, (k + i + 1) % 3) = - sin(theta0[k]);
            lrot((k + i + 1) % 3, (k + i + 1) % 3) =   cos(theta0[k]);
            lrot((k + i + 2) % 3, (k + i    ) % 3) = T(0);
            lrot((k + i + 2) % 3, (k + i + 1) % 3) = T(0);
            lrot((k + i + 2) % 3, (k + i + 2) % 3) = T(1);
            lrot((k + i    ) % 3, (k + i + 2) % 3) = T(0);
            lrot((k + i + 1) % 3, (k + i + 2) % 3) = T(0);
            work.rot = lrot * work.rot;
          }
          Mat2x2 A;
          Vec2   b;
          b[0]       = b[1]       = T(0);
          A(0, 0)    = work.offset.dot(work.offset);// * jmatch.size();
          A(0, 1)    = T(0);
          A(1, 1)    = T(0);
          for(int jj = 0; jj < jmatch.size(); jj ++) {
            A(0, 1) += points[kmatch[jj]].dot(work.rot.transpose() * work.offset);
            A(1, 1) += points[kmatch[jj]].dot(points[kmatch[jj]]);
            b[0]    += shapebase[jmatch[jj]].dot(work.offset);
            b[1]    += shapebase[jmatch[jj]].dot(work.rot * points[kmatch[jj]]);
          }
          A(1, 0) = A(0, 1);
          Eigen::RealSchur<Mat2x2> schur;
          schur.compute(A, true);
          const Mat2x2 U(schur.matrixU());
          b    = U.transpose() * b;
          const Mat2x2 L(schur.matrixT());
          Vec2 bb;
          if(abs(L(0, 0)) < thresh)
            b[0]  = T(0);
          else
            b[0] /= L(0, 0);
          if(abs(L(1, 1)) < thresh)
            b[1]  = T(0);
          else
            b[1] /= L(1, 1);
          b    = U * b;
          work.offset *= b[0];
          work.ratio   = b[1];
          result.push_back(work);
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


template <typename T> class reDig {
public:
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T, 3, 3>                           Mat3x3;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
  typedef Eigen::Matrix<T, 3, 1>              Vec3;
  
  void init();
  Mat deemphasis(const Mat& dst, const Mat& src, const vector<match_t<T> >& match, const T& ratio);
  Mat emphasis(const Mat& dst, const Mat& src, const vector<match_t<T> >& match, const T& ratio);
};


template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> showMatch(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& input, const vector<Eigen::Matrix<T, 3, 1> >& refpoints, const vector<int>& prefpoints) {
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result(input), map(input.rows(), input.cols());
  for(int i = 0; i < map.rows(); i ++)
    for(int j = 0; j < map.cols(); j ++)
      map(i, j) = T(0);
  for(int k = 0; k < prefpoints.size(); k ++) {
    const Eigen::Matrix<T, 3, 1>& lref0(refpoints[prefpoints[k]]);
    for(int l = k + 1; l < prefpoints.size(); l ++) {
      const Eigen::Matrix<T, 3, 1>& lref1(refpoints[prefpoints[l]]);
      if(abs(lref1[1] - lref0[1]) > abs(lref1[0] - lref0[0])) {
        const T tilt((lref1[1] - lref0[1]) / (lref1[0] - lref0[0]));
        T sgndelta(1);
        if(lref1[0] < lref0[0])
          sgndelta = - T(1);
        for(int i = lref0[0], ii = 0; i != int(lref1[0]) && 0 <= lref0[1] + tilt * ii && lref0[1] + tilt * ii < map.cols(); i += sgndelta, ii += sgndelta)
          // XXX magic number:
          map(i, int(lref0[1] + tilt * ii)) = T(.6);
      } else {
        const T tilt((lref1[0] - lref0[0]) / (lref1[1] - lref0[1]));
        T sgndelta(1);
        if(lref1[1] < lref0[1])
          sgndelta = - T(1);
        for(int j = lref0[1], jj = 0; j != int(lref1[1]) && 0 <= lref0[0] + tilt * jj && lref0[0] + tilt * jj < map.rows(); j += sgndelta, jj += sgndelta)
          // XXX magic number:
          map(int(lref0[0] + tilt * jj), j) = T(.6);
      }
    }
  }
  return result + map;
}

#define _SCAN_CONTEXT_
#endif

