#if !defined(_ENLARGE2X_)

#include <Eigen/Core>
#include <Eigen/LU>

using std::complex;
using std::abs;
using std::sqrt;
using std::exp;
using std::vector;
using std::sort;

template <typename T> class enlarger2ex {
public:
  typedef enum {
    ENLARGE_X,
    ENLARGE_Y,
    ENLARGE_FX,
    ENLARGE_FY,
    ENLARGE_BOTH,
    ENLARGE_FBOTH,
    DETECT_X,
    DETECT_Y,
    DETECT_BOTH,
    COLLECT_X,
    COLLECT_Y,
    COLLECT_BOTH } direction_t;
  typedef complex<T> U;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> MatU;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1>              Vec;
  typedef Eigen::Matrix<U, Eigen::Dynamic, 1>              VecU;
  enlarger2ex();
  Mat compute(const Mat& data, const direction_t& dir);
private:
  void seedPattern(const int&, const int&, const int&);
  void initPattern(const int& size, const bool& flag = true);
  U    I;
  T    Pi;
  MatU B;
  MatU G;
  Mat  D;
  Mat  Dop;
  Mat  Iop;
  Vec  F;
  Vec  X;
};

template <typename T> enlarger2ex<T>::enlarger2ex() {
  I  = sqrt(U(- 1.));
  Pi = atan2(T(1), T(1)) * T(4);
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> enlarger2ex<T>::compute(const Mat& data, const direction_t& dir) {
  Mat result;
  switch(dir) {
  case ENLARGE_BOTH:
    {
      Mat rw, rh;
      rw = compute(data, ENLARGE_X);
      rh = compute(data, ENLARGE_Y);
      rw = compute(rw,   ENLARGE_Y);
      rh = compute(rh,   ENLARGE_X);
      result = (rw + rh) / 2.;
    }
    break;
  case ENLARGE_FBOTH:
    {
      Mat rw, rh;
      rw = compute(data, ENLARGE_FX);
      rh = compute(data, ENLARGE_FY);
      rw = compute(rw,   ENLARGE_FY);
      rh = compute(rh,   ENLARGE_FX);
      result = (rw + rh) / 2.;
    }
    break;
  case DETECT_BOTH:
    result = (compute(data, DETECT_X) + compute(data, DETECT_Y)) / 2.;
    break;
  case COLLECT_BOTH:
    result = (compute(data, COLLECT_X) + compute(data, COLLECT_Y)) / 2.;
    break;
  case ENLARGE_X:
    result = compute(data.transpose(), ENLARGE_Y).transpose();
    break;
  case ENLARGE_FX:
    result = compute(data.transpose(), ENLARGE_FY).transpose();
    break;
  case DETECT_X:
    result = compute(data.transpose(), DETECT_Y).transpose();
    break;
  case COLLECT_X:
    result = compute(data.transpose(), COLLECT_Y).transpose();
    break;
  case ENLARGE_Y:
    {
      result = Mat(data.rows() * 2, data.cols());
      cerr << " enlarge_y";
      initPattern(data.rows());
      Mat dd(D * data);
      Vec ff(data.transpose() * F);
      Vec xx(data.transpose() * X);
      Mat co(Dop * data);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int i = 0; i < data.cols(); i ++) {
        ff[i] /= T(2);
        xx[i] /= T(2);
        for(int j = 0; j < data.rows(); j ++) {
          const T delta(xx[i] * co(j, i) + ff[i] * dd(j, i));
          result(j * 2 + 0, i) = data(j, i) - delta;
          result(j * 2 + 1, i) = data(j, i) + delta;
        }
      }
    }
    break;
  case ENLARGE_FY:
    {
      initPattern(data.rows(), false);
      const Mat dcache(Dop * data);
      const Mat work(compute(dcache, ENLARGE_Y));
            Mat work2(data.rows() * 2, data.cols());
      for(int j = 0; j < work2.rows(); j ++)
        work2.row(j) = dcache.row(j / 2);
      initPattern(work.rows(), false);
      result = Iop * (work - work2);
    }
    break;
  case DETECT_Y:
    initPattern(data.rows(), false);
    result = Dop * data;
    break;
  case COLLECT_Y:
    result = compute(data, DETECT_Y);
    for(int i = 0; i < result.rows(); i ++)
      for(int j = 0; j < result.cols(); j ++)
        result(i, j) = abs(result(i, j));
    break;
  default:
    ;
  }
  return result;
}

template <typename T> void enlarger2ex<T>::seedPattern(const int& freqidx, const int& stage, const int& row) {
  for(int i = 0; i < B.cols(); i ++) {
    U   ee(0), eth(0);
    int tt;
    if(stage == 2) {
      ee  = exp(- I * Pi * U(i) * U(freqidx) / U(B.cols()));
      eth = exp(- I * Pi *        U(freqidx) / U(B.cols()));
      tt  = i;
    } else if(stage == i % 2) {
      ee  = exp(- I * Pi * U(i - stage) * U(freqidx) / U(B.cols()));
      eth = exp(- I * Pi *                U(freqidx) / U(B.cols()));
      tt  = (i - stage) / 2 * 2;
    } else
      continue;
    const U c0((U(- 2. * tt) * ee + U(tt) - ee) * eth);
    const U c1( U(  2. * tt) * ee - U(tt));
    B(row, tt) = (c0 + c1) * ee / sqrt(T(2) * T(B.cols()));
    G(row, tt) = (c0 - c1) * ee / sqrt(T(2) * T(B.cols()));
  }
  return;
}

template <typename T> void enlarger2ex<T>::initPattern(const int& size, const bool& flag) {
  if(flag) {
    int wp = size / 4 * 4;
    B = MatU(wp, wp);
    G = MatU(wp, wp);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < wp; i ++)
      for(int j = 0; j < wp; j ++) {
        B(i, j) = U(0);
        G(i, j) = U(0);
      }
    int i = 0;
    for(int j = 0; j < wp / 2; j ++)
      seedPattern(j,         2, i ++);
    for(int j = 0; j < wp / 4; j ++)
      seedPattern(j * 2,     0, i ++);
    for(int j = 0; j < wp / 4; j ++)
      seedPattern(j * 2 + 1, 1, i ++);
    MatU B0(B.inverse() * G);
    D = Mat(size, size);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < size; i ++)
      for(int j = 0; j < size; j ++) {
        const int idx = (j - i + size + wp / 4) % size;
        if(0 <= idx && idx < wp)
          D(i, j) = B0(wp / 4, idx).real() / sqrt(T(size * 2));
        else
          D(i, j) = T(0);
      }
    F = Vec(size);
    X = Vec(size);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < size; i ++) {
      U ee( exp(- I * Pi * U(i * size / 2. / size)));
      U eth(exp(- I * Pi * U(    size / 2. / size)));
      U c0((U(- 2. * i) * ee + U(i) - ee) * eth);
      U c1( U(  2. * i) * ee - U(i));
      F[i] = ((c0 + c1) * ee).real() / sqrt(T(size * 2));
      X[i] = ((c0 - c1) * ee).real() / sqrt(T(size * 2));
    }
    F /= T(2) * Pi * size;
    X /= T(2) * Pi * size;
  }
  MatU DFT( size, size);
  MatU IDFT(size, size);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < DFT.rows(); i ++)
    for(int j = 0; j < DFT.cols(); j ++) {
      DFT( i, j) = exp(U(- 2.) * Pi * I * U(i * j / T(size)));
      IDFT(i, j) = exp(U(  2.) * Pi * I * U(i * j / T(size))) / T(size);
    }
  MatU Dbuf(DFT);
  MatU Ibuf(DFT);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < Dbuf.rows(); i ++)
    Dbuf.row(i) *= U(- 2.) * Pi * I * U(i / T(size));
  Ibuf.row(0) *= T(0);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 1; i < Ibuf.rows(); i ++)
    Ibuf.row(i) /= U(- 2.) * Pi * I * U(i / T(size));
  Dop =   (IDFT * Dbuf).real().template cast<T>();
  Iop = - (IDFT * Ibuf).real().template cast<T>();
  return;
}

#define _ENLARGE2X_
#endif

