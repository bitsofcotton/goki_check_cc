#if !defined(_ENLARGE2X_)

#include <Eigen/Core>
#include <Eigen/LU>

using std::complex;

template <typename T> class enlarger2ex {
public:
  typedef enum {
    ENLARGE_X,
    ENLARGE_Y,
    ENLARGE_BOTH } direction_t;
  typedef complex<T> U;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> MatU;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1>              Vec;
  typedef Eigen::Matrix<U, Eigen::Dynamic, 1>              VecU;
  enlarger2ex();
  Mat enlarge2(const Mat& data, const direction_t& dir);
private:
  void seedPattern(const int&, const int&, const int&);
  void initPattern(const int&);
  U    I;
  T    Pi;
  MatU B;
  MatU G;
  Mat  D;
  Mat  Dop;
  Vec F;
  Vec X;
};

template <typename T> enlarger2ex<T>::enlarger2ex() {
  I  = sqrt(U(- 1.));
  Pi = atan2(T(1.), T(1.)) * T(4.);
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> enlarger2ex<T>::enlarge2(const Mat& data, const direction_t& dir) {
  Mat result;
  switch(dir) {
  case ENLARGE_BOTH:
    {
      Mat rw, rh;
      rw = enlarge2(data, ENLARGE_X);
      rh = enlarge2(data, ENLARGE_Y);
      rw = enlarge2(rw, ENLARGE_Y);
      rh = enlarge2(rh, ENLARGE_X);
      result = (rw + rh) / 2.;
    }
    break;
  case ENLARGE_X:
    result = enlarge2(data.transpose(), ENLARGE_Y).transpose();
    break;
  case ENLARGE_Y:
    {
      result = Mat(data.rows() * 2, data.cols());
      cerr << " enlarge_y";
      initPattern(data.rows());
      Mat dd(D * data);
      Vec ff(data.transpose() * F);
      Vec xx(data.transpose() * X);
      Mat co(Dop * data / T(2));
      for(int i = 0; i < data.cols(); i ++) {
        const T lr(data.rows() * data.rows() * 2. * Pi * 2. * Pi);
        ff[i] /= lr * T(2);
        xx[i] /= lr * T(2);
        for(int j = 0; j < data.rows(); j ++) {
          // XXX fixme:
          // const T delta(- co(j, i) * (xx[i] * co(j, i) + ff[i] * dd(j, i)));
          const T delta(- (xx[i] * co(j, i) + ff[i] * dd(j, i)));
          result(j * 2 + 0, i) = data(j, i) - delta;
          result(j * 2 + 1, i) = data(j, i) + delta;
        }
      }
    }
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
    B(row, tt) = (c0 + c1) * ee;
    G(row, tt) = (c0 - c1) * ee;
  }
  return;
}

template <typename T> void enlarger2ex<T>::initPattern(const int& size) {
  int wp = size / 4 * 4;
  B = MatU(wp, wp);
  G = MatU(wp, wp);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for
#endif
  for(int i = 0; i < wp; i ++)
    for(int j = 0; j < wp; j ++) {
      B(i, j) = U(0);
      G(i, j) = U(0);
    }
  int i = 0;
#if defined(_OPENMP)
#pragma omp for
#endif
  for(int j = 0; j < wp / 2; j ++)
    seedPattern(j,         2, i ++);
#if defined(_OPENMP)
#pragma omp for
#endif
  for(int j = 0; j < wp / 4; j ++)
    seedPattern(j * 2,     0, i ++);
#if defined(_OPENMP)
#pragma omp for
#endif
  for(int j = 0; j < wp / 4; j ++)
    seedPattern(j * 2 + 1, 1, i ++);
  MatU B0(B.inverse() * G);
  D = Mat(size, size);
#if defined(_OPENMP)
#pragma omp for
#endif
  for(int i = 0; i < size; i ++)
    for(int j = 0; j < size; j ++) {
      const int idx = (j - i + size + wp / 4) % size;
      if(0 <= idx && idx < wp)
        D(i, j) = B0(wp / 4, idx).real();
      else
        D(i, j) = T(0);
    }
  F = Vec(size);
  X = Vec(size);
#if defined(_OPENMP)
#pragma omp for
#endif
  for(int i = 0; i < size; i ++) {
    U ee(exp(- I * Pi * U(i * size / 2. / size)));
    U eth(exp(- I * Pi * U(size / 2. / size)));
    U c0((U(- 2. * i) * ee + U(i) - ee) * eth);
    U c1( U(  2. * i) * ee - U(i));
    F[i] = ((c0 + c1) * ee).real();
    X[i] = ((c0 - c1) * ee).real();
  }
  MatU Dbuf(size, size), Iop(size, size);
#if defined(_OPENMP)
#pragma omp for
#endif
  for(int i = 0; i < Dbuf.rows(); i ++)
    for(int j = 0; j < Dbuf.cols(); j ++) {
      Dbuf(i, j) = exp(U(- 2.) * Pi * I * U(i * j / T(size)));
      Iop(i, j)  = exp(U(  2.) * Pi * I * U(i * j / T(size)));
    }
#if defined(_OPENMP)
#pragma omp for
#endif
  for(int i = 0; i < Dbuf.rows(); i ++)
    Dbuf.row(i) *= U(- 2.) * Pi * I * U(i / T(size));
  Dop = (Iop * Dbuf).real().template cast<T>();
  return;
}


template <typename T> class enlarger2exds {
public:
  typedef enum {
    ENLARGE_X,
    ENLARGE_Y,
    ENLARGE_BOTH } direction_t;
  typedef complex<T> U;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> MatU;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1>              Vec;
  typedef Eigen::Matrix<U, Eigen::Dynamic, 1>              VecU;
  Mat enlarge2ds(const Mat& data, const direction_t& dir);
};

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> enlarger2exds<T>::enlarge2ds(const Mat& data, const direction_t& dir) {
  Mat result;
  switch(dir) {
  case ENLARGE_X:
    result = enlarge2ds(data.transpose(), ENLARGE_Y).transpose();
    break;
  case ENLARGE_Y:
    {
      Mat buf(data.rows() - 1, data.cols());
      for(int i = 0; i < data.cols(); i ++)
        for(int j = 1; j < data.rows(); j ++)
          buf(j - 1, i) = data(j, i) - data(j - 1, i);
      enlarger2ex<T> enlarger;
      Mat work(enlarger.enlarge2(buf, enlarger2ex<T>::ENLARGE_Y));
      result = Mat(data.rows() * 2 - 1, data.cols());
      result.row(0) = data.row(0);
      for(int i = 0; i < work.cols(); i ++)
        for(int j = 1; j < work.rows(); j ++)
          result(j, i) = result(j - 1, i) + work(j - 1, i) / 2.;
    }
    break;
  case ENLARGE_BOTH:
    {
      Mat rw, rh;
      rw = enlarge2ds(data, ENLARGE_X);
      rh = enlarge2ds(data, ENLARGE_Y);
      rw = enlarge2ds(rw, ENLARGE_Y);
      rh = enlarge2ds(rh, ENLARGE_X);
      result = (rw + rh) / 2.;
    }
    break;
  }
  return result;
};

#define _ENLARGE2X_
#endif

