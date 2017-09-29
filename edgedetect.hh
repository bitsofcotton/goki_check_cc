#if !defined(_EDGE_DETECT_DFT_)

#include <Eigen/Core>
#include <Eigen/LU>

using std::complex;
using std::abs;

template <typename T> class edgedetect {
public:
  typedef enum {
    DETECT_X,
    DETECT_Y,
    DETECT_BOTH,
    COLLECT_X,
    COLLECT_Y,
    COLLECT_BOTH } direction_t;
  edgedetect();
  typedef complex<T> U;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> MatU;
  Mat detect(const Mat& data, const direction_t& dir);
private:
  void initPattern(const int&);
  U    I;
  T    Pi;
  Mat  Dop;
};

template <typename T> edgedetect<T>::edgedetect() {
  I  = sqrt(U(- 1.));
  Pi = atan2(T(1.), T(1.)) * T(4.);
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> edgedetect<T>::detect(const Mat& data, const direction_t& dir) {
  Mat result(data);
  switch(dir) {
  case DETECT_BOTH:
    {
      Mat w(detect(data, DETECT_X));
      Mat h(detect(data, DETECT_Y));
      result = (w + h) / T(2);
    }
    break;
  case DETECT_X:
    result = detect(data.transpose(), DETECT_Y).transpose();
    break;
  case DETECT_Y:
    initPattern(data.rows());
    result = Dop * data;
    break;
  case COLLECT_BOTH:
    {
      Mat w(detect(data, COLLECT_X));
      Mat h(detect(data, COLLECT_Y));
      result = (w + h) / T(2);
    }
    break;
  case COLLECT_X:
    result = detect(data.transpose(), COLLECT_Y).transpose();
    break;
  case COLLECT_Y:
    {
      result = detect(data, DETECT_Y);
#if defined(_OPENMP)
#pragma omp parallel for
#endif
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          result(i, j) = abs(result(i, j));
    }
    break;
  default:
    ;
  }
  return result;
}

template <typename T> void edgedetect<T>::initPattern(const int& size) {
  MatU Dbuf(size, size), Iop(size, size);
#if defined(_OPENMP)
#pragma omp parallel
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

#define _EDGE_DETECT_DFT_
#endif

