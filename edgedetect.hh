#if !defined(_EDGE_DETECT_DFT_)

#include <Eigen/Core>
#include <Eigen/LU>

using namespace Eigen;

template <typename T, typename U> class edgedetect {
public:
  typedef enum {
    DETECT_X,
    DETECT_Y,
    DETECT_BOTH,
    COLLECT_X,
    COLLECT_Y,
    COLLECT_BOTH } direction_t;
  edgedetect();
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> detect(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>&, const direction_t& dir);
private:
  void initPattern(const int&);
  U    I;
  T    Pi;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Dop;
};

template <typename T, typename U> edgedetect<T,U>::edgedetect() {
  I  = sqrt(U(- 1.));
  Pi = atan2(T(1.), T(1.)) * T(4.);
}

template <typename T, typename U> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> edgedetect<T,U>::detect(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& data, const edgedetect<T,U>::direction_t& dir) {
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result(data);
  switch(dir) {
  case DETECT_BOTH:
    {
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> w(detect(data, DETECT_X));
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> h(detect(data, DETECT_Y));
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
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> w(detect(data, COLLECT_X));
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> h(detect(data, COLLECT_Y));
      result = (w + h) / T(2);
    }
    break;
  case COLLECT_X:
    result = detect(data.transpose(), COLLECT_Y).transpose();
    break;
  case COLLECT_Y:
    {
      result = detect(data, DETECT_Y);
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          result(i, j) = abs(result(i, j));
    }
    break;
  default:
    return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>();
  }
  return result;
}

template <typename T, typename U> void edgedetect<T,U>::initPattern(const int& size) {
  Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> Dbuf(size, size);
  for(int i = 0; i < Dbuf.rows(); i ++)
    for(int j = 0; j < Dbuf.cols(); j ++)
      Dbuf(i, j) = exp(U(- 2.) * Pi * I * U(i * j / T(size)));
  for(int i = 0; i < Dbuf.rows(); i ++)
    Dbuf.row(i) *= U(- 2.) * Pi * I * U(i / T(size));
  Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> Iop(size, size);
  for(int i = 0; i < Dbuf.rows(); i ++)
    for(int j = 0; j < Dbuf.cols(); j ++)
      Iop(i, j) = exp(U(2.) * Pi * I * U(i * j / T(size))) / U(size);
  Dop = (Iop * Dbuf).real().template cast<T>();
  return;
}

#define _EDGE_DETECT_DFT_
#endif

