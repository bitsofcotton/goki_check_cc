#if !defined(_SCAN_CONTEXT_)

#include <Eigen/Core>
#include <Eigen/LU>

using namespace Eigen;

template <typename T> class lowFreq {
public:
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
  typedef complex<T> U;
  typedef Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> MatU;
  lowFreq();
  ~lowFreq();
  void init(const T& thresh);
  
  std::vector<Eigen::Matrix<T, 3, 1> > getLowFreq(const Mat& data);
  Mat getLowFreqImage(const Mat& data);
private:
  Mat prepareCost(const Mat& data);
  T thresh;
  T Pi;
  U I;
};

template <typename T> lowFreq<T>::lowFreq() {
  Pi     = atan2(T(1), T(1)) * T(4);
  I      = std::sqrt(U(- 1));
  thresh = T(.0025);
}

template <typename T> lowFreq<T>::~lowFreq() {
}

template <typename T> void lowFreq<T>::init(const T& thresh) {
  this->thresh = thresh;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> lowFreq<T>::getLowFreqImage(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& data) {
  Mat result(data);
  Mat cost(prepareCost(data));
  for(int i = 1; i < data.rows() - 1; i ++)
    for(int j = 1; j < data.cols() - 1; j ++)
      if(cost(i, j) < thresh) {
        const Mat& work(result);
        result(i, j)  = work(i - 1, j - 1) + work(i - 1, j) + work(i - 1, j + 1);
        result(i, j) += work(i, j - 1)     + work(i, j + 1);
        result(i, j) += work(i + 1, j - 1) + work(i + 1, j) + work(i + 1, j + 1);
        result(i, j) /= T(8);
        result(i, j)  = - T(1);
      }
  return result;
}

template <typename T> std::vector<Eigen::Matrix<T, 3, 1> > lowFreq<T>::getLowFreq(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& data) {
  std::vector<Eigen::Matrix<T, 3, 1> > result;
  Mat cost(prepareCost(data));
  for(int i = 1; i < data.rows() - 1; i ++)
    for(int j = 1; j < data.cols() - 1; j ++)
      if(! (cost(i, j) < thresh)) {
        Eigen::Matrix<T, 3, 1> work;
        work[0] = i;
        work[1] = j;
        work[2] = data(i, j);
        result.push_back(work);
      }
  return result;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> lowFreq<T>::prepareCost(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& data) {
  T   norm2(0);
  Mat cost(data.rows(), data.cols());
  for(int i = 1; i < data.rows() - 1; i ++)
    for(int j = 1; j < data.cols() - 1; j ++) {
      cost(i, j)  = std::abs(data(i - 1, j - 1) - data(i, j));
      cost(i, j) += std::abs(data(i - 1, j)     - data(i, j));
      cost(i, j) += std::abs(data(i - 1, j + 1) - data(i, j));
      cost(i, j) += std::abs(data(i, j - 1) - data(i, j));
      cost(i, j) += std::abs(data(i, j + 1) - data(i, j));
      cost(i, j) += std::abs(data(i + 1, j - 1) - data(i, j));
      cost(i, j) += std::abs(data(i + 1, j)     - data(i, j));
      cost(i, j) += std::abs(data(i + 1, j + 1) - data(i, j));
      cost(i, j) /= T(8);
      norm2      += cost(i, j) * cost(i, j);
    }
  norm2 = std::sqrt(norm2);
  return cost / norm2;
}

template <typename T> class matchPartial {
public:
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T, 3, 3>                           Mat3x3;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
  typedef Eigen::Matrix<T, 3, 1>              Vec3;
  typedef complex<T> U;
  typedef Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> MatU;
  matchPartial();
  ~matchPartial();
  void init(const std::vector<Vec3>& shapebase);
  
  std::vector<Mat3x3> match(const std::vector<std::vector<Vec3> >& points);
private:
  U I;
  T Pi;
  std::vector<std::vector<Vec3> > shapebase;
};

#define _SCAN_CONTEXT_
#endif

