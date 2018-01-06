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
    ENLARGE_3BOTH,
    DETECT_X,
    DETECT_Y,
    DETECT_BOTH,
    COLLECT_X,
    COLLECT_Y,
    COLLECT_BOTH,
    IDETECT_X,
    IDETECT_Y,
    IDETECT_BOTH } direction_t;
  typedef complex<T> U;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> MatU;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1>              Vec;
  typedef Eigen::Matrix<U, Eigen::Dynamic, 1>              VecU;
  enlarger2ex();
  Mat compute(const Mat& data, const direction_t& dir);
private:
  void initPattern(const int& size, const bool& flag = true);
  U    I;
  T    Pi;
  Mat  D;
  Mat  Dop;
  Mat  Iop;
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
  case ENLARGE_3BOTH:
    {
      result = Mat(data.rows() * 2, data.cols() * 2);
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          result(i, j) = data(i / 2, j / 2);
      result += compute(data, ENLARGE_BOTH) +
                compute(data, ENLARGE_FBOTH) / T(2);
    }
    break;
  case DETECT_BOTH:
    result = (compute(data, DETECT_X) + compute(data, DETECT_Y)) / 2.;
    break;
  case COLLECT_BOTH:
    result = (compute(data, COLLECT_X) + compute(data, COLLECT_Y)) / 2.;
    break;
  case IDETECT_BOTH:
    result = (compute(data, IDETECT_X) + compute(data, IDETECT_Y)) / 2.;
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
  case IDETECT_X:
    result = compute(data.transpose(), IDETECT_Y).transpose();
    break;
  case ENLARGE_Y:
    {
      cerr << " enlarge_y";
      initPattern(data.rows());
      result = D * data;
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
      result = work - work2;
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
  case IDETECT_Y:
    {
      initPattern(data.rows(), false);
      Vec avg(data.cols());
      for(int i = 0; i < data.cols(); i ++) {
        avg[i]  = T(0);
        for(int j = 0; j < data.rows(); j ++)
          avg[i] += data(j, i);
        avg[i] /= data.rows();
      }
      result = Iop * data;
      for(int i = 0; i < data.cols(); i ++)
        for(int j = 0; j < data.rows(); j ++)
          result(j, i) += avg[i] * j * j / data.rows() / data.rows();
    }
    break;
  default:
    ;
  }
  return result;
}

template <typename T> void enlarger2ex<T>::initPattern(const int& size, const bool& flag) {
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
  if(flag) {
    MatU DFTc(DFT.rows(), DFT.cols());
    for(int i = 0; i < DFTc.rows(); i ++) {
      // XXX checkme with enlarge.wxm, DFT space plausible one:
      const int ii(DFTc.rows() - i);
      // This can be tricky, this sees IDFT as DFT and both same theta.
      DFTc.row(i) = sin(U(ii * Pi / DFTc.rows())) / (cos(U(ii * Pi / DFTc.rows())) - U(1)) * (DFT.row(i).imag().template cast<U>() + I * DFT.row(i).real().template cast<U>());
    }
    // This also can be tricky, this sees delta and delta must be smaller
    // but the amount isn't known.
    const Mat Dc((IDFT * (DFT - DFTc)).real().template cast<T>() / pow(T(2) * Pi, T(2)));
    D = Mat(Dc.rows() * 2, Dc.cols());
    for(int i = 0; i < D.rows(); i ++) {
      for(int j = 0; j < D.cols(); j ++)
        D(i, j) = i / 2 == j ? T(1) : T(0);
      // XXX select me: i % 2 or ! (i % 2) .
      if(i % 2)
        D.row(i) += Dc.row(i / 2);
      else
        D.row(i) -= Dc.row(i / 2);
    }
  }
  return;
}

#define _ENLARGE2X_
#endif

