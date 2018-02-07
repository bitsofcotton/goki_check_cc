#if !defined(_ENLARGE2X_)

#include <Eigen/Core>
#include <Eigen/LU>

using std::complex;
using std::abs;
using std::sqrt;
using std::exp;
using std::conj;

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
  void initPattern(const int& size);
  U    I;
  T    Pi;
  Mat  D;
  Mat  Dop;
  Mat  Iop;
  Mat  bD;
  Mat  bDop;
  Mat  bIop;
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
    result = compute(data, ENLARGE_BOTH) +
             compute(data, ENLARGE_FBOTH) / T(2);
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
    cerr << " enlarge_y";
    initPattern(data.rows());
    result = D * data;
    break;
  case ENLARGE_FY:
    {
      initPattern(data.rows());
      const Mat dcache(Dop * data);
      const Mat work(compute(dcache, ENLARGE_Y));
            Mat work2(data.rows() * 2, data.cols());
      for(int j = 0; j < work2.rows(); j ++)
        work2.row(j) = dcache.row(j / 2);
      result = work - work2;
    }
    break;
  case DETECT_Y:
    initPattern(data.rows());
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
      initPattern(data.rows());
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

template <typename T> void enlarger2ex<T>::initPattern(const int& size) {
  if(Dop.rows() == size)
    return;
  if(bDop.rows() == size) {
    Mat work(bD);
    bD   = D;
    D    = work;
    work = bDop;
    bDop = Dop;
    Dop  = work;
    work = bIop;
    bIop = Iop;
    Iop  = work;
    return;
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
  {
    MatU DFTa(DFT);
    MatU DFTb(DFT);
    MatU DFTRa(DFT);
    MatU DFTRb(DFT);
    for(int i = 0; i < DFT.rows(); i ++) {
      // XXX checkme with enlarge.wxm, DFT space plausible one:
      const T ii(i + .5);
      const U r(sin(U(ii * Pi / DFTb.rows())) / (cos(U(ii * Pi / DFTb.rows())) - U(1)));
      // This can be tricky, this sees IDFT as DFT and both same theta.
      DFTa.row(i)  = (T(1) - r) * I * DFT.row(i).conjugate();
      DFTb.row(i)  =         r  * I * DFT.row(i).conjugate();
      DFTRa.row(i) = (T(1) - r) * I * DFT.row(DFT.rows() - i - 1).conjugate();
      DFTRb.row(i) =         r  * I * DFT.row(DFT.rows() - i - 1).conjugate();
    }
    MatU DFTRRa(DFTRa), DFTRRb(DFTRb);
    for(int i = 0; i < DFTRa.rows(); i ++) {
      DFTRRa.row(i) = DFTRa.row(DFTRa.rows() - i - 1);
      DFTRRb.row(i) = DFTRb.row(DFTRb.rows() - i - 1);
    }
    // This also can be tricky, this sees delta of b(t).
    const Mat Da((IDFT * (DFTa + DFTRRa)).real().template cast<T>() / (T(2) * Pi * sqrt(T(DFT.rows())) * sqrt(T(8))) / T(2));
    const Mat Db((IDFT * (DFTb + DFTRRb)).real().template cast<T>() / (T(2) * Pi * sqrt(T(DFT.rows())) * sqrt(T(8))) / T(2));
    D = Mat(DFT.rows() * 2, DFT.cols());
    for(int i = 0; i < DFT.rows(); i ++) {
      D.row(i * 2 + 0) = Da.row(i);
      D.row(i * 2 + 1) = Db.row(i);
      D(i * 2 + 0, i) += T(1);
      D(i * 2 + 1, i) += T(1);
    }
  }
  return;
}

#define _ENLARGE2X_
#endif

