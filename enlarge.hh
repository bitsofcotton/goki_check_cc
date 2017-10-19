#if !defined(_ENLARGE2X_)

#include <Eigen/Core>
#include <Eigen/LU>

using std::complex;
using std::abs;
using std::pow;
using std::vector;
using std::sort;

template <typename T> T autoLevel(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>* data, const int& size, int npad = 0) {
  if(npad <= 0)
    npad = (data[0].rows() + data[0].cols()) * 8;
  vector<T> stat;
  for(int i = 0; i < size; i ++)
    for(int j = 0; j < data[i].rows(); j ++)
      for(int k = 0; k < data[i].cols(); k ++)
        stat.push_back(data[i](j, k));
  sort(stat.begin(), stat.end());
  const T mm(stat[npad]);
  T MM(stat[stat.size() - 1 - npad]);
  if(MM == mm)
    MM = mm + 1.;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int k = 0; k < size; k ++) {
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> buf(data[k].rows(), data[k].cols());
    for(int i = 0; i < buf.rows(); i ++)
      for(int j = 0; j < buf.cols(); j ++)
        buf(i, j) = (max(min(data[k](i, j), MM), mm) - mm) / (MM - mm);
    data[k] = buf;
  }
  return stat[stat.size() / 2];
}

template <typename T> class enlarger2ex {
public:
  typedef enum {
    ENLARGE_X,
    ENLARGE_Y,
    ENLARGE_D,
    ENLARGE_DD,
    ENLARGE_BOTH,
    ENLARGE_QUAD } direction_t;
  typedef complex<T> U;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> MatU;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1>              Vec;
  typedef Eigen::Matrix<U, Eigen::Dynamic, 1>              VecU;
  enlarger2ex();
  Mat enlarge2(const Mat& data, const direction_t& dir);
  Mat normQuad(const Mat& rw, const Mat& rh, const Mat& rdr, const Mat& rdl);
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

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> enlarger2ex<T>::normQuad(const Mat& rw, const Mat& rh, const Mat& rdr, const Mat& rdl) {
  Mat result(rdr.rows(), rdr.cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 1; i < rw.rows() - 1; i ++)
    for(int j = 1; j < rh.cols() - 1; j ++) {
      Mat A(4, 4);
      Vec b(4);
      A(0, 0) = rw(i, j * 2)     / T(2);
      A(0, 1) = rw(i, j * 2 + 1) / T(2);
      A(0, 2) = rw(i, j * 2)     / T(2);
      A(0, 3) = rw(i, j * 2 + 1) / T(2);
      A(1, 0) = rh(i * 2,     j) / T(2);
      A(1, 1) = rh(i * 2,     j) / T(2);
      A(1, 2) = rh(i * 2 + 1, j) / T(2);
      A(1, 3) = rh(i * 2 + 1, j) / T(2);
      A(2, 0) = rdr(i * 2,    j * 2);
      A(2, 1) = T(0);
      A(2, 2) = T(0);
      A(2, 3) = rdr(i * 2 + 1, j * 2 + 1);
      A(3, 0) = T(0);
      A(3, 1) = rdl(i * 2, j * 2 + 1);
      A(3, 2) = rdl(i * 2 + 1, j * 2);
      A(3, 3) = T(0);
      A /= T(2);
      const T orig((rw(i, j * 2) + rw(i, j * 2 + 1)) / T(2));
      b[0] = b[1] = b[2] = b[3] = T(1);
      b = A.transpose() * b / T(4);
      result(i * 2 + 0, j * 2 + 0) = b[0];
      result(i * 2 + 0, j * 2 + 1) = b[1];
      result(i * 2 + 1, j * 2 + 0) = b[2];
      result(i * 2 + 1, j * 2 + 1) = b[3];
    }
  return result;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> enlarger2ex<T>::enlarge2(const Mat& data, const direction_t& dir) {
  Mat result;
  switch(dir) {
  case ENLARGE_QUAD:
    result = normQuad(enlarge2(data, ENLARGE_X),
                      enlarge2(data, ENLARGE_Y),
                      enlarge2(data, ENLARGE_D),
                      enlarge2(data, ENLARGE_DD));
    break;
  case ENLARGE_BOTH:
    {
      Mat rw, rh;
      rw = enlarge2(data, ENLARGE_X);
      rh = enlarge2(data, ENLARGE_Y);
      rw = enlarge2(rw,   ENLARGE_Y);
      rh = enlarge2(rh,   ENLARGE_X);
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
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int i = 0; i < data.cols(); i ++) {
        const T lr(pow(T(2) * Pi * data.rows(), T(2)));
        ff[i] /= lr * T(2);
        xx[i] /= lr * T(2);
        for(int j = 0; j < data.rows(); j ++) {
          const T delta(xx[i] * co(j, i) + ff[i] * dd(j, i));
          result(j * 2 + 0, i) = data(j, i) - delta;
          result(j * 2 + 1, i) = data(j, i) + delta;
        }
      }
    }
    break;
  case ENLARGE_D:
    {
      result = Mat(data.rows() * 2, data.cols() * 2);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          result(i, j) = T(0);
      cerr << " enlarge_d";
      const int width0(min(data.rows(), data.cols()));
      initPattern(width0);
/*
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
*/
      for(int i = 0; i < data.cols(); i ++) {
        const int width(min(data.cols() - i, data.rows()));
        Vec dbuf(width0);
        for(int j = width0 / 2 - width / 2, jj = 0; j < dbuf.size(); j ++, jj ++)
          dbuf[j] = data( jj      % width,
                         (jj + i) % width);
        for(int j = width0 / 2 - width / 2, jj = 0; j >= 0; j --, jj --)
          dbuf[j] = data((jj + width * width0)     % width,
                         (jj + width * width0 + i) % width);
        Vec dd(D * dbuf);
        T   ff(dbuf.transpose() * F);
        T   xx(dbuf.transpose() * X);
        Vec co(Dop * dbuf / T(2));
        const T lr(pow(T(2) * Pi * width0, T(2)));
        ff /= lr * T(2);
        xx /= lr * T(2);
        for(int j = 0; j < width; j ++) {
          const int jj(width0 / 2 - width / 2 + j);
          const T delta(xx * co[jj] + ff * dd[jj]);
          if(j * 2 < result.rows() && i * 2 + j * 2 + 0 < result.cols())
            result(j * 2 + 0, i * 2 + j * 2 + 0) = data(j, i + j) - delta;
          if(j * 2 + 1 < result.rows() && i * 2 + j * 2 + 1 < result.cols())
            result(j * 2 + 1, i * 2 + j * 2 + 1) = data(j, i + j) + delta;
        }
      }
/*
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
*/
      for(int i = 0; i < data.rows(); i ++) {
        const int width(min(data.cols(), data.rows() - i));
        Vec dbuf(width0);
        for(int j = width0 / 2 - width / 2, jj = 0; j < dbuf.size(); j ++, jj ++)
          dbuf[j] = data((jj + i) % width,
                          jj      % width);
        for(int j = width0 / 2 - width / 2, jj = 0; j >= 0; j --, jj --)
          dbuf[j] = data((jj + width * width0 + i) % width,
                         (jj + width * width0)     % width);
        Vec dd(D * dbuf);
        T   ff(dbuf.transpose() * F);
        T   xx(dbuf.transpose() * X);
        Vec co(Dop * dbuf / T(2));
        const T lr(pow(T(2) * Pi * width0, T(2)));
        ff /= lr * T(2);
        xx /= lr * T(2);
        for(int j = 0; j < width; j ++) {
          const int jj(width0 / 2 - width / 2 + j);
          const T delta(xx * co[jj] + ff * dd[jj]);
          if(i * 2 + j * 2 < result.rows() && j * 2 + 0 < result.cols())
            result(i * 2 + j * 2 + 0, j * 2 + 0) = data(i + j, j) - delta;
          if(i * 2 + j * 2 + 1 < result.rows() && j * 2 + 1 < result.cols())
            result(i * 2 + j * 2 + 1, j * 2 + 1) = data(i + j, j) + delta;
        }
      }
    }
    break;
  case ENLARGE_DD:
    {
      Mat buf(data.rows(), data.cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int i = 0; i < data.cols(); i ++)
        buf.col(i) = data.col(data.cols() - 1 - i);
      buf = enlarge2(buf, ENLARGE_D);
      result = Mat(buf.rows(), buf.cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int i = 0; i < buf.cols(); i ++)
        result.col(i) = buf.col(buf.cols() - 1 - i);
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
        D(i, j) = B0(wp / 4, idx).real();
      else
        D(i, j) = T(0);
    }
  F = Vec(size);
  X = Vec(size);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
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
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < Dbuf.rows(); i ++)
    for(int j = 0; j < Dbuf.cols(); j ++) {
      Dbuf(i, j) = exp(U(- 2.) * Pi * I * U(i * j / T(size)));
      Iop(i, j)  = exp(U(  2.) * Pi * I * U(i * j / T(size)));
    }
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
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
    ENLARGE_D,
    ENLARGE_DD,
    ENLARGE_BOTH,
    ENLARGE_QUAD } direction_t;
  typedef complex<T> U;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> MatU;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1>              Vec;
  typedef Eigen::Matrix<U, Eigen::Dynamic, 1>              VecU;
  Mat enlarge2ds(const Mat& data, const direction_t& dir, const bool diff = false);
};

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> enlarger2exds<T>::enlarge2ds(const Mat& data, const direction_t& dir, const bool diff) {
  Mat result;
  const T Pi(T(4) * atan2(T(1), T(1)));
  switch(dir) {
  case ENLARGE_X:
    result = enlarge2ds(data.transpose(), ENLARGE_Y, diff).transpose();
    break;
  case ENLARGE_Y:
    {
      Mat buf(data.rows() - 1, data.cols());
      for(int i = 1; i < data.rows(); i ++)
        buf.row(i - 1) = data.row(i) - data.row(i - 1);
      enlarger2ex<T> enlarger;
      Mat work(enlarger.enlarge2(buf, enlarger2ex<T>::ENLARGE_Y));
      if(diff)
        result = work;
      else {
        result = Mat(data.rows() * 2 - 1, data.cols());
        result.row(0) = data.row(0);
        for(int i = 1; i < work.rows() + 1; i ++)
          result.row(i) = result.row(i - 1) + work.row(i - 1) / T(2);
      }
    }
    break;
  case ENLARGE_D:
    {
      Mat buf(data.rows() - 1, data.cols() - 1);
      for(int i = 1; i < data.cols(); i ++)
        for(int j = 1; j < min(data.rows(), data.cols() - i); j ++)
          buf(j - 1, i - 1 + j - 1) = data(j, i + j - 1) - data(j - 1, i + j - 1 - 1);
      for(int i = 1; i < data.rows(); i ++)
        for(int j = 1; j < min(data.cols(), data.rows() - i); j ++)
          buf(i - 1 + j - 1, j - 1) = data(i + j - 1, j) - data(i + j - 1 - 1, j - 1);
      enlarger2ex<T> enlarger;
      buf = enlarger.enlarge2(buf, enlarger2ex<T>::ENLARGE_D);
      if(diff)
        result = buf;
      else {
        result = Mat(data.rows() * 2 - 1, data.cols() * 2 - 1);
        for(int i = 0; i < result.rows(); i ++)
          for(int j = 0; j < result.cols(); j ++)
            result(i, j) = T(0);
        for(int i = 0; i < result.cols(); i ++)
          result(0, i) = data(0, i / 2);
        for(int i = 0; i < result.rows(); i ++)
          result(i, 0) = data(i / 2, 0);
        for(int i = 1; i < buf.cols(); i ++)
          for(int j = 1; j < min(buf.rows(), buf.cols() - i); j ++)
            result(j, i + j - 1) = result(j - 1, i + j - 1 - 1) + buf(j - 1, i + j - 1 - 1) / T(2);
        for(int i = 1; i < buf.rows(); i ++)
          for(int j = 1; j < min(buf.cols(), buf.rows() - i); j ++)
            result(i + j - 1, j) = result(i + j - 1 - 1, j - 1) + buf(i + j - 1 - 1, j - 1) / T(2);
      }
    }
    break;
  case ENLARGE_DD:
    {
      Mat buf(data.rows(), data.cols());
      for(int i = 0; i < data.cols(); i ++)
        buf.col(i) = data.col(data.cols() - 1 - i);
      buf = enlarge2ds(buf, ENLARGE_D, diff);
      result = Mat(buf.rows(), buf.cols());
      for(int i = 0; i < buf.cols() - 1; i ++)
        result.col(i) = buf.col(buf.cols() - 2 - i);
      result.col(buf.cols() - 1) *= T(0);
    }
    break;
  case ENLARGE_BOTH:
    {
      Mat rw, rh;
      rw = enlarge2ds(data, ENLARGE_X, false);
      rh = enlarge2ds(data, ENLARGE_Y, false);
      rw = enlarge2ds(rw,   ENLARGE_Y, false);
      rh = enlarge2ds(rh,   ENLARGE_X, false);
      result = (rw + rh) / 2.;
    }
    break;
  case ENLARGE_QUAD:
    {
      enlarger2ex<T> enlarger;
      Mat buf(enlarger.normQuad(enlarge2ds(data, ENLARGE_X,  true),
                                enlarge2ds(data, ENLARGE_Y,  true),
                                enlarge2ds(data, ENLARGE_D,  true),
                                enlarge2ds(data, ENLARGE_DD, true)));
      result = Mat(data.rows() * 2 - 1, data.cols() * 2 - 1);
      for(int i = 0; i < result.cols(); i ++)
        result(0, i) = data(0, i / 2);
      for(int i = 0; i < result.rows(); i ++)
        result(i, 0) = data(i / 2, 0);
      for(int i = 0; i < data.rows(); i ++)
        for(int j = 0; j < data.cols(); j ++) 
          if(i * 2 < buf.rows() && j * 2 < buf.cols()) {
            result(i * 2 + 0, j * 2 + 0) = data(i, j) + buf(i * 2 + 0, j * 2 + 0);
            if(j * 2 + 1 < buf.cols())
              result(i * 2 + 0, j * 2 + 1) = result(i * 2, j * 2) + buf(i * 2 + 0, j * 2 + 1);
            if(i * 2 + 1 < buf.rows())
              result(i * 2 + 1, j * 2 + 0) = result(i * 2, j * 2) + buf(i * 2 + 1, j * 2 + 0);
            if(i * 2 + 1 < buf.rows() && j * 2 + 1 < buf.cols())
              result(i * 2 + 1, j * 2 + 1) = (result(i * 2, j * 2) + result(i * 2, j * 2 + 1) + result(i * 2 + 1, j * 2)) / T(3) + buf(i * 2 + 1, j * 2 + 1);
        }
    }
    break;
  }
  return result;
};

#define _ENLARGE2X_
#endif

