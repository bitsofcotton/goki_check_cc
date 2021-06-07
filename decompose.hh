/*
BSD 3-Clause License

Copyright (c) 2020-2021, kazunobu watatsu
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#if !defined(_DECOMPOSE_)

template <typename T> class Decompose {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleMatrix<T> Mat;
  inline Decompose();
  inline Decompose(const int& size);
  inline ~Decompose();
  const std::vector<Mat>& A();
  Vec  mimic(const Vec& dst, const Vec& src, const T& intensity = T(1));
  Vec  emphasis(const Vec& dst, const T& intensity = T(1));
  Vec  enlarge(const Vec& in, const int& r = 2);
  Vec  mother(const Vec& in);
  Vec  freq(const Vec& mother, const Vec& in);
  Vec  synth(const Vec& mother, const Vec& freq);
  Mat  represent(const Mat& img, const int& depth = 3);
  Mat  subImage(const Mat& img, const int& x, const int& y, const int& r) const;
private:
  Vec  prepare(const Vec& in, const int& idx = 0) const;
  void apply(Vec& v, const Vec& dst, const Vec& src, const int& idx = 0) const;
  int  flip(const int& x, const int& s) const;
  int  size;
};

template <typename T> inline Decompose<T>::Decompose() {
  size = 0;
}

template <typename T> inline Decompose<T>::Decompose(const int& size) {
  assert(0 < size);
  this->size = size;
}

template <typename T> inline Decompose<T>::~Decompose() {
  ;
}

template <typename T> const std::vector<typename Decompose<T>::Mat>& Decompose<T>::A() {
  static std::vector<std::vector<Mat> > mA;
  if(mA.size() <= size) mA.resize(size + 1, std::vector<Mat>());
  auto& a(mA[size]);
  if(a.size() < size) {
    a.reserve(size);
    for(int i = 0; i < size; i ++) {
      SimpleMatrix<T> AA(size, size);
      for(int j = 0; j < AA.rows(); j ++) {
        const auto jj(T(j) * T(i + 2) / T(size + 1));
        AA.row(j) = taylor<T>(AA.cols(), (jj - floor(jj)) * T(size - 1));
      }
      a.emplace_back(std::move(AA));
    }
    for(int i = 1; i < a.size(); i ++)
      std::swap(a[a.size() - i], a[a.size() - i - 1]);
  }
  return a;
}

template <typename T> typename Decompose<T>::Vec Decompose<T>::mimic(const Vec& dst, const Vec& src, const T& intensity) {
  const auto size2(dst.size() / size);
  const auto size3(src.size() / size);
        auto res(dst);
  for(int i = 0; i < size2; i ++) {
    const auto dd(prepare(dst, i));
    apply(res, synth(mother(prepare(src, i * size3 / size2)),
                     freq(mother(dd), dd)) * intensity +
               dd * (T(1) - intensity), dd, i);
  }
  return res;
}

template <typename T> typename Decompose<T>::Vec Decompose<T>::emphasis(const Vec& dst, const T& intensity) {
  const auto size2(dst.size() / size);
        auto res(dst);
  for(int i = 0; i < size2; i ++) {
    const auto dd(prepare(dst, i));
          auto lfreq(dd);
    for(int j = 0; j < lfreq.size(); j ++)
      lfreq[j] = T(j) / T(lfreq.size());
    apply(res, synth(mother(dd), lfreq) * intensity +
               dd * (T(1) - intensity), dd, i);
  }
  return res;
}

template <typename T> typename Decompose<T>::Vec Decompose<T>::enlarge(const Vec& in, const int& r) {
  assert(0 < r && size == in.size());
  static std::vector<std::vector<Mat> > p;
  if(p.size() <= size)
    p.resize(size + 1, std::vector<Mat>());
  if(p[size].size() <= r)
    p[size].resize(r + 1, Mat());
  auto& pp(p[size][r]);
  if(pp.rows() < size * r) {
    pp.resize(size * r, size);
    for(int i = 0; i < pp.rows(); i ++)
      pp.row(i) = taylor<T>(size, T(i) * T(size - 1) / T(pp.rows() - 1));
  }
  const auto m(mother(in));
  const auto f2(freq(m, in));
  const auto bm(pp * m);
        auto ff(bm);
  for(int i = 0; i < ff.size(); i ++)
    ff[i] = i ? f2[(i % (f2.size() - 1)) + 1] : f2[i];
  Decompose<T> ee(size * r);
  auto result(ee.synth(bm, ff));
  return result *= sqrt(in.dot(in) / result.dot(result) * T(r));
}

template <typename T> typename Decompose<T>::Vec Decompose<T>::prepare(const Vec& in, const int& idx) const {
  const auto cnt(in.size() / size);
  assert(0 < cnt);
  Vec res(size);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < size; i ++) {
    res[i] = T(0);
    const auto ibegin(i * cnt);
    const auto iend(min((i + 1) * cnt, int(in.size()) - idx));
    for(int j = ibegin; j < iend; j ++)
      res[i] += in[j + idx];
    res[i] /= T(iend - ibegin);
  }
  return res;
}

template <typename T> void Decompose<T>::apply(Vec& v, const Vec& dst, const Vec& src, const int& idx) const {
  assert(dst.size() == size && src.size() == size);
  const auto cnt(v.size() / size);
  assert(0 < cnt);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < size; i ++) {
    const auto ratio(dst[i] - src[i]);
    for(int j = i * cnt; j < min((i + 1) * cnt, int(v.size()) - idx); j ++)
      v[j + idx] += ratio;
  }
  return;
}

template <typename T> typename Decompose<T>::Vec Decompose<T>::mother(const Vec& in) {
  std::vector<Mat> A0;
  if(A0.size() <= size) A0.resize(size + 1, Mat());
  auto& a0(A0[size]);
  if(a0.rows() != size || a0.cols() != size) {
    const auto& a(A());
    a0 = a[0];
    for(int i = 1; i < a.size(); i ++)
      a0 += a[i];
  }
  return a0.solve(in);
}

template <typename T> typename Decompose<T>::Vec Decompose<T>::freq(const Vec& mother, const Vec& in) {
  assert(size == mother.size() && size == in.size());
  Mat work(size, size);
  for(int i = 0; i < size; i ++)
    work.setCol(i, A()[i] * mother);
  return work.solve(in);
}

template <typename T> typename Decompose<T>::Vec Decompose<T>::synth(const Vec& mother, const Vec& in) {
  assert(size == mother.size() && size == in.size());
  Vec res(size);
  for(int i = 0; i < size; i ++)
    res[i] = T(0);
  for(int i = 0; i < size; i ++)
    res += A()[i] * mother * in[i];
  return res;
}

template <typename T> typename Decompose<T>::Mat Decompose<T>::represent(const Mat& img, const int& depth) {
  Mat res0(1, size);
  Mat w00(img.rows() - size * 2, size);
  const auto int4(diff<T>(- size));
  const auto int4t(int4.transpose());
  for(int i = size; i < img.rows() - size; i ++) {
    Mat w0(img.cols() - size * 2, size);
    for(int j = size; j < img.cols() - size; j ++) {
      std::vector<Vec> w1;
      for(int r = size;
              r < min(min(i, j),
                    min(img.rows() - i - 1, img.cols() - j - 1));
              r += max(1, min(img.rows() / size, img.cols() / size))) {
        // integrate 4th times because we svd 4 times.
        // svd takes bothside transform, we suppose them as differential op.
        const auto part(int4 * subImage(img, i, j, r) * int4t);
        const auto left(part.SVD() * part);
              Vec  work(left.rows());
        for(int k = 0; k < work.size(); k ++)
          work[k] = sqrt(left.row(k).dot(left.row(k))) + T(1);
        // N.B. normalized singular values on the image with circle region.
        //      If this is flat, the data we have is flat.
        //      If this is edged, the data we have has some data.
        work = mother(work);
        // N.B. enlarging specific bias.
        //      recursive on them.
        w1.emplace_back(std::move(work /= sqrt(work.dot(work))));
      }
      if(! w1.size())
        for(int k = 0; k < w0.cols(); k ++)
          w0(j - size, k) = T(1) / sqrt(T(w0.cols()));
      else if(w1.size() == 1)
        // N.B. line intensity.
        w0.row(j - size) = std::move(w1[0]);
      else {
        Mat w1m(w1.size(), size);
        for(int i = 0; i < w1m.rows(); i ++)
          w1m.row(i) = std::move(w1[i]);
        w1m = w1m.transpose();
        const auto left(w1m.SVD() * w1m);
        for(int k = 0; k < left.rows(); k ++)
          w0(j - size, k) = sqrt(left.row(k).dot(left.row(k))) + T(1);
        w0.row(j - size)  = mother(w0.row(j - size));
        w0.row(j - size) /= sqrt(w0.row(j - size).dot(w0.row(j - size)));
      }
    }
    // N.B. do same on x axis:
    w0 = w0.transpose();
    for(int j = 0; j < w0.rows(); j ++)
      for(int k = 0; k < w0.cols(); k ++)
        assert(isfinite(w0(j, k)) && ! isnan(w0(j, k)));
    const auto left(w0.SVD() * w0);
    for(int k = 0; k < left.rows(); k ++)
      w00(i - size, k) = sqrt(left.row(k).dot(left.row(k))) + T(1);
    w00.row(i - size)  = mother(w00.row(i - size));
    w00.row(i - size) /= sqrt(w00.row(i - size).dot(w00.row(i - size)));
  }
  // N.B. do same on whole image:
  w00 = w00.transpose();
  const auto left(w00.SVD() * w00);
  for(int k = 0; k < left.rows(); k ++)
    res0(0, k) = sqrt(left.row(k).dot(left.row(k))) + T(1);
  res0.row(0)  = mother(res0.row(0));
  res0.row(0) /= sqrt(res0.row(0).dot(res0.row(0)));
  // N.B. recursive on them.
  if(0 < depth && size * 4 <= min(img.rows(), img.cols()) / 2) {
    Mat dimg[5];
    for(int i = 0; i < 5; i ++)
      dimg[i] = Mat(img.rows() / 2, img.cols() / 2);
    for(int i = 0; i < dimg[0].rows(); i ++)
      for(int j = 0; j < dimg[0].cols(); j ++) {
        dimg[0](i, j) = img(i, j);
        dimg[1](i, j) = img(i - dimg[0].rows() + img.rows(), j);
        dimg[2](i, j) = img(i, j - dimg[0].cols() + img.cols());
        dimg[3](i, j) = img(i - dimg[0].rows() + img.rows(),
                            j - dimg[0].cols() + img.cols());
        dimg[4](i, j) = img(i + (img.rows() - dimg[0].rows()) / 2,
                            j + (img.cols() - dimg[0].cols()) / 2);
      }
    Mat dres[5];
    for(int i = 0; i < 5; i ++)
      dres[i] = represent(dimg[i], depth - 1);
    Mat res(1 + dres[0].rows() * 5, size);
    res.row(0) = res0.row(0);
    for(int i = 0; i < 5; i ++)
      for(int j = 0; j < dres[i].rows(); j ++)
        res.row(1 + i * dres[i].rows() + j) = dres[i].row(j);
    return res;
  }
  return res0;
}

template <typename T> typename Decompose<T>::Mat Decompose<T>::subImage(const Mat& img, const int& x, const int& y, const int& r) const {
  Mat res(size, size);
  for(int i = 0; i < res.rows(); i ++)
    for(int j = 0; j < res.cols(); j ++) {
      const auto rr(T(j + 1) / T(res.cols()) * T(r));
      const auto th(T(i) / T(res.rows()) * T(2) * T(4) * atan2(T(1), T(1)));
      res(i, j) = img(flip(x + int(rr * cos(th)), img.rows()),
                      flip(y + int(rr * sin(th)), img.cols()));
    }
  return res;
}

template <typename T> int Decompose<T>::flip(const int& x, const int& s) const {
  const int xx(abs(x % (s * 2)));
  return s <= xx ? s * 2 - xx - 1 : xx;
}

#define _DECOMPOSE_
#endif

