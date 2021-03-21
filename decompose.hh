/*
BSD 3-Clause License

Copyright (c) 2020, kazunobu watatsu
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
         Vec  mimic(const Vec& dst, const Vec& src, const T& intensity = T(1)) const;
         Vec  emphasis(const Vec& dst, const T& intensity = T(1)) const;
         Vec  enlarge(const Vec& in, const int& r = 2) const;
         Vec  mother(const Vec& in) const;
         Vec  freq(const Vec& mother, const Vec& in) const;
         Vec  synth(const Vec& mother, const Vec& freq) const;
         Mat  decompose(const Mat& img, const int& depth = 3) const;
         Mat  subImage(const Mat& img, const int& x, const int& y, const int& r) const;
private:
  std::vector<Mat> A;
         Mat  A0;
         Vec  prepare(const Vec& in, const int& idx = 0) const;
         void apply(Vec& v, const Vec& dst, const Vec& src, const int& idx = 0) const;
         int  flip(const int& x, const int& s) const;
};

template <typename T> inline Decompose<T>::Decompose() {
  ;
}

template <typename T> inline Decompose<T>::Decompose(const int& size) {
  static P0<T> p0;
  for(int i = 0; i < size; i ++) {
    SimpleMatrix<T> AA(size, size);
    for(int j = 0; j < AA.rows(); j ++) {
      const auto jj(T(j) * T(i + 2) / T(size + 1));
      AA.row(j) = p0.taylor(AA.cols(), (jj - floor(jj)) * T(size - 1));
    }
    A.emplace_back(AA);
  }
  for(int i = 1; i < A.size(); i ++)
    std::swap(A[A.size() - i], A[A.size() - i - 1]);
  A0 = A[0];
  for(int i = 1; i < A.size(); i ++)
    A0 += A[i];
}

template <typename T> inline Decompose<T>::~Decompose() {
  ;
}

template <typename T> typename Decompose<T>::Vec Decompose<T>::mimic(const Vec& dst, const Vec& src, const T& intensity) const {
  const int  size(A.size());
  const auto size2(dst.size() / size);
  const auto size3(src.size() / size);
        auto res(dst);
  for(int i = 0; i < size2; i ++) {
    const auto dd(prepare(dst, i));
    apply(res, synth(mother(prepare(src, i * size3 / size2)),
                     freq(mother(dd), dd)) * intensity +
               dd * (T(1) - abs(intensity)), dd, i);
  }
  return res;
}

template <typename T> typename Decompose<T>::Vec Decompose<T>::emphasis(const Vec& dst, const T& intensity) const {
  const int  size(A.size());
  const auto size2(dst.size() / size);
        auto res(dst);
  for(int i = 0; i < size2; i ++) {
    const auto dd(prepare(dst, i));
          auto lfreq(dd);
    lfreq[lfreq.size() - 1] = T(0);
    for(int j = 0; j < lfreq.size() - 1; j ++)
      lfreq[j] = T(j + 1) / T(lfreq.size() - 1);
    apply(res, synth(mother(dd), lfreq) * intensity, dd, i);
  }
  return res;
}

template <typename T> typename Decompose<T>::Vec Decompose<T>::enlarge(const Vec& in, const int& r) const {
  assert(0 < r);
  static P0<T> p0;
  static std::vector<std::vector<Mat> > p;
  static std::vector<Decompose<T> > e;
  if(p.size() <= in.size())
    p.resize(in.size() + 1, std::vector<Mat>());
  if(p[in.size()].size() <= r)
    p[in.size()].resize(r + 1, Mat());
  if(e.size() <= in.size() * r)
    e.resize(in.size() * r + 1, Decompose<T>());
  auto& pp(p[in.size()][r]);
  auto& ee(e[in.size() * r]);
  if(pp.rows() < in.size() * r) {
    pp.resize(in.size() * r, in.size());
    for(int i = 0; i < pp.rows(); i ++)
      pp.row(i) = p0.taylor(in.size(), T(i) * T(in.size() - 1) / T(pp.rows() - 1));
    if(! ee.A.size()) ee = Decompose<T>(in.size() * r);
  }
  const auto m(mother(in));
  const auto f2(freq(m, in));
  const auto bm(pp * m);
        auto ff(bm);
  for(int i = 0; i < ff.size(); i ++)
    ff[i] = i ? f2[(i % (f2.size() - 1)) + 1] : f2[i];
  auto result(ee.synth(bm, ff));
  return result *= sqrt(in.dot(in) / result.dot(result) * T(r));
}

template <typename T> typename Decompose<T>::Vec Decompose<T>::prepare(const Vec& in, const int& idx) const {
  const int  size(A.size());
  const auto cnt(int(in.size() + size - 1) / size - 1);
  assert(0 < cnt);
  Vec res(size);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < res.size(); i ++) {
    res[i] = T(0);
    const auto ibegin(i * cnt);
    const auto iend(i < size - 1 ? min((i + 1) * cnt, int(in.size()) - idx)
                                 : int(in.size()) - idx);
    for(int j = ibegin; j < iend; j ++)
      res[i] += in[j + idx];
    res[i] /= T(iend - ibegin);
  }
  return res;
}

template <typename T> void Decompose<T>::apply(Vec& v, const Vec& dst, const Vec& src, const int& idx) const {
  const int  size(A.size());
  assert(dst.size() == size && src.size() == size);
  const auto cnt(int(v.size() + size - 1) / size - 1);
  assert(0 < cnt);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < size; i ++) {
    const auto ratio(dst[i] - src[i]);
    for(int j = i * cnt;
            j < (i < size - 1 ? min((i + 1) * cnt, int(v.size()) - idx)
                              : int(v.size()) - idx);
            j ++)
      v[j + idx] += ratio;
  }
  return;
}

template <typename T> typename Decompose<T>::Vec Decompose<T>::mother(const Vec& in) const {
  assert(in.size() && A.size() == in.size() &&
         A[0].rows() == in.size() && A[0].cols() == in.size());
  return A0.solve(in);
}

template <typename T> typename Decompose<T>::Vec Decompose<T>::freq(const Vec& mother, const Vec& in) const {
  assert(mother.size() && A.size() == mother.size() &&
         A[0].cols() == mother.size() && mother.size() == in.size());
  Mat work(in.size(), in.size());
  for(int i = 0; i < mother.size(); i ++)
    work.setCol(i, A[i] * mother);
  return work.solve(in);
}

template <typename T> typename Decompose<T>::Vec Decompose<T>::synth(const Vec& mother, const Vec& in) const {
  assert(mother.size() && A.size() == mother.size() &&
         A[0].cols() == mother.size() && mother.size() == in.size());
  Vec res(mother.size());
  for(int i = 0; i < res.size(); i ++)
    res[i] = T(0);
  for(int i = 0; i < A.size(); i ++)
    res += A[i] * mother * in[i];
  return res;
}

template <typename T> typename Decompose<T>::Mat Decompose<T>::decompose(const Mat& img, const int& depth) const {
  static P0<T> p;
  Mat res0(1, A.size());
  Mat w00(img.rows() - A.size() * 2, A.size());
  for(int i = A.size(); i < img.rows() - A.size(); i ++) {
    Mat w0(img.cols() - A.size() * 2, A.size());
    for(int j = A.size(); j < img.cols() - A.size(); j ++) {
      std::vector<Vec> w1;
      for(int r = A.size();
              r < min(min(i, j),
                    min(img.rows() - i - 1, img.cols() - j - 1));
              r += min(img.rows() / A.size(), img.cols() / A.size())) {
        const auto part(subImage(img, i, j, r));
        const auto left(part.LSVD().transpose() * part);
              Vec  work(left.rows());
        for(int k = 0; k < work.size(); k ++)
          work[k] = sqrt(left.row(k).dot(left.row(k))) + T(1);
        w1.emplace_back(std::move(work /= sqrt(work.dot(work))));
      }
      if(! w1.size())
        for(int k = 0; k < w0.cols(); k ++)
          w0(j - A.size(), k) = T(1) / sqrt(T(w0.cols()));
      else if(w1.size() == 1)
        w0.row(j - A.size()) = std::move(w1[0]);
      else {
        Mat w1m(w1.size(), A.size());
        for(int i = 0; i < w1m.rows(); i ++)
          w1m.row(i) = std::move(w1[i]);
        w1m = w1m.transpose();
        const auto left(w1m.LSVD().transpose() * w1m);
        for(int k = 0; k < left.rows(); k ++)
          w0(j - A.size(), k) = sqrt(left.row(k).dot(left.row(k))) + T(1);
        w0.row(j - A.size())  = p.diff(- w0.cols()) * mother(w0.row(j - A.size()));
        w0.row(j - A.size()) /= sqrt(w0.row(j - A.size()).dot(w0.row(j - A.size())));
      }
    }
    w0 = w0.transpose();
    for(int j = 0; j < w0.rows(); j ++)
      for(int k = 0; k < w0.cols(); k ++)
        assert(isfinite(w0(j, k)) && ! isnan(w0(j, k)));
    const auto left(w0.LSVD().transpose() * w0);
    for(int k = 0; k < left.rows(); k ++)
      w00(i - A.size(), k) = sqrt(left.row(k).dot(left.row(k))) + T(1);
    w00.row(i - A.size())  = p.diff(- w00.cols()) * mother(w00.row(i - A.size()));
    w00.row(i - A.size()) /= sqrt(w00.row(i - A.size()).dot(w00.row(i - A.size())));
  }
  w00 = w00.transpose();
  const auto left(w00.LSVD().transpose() * w00);
  for(int k = 0; k < left.rows(); k ++)
    res0(0, k) = sqrt(left.row(k).dot(left.row(k))) + T(1);
  res0.row(0)  = p.diff(- res0.cols()) * mother(res0.row(0));
  res0.row(0) /= sqrt(res0.row(0).dot(res0.row(0)));
  // N.B. recursive on them.
  if(0 < depth && A.size() * 4 <= min(img.rows(), img.cols()) / 2) {
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
      dres[i] = decompose(dimg[i], depth - 1);
    Mat res(1 + dres[0].rows() * 5, A.size());
    res.row(0) = res0.row(0);
    for(int i = 0; i < 5; i ++)
      for(int j = 0; j < dres[i].rows(); j ++)
        res.row(1 + i * dres[i].rows() + j) = dres[i].row(j);
    return res;
  }
  return res0;
}

template <typename T> typename Decompose<T>::Mat Decompose<T>::subImage(const Mat& img, const int& x, const int& y, const int& r) const {
  Mat res(A.size(), A.size());
  for(int i = 0; i < res.rows(); i ++)
    for(int j = 0; j < res.cols(); j ++) {
      const auto rr(T(j - res.cols() / 2) / T(res.cols() / 2) * r);
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

