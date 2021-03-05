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
private:
  std::vector<Mat> A;
         Vec  prepare(const Vec& in, const int& idx = 0) const;
         void apply(Vec& v, const Vec& dst, const Vec& src, const int& idx = 0) const;
};

template <typename T> inline Decompose<T>::Decompose() {
  ;
}

template <typename T> inline Decompose<T>::Decompose(const int& size) {
  static P0<T> p0;
  for(int i = 0; i < size; i ++) {
    SimpleMatrix<T> AA(size, size);
    for(int j = 0; j < AA.rows(); j ++) {
      // XXX: we avoid const. function.
      const auto jj(T(j) * T(i + 1) / T(size + 1));
      AA.row(j) = p0.taylor(AA.cols(), (jj - floor(jj)) * T(size - 1));
    }
    A.emplace_back(AA);
  }
}

template <typename T> inline Decompose<T>::~Decompose() {
  ;
}

template <typename T> typename Decompose<T>::Vec Decompose<T>::mimic(const Vec& dst, const Vec& src, const T& intensity) const {
  const int  size(A.size());
  const auto size2(dst.size() / size);
  const auto size3(src.size() / size);
        auto res(dst);
  Vec ddst;
  Vec dsrc;
  for(int i = 0; i < size2; i ++) {
    const auto dd(prepare(dst, i));
    if(i) {
      dsrc += dd;
      ddst += synth(mother(prepare(src, i * size3 / size2)),
                    freq(mother(dd), dd)) * intensity +
              dd * (T(1) - abs(intensity));
   } else {
      dsrc  = dd;
      ddst  = synth(mother(prepare(src, i * size3 / size2)),
                    freq(mother(dd), dd)) * intensity +
              dd * (T(1) - abs(intensity));
   }
  }
  dsrc /= T(size2);
  ddst /= T(size2);
  for(int i = 0; i < size2; i ++)
    apply(res, ddst, dsrc, i);
  return res;
}

template <typename T> typename Decompose<T>::Vec Decompose<T>::emphasis(const Vec& dst, const T& intensity) const {
  const int  size(A.size());
  const auto size2(dst.size() / size);
        auto res(dst);
  Vec   ddst;
  Vec   dsrc;
  for(int i = 0; i < size2; i ++) {
    const auto dd(prepare(dst, i));
          auto lfreq(dd);
    lfreq[lfreq.size() - 1] = T(0);
    for(int j = 0; j < lfreq.size() - 1; j ++)
      lfreq[j] = T(j + 1) / T(lfreq.size() - 1);
    if(i) {
      dsrc += dd;
      ddst += synth(mother(dd), lfreq) * intensity;
    } else {
      dsrc  = dd;
      ddst  = synth(mother(dd), lfreq) * intensity;
    }
  }
  dsrc /= T(size2);
  ddst /= T(size2);
  for(int i = 0; i < size2; i ++)
    apply(res, ddst, dsrc, i);
  return res;
}

template <typename T> typename Decompose<T>::Vec Decompose<T>::enlarge(const Vec& in, const int& r) const {
  assert(1 < r);
  static P0<T> p0;
  static std::vector<std::vector<Mat> > p;
  static std::vector<std::vector<Vec> > f;
  static std::vector<Decompose<T> > e;
  if(p.size() < in.size()) {
    p.resize(in.size() + 1, std::vector<Mat>());
    f.resize(in.size() + 1, std::vector<Vec>());
  }
  if(p[in.size()].size() < r) {
    p[in.size()].resize(r + 1, Mat());
    f[in.size()].resize(r + 1, Vec());
  }
  if(e.size() < in.size() * r)
    e.resize(in.size() * r + 1, Decompose<T>());
  auto& pp(p[in.size()][r]);
  auto& ff(f[in.size()][r]);
  auto& ee(e[in.size() * r]);
  if(pp.rows() < in.size() * r) {
    pp.resize(in.size() * r, in.size());
    ff.resize(pp.rows());
    for(int i = 0; i < pp.rows(); i ++) {
      pp.row(i) = p0.taylor(in.size(), T(i) / T(r));
      ff[i] = T(1);
    }
    ff /= sqrt(ff.dot(ff));
    ee  = Decompose<T>(in.size() * r);
  }
  return ee.synth(pp * mother(in), ff) * sqrt(in.dot(in));
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
  Mat work(in.size(), in.size());
  Vec one(in.size());
  for(int i = 0; i < A.size(); i ++) {
    work.setCol(i, A[i] * in);
    one[i] = T(1);
  }
  auto f(work.solve(one));
  return f /= sqrt(f.dot(f));
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

#define _DECOMPOSE_
#endif

