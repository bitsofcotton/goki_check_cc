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
#if !defined(_P1_)

// Get invariant structure that
// \[- &alpha, &alpha;\[ register computer with deterministic calculation
// Please refer bitsofcotton/randtools .
// with noised ones, we can use catgp instead of this.
template <typename T, typename feeder, bool tanspace = true> class P1I {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleMatrix<T> Mat;
  inline P1I();
  inline P1I(const int& stat, const int& var, const int& comp = 0);
  inline ~P1I();
  inline T next(const T& in);
private:
  feeder f;
  Mat pp;
  int statlen;
};

template <typename T, typename feeder, bool tanspace> inline P1I<T,feeder,tanspace>::P1I() {
  statlen = 0;
}

template <typename T, typename feeder, bool tanspace> inline P1I<T,feeder,tanspace>::P1I(const int& stat, const int& var, const int& comp0) {
  assert(0 <= comp0);
  auto comp(comp0);
  T    work(var + 2);
  while(true) {
    work /= T(2);
    if(work < T(1)) break;
    work  = ceil(work);
    comp += int(work);
  }
  assert(var <= comp);
  assert(0 <= stat && 1 <= var);
  f = feeder((statlen = stat + comp + 2) + comp - 1);
  pp = Mat(comp, var);
  for(int i = 0; i < pp.rows() - 1; i ++) {
    const auto work(taylor(pp.cols() - 1, T(1) / T(pp.rows() - 2) * T(pp.cols() - 2)));
    for(int j = 0; j < work.size(); j ++)
      pp(i, j) = work[j];
    pp(i, work.size()) = T(0);
  }
  for(int i = 0; i < pp.cols(); i ++)
    pp(pp.rows() - 1, i) = T(i == pp.cols() - 1 ? 1 : 0);
  T M(0);
  for(int i = 0; i < pp.rows(); i ++)
    for(int j = 0; j < pp.cols(); j ++)
      M = max(M, abs(pp(i, j)));
  pp /= M * T(4);
}

template <typename T, typename feeder, bool tanspace> inline P1I<T,feeder,tanspace>::~P1I() {
  ;
}

template <typename T, typename feeder, bool tanspace> T P1I<T,feeder,tanspace>::next(const T& in) {
  const auto& buf(f.next(in));
  if(! f.full) return T(0);
  // N.B. to compete with noise, we calculate each.
  //      we can use catgp on worse noised ones.
  const auto nin(sqrt(buf.dot(buf)));
  vector<SimpleVector<T>> toeplitz;
  toeplitz.reserve(statlen);
  for(int i = 0; i < statlen; i ++) {
    SimpleVector<T> work(pp.cols());
    for(int j = 0; j < work.size(); j ++)
      work[j] = buf[i + j] / nin;
    toeplitz.emplace_back(tanspace
      ? makeProgramInvariant<T>(pp * work, T(i + 1) / T(statlen + 1))
      : pp * work);
  }
  const auto invariant(linearInvariant<T>(toeplitz));
  SimpleVector<T> work(pp.cols());
  for(int i = 1; i < work.size(); i ++)
    work[i - 1] = buf[i - work.size() + buf.size()] / nin;
  work[work.size() - 1] = work[work.size() - 2];
  work = pp * work;
  if(invariant[pp.rows() - 1] == T(0)) return T(0);
  if(tanspace)
    work = makeProgramInvariant<T>(work, T(1));
  const auto p0((invariant.dot(work) - invariant[pp.rows() - 1] * work[pp.rows() - 1]) / invariant[pp.rows() - 1]);
  return (tanspace ? revertProgramInvariant<T>(p0) : p0) * nin;
}

#define _P1_
#endif

