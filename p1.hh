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
template <typename T> class P1I {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleMatrix<T> Mat;
  inline P1I();
  inline P1I(const int& stat, const int& var);
  inline ~P1I();
  inline T next(const T& in, const int& ratio = 0);
private:
  Vec buf;
  int statlen;
  int varlen;
  int t;
};

template <typename T> inline P1I<T>::P1I() {
  varlen = t = 0;
}

template <typename T> inline P1I<T>::P1I(const int& stat, const int& var) {
  assert(0 <= stat && 1 <= var);
  buf.resize((statlen = stat + var + 2) + (varlen = var) - 1);
  for(int i = 0; i < buf.size(); i ++)
    buf[i] = T(0);
  t = 0;
}

template <typename T> inline P1I<T>::~P1I() {
  ;
}

template <typename T> T P1I<T>::next(const T& in, const int& ratio) {
  for(int i = 0; i < buf.size() - 1; i ++)
    buf[i]  = buf[i + 1];
  buf[buf.size() - 1] = in;
  if(t ++ < buf.size()) return T(0);
  // N.B. to compete with noise, we calculate each.
  //      we can use catgp on worse noised ones.
  const auto nin(sqrt(buf.dot(buf)));
  vector<SimpleVector<T>> toeplitz;
  toeplitz.reserve(statlen);
  for(int i = 0; i < statlen; i ++) {
    SimpleVector<T> work(varlen);
    for(int j = 0; j < varlen; j ++)
      work[j] = buf[i + 1] / nin;
    toeplitz.emplace_back(makeProgramInvariant<T>(work, ratio, T(i + 1) / T(statlen + 1)));
  }
  const auto invariant(linearInvariant<T>(toeplitz));
        auto work(invariant);
  for(int i = 1; i < varlen; i ++)
    work[i - 1] = toeplitz[toeplitz.size() - 1][i];
  work[varlen - 1] = work[varlen - 2];
  work[varlen + 1] = work[varlen] = T(1);
  return (atan((invariant.dot(work) - invariant[varlen - 1] * work[varlen - 1]) / invariant[varlen - 1]) * T(2) / atan2(T(1), T(1)) - T(1)) * nin;
}

#define _P1_
#endif

