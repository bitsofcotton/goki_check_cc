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
template <typename T, typename feeder> class P1I {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleMatrix<T> Mat;
  inline P1I();
  inline P1I(const int& stat, const int& var);
  inline ~P1I();
  inline T next(const T& in);
  feeder f;
private:
  int varlen;
};

template <typename T, typename feeder> inline P1I<T,feeder>::P1I() {
  varlen = 0;
}

template <typename T, typename feeder> inline P1I<T,feeder>::P1I(const int& stat, const int& var) {
  assert(0 <= stat && 1 <= var);
  f = feeder(stat + (varlen = var) - 1);
}

template <typename T, typename feeder> inline P1I<T,feeder>::~P1I() {
  ;
}

template <typename T, typename feeder> T P1I<T,feeder>::next(const T& in) {
  const auto& buf(f.next(in));
  if(! f.full) return in;
  // N.B. to compete with noise, we calculate each.
  //      we can use catgp on worse noised ones.
  const auto nin(sqrt(buf.dot(buf)));
  if(nin == T(0)) return in;
  SimpleMatrix<T> toeplitz(buf.size() - varlen + 1, varlen + 3);
  for(int i = 0; i < toeplitz.rows(); i ++)
    toeplitz.row(i) =
      makeProgramInvariant<T>(
        buf.subVector(i, varlen) / nin, T(i + 1) / T(buf.size() - varlen + 2)
      ).first;
  const auto invariant(linearInvariant<T>(toeplitz));
  if(invariant[varlen - 1] == T(0)) return in;
  SimpleVector<T> work(varlen);
  for(int i = 1; i < work.size(); i ++)
    work[i - 1] = buf[i - work.size() + buf.size()] / nin;
  work[work.size() - 1] = work[work.size() - 2];
  auto work2(makeProgramInvariant<T>(work, T(1)));
  work = move(work2.first);
  return revertProgramInvariant<T>(make_pair((invariant.dot(work) - invariant[varlen - 1] * work[varlen - 1]) / invariant[varlen - 1], work2.second)) * nin;
}

#define _P1_
#endif

