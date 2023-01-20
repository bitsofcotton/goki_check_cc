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
// \[- &alpha, &alpha;\[ register computer with deterministic calculation.
// cf. bitsofcotton/randtools .
template <typename T, typename feeder> class P1I {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleMatrix<T> Mat;
  inline P1I() { varlen = 0; }
  inline P1I(const int& stat, const int& var, const int& step = 1) {
    assert(0 < stat && 0 < var && 0 < step);
    f = feeder(stat + (this->step = step) - 1 + (varlen = var) - 1);
  }
  inline ~P1I() { ; }
  inline T next(const T& in) {
    static const T zero(0);
    const auto& buf(f.next(in));
    if(! f.full) return T(0);
    // N.B. please use catgp to compete with over learning.
    // XXX: division accuracy glitch.
    const auto nin(sqrt(buf.dot(buf)) * T(int(2)));
    if(! isfinite(nin) || nin == zero) return zero;
    SimpleMatrix<T> toeplitz(buf.size() - varlen - step + 2, varlen + 2);
    for(int i = 0; i < toeplitz.rows(); i ++) {
      auto work(buf.subVector(i, varlen) / nin);
      work[work.size() - 1] = buf[i + varlen + step - 2] / nin;
      auto mp(makeProgramInvariant<T>(move(work),
                T(i + 1) / T(toeplitz.rows() + 1)));
      toeplitz.row(i)  = move(mp.first);
      toeplitz.row(i) *= pow(mp.second, ceil(- log(toeplitz.epsilon()) ));
    }
    const auto invariant(linearInvariant<T>(toeplitz));
    if(invariant[varlen - 1] == zero) return zero;
    SimpleVector<T> work(varlen);
    for(int i = 1; i < work.size(); i ++)
      work[i - 1] = buf[i - work.size() + buf.size()] / nin;
    work[work.size() - 1] = zero;
    const auto work2(makeProgramInvariant<T>(work, T(1)));
    return revertProgramInvariant<T>(make_pair(
             - (invariant.dot(work2.first) -
                    invariant[varlen - 1] * work2.first[varlen - 1]) /
               invariant[varlen - 1], work2.second)) * nin;
  }
  feeder f;
private:
  int varlen;
  int step;
};

#define _P1_
#endif

