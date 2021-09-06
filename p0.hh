/*
 BSD 3-Clause License

Copyright (c) 2019-2021, bitsofcotton (kazunobu watatsu)
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
#if !defined(_P0_)

template <typename T, typename feeder> class P0 {
public:
  typedef SimpleVector<T> Vec;
  inline P0() { ; }
  inline P0(const int& size, const int& step = 1) {
    assert(2 < size);
    f = feeder(size);
    p = taylor(size, T(size + step - 1));
    for(int i = 3; i < p.size(); i ++) {
      const auto pp(taylor(i, T(i + step - 1)));
      for(int j = 0; j < pp.size(); j ++)
        p[j - pp.size() + p.size()] += pp[j];
    }
    p /= T(p.size() + 1 - 3);
  }
  inline ~P0() { ; };
  inline T next(const T& in) {
    const auto& ff(f.next(in));
    return f.full ? p.dot(ff) - ff[ff.size() - 1] : T(0);
  }
  Vec p;
  feeder f;
};

#define _P0_
#endif

