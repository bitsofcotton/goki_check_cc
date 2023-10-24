/* BSD 3-Clause License:
 * Copyright (c) 2018-2021, bitsofcotton.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 *    Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *    Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation or other materials provided with the distribution.
 *    Neither the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#if !defined(_GOKICHECK_)

using std::max;
using std::min;
using std::cerr;
using std::endl;
using std::string;
using std::to_string;
using std::ifstream;
using std::getline;
using std::istringstream;
using std::stringstream;
using std::ofstream;
using std::vector;
using std::pair;
using std::make_pair;
using std::move;
using std::flush;
using std::sort;
using std::abs;
using std::lower_bound;
using std::upper_bound;
using std::distance;
using std::unique;
using std::istream;
using std::ostream;
using std::swap;
using std::binary_search;

template <typename T> static inline bool less0(const T& x, const T& y) {
  return x.first[0] < y.first[0] || (x.first[0] == y.first[0] && x.first[1] < y.first[1]);
}

template <typename T> static inline bool lessf(const T& x, const T& y) {
  return x.first < y.first;
}

template <typename T> using triangles_t = pair<SimpleMatrix<T>, T>;

typedef enum {
  SHARPEN_X,
  SHARPEN_Y,
  SHARPEN_BOTH,
  ENLARGE_X,
  ENLARGE_Y,
  ENLARGE_BOTH,
  FLARGE_X,
  FLARGE_Y,
  FLARGE_BOTH,
  BLUR_X,
  BLUR_Y,
  BLUR_BOTH,
  COLLECT_BOTH,
  BUMP_BOTH,
  BLINK_X,
  BLINK_Y,
  BLINK_BOTH,
  REPRESENT,
  CLIP } direction_t;

template <typename T> static inline T getImgPt(const T& y, const T& h) {
  auto yy(y % (2 * h));
  if(yy < 0)
    yy = - yy;
  if(yy >= h)
    yy = h - (yy - h);
  return yy % h;
}

template <typename T> bool saveobj(const vector<SimpleVector<T> >& data, const T& Mw0, const T& Mh0, const vector<SimpleVector<int> >& polys, const char* filename) {
  ofstream output;
  output.open(filename, std::ios::out);
  if(output.is_open()) {
    int lfs(0);
    const T Mh(Mh0 / T(2));
    const T Mw(Mw0 / T(2));
    for(int fslash(0) ; filename[fslash]; fslash ++)
      if(filename[fslash] == '/')
        lfs = fslash;
    if(lfs) lfs ++;
    output << "mtllib " << &filename[lfs] << ".mtl" << endl;
    output << "usemtl material0" << endl;
    for(int i = 0; i < data.size(); i ++)
      output << "v " << data[i][1] << " " << - data[i][0] << " " << data[i][2] << endl;
    for(int i = 0; i < data.size(); i ++)
      output << "vt " << data[i][1] / T(Mh) / T(2) << " " << T(1) - data[i][0] / T(Mw) / T(2) << endl;
    // xchg with clockwise/counter clockwise.
    for(int i = 0; i < polys.size(); i ++) {
      const int i0(polys[i][0] + 1);
      const int i1(polys[i][1] + 1);
      const int i2(polys[i][2] + 1);
      if(i0 != i1 && i1 != i2 && i2 != i0) {
        output << "f " << i0 << "/" << i0 << "/" << i0;
        output << " "  << i1 << "/" << i1 << "/" << i1;
        output << " "  << i2 << "/" << i2 << "/" << i2 << endl;
      }
    }
    output.close();
  } else {
    cerr << "Unable to open file for write: " << filename << endl;
    return false;
  }
  return true;
}

template <typename T> bool loadobj(vector<SimpleVector<T> >& data, vector<SimpleVector<int> >& polys, const char* filename) {
  ifstream input;
  input.open(filename);
  if(input.is_open()) {
    string work;
    while(getline(input, work) && !input.eof() && !input.bad()) {
      int i = 0;
      for( ; i < work.size() && work[i] == ' '; i ++);
      if(i + 1 < work.size() && work[i] == 'v' && work[i + 1] == ' ') {
        stringstream sub(work.substr(i + 2, work.size() - (i + 2)));
        SimpleVector<T> buf(3);
        sub >> buf[1];
        sub >> buf[0];
        sub >> buf[2];
        buf[0] = - buf[0];
        data.emplace_back(move(buf));
      } else if(i + 1 < work.size() && work[i] == 'f' && work[i + 1] == ' ') {
        stringstream sub(work.substr(i + 2, work.size() - (i + 2)));
        SimpleVector<int> wbuf(3);
        int  widx(0);
        bool flag(false);
        while(!sub.eof() && !sub.bad()) {
          sub >> wbuf[widx];
          if(wbuf[widx] >= 0)
            wbuf[widx] --;
          widx ++;
          if(widx > 2)
            flag = true;
          if(flag)
            polys.emplace_back(wbuf);
          widx %= 3;
          if(sub.eof() || sub.bad())
            break;
          sub.ignore(20, ' ');
        }
      }
    }
    for(int i = 0; i < polys.size(); i ++)
      for(int j = 0; j < polys[i].size(); j ++) {
        while(polys[i][j] < 0) polys[i][j] += data.size();
        polys[i][j] = polys[i][j] % data.size();
      }
    input.close();
  } else {
    cerr << "Unable to open file for read: " << filename << endl;
    return false;
  }
  return true;
}
  
template <typename T> bool saveMTL(const char* photo, const char* filename) {
  ofstream output;
  output.open(filename, std::ios::out);
  if(output.is_open()) {
    string pstr(photo);
    for(int i = 0; i < pstr.size(); i ++)
      if(pstr[pstr.size() - i - 1] == '.') {
        pstr = pstr.substr(0, max(int(0), int(pstr.size()) - i - 1)) + string(".ppm");
        break;
      }
    output << "newmtl material0" << endl;
    output << "Ka 1.000000 1.000000 1.000000" << endl;
    output << "Kd 1.000000 1.000000 1.000000" << endl;
    output << "Ks 0.000000 0.000000 0.000000" << endl;
    output << "illum 1" << endl;
    output << "map_Ka " << pstr << endl;
    output << "map_Kd " << pstr << endl << endl;
    output.close();
  } else {
    cerr << "Unable to open file for write: " << filename << endl;
    return false;
  }
  return true;
}

template <typename T> bool loaddat(const char* filename, string& header, vector<vector<T> >& data) {
  ifstream input;
  input.open(filename);
  if(input.is_open()) {
    string work;
    header = string("");
    data   = vector<vector<T> >();
    while(getline(input, work) && !input.eof() && !input.bad())
      if(whiteline(work))
        continue;
      else if(work[0] == ';')
        header += work + string("\n");
      else {
        stringstream ss(work);
        for(int i = 0, j = 0; ss.tellg() <= work.size(); j ++) {
          if(data.size() <= j)
            data.resize(j + 1, vector<T>());
          data[j].emplace_back(T(0));
          ss >> data[j][data[j].size() - 1];
        }
      }
    input.close();
  } else {
    cerr << "Unable to open file for read: " << filename << endl;
    return false;
  }
  return true;
}
  
template <typename T> bool savedat(const char* filename, string& header, vector<vector<T> >& data) {
  ofstream output;
  output.open(filename, std::ios::out);
  if(output.is_open()) {
    output << header;
    for(int i = 0; i < data[0].size(); i ++) {
      for(int j = 0; j < data.size(); j ++)
        output << (i < data[j].size() ? data[j][i] : T(0)) << " ";
      output << endl;
    }
    output.close();
  } else {
    cerr << "Unable to open file for write: " << filename << endl;
    return false;
  }
  return true;
}
  
template <typename T> bool loadcenterr(vector<SimpleVector<T> >& center, vector<T>& r, const char* filename) {
  center = vector<SimpleVector<T> >();
  r      = vector<T>();
  ifstream input;
  try {
    input.open(filename);
    string buf;
    while(getline(input, buf) && !input.eof() && !input.bad()) {
      stringstream sbuf(buf);
      SimpleVector<T> work(3);
      sbuf >> work[0];
      sbuf >> work[1];
      sbuf >> work[2];
      center.emplace_back(move(work));
      T workr;
      sbuf >> workr;
      r.emplace_back(workr);
    }
    input.close();
  } catch(...) {
    cerr << "Something had occured when reading center - r txt." << endl;
    return false;
  }
  return center.size() == r.size();
}

template <typename T> bool savecenterr(const char* filename, const vector<SimpleVector<T> >& center, const vector<T>& r) {
  ofstream output;
  output.open(filename, std::ios::out);
  if(output.is_open()) {
    assert(center.size() == r.size());
    for(int i = 0; i < center.size(); i ++) {
      assert(center[i].size() == 3);
      output << center[i][0] << " " << center[i][1] << " " << center[i][2] << " " << r[i] << endl;
    }
    output.close();
  } else {
    cerr << "Unable to open file for write: " << filename << endl;
    return false;
  }
  return true;
}

template <typename T> SimpleMatrix<T> sharpen(const int& size) {
  assert(0 < size);
  SimpleMatrix<T> s;
  const auto file(string("./.cache/lieonn/sharpen-") + to_string(size) +
#if defined(_FLOAT_BITS_)
    string("-") + to_string(_FLOAT_BITS_)
#else
    string("-ld")
#endif
  );
  ifstream cache(file.c_str());
  if(cache.is_open()) {
    cache >> s;
    cache.close();
  } else {
/*
    if(2 < size) {
      const auto s0(sharpen<T>(size - 1) * T(size - 1));
      s = SimpleMatrix<T>(size, size).O().setMatrix(0, 0, s0);
      s.setMatrix(1, 1, s.subMatrix(1, 1, size - 1, size - 1) + s0);
      s.row(0) *= T(2);
      s.row(s.rows() - 1) *= T(2);
      s /= T(2);
    } else
      s  = SimpleMatrix<T>(size, size).O();
*/
    auto dfts(dft<T>(size));
    static const auto Pi(atan(T(1)) * T(4));
    dfts.row(0) *= complex<T>(T(0));
    for(int i = 1; i < dfts.rows(); i ++) {
      // N.B. d/dt((d^(t)/dy^(t)) f), differential-integral space tilt on f.
      // DFTH.row(i) *= log(phase);
      // N.B. please refer enlarge.wxm, half freq space refer and uses each.
      //   -> This is sharpen operation at all because this is same as original
      //      picture when {x0 + x0.5, x0.5 + x1, x1 + x1.5, x1.5 + x2, ...}
      //      series, and both picture of dft is same, them, pick {x0, x1, ...}.
      dfts.row(i) /= exp(complex<T>(T(0), Pi * T(i) / T(dfts.rows()))) - complex<T>(T(1));
    }
    s = (dft<T>(- size) * dfts).template real<T>() / T(size - 1);
/*
    if(2 < size)
      s /= T(size);
*/
    ofstream ocache(file.c_str());
    ocache << s;
    ocache.close();
    cerr << "." << flush;
  }
  return s;
}



template <typename T> SimpleMatrix<T> rotate(const SimpleMatrix<T>& d, const T& theta) {
  assert(abs(theta) < atan(T(1)));
  const auto c(cos(theta));
  const auto s(sin(theta));
  const auto h0(abs(int(c * T(d.rows()) - s * T(d.cols()))));
  const auto h1(h0 + abs(int(s * T(d.cols()))) * 2);
  const auto w0(abs(int(s * T(d.rows()) + c * T(d.cols()))));
  const auto w1(w0 + abs(int(s * T(d.rows()))) * 2);
  SimpleMatrix<T> res(h0 < d.rows() ? h1 : h0,
                      w0 < d.cols() ? w1 : w0);
  const T offy(h0 < d.rows() ? abs(int(s * T(d.cols()))) : 0);
  const T offx(w0 < d.cols() ? abs(int(s * T(d.rows()))) : 0);
  res.O();
  const auto diag(int(sqrt(res.rows() * res.rows() +
                           res.cols() * res.cols())) + 1);
  for(int j = - diag; j < diag; j ++)
    for(int k = - diag; k < diag; k ++) {
      const int yy(c * T(j) - s * T(k) + offy);
      const int xx(s * T(j) + c * T(k) + offx);
      if(0 <= yy && yy < res.rows() &&
         0 <= xx && xx < res.cols()) {
        const auto dyy(getImgPt<int>(j, d.rows()));
        const auto dxx(getImgPt<int>(k, d.cols()));
        {
          res(yy, xx) = res(min(yy + 1, int(res.rows()) - 1), xx) =
            res(yy, min(xx + 1, int(res.cols()) - 1)) =
            res(min(yy + 1, int(res.rows()) - 1),
                min(xx + 1, int(res.cols()) - 1)) =
              d(dyy, dxx);
        }
      }
    }
  return res;
}

template <typename T> static inline SimpleMatrix<T> center(const SimpleMatrix<T>& dr, const SimpleMatrix<T>& d) {
  SimpleMatrix<T> res(d.rows(), d.cols());
  for(int i = 0; i < res.rows(); i ++)
    for(int j = 0; j < res.cols(); j ++)
      res(i, j) = dr(max(0, min(i + (dr.rows() - d.rows()) / 2, dr.rows() - 1)),
                     max(0, min(j + (dr.cols() - d.cols()) / 2, dr.cols() - 1)));
  return res;
}

template <typename T> static inline SimpleMatrix<T> flip(const SimpleMatrix<T>& d) {
  auto res(d);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < d.rows(); i ++)
    res.row(res.rows() - 1 - i) = d.row(i);
  return res;
}

template <typename T> static inline SimpleMatrix<T> flop(const SimpleMatrix<T>& d) {
  auto res(d);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < d.cols(); i ++)
    res.setCol(res.cols() - 1 - i, d.col(i));
  return res;
}

template <typename T> static inline SimpleMatrix<T> normalize(const SimpleMatrix<T>& data, const T& upper = T(1)) {
  vector<SimpleMatrix<T> > work;
  work.emplace_back(data);
  return normalize<T>(work, upper)[0];
}

template <typename T> static inline SimpleMatrix<T> autoLevel(const SimpleMatrix<T>& data, const int& count = 0) {
  vector<SimpleMatrix<T> > work;
  work.emplace_back(data);
  return autoLevel(work, count)[0];
}

// N.B. this function is NOT thread safe.
template <typename T> SimpleMatrix<T> filter(const SimpleMatrix<T>& data, const direction_t& dir, const int& recur = 2, const int& rot = 0) {
  assert(0 <= rot);
  if(0 < rot) {
    auto res(filter<T>(data, dir, recur));
    if(rot <= 1) return res;
    for(int i = 0; i < rot; i ++) {
      cerr << "r" << flush;
      const auto theta((T(i) - T(rot - 1) / T(2)) * atan(T(1)) / (T(rot) / T(2)));
      res += center<T>(rotate<T>(filter<T>(rotate<T>(data, theta),
                       dir, recur), - theta), res);
    }
    return res /= T(rot + 1);
  }
  SimpleMatrix<T> result;
  static const auto Pi(atan2(T(1), T(1)) * T(4));
  switch(dir) {
  case SHARPEN_BOTH:
    result = filter<T>(filter<T>(data, SHARPEN_X, recur), SHARPEN_Y, recur);
    break;
  case ENLARGE_BOTH:
    result = filter<T>(filter<T>(data, ENLARGE_X, recur), ENLARGE_Y, recur);
    break;
  case FLARGE_BOTH:
    result = filter<T>(filter<T>(data, FLARGE_X, recur), FLARGE_Y, recur);
    break;
  case BLUR_BOTH:
    result = filter<T>(filter<T>(data, BLUR_Y, recur), BLUR_X, recur);
    break;
  case BLINK_BOTH:
    result = filter<T>(filter<T>(data, BLINK_X, recur), BLINK_Y, recur);
    break;
  case SHARPEN_X:
    result = filter<T>(data.transpose(), SHARPEN_Y, recur).transpose();
    break;
  case ENLARGE_X:
    result = filter<T>(data.transpose(), ENLARGE_Y, recur).transpose();
    break;
  case FLARGE_X:
    result = filter<T>(data.transpose(), FLARGE_Y, recur).transpose();
    break;
  case BLUR_X:
    result = filter<T>(data.transpose(), BLUR_Y, recur).transpose();
    break;
  case BLINK_X:
    result = filter<T>(data.transpose(), BLINK_Y, recur).transpose();
    break;
  case BLUR_Y:
    {
      // N.B. d^exp(t)/dx^exp(t) f(x) == f(x + t dx), t != 0.
      //   so d^exp(- t)/dx^exp(- t)
      //   == d^(- exp(t))/dx^(- exp(t)) == f(x - t dx).
      //   we normalize with f(x - dx) + f(x + dx) in weak differential meaning.
      // Cor: d/dt d^t/d(x^t) f(x)
      //   == d/dt f(x + log(t) dx) == f'(x + log(t) dx) / g(x, t).
      // N.B. d^0/dx^0 f(x) == f(x)
      //   == d^exp(- inf)/dx^exp(- inf) f(x) == f(x - inf dx)
      //   == d^(- exp(- inf))/dx^(- exp(- inf dx))
      //   == d^(exp(inf))/dx^exp(inf) == f(x + inf dx)
      auto normalize(dft<T>(data.rows()) * data.template cast<complex<T> >());
      for(int i = 0; i < normalize.rows(); i ++) {
        const auto n(complex<T>(T(0), - T(2) * Pi * T(i) / T(normalize.rows())));
        normalize.row(i) *= exp(n) + exp(- n);
      }
      result = (dft<T>(- data.rows()) * normalize).template real<T>();
    }
    break;
  case SHARPEN_Y:
    result = filter(data - sharpen<T>(int(data.rows())) * data, CLIP);
    break;
  case ENLARGE_Y:
    {
      assert(2 <= data.rows());
      static vector<vector<SimpleMatrix<T> > > Eop;
      const auto& size(data.rows());
      if(Eop.size() <= size)
        Eop.resize(size + 1, vector<SimpleMatrix<T> >());
      else if(recur < Eop[size].size())
        goto eopi;
      if(Eop[size].size() <= recur)
        Eop[size].resize(recur + 1, SimpleMatrix<T>());
      {
        auto& eop(Eop[size][recur]);
        if(eop.cols() == size)
          goto eopi;
        cerr << "e" << flush;
        eop.resize((size - 1) * recur + 1, size);
        for(int j = 0; j < eop.rows(); j ++)
          // N.B. sampling th. hack isn't work well.
          eop.row(j) = taylor<T>(eop.cols(), T(j) / T(eop.rows() - 1) * T(eop.cols() - 1));
      }
     eopi:
      result = Eop[size][recur] * data;
    }
    break;
  case FLARGE_Y:
    {
      SimpleMatrix<T> work(data);
      for(int i = 0; i < work.rows(); i ++)
        for(int j = 0; j < work.cols(); j ++)
          work(i, j) += T(1) / T(256);
      Decompose<T> e(work.rows());
      result = SimpleMatrix<T>(work.rows() * recur, work.cols());
      for(int i = 0; i < work.cols(); i ++)
        result.setCol(i, e.enlarge(work.col(i), recur));
    }
    break;
  case COLLECT_BOTH:
    {
      const auto zy(diff<T>(data.rows()) * data);
      const auto zx(data * diff<T>(data.cols()).transpose());
      const auto zxx(diff<T>(data.rows()) * diff<T>(data.rows()) * data);
      const auto zxy(diff<T>(data.rows()) * data * diff<T>(data.cols()).transpose());
      const auto zyy(data * diff<T>(data.cols()).transpose() * diff<T>(data.cols()).transpose());
      result = SimpleMatrix<T>(data.rows(), data.cols()).O();
      for(int i = 0; i < result.rows(); i ++)
        for(int j = 0; j < result.cols(); j ++)
          // N.B. thanks to https://en.wikipedia.org/wiki/Gaussian_curvature .
          result(i, j) = abs(
            (zxx(i, j) * zyy(i, j) - zxy(i, j) * zxy(i, j)) /
            (zx(i, j) * zx(i, j) + zy(i, j) * zy(i, j) + T(int(1))) /
            (zx(i, j) * zx(i, j) + zy(i, j) * zy(i, j) + T(int(1))) );
    }
    break;
  case BUMP_BOTH:
    {
      SimpleMatrix<T> zscore(data.rows(), data.cols());
      result.resize(data.rows(), data.cols());
      result.O();
      zscore.O(- T(1));
      const auto rxy(T(min(data.rows(), data.cols())));
      const int  dratio(sqrt(sqrt(rxy)));
            SimpleVector<T> camera(2);
            SimpleVector<T> cpoint(2);
      camera[0] = T(0);
      camera[1] = T(1);
      cpoint[0] = T(1) / T(2 * dratio);
      for(int zi = 0; zi < dratio; zi ++) {
        cerr << "z" << flush;
        // N.B. projection scale is linear.
        cpoint[1] = T(zi) / T(dratio);
        // x-z plane projection of point p with camera geometry c to z=0.
        // c := camera, p := cpoint.
        // <c + (p - c) * t, [0, 1]> = 0
        const auto t(- camera[1] / (cpoint[1] - camera[1]));
        const auto y0((camera + (cpoint - camera) * t)[0] * rxy);
        if(abs(int(y0)) < 3 || rxy < abs(y0) * T(2)) continue;
        const auto Dop(diff<T>(abs(int(y0) & ~ int(1))));
        const auto Dop0((Dop.row(int(y0) / 2) + Dop.row(int(y0) / 2 + 1)) / T(2));
        const auto DDop(Dop * Dop);
        const auto DDop0((DDop.row(int(y0) / 2) + DDop.row(int(y0) / 2 + 1)) / T(2));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
        for(int i = 0; i < data.rows(); i ++)
          for(int j = 0; j < data.cols(); j ++) {
            T fu(0), fv(0), L(0), M(0), N(0);
            for(int kk = 0; kk < Dop0.size(); kk ++) {
              fu += data(getImgPt<int>(i + kk - Dop0.size() / 2,
                data.rows()), j) * Dop0[kk];
              fv += data(i, getImgPt<int>(j + kk - Dop0.size() / 2,
                data.cols()) ) * Dop0[kk];
              L += data(getImgPt<int>(i + kk - Dop0.size() / 2,
                data.rows()), j) * DDop0[kk];
              N += data(i, getImgPt<int>(j + kk - Dop0.size() / 2,
                data.cols()) ) * DDop0[kk];
              for(int ll = 0; ll < Dop0.size(); ll ++)
                M += Dop0[kk] * Dop0[ll] *
                  data(getImgPt<int>(i + kk - Dop0.size() / 2, data.rows()),
                       getImgPt<int>(j + ll - Dop0.size() / 2, data.cols()));
            }
            // N.B. thanks to https://en.wikipedia.org/wiki/Gaussian_curvature .
            const auto lscore(abs((L * N - M * M) / (fu * fu + fv * fv + T(int(1))) / (fu * fu + fv * fv + T(int(1))) ));
            if(zscore(i, j) <= lscore) {
              result(i, j) = T(zi + 1) / T(dratio);
              zscore(i, j) = lscore;
            }
          }
      }
      // N.B. we don't need local to global with correct gaussian curvature.
      assert(result.rows() == data.rows() && result.cols() == data.cols());
    }
    break;
  case BLINK_Y:
    {
      auto dif(dft<T>(data.rows()) * data.template cast<complex<T> >());
      for(int i = 1; i < data.rows(); i ++)
        dif.row(i) *= - complex<T>(T(0), T(2)) * T(i) / T(data.rows());
      dif = dft<T>(- data.rows()) * dif;
      for(int i = 1; i < dif.rows() - 1; i ++)
        dif.row(i) += (dif.row(i - 1) + dif.row(i + 1)) * complex<T>(T(recur) / T(256));
      dif = dft<T>(data.rows()) * dif;
      for(int i = 1; i < data.rows(); i ++)
        dif.row(i) /= - complex<T>(T(0), T(2)) * T(i) / T(data.rows());
      result = (dft<T>(- data.rows()) * dif).template real<T>();
    }
    break;
  case REPRESENT:
    result = Decompose<T>(recur).represent(data, 2);
    break;
  case CLIP:
    result.resize(data.rows(), data.cols());
    for(int i = 0; i < result.rows(); i ++)
      for(int j = 0; j < result.cols(); j ++)
        result(i, j) = max(T(0), min(T(1), data(i, j)));
    break;
  default:
    assert(0 && "unknown command in filter (should not be reached.)");
  }
  return result;
}

template <typename T> class match_t {
public:
  typedef SimpleMatrix<T>   Mat;
  typedef SimpleVector<T>   Vec;
  typedef SimpleVector<int> Veci;
  Mat         rot;
  Vec         offset;
  T           ratio;
  T           rdepth;
  vector<int> dst;
  vector<int> src;
  T           thresh;
  T           othresh;
  inline match_t() {
    thresh  = T(0);
    othresh = T(0);
    initId();
  }
  inline match_t(const T& thresh, const T& d) {
    this->thresh  = thresh;
    this->othresh = thresh * d;
    initId();
  }
  inline match_t(const match_t<T>& s) {
    *this = s;
  }
  inline void initId() {
    rot       = Mat(3, 3);
    rot(0, 0) = rot(1, 1) = rot(2, 2) = T(1);
    rot(1, 0) = rot(2, 0) = rot(0, 1) = rot(2, 1)
              = rot(0, 2) = rot(1, 2) = T(0);
    offset    = Vec(3);
    offset[0] = offset[1] = offset[2] = T(0);
    ratio     = T(1);
    rdepth    = T(0);
  }
  inline match_t<T>  operator ~ () const {
    match_t<T> result(*this);
    result.rot    = rot.transpose();
    result.ratio  = T(1) / ratio;
    result.offset = - result.rot * offset * result.ratio;
    swap(result.src, result.dst);
    return result;
  }
  inline match_t<T>  operator / (const match_t<T>& s) const {
    match_t<T> result;
    result.rot    = rot   * s.rot.transpose();
    result.ratio  = ratio / s.ratio;
    result.offset = offset - result.rot * s.offset * result.ratio;
    result.rdepth = rdepth + s.rdepth;
    result.dst = vector<int>();
    result.src = vector<int>();
    for(int i = 0; i < src.size(); i ++)
      for(int j = 0; j < s.src.size(); j ++)
        if(src[i] == s.src[j]) {
          result.dst.emplace_back(  dst[i]);
          result.src.emplace_back(s.dst[j]);
        }
    return result;
  }
  inline match_t<T>& operator = (const match_t<T>& other) {
    rot        = other.rot;
    offset     = other.offset;
    ratio      = other.ratio;
    rdepth     = other.rdepth;
    dst  = other.dst;
    src  = other.src;
    thresh     = other.thresh;
    othresh    = other.othresh;
    return *this;
  }
  inline T distance(const match_t<T>& other, const Vec& p) {
    const auto d(transform(p) - other.transform(p));
    return sqrt(d.dot(d));
  }
  inline vector<Veci> hullConv(const vector<Veci>& srchull) const {
    assert(src.size() == dst.size());
    vector<Veci> res;
    res.reserve(srchull.size());
    for(int i = 0; i < srchull.size(); i ++) {
      Veci tmp(3);
      tmp[0] = tmp[1] = tmp[2] = - 1;
      for(int j = 0; j < srchull[i].size(); j ++)
        for(int k = 0; k < src.size(); k ++)
          if(src[k] == srchull[i][j]) {
            tmp[j] = dst[k];
            break;
          }
      assert(0 <= tmp[0] && 0 <= tmp[1] && 0 <= tmp[2]);
      res.emplace_back(tmp);
    }
    return res;
  }
  inline Vec transform(const Vec& x) const {
    return rot * x * ratio + offset;
  }
  inline vector<Vec> transform(const vector<Vec>& x) const {
    vector<Vec> result(x);
    for(int i = 0; i < result.size(); i ++)
      result[i] = transform(result[i]);
    return result;
  }
  inline bool operator < (const match_t<T>& x1) const {
    const T rratio(max(abs(   ratio), T(1) / abs(   ratio)));
    const T xratio(max(abs(x1.ratio), T(1) / abs(x1.ratio)));
    return dst.size() > x1.dst.size() || (dst.size() == x1.dst.size() && (rdepth < x1.rdepth || (rdepth == x1.rdepth && rratio < xratio)));
  }
  inline bool operator != (const match_t<T>& x) const {
    const auto test(offset - x.offset);
    const auto roterr(rot * x.rot.transpose());
    return !(abs(T(1) - roterr(0, 0)) <= thresh) ||
           !(abs(T(1) - roterr(1, 1)) <= thresh) ||
           !(abs(T(1) - roterr(2, 2)) <= thresh) ||
           !(sqrt(test.dot(test) / 
               (offset.dot(offset) + x.offset.dot(x.offset))) <= othresh) ||
           ratio * x.ratio < T(0) ||
           !(abs(ratio - x.ratio) / sqrt(ratio * x.ratio) <= thresh);
  }
  inline bool operator == (const match_t<T>& x) const {
    return ! (*this != x);
  }
  friend ostream& operator << (ostream& os, const match_t<T>& x) {
    os << x.rot;
    os << x.offset;
    os << x.ratio  << endl;
    os << x.rdepth << endl;
    assert(x.dst.size() == x.src.size());
    os << x.dst.size() << endl;
    for(int i = 0; i < x.dst.size(); i ++)
      os << x.dst[i] << " ";
    os << endl;
    for(int i = 0; i < x.src.size(); i ++)
      os << x.src[i] << " ";
    os << endl;
    os << x.thresh  << endl;
    return os;
  }
  friend istream& operator >> (istream& is, match_t<T>& x) {
    try {
      is >> x.rot;
      is >> x.offset;
      is >> x.ratio;
      is >> x.rdepth;
      int size(0);
      is >> size;
      assert(size > 0);
      x.dst.resize(size);
      x.src.resize(size);
      for(int i = 0; i < size; i ++)
        is >> x.dst[i];
      for(int i = 0; i < size; i ++)
        is >> x.src[i];
      is >> x.thresh;
    } catch(...) {
      assert(0 && "match_t input failed.");
    }
    return is;
  }
};

template <typename T> static inline pair<SimpleVector<T>, vector<SimpleVector<T> > > makeG(const vector<SimpleVector<T> >& in) {
  pair<SimpleVector<T>, vector<SimpleVector<T> > > res;
  res.second.reserve(in.size());
  assert(in.size() && in[0].size() == 3);
  res.first = in[0];
  for(int i = 1; i < in.size(); i ++)
    res.first += in[i];
  res.first /= T(in.size());
  for(int i = 0; i < in.size(); i ++)
    res.second.emplace_back(in[i] - res.first);
  assert(res.second.size() == in.size());
  return move(res);
}

template <typename T> static inline pair<T, vector<SimpleVector<T> > > normalizeG(const vector<SimpleVector<T> >& s) {
  pair<T, vector<SimpleVector<T> > > res;
  res.first = T(0);
  res.second.reserve(s.size());
  for(int i = 0; i < s.size(); i ++)
    res.first = max(max(res.first, abs(s[i][0])), max(abs(s[i][1]), abs(s[i][2])));
  for(int i = 0; i < s.size(); i ++)
    res.second.emplace_back(s[i] / res.first);
  assert(res.second.size() == s.size());
  return move(res);
}

template <typename T> static inline SimpleVector<T> toQuarterNormalize(const SimpleVector<T>& xyz) {
  assert(xyz.size() == 3);
  static const T    zero(int(0));
  static const auto twoPi(atan(T(int(1))) * T(int(4)));
  SimpleVector<T> quat(4);
  quat[0] = sqrt(xyz.dot(xyz) / T(int(6)));
  // y-z plane
  quat[1] = xyz[1] == zero && xyz[2] == zero ? zero : atan2(xyz[1], xyz[2]) / twoPi;
  // same for z-x.
  quat[2] = xyz[2] == zero && xyz[0] == zero ? zero : atan2(xyz[2], xyz[0]) / twoPi;
  // same for x-y.
  quat[3] = xyz[0] == zero && xyz[1] == zero ? zero : atan2(xyz[0], xyz[1]) / twoPi;
  return quat;
}

template <typename T> match_t<T> reconfigureMatch(match_t<T>& m, const vector<SimpleVector<T> >& dst0, const vector<SimpleVector<T> >& src0) {
  SimpleVector<T> off(3);
  off.O();
  for(int k = 0; k < m.dst.size(); k ++)
    off += dst0[m.dst[k]] - m.transform(src0[m.src[k]]);
  m.offset += (off /= T(int(m.dst.size())));
  if(m.dst.size() < 4) return m;
  SimpleMatrix<T> rot(3, 3);
  rot.I();
  for(int k = 0; k < m.dst.size() - 3; k ++) {
    SimpleMatrix<T> rotl(3, 3);
    SimpleMatrix<T> rotr(3, 3);
    rotl.O();
    rotr.O();
    for(int kk = 0; kk < 3; kk ++) {
      rotl.setCol(kk, dst0[m.dst[k + kk]]);
      rotr.setCol(kk, m.transform(src0[m.src[k + kk]]));
    }
    // dst == Q R == P R' == src, dst !~ avg(Q P^t) src'.
    rot *= rotl.QR().transpose() * rotr.QR();
  }
  m.rot  = pow(rot, T(1) / T(int(m.dst.size() - 3)));
  m.rot /= pow(abs(m.rot.determinant()), T(1) / T(3));
  off.O();
  for(int k = 0; k < m.dst.size(); k ++)
    off += dst0[m.dst[k]] - m.transform(src0[m.src[k]]);
  m.offset += (off /= T(int(m.dst.size())));
  T rlog0(int(0));
  for(int k = 0; k < m.dst.size(); k ++) {
    const auto& dstk(dst0[m.dst[k]]);
    const auto  srck(m.transform(src0[m.src[k]]));
    const auto  r(abs(dstk.dot(srck) / srck.dot(srck)));
    if(r != T(0)) rlog0 += log(r);
  }
  m.ratio *= exp(rlog0 / T(int(m.dst.size())));
  off.O();
  for(int k = 0; k < m.dst.size(); k ++)
    off += dst0[m.dst[k]] - m.transform(src0[m.src[k]]);
  m.offset += (off /= T(int(m.dst.size())));
  for(int k = 0; k < m.dst.size(); k ++) {
    const auto& dstk(dst0[m.dst[k]]);
    const auto  srck(m.transform(src0[m.src[k]]));
    const auto  err(dstk - srck);
    m.rdepth += sqrt(err.dot(err) / sqrt(dstk.dot(dstk) * srck.dot(srck)));
  }
  m.rdepth /= m.dst.size() * m.dst.size();
  return m;
}

template <typename T> vector<match_t<T> > matchPartialR(const vector<SimpleVector<T> >& dst0, const vector<SimpleVector<T> >& src0, const int& n = 1) {
  const auto  gs(normalizeG((makeG(dst0).second)));
  const auto  gp(normalizeG((makeG(src0).second)));
  const auto& dst(gs.second);
  const auto& src(gp.second);
  assert(dst.size() == dst0.size());
  assert(src.size() == src0.size());
  cerr << "match(" << dst.size() << ", " << src.size() << ")" << endl;
  SimpleMatrix<T> qdst(dst.size(), 4);
  SimpleMatrix<T> qsrc(src.size(), 4);
  for(int i = 0; i < qdst.rows(); i ++)
    qdst.row(i) = toQuarterNormalize<T>(dst[i]);
  for(int i = 0; i < qsrc.rows(); i ++)
    qsrc.row(i) = toQuarterNormalize<T>(src[i]);
  vector<SimpleVector<T> > test;
  vector<pair<int, int> > idx;
  test.reserve(qdst.rows() * qsrc.rows());
  idx.reserve(qdst.rows() * qsrc.rows());
  for(int i = 0; i < qdst.rows(); i ++)
    for(int j = 0; j < qsrc.rows(); j ++) {
      idx.emplace_back(make_pair(i, j));
      test.emplace_back(qdst.row(i) - qsrc.row(j));
    }
  const auto cr(crush<T>(test, test[0].size(), n));
  vector<match_t<T> > mm;
  mm.reserve(cr.size());
  for(int i = 0; i < cr.size(); i ++) {
    if(! cr[i].first.size()) continue;
    match_t<T> m(T(int(1)) / T(int(100)), max(gs.first, gp.first));
    SimpleVector<int> dfix, sfix;
    dfix.resize(dst.size());
    sfix.resize(src.size());
    dfix.I(false);
    sfix.I(false);
    m.dst.reserve(min(dst.size(), src.size()));
    m.src.reserve(min(dst.size(), src.size()));
    auto avg(cr[i].first[0]);
    for(int j = 1; j < cr[i].first.size(); j ++)
      avg += cr[i].first[j];
    avg /= T(int(cr[i].first.size()));
    vector<pair<T, int> > pp;
    pp.reserve(cr[i].first.size());
    for(int j = 0; j < cr[i].first.size(); j ++) {
      const auto err(cr[i].first[j] - avg);
      pp.emplace_back(make_pair(err.dot(err), j));
    }
    sort(pp.begin(), pp.end());
    for(int j = 0; j < pp.size(); j ++) {
      const auto& lidx(idx[cr[i].second[pp[j].second]]);
      if(dfix[lidx.first] || sfix[lidx.second]) continue;
      dfix[lidx.first] = sfix[lidx.second] = true;
      m.dst.emplace_back(lidx.first);
      m.src.emplace_back(lidx.second);
    }
    if(m.dst.size() < 4) continue;
    m.dst.reserve(m.dst.size());
    m.src.reserve(m.src.size());
    m = reconfigureMatch<T>(m, dst0, src0);
    for(int i = 0; i < m.rot.rows(); i ++)
      for(int j = 0; j < m.rot.cols(); j ++)
        if(! isfinite(m.rot(i, j))) goto nofix;
    for(int i = 0; i < m.offset.size(); i ++)
      if(! isfinite(m.offset[i])) goto nofix;
    if(! isfinite(m.ratio)) goto nofix;
    mm.emplace_back(move(m));
    cerr << mm[mm.size() - 1] << endl;
   nofix:
    ;
  }
  sort(mm.begin(), mm.end());
  return mm;
}

template <typename T> static inline vector<match_t<T> > matchPartial(const vector<SimpleVector<T> >& dst0, const vector<SimpleVector<T> >& src0, const int& n = 1) {
  auto m(matchPartialR<T>(src0, dst0, n));
  for(int i = 0; i < m.size(); i ++)
    m[i] = ~ m[i];
  return m;
}


template <typename T> void drawMatchLine(SimpleMatrix<T>& map, const SimpleVector<T>& lref0, const SimpleVector<T>& lref1, const T& c) {
  int idxm(0);
  int idxM(1);
  if(abs(lref1[idxM] - lref0[idxM]) < abs(lref1[idxm] - lref0[idxm]))
    swap(idxm, idxM);
  const auto d10(lref1 - lref0);
  const auto dlt(abs(lref0[idxM] - lref1[idxM]));
  if(dlt == T(0)) return;
  const auto denom(T(1) / dlt);
  for(int i = 0; i <= int(ceil(dlt)); i ++) {
    const auto gidx(lref0 + d10 * T(i) * denom);
    map(max(int(0), min(int(gidx[0]), int(map.rows() - 1))),
        max(int(0), min(int(gidx[1]), int(map.cols() - 1)))) = c;
  }
  return;
}

template <typename T> void drawMatchTriangle(SimpleMatrix<T>& map, SimpleVector<T> lref0, SimpleVector<T> lref1, SimpleVector<T> lref2, const T& c) {
  // make middle point to lref2 on index 0.
  if((lref0[0] <= lref1[0] && lref1[0] <= lref2[0]) ||
     (lref2[0] <= lref1[0] && lref1[0] <= lref0[0]))
    swap(lref1, lref2);
  if((lref1[0] <= lref0[0] && lref0[0] <= lref2[0]) ||
     (lref2[0] <= lref0[0] && lref0[0] <= lref1[0]) )
    swap(lref0, lref2);
  const auto d0(lref0 - lref2);
  const auto d1(lref1 - lref2);
  const auto d2(lref1 - lref0);
  const auto idx(abs(d2[0]) < abs(d2[1]) ? 1 : 0);
  if(abs(d0[idx]) != T(int(0)) && abs(d2[idx]) != T(int(0)))
    try {
      for(int i = 0; i <= int(abs(d0[idx])); i ++)
        drawMatchLine<T>(map, d0 * (abs(d0[idx]) - T(i)) / abs(d0[idx]) + lref2,
                              d2 * T(i) / abs(d2[idx]) + lref0, c);
    } catch (const char* e) {
      ; /* fall through */
    }
  if(abs(d1[idx]) != T(int(0)) && abs(d2[idx]) != T(int(0)))
    try {
      for(int i = 0; i <= int(abs(d1[idx])); i ++)
        drawMatchLine<T>(map, d1 * (abs(d1[idx]) - T(i)) / abs(d1[idx]) + lref2,
                            - d2 * T(i) / abs(d2[idx]) + lref1, c);
    } catch (const char* e) {
      ; /* fall through */
    }
  return;
}

template <typename T> void addMeshTri(vector<SimpleVector<int> >& res, vector<pair<SimpleVector<T>, int> >& scan, const vector<SimpleVector<T> >& p, const int& idx) {
  assert(0 <= idx && idx < scan.size());
  vector<int> elim;
  if(0 <= idx - 1 &&
     scan[idx].first[0] < scan[idx - 1].first[0] &&
     scan[idx].first[0] < scan[idx + 1].first[0]) {
    elim.emplace_back(idx);
    SimpleVector<int> lres(3);
    lres[0] = scan[idx - 1].second;
    lres[1] = scan[idx].second;
    lres[2] = scan[idx + 1].second;
    bool psize(false);
    for(int k = 0; k < 3; k ++)
      psize = psize || p.size() <= lres[k];
    if(! psize)
      res.emplace_back(move(lres));
  }
  if(idx + 3 < scan.size() &&
     scan[idx + 2].first[0] < scan[idx + 1].first[0] &&
     scan[idx + 2].first[0] < scan[idx + 3].first[0]) {
    elim.emplace_back(idx + 2);
    SimpleVector<int> lres(3);
    lres[0] = scan[idx + 1].second;
    lres[1] = scan[idx + 2].second;
    lres[2] = scan[idx + 3].second;
    bool psize(false);
    for(int k = 0; k < 3; k ++)
      psize = psize || p.size() <= lres[k];
    if(! psize)
      res.emplace_back(move(lres));
  }
  {
    SimpleVector<int> lres(3);
    lres[0] = scan[idx].second;
    lres[1] = scan[idx + 1].second;
    lres[2] = scan[idx + 2].second;
    bool psize(false);
    for(int k = 0; k < 3; k ++)
      psize = psize || p.size() <= lres[k];
    if(! psize)
      res.emplace_back(move(lres));
  }
  sort(elim.begin(), elim.end());
  for(int j = 0; j < elim.size(); j ++)
    scan.erase(scan.begin() + elim[j] - j);
  return;
}

template <typename T> vector<SimpleVector<int> > mesh2(const vector<SimpleVector<T> >& p, const vector<int>& pp) {
  vector<pair<SimpleVector<T>, int> > sp;
  sp.reserve(pp.size());
  T Mxy(int(0));
  for(int i = 0; i < p.size(); i ++)
    Mxy = max(Mxy, max(abs(p[i][0]), abs(p[i][1])));
  Mxy *= T(int(2));
  SimpleMatrix<T> lrot(3, 3);
  lrot.I();
  lrot(0, 0) =    lrot(1, 1) = cos(T(int(1)) / max(Mxy, T(pp.size())));
  lrot(0, 1) = - (lrot(1, 0) = sin(T(int(1)) / max(Mxy, T(pp.size()))));
  T    m1((lrot * p[pp[0]])[1]);
  auto M1(m1);
  for(int i = 0; i < pp.size(); i ++) {
    sp.emplace_back(make_pair(lrot * p[pp[i]], pp[i]));
    sp[i].first[2] = T(0);
    m1 = min(m1, sp[i].first[1]);
    M1 = max(M1, sp[i].first[1]);
  }
  sort(sp.begin(), sp.end(), less0<pair<SimpleVector<T>, int> >);
  vector<pair<SimpleVector<T>, int> > scan;
  scan.reserve(sp.size() + 2);
  scan.emplace_back(sp[0]);
  scan[scan.size() - 1].first[0] -= T(2);
  scan[scan.size() - 1].first[1]  = m1 - T(1);
  scan.emplace_back(sp[0]);
  scan[scan.size() - 1].first[0] -= T(1);
  scan[scan.size() - 1].first[1]  = M1 + T(1);
  vector<SimpleVector<int> > res;
  res.reserve(sp.size());
  for(int i = 0; i < sp.size(); i ++) {
    // N.B. lrot support this on lattice.
    assert(! i || (sp[i].first[0] != sp[i - 1].first[0] &&
                   sp[i].first[1] != sp[i - 1].first[1]) );
    // scanline update
    int idx;
    for(idx = 0; idx < scan.size(); idx ++)
      if(sp[i].first[1] < scan[idx].first[1]) break;
    idx = max(0, min(int(scan.size()) - 2, idx - 1));
    assert(scan[idx].first[1] < sp[i].first[1]);
    assert(sp[i].first[1] < scan[idx + 1].first[1]);
    assert(scan[idx].first[0] < sp[i].first[0]);
    assert(scan[idx + 1].first[0] < sp[i].first[0]);
    scan.insert(scan.begin() + idx + 1, pair<SimpleVector<T>, int>(sp[i]));
    assert(scan[idx].first[1] < scan[idx + 1].first[1] &&
           scan[idx + 1].first[1] < scan[idx + 2].first[1]);
    addMeshTri<T>(res, scan, p, idx);
  }
  while(6 < scan.size()) {
    const auto before(scan.size());
    for(int i = 2; i < scan.size() - 3; i ++)
      addMeshTri<T>(res, scan, p, i);
    if(before == scan.size()) break;
  }
  res.reserve(res.size());
  for(int i = 0; i < res.size(); i ++)
    if(p[res[i][0]][0] * p[res[i][1]][1]
     + p[res[i][1]][0] * p[res[i][2]][1]
     + p[res[i][2]][0] * p[res[i][0]][1]
     - p[res[i][0]][1] * p[res[i][1]][0]
     - p[res[i][1]][1] * p[res[i][2]][0]
     - p[res[i][2]][1] * p[res[i][0]][0] < T(0))
      swap(res[i][0], res[i][1]);
  return res;
}

template <typename T> static inline vector<SimpleVector<int> > mesh2(const vector<SimpleVector<T> >& p) {
  vector<int> pp;
  pp.reserve(p.size());
  for(int i = 0; i < p.size(); i ++) pp.emplace_back(i);
  return mesh2<T>(p, pp);
}

// get bump with multiple scale and vectorized result.
template <typename T> vector<SimpleVector<T> > getTileVec(const SimpleMatrix<T>& in, const int& vbox = 1) {
  vector<SimpleVector<T> > geoms;
  geoms.reserve((in.rows() / vbox + 1) * (in.cols() / vbox + 1));
  // N.B. align with BUMP_BOTH z-axis rxy.
  const auto diag(sqrt(sqrt(T(min(in.rows(), in.cols())) )) );
  for(int i = 0; i < in.rows() / vbox + 1; i ++)
    for(int j = 0; j < in.cols() / vbox + 1; j ++) {
      if(in.rows() < (i + 1) * vbox ||
         in.cols() < (j + 1) * vbox) {
        SimpleVector<T> gbuf(3);
        gbuf[0] = T(i * vbox);
        gbuf[1] = T(j * vbox);
        gbuf[2] = geoms[geoms.size() - 1][2];
        geoms.emplace_back(gbuf);
      } else {
        SimpleVector<T> work(3);
        work[0] = T(i * vbox);
        work[1] = T(j * vbox);
        work[2] = diag * in(i * vbox, j * vbox);
        geoms.emplace_back(work);
      }
    }
  T avg(int(0));
  for(int i = 0; i < geoms.size(); i ++)
    avg += geoms[i][2];
  avg /= T(geoms.size());
  for(int i = 0; i < geoms.size(); i ++)
    geoms[i][2] -= avg;
  geoms.reserve(geoms.size());
  return geoms;
}

template <typename T> vector<SimpleVector<T> > getHesseVec(const SimpleMatrix<T>& in, const int& vbox = 300) {
  const auto guard(max(1, int(sqrt(T(in.rows() * in.cols() / vbox)))));
  vector<SimpleVector<T> > geoms;
  geoms.reserve(vbox + 4);
  const auto x(in * diff<T>(in.cols()).transpose());
  const auto y(diff<T>(in.rows()) * in);
  const auto xx(in * diff<T>(in.cols()).transpose() * diff<T>(in.cols()).transpose());
  const auto xy(diff<T>(in.rows()) * in * diff<T>(in.cols()).transpose());
  const auto yy(diff<T>(in.rows()) * diff<T>(in.rows()) * in);
  // N.B. align with BUMP_BOTH z-axis rxy.
  const auto diag(sqrt(sqrt(T(min(in.rows(), in.cols())) )) );
  vector<pair<T, pair<int, int> > > score;
  score.reserve(in.rows() * in.cols());
  for(int i = 0; i < in.rows(); i ++)
    for(int j = 0; j < in.cols(); j ++)
      score.emplace_back(make_pair(abs((xx(i, j) * yy(i, j) - xy(i, j) * xy(i, j)) / ((T(int(1)) + x(i, j)) * (T(int(1)) + y(i, j)) - x(i, j) * y(i, j))), make_pair(i, j)));
  sort(score.begin(), score.end());
  vector<pair<int, int> > cache;
  cache.reserve(score.size());
  for(int i = score.size() - 1; 0 <= i && geoms.size() < abs(vbox); i --)
    if(! binary_search(cache.begin(), cache.end(),
           make_pair(score[i].second.first / guard,
                     score[i].second.second / guard)) ) {
      SimpleVector<T> g(3);
      g[0] = T(int(score[i].second.first));
      g[1] = T(int(score[i].second.second));
      g[2] = diag *
        in(score[i].second.first, score[i].second.second);
      geoms.emplace_back(move(g));
      cache.emplace_back(make_pair(score[i].second.first / guard,
                                   score[i].second.second / guard));
      sort(cache.begin(), cache.end());
    }
  SimpleVector<T> g(3);
  g[0] = T(int(0));
  g[1] = T(int(0));
  g[2] = diag * in(0, 0);
  geoms.emplace_back(g);
  g[0] = T(int(in.rows() - 1));
  g[1] = T(int(0));
  g[2] = diag * in(in.rows() - 1, 0);
  geoms.emplace_back(g);
  g[0] = T(int(0));
  g[1] = T(int(in.cols() - 1));
  g[2] = diag * in(0, in.cols() - 1);
  geoms.emplace_back(g);
  g[0] = T(int(in.rows() - 1));
  g[1] = T(int(in.cols() - 1));
  g[2] = diag * in(in.rows() - 1, in.cols() - 1);
  geoms.emplace_back(g);
  T avg(int(0));
  for(int i = 0; i < geoms.size(); i ++)
    avg += geoms[i][2];
  avg /= T(geoms.size());
  for(int i = 0; i < geoms.size(); i ++)
    geoms[i][2] -= avg;
  geoms.reserve(geoms.size());
  return geoms;
}

template <typename T> SimpleMatrix<T> tilt(const SimpleMatrix<T>& in, vector<triangles_t<T> >& triangles, const T& depth = - T(10000)) {
  cerr << "t" << flush;
  SimpleMatrix<T> result(in.rows(), in.cols());
  result.O();
  vector<pair<T, triangles_t<T>> > zbuf;
  zbuf.resize(triangles.size(), make_pair(T(0), triangles_t<T>()));
  assert(zbuf.size() == triangles.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int j = 0; j < triangles.size(); j ++) {
    auto& tri(triangles[j]);
    assert(tri.first.rows() == tri.first.cols() && tri.first.rows() == 3);
    // N.B. /= 3 isn't needed because only the order is the matter.
    zbuf[j].first  = - (tri.first(0, 2) + tri.first(1, 2) + tri.first(2, 2));
    zbuf[j].second = move(tri);
  }
  sort(zbuf.begin(), zbuf.end(), lessf<pair<T, triangles_t<T> > >);
  int i;
  // XXX: patent???
  // N.B. we could avoid with this because no z-buffer matrix on them,
  //      but this is obscure.
  for(i = 0; i < zbuf.size() && zbuf[i].first < depth; i ++) ;
  for( ; i < zbuf.size(); i ++) {
    const auto& zbi(zbuf[i].second);
    drawMatchTriangle<T>(result, zbi.first.row(0), zbi.first.row(1), zbi.first.row(2), zbi.second);
  }
  return result;
}

template <typename T> static inline SimpleMatrix<T> tilt(const SimpleMatrix<T>& in, const vector<triangles_t<T> >& triangles, const T& depth = - T(10000)) {
  auto tris(triangles);
  return tilt<T>(in, tris, depth);
}

template <typename T> vector<triangles_t<T> > triangles(const SimpleMatrix<T>& in, const SimpleMatrix<T>& bump, const match_t<T>& m) {
  assert(in.rows() == bump.rows() && in.cols() == bump.cols());
  auto points(getTileVec<T>(bump));
  auto facets(mesh2<T>(points));
  vector<triangles_t<T> > triangles;
  triangles.resize(facets.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < facets.size(); i ++) {
    triangles_t<T> work;
    work.first = SimpleMatrix<T>(3, 3);
    for(int j = 0; j < 3; j ++) {
      assert(0 <= facets[i][j] && facets[i][j] < points.size());
      work.first.row(j) = m.transform(points[facets[i][j]]);
    }
    for(int j = 0; j < 2; j ++) {
      if(work.first(0, j) <= work.first(1, j) && work.first(0, j) <= work.first(2, j))
        work.first(0, j) = floor(work.first(0, j));
      else if(work.first(1, j) <= work.first(0, j) && work.first(1, j) <= work.first(2, j))
        work.first(1, j) = floor(work.first(1, j));
      else if(work.first(2, j) <= work.first(0, j) && work.first(2, j) <= work.first(1, j))
        work.first(2, j) = floor(work.first(2, j));
      if(work.first(1, j) <= work.first(0, j) && work.first(2, j) <= work.first(0, j))
        work.first(0, j) = ceil(work.first(0, j));
      else if(work.first(0, j) <= work.first(1, j) && work.first(2, j) <= work.first(1, j))
        work.first(1, j) = ceil(work.first(1, j));
      else if(work.first(0, j) <= work.first(2, j) && work.first(1, j) <= work.first(2, j))
        work.first(2, j) = ceil(work.first(2, j));
    }
    if(T(0) <= points[facets[i][0]][0] && points[facets[i][0]][0] < T(in.rows()) &&
       T(0) <= points[facets[i][0]][1] && points[facets[i][0]][1] < T(in.cols()))
      work.second = in(int(points[facets[i][0]][0]),
                       int(points[facets[i][0]][1]));
    else
      work.second = T(0);
    triangles[i] = move(work);
  }
  return triangles;
}

template <typename T> SimpleMatrix<T> draw(const SimpleMatrix<T>& img, const vector<SimpleVector<T> >& shape, const vector<SimpleVector<T> >& emph, const vector<SimpleVector<int> >& hull) {
  assert(shape.size() == emph.size());
  vector<triangles_t<T> > tris;
  tris.resize(hull.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < hull.size(); i ++) {
    assert(hull[i].size() == 3);
    assert(0 <= hull[i][0] && hull[i][0] < shape.size());
    assert(0 <= hull[i][1] && hull[i][1] < shape.size());
    assert(0 <= hull[i][2] && hull[i][2] < shape.size());
    triangles_t<T> work;
    work.first = SimpleMatrix<T>(3, 3);
    for(int j = 0; j < 3; j ++)
      work.first.row(j) = emph[hull[i][j]];
    work.second = img(max(int(0), min(int(img.rows() - 1),
                          int(shape[hull[i][0]][0]))),
                      max(int(0), min(int(img.cols() - 1),
                          int(shape[hull[i][0]][1]))));
    tris[i] = move(work);
  }
  return tilt<T>(img * T(0), tris);
}

template <typename T> SimpleMatrix<T> draw(const SimpleMatrix<T>& img, const vector<SimpleVector<T> >& shape, const vector<SimpleVector<int> >& hull, const bool& elim = false) {
  auto result(img);
  T    M(0);
  T    m(0);
  for(int i = 0; i < shape.size(); i ++) {
    if(i) {
      M = max(M, shape[i][2]);
      m = min(m, shape[i][2]);
    } else
      M = m = shape[i][2];
  }
  auto tsrc(shape);
  if(M - m != T(0))
    for(int i = 0; i < tsrc.size(); i ++)
      tsrc[i][2] = elim ? T(0) : (tsrc[i][2] - m) / (M - m);
  for(int ii = 0; ii < hull.size(); ii ++)
    drawMatchTriangle<T>(result, tsrc[hull[ii][0]],
                                 tsrc[hull[ii][1]],
                                 tsrc[hull[ii][2]],
                                 (tsrc[hull[ii][0]][2] + tsrc[hull[ii][1]][2] +
                                  tsrc[hull[ii][2]][2]) / T(3));
  return result;
}

template <typename T> static inline vector<SimpleVector<T> > takeShape(const vector<SimpleVector<T> >& dst, const vector<SimpleVector<T> >& src, const match_t<T>& match, const T& ratio) {
  auto result(dst);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < match.src.size(); i ++)
    result[match.dst[i]] += (match.transform(src[match.src[i]]) - dst[match.dst[i]]) * ratio;
  return result;
}

template <typename T> static inline SimpleMatrix<T> showMatch(const SimpleMatrix<T>& dstimg, const vector<SimpleVector<T> >& dst, const vector<SimpleVector<int> >& hull, const T& emph = T(1)) {
  auto map(dstimg);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int k = 0; k < hull.size(); k ++) {
    drawMatchLine<T>(map, dst[hull[k][0]], dst[hull[k][1]], emph);
    drawMatchLine<T>(map, dst[hull[k][1]], dst[hull[k][2]], emph);
    drawMatchLine<T>(map, dst[hull[k][2]], dst[hull[k][0]], emph);
  }
  return map;
}

template <typename T> static inline SimpleMatrix<T> makeRefMatrix(const SimpleMatrix<T>& orig, const int& start) {
  SimpleMatrix<T> result(orig.rows(), orig.cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < orig.rows() * orig.cols(); i ++)
    result(i % orig.rows(), i / orig.rows()) = i + start;
  return result;
}

template <typename T> static inline SimpleMatrix<T> pullRefMatrix(const SimpleMatrix<T>& ref, const int& start, const SimpleMatrix<T>& orig) {
  assert(orig.rows() == ref.rows() && orig.cols() == ref.cols());
  SimpleMatrix<T> result(ref.rows(), ref.cols());
  result.O();
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < ref.rows() * ref.cols(); i ++) {
    const int ly(i % ref.rows());
    const int lx(i / ref.rows());
    const int v(int(ref(ly, lx)) - start);
    if(0 <= v && v < orig.rows() * orig.cols())
      result(ly, lx) = orig(v % orig.rows(), v / orig.rows());
    else
      result(ly, lx) = T(0);
  }
  return result;
}

template <typename T> SimpleMatrix<T> reShape(const SimpleMatrix<T>& cbase, const SimpleMatrix<T>& vbase, const int& count, const T& thresh) {
  assert(cbase.rows() && cbase.cols() && vbase.rows() && vbase.cols());
  assert(cbase.rows() == vbase.rows() && cbase.cols() == vbase.cols());
  vector<pair<T, pair<int, int> > > vpoints;
  vpoints.reserve(vbase.rows() * vbase.cols());
  for(int i = 0; i < vbase.rows(); i ++)
    for(int j = 0; j < vbase.cols(); j ++)
      vpoints.emplace_back(make_pair(vbase(i, j), make_pair(i, j)));
  sort(vpoints.begin(), vpoints.end());
  SimpleMatrix<T> res(vbase.rows(), vbase.cols());
  T avg(0);
  for(int i = 0, ii = 0; i < vpoints.size(); i ++) {
    if(abs(vpoints[i].first - vpoints[ii].first) < abs(vpoints[vpoints.size() - 1].first - vpoints[0].first) / T(count) && i < vpoints.size() - 1)
      avg += cbase(vpoints[i].second.first, vpoints[i].second.second);
    else if(i != ii) {
      if(i == vpoints.size() - 1) {
        avg += cbase(vpoints[i].second.first, vpoints[i].second.second);
        i ++;
      }
      avg /= i - ii;
      for(int j = ii; j < i; j ++)
        res(vpoints[j].second.first, vpoints[j].second.second) = avg;
      avg  = T(0);
      ii   = i;
    }
  }
  SimpleMatrix<T> masker(res);
  SimpleMatrix<T> mask(res.rows(), res.cols());
  mask.O(false);
  for(int i = 0; i < mask.rows(); i ++)
    for(int j = 0; j < mask.cols(); j ++)
      if(! mask(i, j)) {
        vector<pair<int, int> > store;
        vector<pair<int, int> > tries;
        tries.emplace_back(make_pair(+ 1,   0));
        tries.emplace_back(make_pair(  0, + 1));
        tries.emplace_back(make_pair(- 1,   0));
        tries.emplace_back(make_pair(  0, - 1));
        vector<pair<int, int> > stack;
        stack.emplace_back(make_pair(i, j));
        T   avg(0);
        int cnt(0);
        while(stack.size()) {
          const auto pop(stack[stack.size() - 1]);
          stack.pop_back();
          const int& yy(pop.first);
          const int& xx(pop.second);
          if(0 <= yy && yy < mask.rows() && 0 <= xx && xx < mask.cols() &&
             masker(i, j) == masker(yy, xx) &&
             abs(cbase(i, j) - cbase(yy, xx)) <= thresh &&
             ! mask(yy, xx)) {
            mask(yy, xx) = true;
            store.emplace_back(make_pair(yy, xx));
            for(int ii = 0; ii < tries.size(); ii ++)
              stack.emplace_back(make_pair(yy + tries[ii].first, xx + tries[ii].second));
            avg += cbase(yy, xx);
            cnt ++;
          }
        }
        avg /= T(cnt);
        for(int i = 0; i < store.size(); i ++)
          res(store[i].first, store[i].second) = avg;
      }
  return res;
}

template <typename T> SimpleMatrix<T> reColor(const SimpleMatrix<T>& cbase, const SimpleMatrix<T>& vbase, const int& count, const T& intensity) {
  assert(cbase.rows() && cbase.cols() && vbase.rows() && vbase.cols());
  vector<pair<T, pair<int, int> > > vpoints;
  vector<pair<T, pair<int, int> > > cpoints;
  vpoints.reserve(vbase.rows() * vbase.cols());
  cpoints.reserve(cbase.rows() * cbase.cols());
  for(int i = 0; i < vbase.rows(); i ++)
    for(int j = 0; j < vbase.cols(); j ++)
      vpoints.emplace_back(make_pair(vbase(i, j), make_pair(i, j)));
  for(int i = 0; i < cbase.rows(); i ++)
    for(int j = 0; j < cbase.cols(); j ++)
      cpoints.emplace_back(make_pair(cbase(i, j), make_pair(i, j)));
  sort(vpoints.begin(), vpoints.end());
  sort(cpoints.begin(), cpoints.end());
  SimpleVector<T> vv(vpoints.size());
  SimpleVector<T> cc(cpoints.size());
  for(int i = 0; i < vv.size(); i ++)
    vv[i] = vpoints[i].first;
  for(int i = 0; i < cc.size(); i ++)
    cc[i] = cpoints[i].first;
  const auto ccc(Decompose<T>(count).mimic(cc, vv, intensity));
        auto res(cbase);
  for(int i = 0; i < ccc.size(); i ++)
    res(cpoints[i].second.first, cpoints[i].second.second) = ccc[i];
  return res;
}

template <typename T> SimpleMatrix<T> reColor3(const SimpleMatrix<T>& ccbase, const SimpleMatrix<T>& vbase, const int& count) {
  assert(ccbase.rows() && ccbase.cols() && vbase.rows() && vbase.cols());
  auto cbase(ccbase);
  for(int i = 0; i < cbase.rows(); i ++)
    for(int j = 0; j < cbase.cols(); j ++)
      cbase(i, j) += T(int(1)) / T(int(256));
  vector<pair<T, pair<int, int> > > vpoints;
  vector<pair<T, pair<int, int> > > cpoints;
  vpoints.reserve(vbase.rows() * vbase.cols());
  cpoints.reserve(cbase.rows() * cbase.cols());
  for(int i = 0; i < vbase.rows(); i ++)
    for(int j = 0; j < vbase.cols(); j ++)
      vpoints.emplace_back(make_pair(vbase(i, j), make_pair(i, j)));
  for(int i = 0; i < cbase.rows(); i ++)
    for(int j = 0; j < cbase.cols(); j ++)
      cpoints.emplace_back(make_pair(cbase(i, j), make_pair(i, j)));
  sort(vpoints.begin(), vpoints.end());
  sort(cpoints.begin(), cpoints.end());
  SimpleMatrix<T> res(cbase.rows(), cbase.cols());
  for(int i = 0; i < count; i ++) {
    T scorev(0);
    T scorec(0);
    const auto vstart(i * int(vpoints.size() / count));
    const auto vend(min((i + 1) * int(vpoints.size() / count), int(vpoints.size())));
    const auto cstart(i * int(cpoints.size() / count));
    const auto cend(min((i + 1) * int(cpoints.size() / count), int(cpoints.size())));
    for(int j = vstart; j < vend; j ++)
      scorev += vpoints[j].first;
    for(int j = cstart; j < cend; j ++)
      scorec += cpoints[j].first;
    scorev /= T(vend - vstart);
    scorec /= T(cend - cstart);
    for(int j = cstart; j < cend; j ++)
      res(cpoints[j].second.first, cpoints[j].second.second) =
        cpoints[j].first * scorev / scorec;
  }
  return res;
}

template <typename T> SimpleMatrix<T> reColor(const SimpleMatrix<T>& cbase, const int& count, const T& intensity) {
  assert(cbase.rows() && cbase.cols());
  vector<pair<T, pair<int, int> > > cpoints;
  cpoints.reserve(cbase.rows() * cbase.cols());
  for(int i = 0; i < cbase.rows(); i ++)
    for(int j = 0; j < cbase.cols(); j ++)
      cpoints.emplace_back(make_pair(cbase(i, j), make_pair(i, j)));
  sort(cpoints.begin(), cpoints.end());
  SimpleVector<T> cc(cpoints.size());
  for(int i = 0; i < cc.size(); i ++)
    cc[i] = cpoints[i].first;
  const auto ccc(Decompose<T>(count).emphasis(cc, intensity));
        auto res(cbase);
  for(int i = 0; i < ccc.size(); i ++)
    res(cpoints[i].second.first, cpoints[i].second.second) = ccc[i];
  return res;
}

template <typename T> vector<vector<int> > catImage(const vector<SimpleMatrix<T> >& imgs, const int& cs = 40) {
  for(int i = 1; i < imgs.size(); i ++) {
    assert(imgs[i].rows() == imgs[0].rows());
    assert(imgs[i].cols() == imgs[0].cols());
  }
  vector<SimpleVector<T> > work;
  work.reserve(imgs.size());
  for(int i = 0; i < imgs.size(); i ++) {
    work.emplace_back(SimpleVector<T>(imgs[i].rows() * imgs[i].cols()).O());
    for(int j = 0; j < imgs[i].rows(); j ++)
      work[i].setVector(imgs[i].cols() * j, imgs[i].row(j));
  }
  auto cg(crush<T>(work, cs, 0));
  vector<vector<int> > res;
  res.reserve(cg.size());
  for(int i = 0; i < cg.size(); i ++)
    res.emplace_back(move(cg[i].second));
  return res;
}

template <typename T> static inline vector<SimpleMatrix<T> > rgb2xyz(const vector<SimpleMatrix<T> >& rgb) {
  // CIE 1931 XYZ from wikipedia.org
  SimpleMatrix<T> mRGB2XYZ(3, 3);
  mRGB2XYZ(0, 0) = T(49000);
  mRGB2XYZ(0, 1) = T(31000);
  mRGB2XYZ(0, 2) = T(20000);
  mRGB2XYZ(1, 0) = T(17697);
  mRGB2XYZ(1, 1) = T(81240);
  mRGB2XYZ(1, 2) = T( 1063);
  mRGB2XYZ(2, 0) = T(0);
  mRGB2XYZ(2, 1) = T( 1000);
  mRGB2XYZ(2, 2) = T(99000);
  mRGB2XYZ /= T(17697);
  assert(rgb.size() == 3);
  assert(rgb[0].rows() == rgb[1].rows() && rgb[1].rows() == rgb[2].rows());
  assert(rgb[0].cols() == rgb[1].cols() && rgb[1].cols() == rgb[2].cols());
  auto xyz(rgb);
  xyz[0] = rgb[0] * mRGB2XYZ(0, 0) + rgb[1] * mRGB2XYZ(0, 1) + rgb[2] * mRGB2XYZ(0, 2);
  xyz[1] = rgb[0] * mRGB2XYZ(1, 0) + rgb[1] * mRGB2XYZ(1, 1) + rgb[2] * mRGB2XYZ(1, 2);
  xyz[2] = rgb[0] * mRGB2XYZ(2, 0) + rgb[1] * mRGB2XYZ(2, 1) + rgb[2] * mRGB2XYZ(2, 2);
  assert(xyz.size() == 3);
  assert(xyz[0].rows() == xyz[1].rows() && xyz[1].rows() == xyz[2].rows());
  assert(xyz[0].cols() == xyz[1].cols() && xyz[1].cols() == xyz[2].cols());
  return xyz;
}

template <typename T> static inline vector<SimpleMatrix<T> > xyz2rgb(const vector<SimpleMatrix<T> >& xyz) {
  // CIE 1931 XYZ from wikipedia.org
  SimpleMatrix<T> mRGB2XYZ(3, 3);
  mRGB2XYZ(0, 0) = T(49000);
  mRGB2XYZ(0, 1) = T(31000);
  mRGB2XYZ(0, 2) = T(20000);
  mRGB2XYZ(1, 0) = T(17697);
  mRGB2XYZ(1, 1) = T(81240);
  mRGB2XYZ(1, 2) = T( 1063);
  mRGB2XYZ(2, 0) = T(0);
  mRGB2XYZ(2, 1) = T( 1000);
  mRGB2XYZ(2, 2) = T(99000);
  mRGB2XYZ /= T(17697);
  const auto mXYZ2RGB(mRGB2XYZ.inverse());
  assert(xyz.size() == 3);
  assert(xyz[0].rows() == xyz[1].rows() && xyz[1].rows() == xyz[2].rows());
  assert(xyz[0].cols() == xyz[1].cols() && xyz[1].cols() == xyz[2].cols());
  auto rgb(xyz);
  rgb[0] = xyz[0] * mXYZ2RGB(0, 0) + xyz[1] * mXYZ2RGB(0, 1) + xyz[2] * mXYZ2RGB(0, 2);
  rgb[1] = xyz[0] * mXYZ2RGB(1, 0) + xyz[1] * mXYZ2RGB(1, 1) + xyz[2] * mXYZ2RGB(1, 2);
  rgb[2] = xyz[0] * mXYZ2RGB(2, 0) + xyz[1] * mXYZ2RGB(2, 1) + xyz[2] * mXYZ2RGB(2, 2);
  assert(rgb.size() == 3);
  assert(rgb[0].rows() == rgb[1].rows() && rgb[1].rows() == rgb[2].rows());
  assert(rgb[0].cols() == rgb[1].cols() && rgb[1].cols() == rgb[2].cols());
  return rgb;
}

template <typename T> static inline SimpleMatrix<T> rgb2d(const vector<SimpleMatrix<T> > rgb) {
  auto xyz(rgb2xyz<T>(rgb));
  SimpleMatrix<T> result(rgb[0].rows(), rgb[0].cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int j = 0; j < rgb[0].rows(); j ++) {
    for(int k = 0; k < rgb[0].cols(); k ++)
      result(j, k) = sqrt(xyz[0](j, k) * xyz[0](j, k) + xyz[1](j, k) * xyz[1](j, k) + xyz[2](j, k) * xyz[2](j, k));
  }
  return result;
}

template <typename T> static inline match_t<T> tiltprep(const SimpleMatrix<T>& in, const int& idx, const int& samples, const T& psi, const T& z0 = T(int(0))) {
  const auto Pi(atan2(T(int(1)), T(int(1))) * T(int(4)));
  const auto theta(T(2) * Pi * T(idx) / T(samples));
  const auto lpsi(Pi * psi);
  SimpleMatrix<T> R0(3, 3);
  SimpleMatrix<T> R1(3, 3);
  R0(0, 0) =   cos(theta);
  R0(0, 1) =   sin(theta);
  R0(0, 2) = 0.;
  R0(1, 0) = - sin(theta);
  R0(1, 1) =   cos(theta);
  R0(1, 2) = 0.;
  R0(2, 0) = 0.;
  R0(2, 1) = 0.;
  R0(2, 2) = 1.;
  R1(0, 0) = 1.;
  R1(0, 1) = 0.;
  R1(0, 2) = 0.;
  R1(1, 0) = 0.;
  R1(1, 1) = 1.;
  R1(1, 2) = 0.;
  R1(2, 0) = 0.;
  R1(2, 1) = 0.;
  R1(2, 2) = 1.;
  R1(0, 0) =   cos(lpsi);
  R1(0, 2) = - sin(lpsi);
  R1(2, 0) =   sin(lpsi);
  R1(2, 2) =   cos(lpsi);
  match_t<T> m;
  m.rot    = R0.transpose() * R1 * R0;
  SimpleVector<T> pcenter(3);
  pcenter[0] = T(in.rows() - 1) / T(2);
  pcenter[1] = T(in.cols() - 1) / T(2);
  pcenter[2] = z0;
  // x -> m.rot * x, same center
  // x - origin -> m.rot * (x - origin)
  // x -> m.rot * x - m.rot * origin + origin.
  m.offset = pcenter - m.rot * pcenter;
  m.ratio  = T(1);
  return m;
}

#define _GOKICHECK_
#endif

