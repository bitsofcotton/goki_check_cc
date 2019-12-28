#if !defined(_P0_)

template <typename T> class P0 {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleVector<complex<T> > VecU;
  typedef SimpleMatrix<T> Mat;
  typedef SimpleMatrix<complex<T> > MatU;
  inline P0();
  inline P0(const int& range, const int& look = 1);
  inline ~P0();
  inline void nextVoid(const T& in);
  inline T    next(const T& in);
private:
  int look;
  Vec buf;
  const MatU& seed(const int& size, const bool& idft);
  const Mat&  diff(const int& size, const bool& integrate);
  const Mat&  lhpf(const int& size, const bool& hpf);
  const Vec&  nextTaylor(const int& size, const int& step);
  const T& Pi() const;
  const complex<T>& J() const;
};

template <typename T> inline P0<T>::P0() {
  look = 0;
}

template <typename T> inline P0<T>::P0(const int& range, const int& look) {
  assert(1 < range && 0 < look);
  buf.resize(range);
  for(int i = 0; i < buf.size(); i ++)
    buf[i] = T(0);
  this->look = look;
}

template <typename T> inline P0<T>::~P0() {
  ;
}

template <typename T> inline void P0<T>::nextVoid(const T& in) {
  if(in == buf[buf.size() - 1] || in == T(0))
    return;
  for(int i = 1; i < buf.size(); i ++)
    buf[i - 1] = buf[i];
  buf[buf.size() - 1] = in;
}

template <typename T> inline T P0<T>::next(const T& in) {
  nextVoid(in);
  return nextTaylor(buf.size(), look).dot(lhpf(buf.size(), false) * buf);
}

template <typename T> const T& P0<T>::Pi() const {
  const static auto pi(atan2(T(1), T(1)) * T(4));
  return pi;
}

template <typename T> const complex<T>& P0<T>::J() const {
  const static auto i(complex<T>(T(0), T(1)));
  return i;
}

template <typename T> const typename P0<T>::MatU& P0<T>::seed(const int& size, const bool& f_idft) {
  assert(0 < size);
  static vector<MatU> dft;
  static vector<MatU> idft;
  if(dft.size() <= size)
    dft.resize(size + 1, MatU());
  if(idft.size() <= size)
    idft.resize(size + 1, MatU());
  if((!f_idft) &&  dft[size].rows() == size &&  dft[size].cols() == size)
    return dft[size];
  if(  f_idft  && idft[size].rows() == size && idft[size].cols() == size)
    return idft[size];
  auto& edft( dft[size]);
  auto& eidft(idft[size]);
  edft.resize(size, size);
  eidft.resize(size, size);
  for(int i = 0; i < edft.rows(); i ++)
    for(int j = 0; j < edft.cols(); j ++) {
      const auto theta(- T(2) * Pi() * T(i) * T(j) / T(edft.rows()));
      edft(i, j)  = complex<T>(cos(  theta), sin(  theta));
      eidft(i, j) = complex<T>(cos(- theta), sin(- theta)) / complex<T>(T(size));
    }
  if(f_idft)
    return eidft;
  return edft;
}

template <typename T> const typename P0<T>::Mat& P0<T>::diff(const int& size, const bool& integrate) {
  assert(0 < size);
  static vector<Mat> D;
  static vector<Mat> I;
  if(D.size() <= size)
    D.resize(size + 1, Mat());
  if(I.size() <= size)
    I.resize(size + 1, Mat());
  if((!integrate) && D[size].rows() == size && D[size].cols() == size)
    return D[size];
  if(  integrate  && I[size].rows() == size && I[size].cols() == size)
    return I[size];
  auto& d(D[size]);
  auto& i(I[size]);
  d.resize(size, size);
  i.resize(size, size);
  auto  DD(seed(size, false));
  DD.row(0) *= complex<T>(T(0));
  auto  II(DD);
  T nd(0);
  T ni(0);
  for(int i = 1; i < DD.rows(); i ++) {
    const auto phase(- J() * T(2) * Pi() * T(i) / T(DD.rows()));
    const auto phase2(complex<T>(T(1)) / phase);
    DD.row(i) *= phase;
    II.row(i) /= phase;
    nd += abs(phase)  * abs(phase);
    ni += abs(phase2) * abs(phase2);
  }
  d = (seed(size, true) * DD).template real<T>() * sqrt(sqrt(T(DD.rows() - 1) / (nd * ni)));
  i = (seed(size, true) * II).template real<T>() * sqrt(sqrt(T(II.rows() - 1) / (nd * ni)));
  if(integrate)
    return i;
  return d;
}

template <typename T> const typename P0<T>::Mat& P0<T>::lhpf(const int& size, const bool& hpf) {
  assert(0 < size);
  static vector<Mat> L;
  static vector<Mat> H;
  if(L.size() <= size)
    L.resize(size + 1, Mat());
  if(H.size() <= size)
    H.resize(size + 1, Mat());
  if((!hpf) && L[size].rows() == size && L[size].cols() == size)
    return L[size];
  if(  hpf  && H[size].rows() == size && H[size].cols() == size)
    return H[size];
  auto& l(L[size]);
  auto& h(H[size]);
  auto  ll(seed(size, false));
  auto  hh(seed(size, false));
  for(int i = 0; i < size / 2; i ++)
    hh.row(i) *= complex<T>(T(0));
  for(int i = size / 2; i < size; i ++)
    ll.row(i) *= complex<T>(T(0));
  l = (seed(size, true) * ll).template real<T>();
  h = (seed(size, true) * hh).template real<T>();
  if(hpf)
    return h;
  return l;
}

template <typename T> const typename P0<T>::Vec& P0<T>::nextTaylor(const int& size, const int& step) {
  assert(0 < size);
  static vector<vector<Vec> > ntayl;
  if(ntayl.size() <= size)
    ntayl.resize(size + 1, vector<Vec>());
  if(ntayl[size].size() <= step)
    ntayl[size].resize(step + 1, Vec());
  if(ntayl[size][step].size() == size)
    return ntayl[size][step];
        auto& v(ntayl[size][step]);
  const auto& D(diff(size, false));
        auto  dd(D * T(step));
  v.resize(size);
  for(int i = 0; i < v.size() - 1; i ++)
    v[i] = T(0);
  v[v.size() - 1] = T(1);
  for(int i = 2; 0 <= i; i ++) {
    const auto bv(v);
    v += dd.row(dd.rows() - 1);
    if(bv == v)
      break;
    dd = (D * dd) * (T(step) / T(i));
  }
  return v;
}

#define _P0_
#endif

