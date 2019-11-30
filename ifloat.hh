#if !defined(_INTEGER_FLOAT_)

using std::move;
using std::max;
using std::min;
using std::vector;

// Double int to new int class.
template <typename T, int bits> class DUInt {
public:
  inline DUInt();
  inline DUInt(const int& src);
  inline DUInt(const T& src);
  inline DUInt(const DUInt<T,bits>& src);
  inline DUInt(const DUInt<DUInt<T,bits>,bits*2>& src);
  inline DUInt(DUInt<T,bits>&& src);
  inline ~DUInt();
  
  inline DUInt<T,bits>& operator ++ ();
  inline DUInt<T,bits>  operator ++ (int);
  inline DUInt<T,bits>& operator -- ();
  inline DUInt<T,bits>  operator -- (int);
  inline DUInt<T,bits>  operator -  () const;
  inline DUInt<T,bits>  operator +  (const DUInt<T,bits>& src) const;
  inline DUInt<T,bits>& operator += (const DUInt<T,bits>& src);
  inline DUInt<T,bits>  operator -  (const DUInt<T,bits>& src) const;
  inline DUInt<T,bits>& operator -= (const DUInt<T,bits>& src);
  inline DUInt<T,bits>  operator *  (const DUInt<T,bits>& src) const;
  inline DUInt<T,bits>& operator *= (const DUInt<T,bits>& src);
  inline DUInt<T,bits>  operator /  (const DUInt<T,bits>& src) const;
  inline DUInt<T,bits>& operator /= (const DUInt<T,bits>& src);
  inline DUInt<T,bits>  operator %  (const DUInt<T,bits>& src) const;
  inline DUInt<T,bits>& operator %= (const DUInt<T,bits>& src);
  inline DUInt<T,bits>  operator << ( const int& b)            const;
  inline DUInt<T,bits>& operator <<= (const int& b);
  inline DUInt<T,bits>  operator >> ( const int& b)            const;
  inline DUInt<T,bits>& operator >>= (const int& b);
  inline DUInt<T,bits>  operator &  (const DUInt<T,bits>& src) const;
  inline DUInt<T,bits>& operator &= (const DUInt<T,bits>& src);
  inline DUInt<T,bits>  operator |  (const DUInt<T,bits>& src) const;
  inline DUInt<T,bits>& operator |= (const DUInt<T,bits>& src);
  inline DUInt<T,bits>  operator ^  (const DUInt<T,bits>& src) const;
  inline DUInt<T,bits>& operator ^= (const DUInt<T,bits>& src);
  inline DUInt<T,bits>  operator ~  ()                         const;
  inline DUInt<T,bits>& operator =  (const DUInt<T,bits>& src);
  inline DUInt<T,bits>& operator =  (const DUInt<DUInt<T,bits>,bits*2>& src);
  inline DUInt<T,bits>& operator =  (const int& src);
  inline DUInt<T,bits>& operator =  (DUInt<T,bits>&& src);
  inline bool           operator <  (const DUInt<T,bits>& src) const;
  inline bool           operator <= (const DUInt<T,bits>& src) const;
  inline bool           operator >  (const DUInt<T,bits>& src) const;
  inline bool           operator >= (const DUInt<T,bits>& src) const;
  inline bool           operator == (const DUInt<T,bits>& src) const;
  inline bool           operator != (const DUInt<T,bits>& src) const;
  inline bool           operator && (const DUInt<T,bits>& src) const;
  inline bool           operator || (const DUInt<T,bits>& src) const;
  inline bool           operator !    () const;
  inline                operator bool () const;
  inline                operator int  () const;
  inline                operator T    () const;
  inline                operator DUInt<T,bits> () const;

/*
  friend std::ostream&  operator << (std::ostream& os, DUInt<T,bits>  v);
  friend std::istream&  operator >> (std::istream& is, DUInt<T,bits>& v);
*/

  T e[2];
private:
  inline int getMSB() const;
};

template <typename T, int bits> inline DUInt<T,bits>::DUInt() {
  assert(0 < bits && ! (bits & 3));
}

template <typename T, int bits> inline DUInt<T,bits>::DUInt(const int& src) {
  e[0]   = src;
  e[1]  ^= e[1];
}

template <typename T, int bits> inline DUInt<T,bits>::DUInt(const T& src) {
  e[0]   = src;
  e[1]  ^= e[1];
}

template <typename T, int bits> inline DUInt<T,bits>::DUInt(const DUInt<T,bits>& src) {
  *this = src;
}

template <typename T, int bits> inline DUInt<T,bits>::DUInt(const DUInt<DUInt<T,bits>,bits*2>& src) {
  *this = src;
}

template <typename T, int bits> inline DUInt<T,bits>::DUInt(DUInt<T,bits>&& src) {
  *this = src;
}

template <typename T, int bits> inline DUInt<T,bits>::~DUInt() {
  ;
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator ++ () {
  ++ e[0];
  if(!e[0])
    ++ e[1];
  return *this;
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator ++ (int) {
  const auto work(*this);
  ++ *this;
  return work;
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator -- () {
  if(!e[0])
    -- e[1];
  -- e[0];
  return *this;
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator -- (int) {
  const auto work(*this);
  -- *this;
  return work;
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator -  () const {
  auto work(~ *this);
  return ++ work;
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator +  (const DUInt<T,bits>& src) const {
  auto work(*this);
  return work += src;
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator += (const DUInt<T,bits>& src) {
  // N.B. assembler can boost dramatically this code. but not here.
  const auto e0(max(e[0], src.e[0]));
  e[0] += src.e[0];
  if(e[0] < e0)
    e[1] ++;
  e[1] += src.e[1];
  return *this;
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator -  (const DUInt<T,bits>& src) const {
  auto work(*this);
  return work -= src;
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator -= (const DUInt<T,bits>& src) {
  return *this += - src;
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator *  (const DUInt<T,bits>& src) const {
  const static auto hbits(bits >> 1);
  const static auto lmask((T(1) << hbits) - T(1));
  const auto d0(e[0] & lmask);
  const auto d1(e[0] >> hbits);
  const auto d2(e[1] & lmask);
  const auto d3(e[1] >> hbits);
  const auto s0(src.e[0] & lmask);
  const auto s1(src.e[0] >> hbits);
  const auto s2(src.e[1] & lmask);
  const auto s3(src.e[1] >> hbits);
  // (d0 + d1 * p1 + d2 * p2 + d3 * p3) * (s0 + s1 * p1 + s2 * p2 + s3 * p3) ==
  // ... ==
  // d0 * s0 + (d0 * s1 + s0 * d1) * p1 +
  //   (d0 * s2 + d2 * s0 + d1 * s1) * p2 +
  //   (d0 * s3 + d2 * s1 + d1 * s2 + d3 * s0) * p3
  // N.B. not used:
  //   dk * sl + dl * sk == dk * sk + sl * dl - (dk - dl) * (sk - sl)
  return DUInt<T,bits>(d0 * s0) +
       ((DUInt<T,bits>(s0 * d1) + DUInt<T,bits>(s1 * d0)) << hbits) +
       (DUInt<T,bits>(s0 * d2 + s2 * d0 + s1 * d1) << bits) +
       (DUInt<T,bits>(s0 * d3 + s3 * d0 + s1 * d2 + s2 * d1) << (bits + hbits));
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator *= (const DUInt<T,bits>& src) {
  return *this = *this * src;
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator /  (const DUInt<T,bits>& src) const {
  auto work(*this);
  return work /= src;
}

template <typename T, int bits> inline int DUInt<T,bits>::getMSB() const {
  int shift(- 1);
  if(! e[1]) {
    if(! e[0])
      return shift;
    for(int i = bits - 1; 0 <= i; i --)
      if(int(e[0] >> i) & 1) {
        shift = i;
        break;
      }
    return shift;
  }
  for(int i = bits - 1; 0 <= i; i --)
    if(int(e[1] >> i) & 1) {
      shift = i + bits;
      break;
    }
  return shift;
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator /= (const DUInt<T,bits>& src) {
  const static auto hbits(bits >> 1);
  const static auto lmask((T(1) << hbits) - T(1));
  if(! *this)
    return *this;
  const auto shift(src.getMSB());
  if(shift < 0)
    throw "Zero division";
  const auto dblocks(shift / hbits + 1);
  const auto lshift(dblocks * hbits - shift - 1);
  assert(0 <= lshift && lshift < hbits && !((shift + lshift + 1) % hbits));
  const auto dd(src << lshift);
  const auto dlast((dblocks - 1) & 1 ? dd.e[(dblocks - 1) >> 1] >> hbits
                                     : dd.e[(dblocks - 1) >> 1] &  lmask);
  assert(dlast);
        auto ltshift(getMSB());
  assert(0 <= ltshift);
  ltshift = bits * 2 - 1 - ltshift;
  *this <<= ltshift;
  assert(*this);
  // N.B.
  //   block division with better condition.
  //   de[0] ... de[n], de[n] = 0...1..., each de uses half of space.
  //                                ^ hbits - 1
  //   dlast := de[n].
  //   res = *this / (src == dd >> lshift ~= (dlast << ...))
  auto res(src ^ src);
  auto div(res);
  auto d(e[0] ^ e[0]);
  for(int i = 2; - 1 <= i; i --) {
    switch(i) {
    case - 1:
      d =  (e[0] << hbits) / dlast;
      break;
    case 0:
      d =   e[0] / dlast;
      break;
    case 1:
      d = ((e[0] >> hbits) | (e[1] << hbits)) / dlast;
      break;
    case 2:
      d =   e[1] / dlast;
      break;
    default:
      assert(0 && "Should not be reached.");
    }
    // N.B.
    //   d(0)  := original.
    //   d(k+1) = (d << tshift) * src + r + d(k) == orig.
    //   ( |r| < |src|, ((d << tshift) * src) <= orig - d(k) )
    const auto tshift(hbits * i + lshift - (dblocks - 1) * hbits);
    div = (DUInt<T,bits>(d) * src) << tshift;
    for(int j = 0; j < 4 && ! (div <= *this); j ++) {
      -- d;
      div -= src << tshift;
    }
    assert(div <= *this);
    *this -= div;
    res   += DUInt<T,bits>(d) << tshift;
  }
  return *this = (res >>= ltshift);
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator %  (const DUInt<T,bits>& src) const {
  return *this - ((*this / src) * src);
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator %= (const DUInt<T,bits>& src) {
  return *this = *this % src;
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator << (const int& b)             const {
  auto work(*this);
  return work <<= b;
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator <<= (const int& b) {
  if(! b)
    return *this;
  else if(b < 0)
    return *this >>= (- b);
  else if(b > bits * 2)
    return *this ^= *this;
  else if(b > bits) {
    e[1]  = e[0] << (b - bits);
    e[0] ^= e[0];
  } else if(b == bits) {
    e[1]  = e[0];
    e[0] ^= e[0];
  } else {
    e[1] <<= b;
    e[1]  |= e[0] >> (bits - b);
    e[0] <<= b;
  }
  return *this;
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator >> (const int& b)             const {
  auto work(*this);
  return work >>= b;
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator >>= (const int& b) {
  if(! b)
    return *this;
  else if(b < 0)
    return *this <<= (- b);
  else if(b > bits * 2)
    return *this ^= *this;
  else if(b > bits) {
    e[0]  = e[1] >> (b - bits);
    e[1] ^= e[1];
  } else if(b == bits) {
    e[0]  = e[1];
    e[1] ^= e[1];
  } else {
    e[0] >>= b;
    e[0]  |= e[1] << (bits - b);
    e[1] >>= b;
  }
  return *this;
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator &  (const DUInt<T,bits>& src) const {
  auto work(*this);
  return work &= src;
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator &= (const DUInt<T,bits>& src) {
  e[0] &= src.e[0];
  e[1] &= src.e[1];
  return *this;
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator |  (const DUInt<T,bits>& src) const {
  auto work(*this);
  return work |= src;
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator |= (const DUInt<T,bits>& src) {
  e[0] |= src.e[0];
  e[1] |= src.e[1];
  return *this;
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator ^  (const DUInt<T,bits>& src) const {
  auto work(*this);
  return work ^= src;
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator ^= (const DUInt<T,bits>& src) {
  e[0] ^= src.e[0];
  e[1] ^= src.e[1];
  return *this;
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator ~  () const {
  DUInt<T,bits> work;
  work.e[0] = ~ e[0];
  work.e[1] = ~ e[1];
  return work;
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator =  (const DUInt<T,bits>& src) {
  e[0] = src.e[0];
  e[1] = src.e[1];
  return *this;
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator =  (const DUInt<DUInt<T,bits>,bits*2>& src) {
  return *this = src.e[0];
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator =  (const int& src) {
  e[0]  = src;
  e[1] ^= e[1];
  return *this;
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator =  (DUInt<T,bits>&& src) {
  e[0] = move(src.e[0]);
  e[1] = move(src.e[1]);
  return *this;
}

template <typename T, int bits> inline bool      DUInt<T,bits>::operator <  (const DUInt<T,bits>& src) const {
  if(e[1])
    return e[1] == src.e[1] ? e[0] < src.e[0] : e[1] < src.e[1];
  return src.e[1] ? true : e[0] < src.e[0];
}

template <typename T, int bits> inline bool      DUInt<T,bits>::operator <= (const DUInt<T,bits>& src) const {
  return *this < src || *this == src;
}

template <typename T, int bits> inline bool      DUInt<T,bits>::operator >  (const DUInt<T,bits>& src) const {
  return ! (*this <= src);
}

template <typename T, int bits> inline bool      DUInt<T,bits>::operator >= (const DUInt<T,bits>& src) const {
  return ! (*this < src);
}

template <typename T, int bits> inline bool      DUInt<T,bits>::operator == (const DUInt<T,bits>& src) const {
  return ! (*this - src);
}

template <typename T, int bits> inline bool      DUInt<T,bits>::operator != (const DUInt<T,bits>& src) const {
  return ! (*this == src);
}

template <typename T, int bits> inline bool      DUInt<T,bits>::operator && (const DUInt<T,bits>& src) const {
  return ! ( (! *this) || (! src) );
}

template <typename T, int bits> inline bool      DUInt<T,bits>::operator || (const DUInt<T,bits>& src) const {
  return ! ( (! *this) && (! src) );
}

template <typename T, int bits> inline bool      DUInt<T,bits>::operator !    () const {
  return (! e[0]) && (! e[1]);
}

template <typename T, int bits> inline           DUInt<T,bits>::operator bool () const {
  return ! (! *this);
}

template <typename T, int bits> inline           DUInt<T,bits>::operator int () const {
  return int(e[0]);
}

template <typename T, int bits> inline           DUInt<T,bits>::operator T   () const {
  return e[0];
}

template <typename T, int bits> inline           DUInt<T,bits>::operator DUInt<T,bits> () const {
  return *this;
}

template <typename T, int bits> std::ostream&  operator << (std::ostream& os, DUInt<T,bits> v) {
  const static DUInt<T,bits> ten(10);
  vector<char> buf;
  while(v) {
    const auto div(v / ten);
    buf.push_back(int(v - div * ten));
    v = div;
  }
  if(buf.size()) {
    for(int i = 0; 0 <= i && i < buf.size(); i ++)
      os << int(buf[buf.size() - 1 - i]);
    return os;
  }
  return os << '0';
}

template <typename T, int bits> std::istream&  operator >> (std::istream& is, DUInt<T,bits>& v) {
  const static DUInt<T,bits> ten(10);
  v = DUInt<T,bits>(0);
  // skip white spaces.
  while(! is.eof()) {
    const auto buf(is.get());
    if(buf != ' ' && buf != '\t') {
      is.unget();
      break;
    }
  }
  while(! is.eof() ) {
    const auto buf(is.get());
    if('0' <= buf && buf <= '9') {
      v *= ten;
      v += DUInt<T,bits>(int(buf - '0'));
    } else
      goto ensure;
  }
 ensure:
  return is;
}


// integer to integer float part.
template <typename T, typename W, int bits, typename U> class SimpleFloat {
public:
  inline SimpleFloat();
  inline SimpleFloat(const int& src);
  inline SimpleFloat(const SimpleFloat<T,W,bits,U>& src);
  inline SimpleFloat(SimpleFloat<T,W,bits,U>&& src);
  inline ~SimpleFloat();
  
  inline SimpleFloat<T,W,bits,U>  operator -  () const;
  inline SimpleFloat<T,W,bits,U>  operator +  (const SimpleFloat<T,W,bits,U>& src) const;
         SimpleFloat<T,W,bits,U>& operator += (const SimpleFloat<T,W,bits,U>& src);
  inline SimpleFloat<T,W,bits,U>  operator -  (const SimpleFloat<T,W,bits,U>& src) const;
  inline SimpleFloat<T,W,bits,U>& operator -= (const SimpleFloat<T,W,bits,U>& src);
  inline SimpleFloat<T,W,bits,U>  operator *  (const SimpleFloat<T,W,bits,U>& src) const;
         SimpleFloat<T,W,bits,U>& operator *= (const SimpleFloat<T,W,bits,U>& src);
  inline SimpleFloat<T,W,bits,U>  operator /  (const SimpleFloat<T,W,bits,U>& src) const;
         SimpleFloat<T,W,bits,U>& operator /= (const SimpleFloat<T,W,bits,U>& src);
  inline SimpleFloat<T,W,bits,U>  operator %  (const SimpleFloat<T,W,bits,U>& src) const;
  inline SimpleFloat<T,W,bits,U>& operator %= (const SimpleFloat<T,W,bits,U>& src);
  inline SimpleFloat<T,W,bits,U>  operator <<  (const U& b) const;
  inline SimpleFloat<T,W,bits,U>& operator <<= (const U& b);
  inline SimpleFloat<T,W,bits,U>  operator >>  (const U& b) const;
  inline SimpleFloat<T,W,bits,U>& operator >>= (const U& b);
  inline SimpleFloat<T,W,bits,U>& operator =  (const SimpleFloat<T,W,bits,U>& src);
  inline SimpleFloat<T,W,bits,U>& operator =  (SimpleFloat<T,W,bits,U>&& src);
  inline bool             operator == (const SimpleFloat<T,W,bits,U>& src) const;
  inline bool             operator != (const SimpleFloat<T,W,bits,U>& src) const;
  inline bool             operator <  (const SimpleFloat<T,W,bits,U>& src) const;
  inline bool             operator <= (const SimpleFloat<T,W,bits,U>& src) const;
  inline bool             operator >  (const SimpleFloat<T,W,bits,U>& src) const;
  inline bool             operator >= (const SimpleFloat<T,W,bits,U>& src) const;
  inline bool             operator !  () const;
  inline                  operator bool () const;
  inline                  operator int  () const;
  inline                  operator T    () const;
  inline                  operator SimpleFloat<T,W,bits,U> () const;
  inline SimpleFloat<T,W,bits,U>  floor() const;
  inline SimpleFloat<T,W,bits,U>  ceil() const;
  inline SimpleFloat<T,W,bits,U>  abs()  const;
         SimpleFloat<T,W,bits,U>  log()  const;
         SimpleFloat<T,W,bits,U>  exp()  const;
         SimpleFloat<T,W,bits,U>  sin()  const;
         SimpleFloat<T,W,bits,U>  cos()  const;
         SimpleFloat<T,W,bits,U>  atan() const;
  inline SimpleFloat<T,W,bits,U>  sqrt() const;
  
/*
  friend std::ostream&    operator << (std::ostream& os, const SimpleFloat<T,W,bits,U>& v);
  friend std::istream&    operator >> (std::istream& is, SimpleFloat<T,W,bits,U>& v);
*/
  
  unsigned char s;
  typedef enum {
    INF = 0,
    NaN = 1,
    SIGN = 2,
    DWRK = 3
  } state_t;
  T m;
  U e;
  const SimpleFloat<T,W,bits,U>& zero()   const;
  const SimpleFloat<T,W,bits,U>& one()    const;
  const SimpleFloat<T,W,bits,U>& two()    const;
  const SimpleFloat<T,W,bits,U>& pi()     const;
  const SimpleFloat<T,W,bits,U>& halfpi() const;
  const SimpleFloat<T,W,bits,U>& quatpi() const;
  const SimpleFloat<T,W,bits,U>& twopi()  const;
  const SimpleFloat<T,W,bits,U>& sqrt2()  const;
private:
  template <typename V> inline int normalize(V& src) const;
  inline SimpleFloat<T,W,bits,U>& ensureFlag();
  inline unsigned char safeAdd(U& dst, const U& src);
  inline char residue2() const;

  const vector<SimpleFloat<T,W,bits,U> >& exparray()    const;
  const vector<SimpleFloat<T,W,bits,U> >& invexparray() const;
};

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>::SimpleFloat() {
  assert(0 < bits && ! (bits & 1));
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>::SimpleFloat(const int& src) {
  s ^= s;
  m  = std::abs(src);
  e ^= e;
  s |= safeAdd(e, normalize(m));
  if(src < 0)
    s |= 1 << SIGN;
  ensureFlag();
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>::SimpleFloat(const SimpleFloat<T,W,bits,U>& src) {
  *this = src;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>::SimpleFloat(SimpleFloat<T,W,bits,U>&& src) {
  *this = src;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>::~SimpleFloat() {
  ;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::operator -  () const {
  auto work(*this);
  work.s ^= 1 << SIGN;
  return work;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::operator +  (const SimpleFloat<T,W,bits,U>& src) const {
  auto work(*this);
  return work += src;
}

template <typename T, typename W, int bits, typename U>        SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::operator += (const SimpleFloat<T,W,bits,U>& src) {
  if((s |= src.s & (1 << NaN)) & (1 << NaN))
    return *this;
  if(s & (1 << INF)) {
    if(src.s & (1 << INF))
      s |= 1 << NaN;
    return *this;
  }
  if(src.s & (1 << INF))
    return *this = src;
  if(! m)
    return *this = src;
  if(! src.m)
    return *this;
  if(! ((s & (1 << SIGN)) ^ (src.s & (1 << SIGN)))) {
    if(e >= src.e) {
      m >>= 1;
      s |= safeAdd(e, 1);
      U se(e);
      if(! safeAdd(se, - src.e) && se < bits)
        m += src.m >> se;
    } else
      return *this = src + *this;
  } else {
    if(e > src.e) {
      U se(e);
      if(! safeAdd(se, - src.e) && se < bits)
        m -= src.m >> se;
    } else if(e == src.e) {
      if(m >= src.m)
        m -= src.m;
      else
        return *this = src + *this;
    } else
      return *this = src + *this;
  }
  s |= safeAdd(e, normalize(m));
  return ensureFlag();
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::operator -  (const SimpleFloat<T,W,bits,U>& src) const {
  auto work(*this);
  return work -= src;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::operator -= (const SimpleFloat<T,W,bits,U>& src) {
  return *this += - src;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::operator *  (const SimpleFloat<T,W,bits,U>& src) const {
  auto work(*this);
  return work *= src;
}

template <typename T, typename W, int bits, typename U>        SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::operator *= (const SimpleFloat<T,W,bits,U>& src) {
  if((s |= src.s & (1 << NaN)) & (1 << NaN))
    return *this;
  s ^= src.s & (1 << SIGN);
  if((s |= src.s & (1 << INF)) & (1 << INF))
    return *this;
  if((! m) || (! src.m)) {
    m ^= m;
    e ^= e;
    return *this;
  }
  auto mm(W(m) * W(src.m));
  s |= safeAdd(e, src.e);
  s |= safeAdd(e, normalize(mm));
  s |= safeAdd(e, bits);
  m  = T(mm >> U(bits));
  return ensureFlag();
}

template <typename T, typename W, int bits, typename U> inline char SimpleFloat<T,W,bits,U>::residue2() const {
  if(0 < e || bits <= - e)
    return 0;
  if(! e)
    return char(int(m) & 1);
  return char(int(m >> - e) & 1);
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::operator /  (const SimpleFloat<T,W,bits,U>& src) const {
  auto work(*this);
  return work /= src;
}

template <typename T, typename W, int bits, typename U>        SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::operator /= (const SimpleFloat<T,W,bits,U>& src) {
  if((s |= src.s & (1 << NaN)) & (1 << NaN))
    return *this;
  s ^= src.s & (1 << SIGN);
  if(! (s & (1 << INF)) && (src.s & (1 << INF))) {
    s |= 1 << DWRK;
    return ensureFlag();
  }
  if(s & (1 << INF))
    return *this;
  if(! src.m) {
    throw "Zero division";
    s |= 1 << NaN;
    return *this;
  }
  if(! m)
    return *this;
  auto mm((W(m) << bits) / W(src.m));
  s |= safeAdd(e, - src.e);
  s |= safeAdd(e, normalize(mm));
  m  = T(mm >> U(bits));
  return ensureFlag();
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>  SimpleFloat<T,W,bits,U>::operator %  (const SimpleFloat<T,W,bits,U>& src) const {
  return *this - (*this / src).floor() * src;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::operator %= (const SimpleFloat<T,W,bits,U>& src) {
  return *this = *this % src;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>  SimpleFloat<T,W,bits,U>::operator <<  (const U& b) const {
  auto work(*this);
  return work <<= b;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::operator <<= (const U& b) {
  s |= safeAdd(e, b);
  return ensureFlag();
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>  SimpleFloat<T,W,bits,U>::operator >>  (const U& b) const {
  auto work(*this);
  return work >>= b;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::operator >>= (const U& b) {
  s |= safeAdd(e, - b);
  return ensureFlag();
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::operator =  (const SimpleFloat<T,W,bits,U>& src) {
  s = src.s;
  e = src.e;
  m = src.m;
  return *this;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::operator =  (SimpleFloat<T,W,bits,U>&& src) {
  s = move(src.s);
  e = move(src.e);
  m = move(src.m);
  return *this;
}

template <typename T, typename W, int bits, typename U> inline bool             SimpleFloat<T,W,bits,U>::operator == (const SimpleFloat<T,W,bits,U>& src) const {
  return ! (*this != src);
}

template <typename T, typename W, int bits, typename U> inline bool             SimpleFloat<T,W,bits,U>::operator != (const SimpleFloat<T,W,bits,U>& src) const {
  return ((s | src.s) & ((1 << INF) | (1 << NaN))) ||
           (s != src.s || e != src.e || m != src.m);
}

template <typename T, typename W, int bits, typename U> inline bool             SimpleFloat<T,W,bits,U>::operator <  (const SimpleFloat<T,W,bits,U>& src) const {
  if((s | src.s) & (1 << NaN))
    throw "compair NaN";
  const auto s_is_minus(s & (1 << SIGN));
  if(s_is_minus ^ (src.s & (1 << SIGN)))
    return s_is_minus;
  if(s & (1 << INF)) {
    if(src.s & (1 << INF))
      throw "compair INF";
    return s_is_minus;
  }
  if(src.s & (1 << INF))
    return ! s_is_minus;
  if(m && src.m) {
    if(e < src.e)
      return ! s_is_minus;
    if(e == src.e)
      return s_is_minus ? src.m < m : m < src.m;
    return s_is_minus;
  }
  return !m ? (!src.m ? false : ! s_is_minus) : s_is_minus;
}

template <typename T, typename W, int bits, typename U> inline bool             SimpleFloat<T,W,bits,U>::operator <= (const SimpleFloat<T,W,bits,U>& src) const {
  return *this < src || *this == src;
}

template <typename T, typename W, int bits, typename U> inline bool             SimpleFloat<T,W,bits,U>::operator >  (const SimpleFloat<T,W,bits,U>& src) const {
  return ! (*this <= src);
}

template <typename T, typename W, int bits, typename U> inline bool             SimpleFloat<T,W,bits,U>::operator >= (const SimpleFloat<T,W,bits,U>& src) const {
  return ! (*this < src);
}

template <typename T, typename W, int bits, typename U> inline bool             SimpleFloat<T,W,bits,U>::operator !  () const {
  return ! m && isfinite(*this);
}

template <typename T, typename W, int bits, typename U> inline                  SimpleFloat<T,W,bits,U>::operator bool () const {
  return ! (!*this);
}

template <typename T, typename W, int bits, typename U> inline                  SimpleFloat<T,W,bits,U>::operator int  () const {
  return int(T(*this));
}

template <typename T, typename W, int bits, typename U> inline                  SimpleFloat<T,W,bits,U>::operator T    () const {
  auto deci(*this);
  if(deci.s & (1 << INF))
    throw "Inf to convert int";
  if(deci.s & (1 << NaN))
    throw "NaN to convert int";
  if(! deci.m)
    return T(0);
  if(bits <= deci.e || (U(0) < deci.e && ! (deci.m <<= deci.e)))
    throw "Overflow to convert int.";
  if(deci.e <= - bits)
    return T(0);
  if(deci.e <  U(0))
    deci.m >>= - deci.e;
  return deci.m;
}

template <typename T, typename W, int bits, typename U> inline                  SimpleFloat<T,W,bits,U>::operator SimpleFloat<T,W,bits,U> () const {
  return *this;
}

template <typename T, typename W, int bits, typename U> template <typename V> inline int SimpleFloat<T,W,bits,U>::normalize(V& src) const {
  V   bt(1);
  int b(0);
  int tb;
  for(tb = 0; bt; tb ++) {
    if(src & bt)
      b = tb;
    bt <<= 1;
  }
  const auto shift(tb - b - 1);
  assert(0 <= shift);
  src <<= shift;
  return - shift;
}

template <typename T, typename W, int bits, typename U> SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::ensureFlag() {
  if(! m || (s & (1 << DWRK))) {
    e ^= e;
    m ^= m;
    s &= ~ (1 << DWRK);
  }
  return * this;
}

template <typename T, typename W, int bits, typename U> inline unsigned char SimpleFloat<T,W,bits,U>::safeAdd(U& dst, const U& src) {
  const auto dst0(dst);
  unsigned char ss(0);
  dst += src;
  if((dst0 > 0 && src > 0 && dst < 0) ||
     (dst0 < 0 && src < 0 && dst > 0)) {
    if(dst < 0)
      ss |= 1 << INF;
    else
      ss |= 1 << DWRK;
  }
  return ss;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::floor() const {
  if(0 <= e)
    return *this;
  if(e <= - bits)
    return zero();
  auto deci(*this);
  deci.m >>= - deci.e;
  if(! deci.m)
    return zero();
  deci.m <<= - deci.e;
  return deci;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::ceil() const {
  if(*this - this->floor()) {
    auto pmone(one());
    pmone.s |= s & (1 << SIGN);
    return this->floor() + pmone;
  }
  return this->floor();
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::abs() const {
  auto work(*this);
  work.s &= ~ (1 << SIGN);
  return work;
}

template <typename T, typename W, int bits, typename U> SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::log() const {
  const static auto einv(one() / one().exp());
  const static auto one_einv(one() + einv);
  if((s & (1 << SIGN)) && this->m)
    throw "Negative log";
  if(s & ((1 << INF) | (1 << NaN)))
    return *this;
  if(! this->m) {
    auto work(*this);
    work.s |= (1 << INF) | (1 << SIGN);
    return work;
  }
  if(einv <= *this && *this <= one_einv) {
    // ln(x) = (x - 1) - (x - 1)^2/2 + (x-1)^3/3- ...
    const auto dx(*this - one());
          auto x(dx);
          auto before(one());
          auto res(zero());
    for(int t = 1; (res - before).m; t ++, x *= dx) {
      const auto abst(x / SimpleFloat<T,W,bits,U>(t));
      before = res;
      res   += (t % 2 ? abst : - abst);
    }
    return res;
  }
  const auto& ea(exparray());
  const auto& iea(invexparray());
        auto  result(zero());
        auto  work(*this);
  if(one_einv < work) {
    for(int i = min(ea.size(), iea.size()) - 1; 0 < i; i --)
      if(ea[i] <= work) {
        result += one() << U(i - 1);
        work   *= iea[i];
      }
    if(! (work <= one_einv)) {
      result += one();
      work   *= iea[1];
    }
  } else {
    for(int i = min(ea.size(), iea.size()) - 1; 0 < i; i --)
      if(work <= iea[i]) {
        result -= one() << U(i - 1);
        work   *= ea[i];
      }
    if(! (einv <= work)) {
      result -= one();
      work   *= ea[1];
    }
  }
  assert(einv <= work && work <= one_einv);
  return result += work.log();
}

template <typename T, typename W, int bits, typename U> SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::exp() const {
  if(s & ((1 << INF) | (1 << NaN))) {
    if(! (s & (1 << NaN)) && (s & (1 << SIGN)))
      return one();
    return *this;
  }
  if(this->abs() <= one()) {
    // exp(x) = 1 + x/1! + x^2/2! + ...
    auto denom(one());
    auto x(*this);
    auto before(zero());
    auto res(one());
    for(int t = 1; (res - before).m; t ++, x *= *this) {
      before = res;
      denom *= SimpleFloat<T,W,bits,U>(t);
      res   += x / denom;
    }
    return res;
  }
  const auto& en(exparray());
  const auto& ien(invexparray());
        auto  work(this->abs());
        auto  result(one());
  for(int i = 1; 0 <= i && i < min(en.size(), ien.size()) && work.floor(); i ++, work >>= U(1))
    if(work.residue2()) {
      if(s & (1 << SIGN))
        result *= ien[i];
      else
        result *= en[i];
    }
  if(work.floor()) {
    work.s |= 1 << INF;
    return work;
  }
  const auto residue(*this - this->floor());
  assert(residue.abs() <= one());
  return result *= residue.exp();
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::sin() const {
  if(s & ((1 << INF) | (1 << NaN)))
    return *this;
  if(- quatpi() <= *this && *this <= quatpi()) {
    // sin(x) = x - x^3/3! + x^5/5! - ...
    const auto sqx(*this * *this);
          auto denom(one());
          auto x(sqx * *this);
          auto before(zero());
          auto res(*this);
    for(int t = 1; (res - before).m; t ++, x *= sqx) {
      SimpleFloat<T,W,bits,U> tt(t);
      tt   <<= U(1);
      before = res;
      denom *= - tt * (tt + one());
      res   += x / denom;
    }
    return res;
  }
  if(- halfpi() <= *this && *this <= halfpi())
    return ((*this - quatpi()).cos() + (*this - quatpi()).sin()) / sqrt2();
  if(- pi() <= *this && *this <= pi())
    return (halfpi() - *this).cos();
  if(- twopi() <= *this && *this <= twopi())
    return - (*this + pi()).sin();
  return (*this % twopi()).sin();
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::cos() const {
  if(s & ((1 << INF) | (1 << NaN)))
    return *this;
  if(- quatpi() <= *this && *this <= quatpi()) {
    // cos(x) = 1 - x^2/2! + x^4/4! - ...
    const auto sqx(*this * *this);
          auto denom(one());
          auto x(sqx);
          auto before(zero());
          auto res(one());
    for(int t = 1; (res - before).m; t ++, x *= sqx) {
      SimpleFloat<T,W,bits,U> tt(t);
      tt   <<= U(1);
      before = res;
      denom *= - tt * (tt - one());
      res   += x / denom;
    }
    return res;
  }
  if(- halfpi() <= *this && *this <= halfpi())
    return ((*this - quatpi()).cos() - (*this - quatpi()).sin()) / sqrt2();
  if(- pi() <= *this && *this <= pi())
    return (halfpi() - *this).sin();
  if(- twopi() <= *this && *this <= twopi())
    return - (*this + pi()).cos();
  return (*this % twopi()).cos();
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::atan() const {
  if(s & ((1 << INF) | (1 << NaN)))
    return *this;
  static const auto half(one() >> U(1));
  static const auto four(one() << U(2));
  static const auto five((one() << U(2)) + one());
  if(- half <= *this && *this <= half) {
    // arctan(x) = x - x^3/3 + x^5/5 - ...
    const auto sqx(*this * *this);
          auto x(sqx * *this);
          auto before(zero());
          auto res(*this);
    for(int t = 1; (res - before).m; t ++, x *= sqx) {
      const auto abst(x / ((SimpleFloat<T,W,bits,U>(t) << U(1)) + one()));
      before = res;
      res   += (t % 2 ? - abst : abst);
    }
    return res;
  }
  // N.B.
  //  atan(u) + atan(v) = atan((u + v) / (1 - uv)) mod pi, uv != 1.
  //    in u = 0.5, v = x - 0.5 case,
  //  atan(x / (1 - x / 2 + 1 / 4)) = atan(.5) + atan(x - .5) =
  //  atan(x / (1.25 - .5 * x)) 
  //  y := x / (1.25 - .5 * x) then,
  //  (1.25 - .5 * x) * y = x,
  //  (5 - 2x) * y = 4 x
  //  x = 5y / (4 + 2y),
  //     y - x = ((4 + 2y) * y - 5y) / (4 + 2y)
  //           = y * (2y - 1) / (4 + 2y)
  //     so 0 <= y and 0 < y case, this makes decreasing function.
  //       (v = x - .5 and 0 <= 2y - 1)
  if(- two() <= *this && *this <= two()) {
    static const auto atanhalf(half.atan());
    if(s & (1 << SIGN))
      return - (- *this).atan();
    const auto v(five * *this / (four + (*this << U(1))) - half);
    assert(v < *this);
    return atanhalf + v.atan();
  }
  // N.B.
  //    in u = v case,
  //  2 atan(u) = atan(2 * u / (1 - u * u))
  //    in u := x + 1 case,
  //  2 atan(1 + x) = atan(2 * (1 + x) / (x + x * x))
  //                = atan(2 / x)
  //    in Y := 2 / x case,
  //  atan(Y) = 2 atan(1 + 2 / Y)
  const auto y(one() + two() / (*this));
  assert(- two() <= y && y <= two());
  return y.atan() << U(1);
}

template <typename T, typename W, int bits, typename U> const vector<SimpleFloat<T,W,bits,U> >& SimpleFloat<T,W,bits,U>::exparray() const {
  static vector<SimpleFloat<T,W,bits,U> > ebuf;
  if(ebuf.size())
    return ebuf;
  ebuf.push_back(one());
  ebuf.push_back(ebuf[0].exp());
  for(int i = 1; 0 <= i; i ++) {
    const auto en(ebuf[i] * ebuf[i]);
    if(en && isfinite(en))
      ebuf.push_back(en);
    else
      break;
  }
  return ebuf;
}

template <typename T, typename W, int bits, typename U> const vector<SimpleFloat<T,W,bits,U> >& SimpleFloat<T,W,bits,U>::invexparray() const {
  static vector<SimpleFloat<T,W,bits,U> > iebuf;
  if(iebuf.size())
    return iebuf;
  const auto& ea(exparray());
  for(int i = 0; 0 <= i && i < ea.size(); i ++) {
    const auto ien(one() / ea[i]);
    if(ien && isfinite(ien))
      iebuf.push_back(ien);
    else
      break;
  }
  return iebuf;
}

template <typename T, typename W, int bits, typename U> const SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::zero() const {
  const static SimpleFloat<T,W,bits,U> vzero(0);
  return vzero;
}

template <typename T, typename W, int bits, typename U> const SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::one() const {
  const static SimpleFloat<T,W,bits,U> vone(1);
  return vone;
}

template <typename T, typename W, int bits, typename U> const SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::two() const {
  const static auto vtwo(one() << U(1));
  return vtwo;
}

template <typename T, typename W, int bits, typename U> const SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::pi() const {
  const static auto vpi(quatpi() << U(2));
  return vpi;
}

template <typename T, typename W, int bits, typename U> const SimpleFloat<T,W,bits,U>
& SimpleFloat<T,W,bits,U>::halfpi() const {
  const static auto vhalfpi(quatpi() << U(1));
  return vhalfpi;
}

template <typename T, typename W, int bits, typename U> const SimpleFloat<T,W,bits,U>
& SimpleFloat<T,W,bits,U>::quatpi() const {
  const static auto vquatpi(one().atan());
  return vquatpi;
}

template <typename T, typename W, int bits, typename U> const SimpleFloat<T,W,bits,U>
& SimpleFloat<T,W,bits,U>::twopi() const {
  const static auto vtwopi(quatpi() << U(3));
  return vtwopi;
}

template <typename T, typename W, int bits, typename U> const SimpleFloat<T,W,bits,U>
& SimpleFloat<T,W,bits,U>::sqrt2() const {
  const static auto vsqrt2((one() << U(1)).sqrt());
  return vsqrt2;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::sqrt() const {
  if(s & ((1 << INF) | (1 << NaN)))
    return *this;
  return (this->log() >> U(1)).exp();
}

template <typename T, typename W, int bits, typename U> std::ostream& operator << (std::ostream& os, const SimpleFloat<T,W,bits,U>& v) {
  if(isnan(v))
    return os << "NaN";
  if(isinf(v))
    return os << (const char*)(v.s & (1 << v.SIGN) ? "-" : "") << "Inf";
  return os << (const char*)(v.s & (1 << v.SIGN) ? "-" : "") << T(v.m) << "*2^" << int(v.e);
}

template <typename T, typename W, int bits, typename U> std::istream& operator >> (std::istream& is, SimpleFloat<T,W,bits,U>& v) {
  const static SimpleFloat<T,W,bits,U> ten(10);
               SimpleFloat<T,W,bits,U> e(0);
  bool mode(false);
  bool sign(false);
  bool fsign(false);
  v = SimpleFloat<T,W,bits,U>(0);
  // skip white spaces.
  while(! is.eof()) {
    const auto buf(is.get());
    if(buf != ' ' && buf != '\t') {
      is.unget();
      break;
    }
  }
  while(! is.eof() ) {
    const auto buf(is.get());
    switch(buf) {
    case '-':
      sign  = true;
    case '+':
      if(fsign)
        throw "Wrong input";
      fsign = true;
      break;
    case 'e':
      if(mode)
        goto ensure;
      if(sign)
        v   = - v;
      mode  = true;
      sign  = false;
      fsign = false;
      break;
    case '.':
      throw "not implemented now";
      break;
    case '0': case '1': case '2': case '3': case '4':
    case '5': case '6': case '7': case '8': case '9':
      if(mode) {
        e *= ten;
        e += SimpleFloat<T,W,bits,U>(int(buf - '0'));
      } else {
        v *= ten;
        v += SimpleFloat<T,W,bits,U>(int(buf - '0'));
      }
      fsign = true;
      break;
    default:
      goto ensure;
    }
  }
 ensure:
  if(sign) {
    if(mode)
      e = - e;
    else
      v = - v;
  }
  v *= pow(ten, e);
  return is;
}

template <typename T, typename W, int bits, typename U> inline bool isinf(const SimpleFloat<T,W,bits,U>& src) {
  return src.s & (1 << src.INF);
}

template <typename T, typename W, int bits, typename U> inline bool isnan(const SimpleFloat<T,W,bits,U>& src) {
  return src.s & (1 << src.NaN);
}

template <typename T, typename W, int bits, typename U> inline bool isfinite(const SimpleFloat<T,W,bits,U>& src) {
  return ! (src.s & ((1 << src.INF) | (1 << src.NaN)));
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> floor(const SimpleFloat<T,W,bits,U>& src) {
  return src.floor();
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> ceil(const SimpleFloat<T,W,bits,U>& src) {
  return src.ceil();
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> abs(const SimpleFloat<T,W,bits,U>& src) {
  return src.abs();
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> sqrt(const SimpleFloat<T,W,bits,U>& src) {
  return src.sqrt();
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> exp(const SimpleFloat<T,W,bits,U>& src) {
  return src.exp();
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> log(const SimpleFloat<T,W,bits,U>& src) {
  return src.log();
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> sin(const SimpleFloat<T,W,bits,U>& src) {
  return src.sin();
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> cos(const SimpleFloat<T,W,bits,U>& src) {
  return src.cos();
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> tan(const SimpleFloat<T,W,bits,U>& src) {
  return src.sin() / src.cos();
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> atan2(const SimpleFloat<T,W,bits,U>& y, const SimpleFloat<T,W,bits,U>& x) {
  auto atan0(y.halfpi());
  if(! x && ! y)
    return x / y;
  else if(isfinite(x)) {
    if(! isfinite(y) )
      goto ensure;
    const auto yoverx((y / x).abs());
    if(! isfinite(yoverx) )
      goto ensure;
    const auto atan00(yoverx.atan());
    if(! isfinite(atan00) )
      goto ensure;
    atan0 = atan00;
    goto ensure;
  } else if(isfinite(y)) {
    atan0 = x.zero();
    goto ensure;
  }
  return y;
 ensure:
  if(y.s & (1 << y.SIGN)) {
    if(x.s & (1 << x.SIGN))
      atan0 = - (x.pi() - atan0);
    else
      atan0 = - atan0;
  } else if(x.s & (1 << x.SIGN))
    atan0  = x.pi() - atan0;
  return atan0;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> pow(const SimpleFloat<T,W,bits,U>& src, const SimpleFloat<T,W,bits,U>& dst) {
  if(! dst) {
    if(! src)
      throw "0^0";
    return dst.one();
  }
  return exp(log(src) * dst);
}


// class complex part:
template <typename T> class Complex {
public:
  inline Complex();
  inline Complex(const Complex<T>& s);
  inline Complex(Complex<T>&& s);
  inline Complex(const T& real, const T& imag = T(0));
  inline Complex(T&& real);
  inline Complex(T&& real, T&& imag);
  inline ~Complex();

  inline Complex<T>  operator ~  ()                    const;
  inline Complex<T>  operator -  ()                    const;
  inline Complex<T>  operator +  (const Complex<T>& s) const;
  inline Complex<T>& operator += (const Complex<T>& s);
  inline Complex<T>  operator -  (const Complex<T>& s) const;
  inline Complex<T>& operator -= (const Complex<T>& s);
  inline Complex<T>  operator *  (const T& s)          const;
  inline Complex<T>& operator *= (const T& s);
  inline Complex<T>  operator *  (const Complex<T>& s) const;
  inline Complex<T>& operator *= (const Complex<T>& s);
  inline Complex<T>  operator /  (const T& s)          const;
  inline Complex<T>& operator /= (const T& s);
  inline Complex<T>  operator /  (const Complex<T>& s) const;
  inline Complex<T>& operator /= (const Complex<T>& s);
  inline bool        operator == (const Complex<T>& s) const;
  inline bool        operator != (const Complex<T>& s) const;
  inline bool        operator !  ()                    const;
  inline Complex<T>  operator &  (const Complex<T>& s) const;
  inline Complex<T>& operator &= (const Complex<T>& s);
  inline Complex<T>  operator |  (const Complex<T>& s) const;
  inline Complex<T>& operator |= (const Complex<T>& s);
  inline Complex<T>  operator ^  (const Complex<T>& s) const;
  inline Complex<T>& operator ^= (const Complex<T>& s);
  inline bool        operator && (const Complex<T>& s) const;
  inline bool        operator || (const Complex<T>& s) const;
  inline Complex<T>& operator =  (const Complex<T>& s);
  inline Complex<T>& operator =  (Complex<T>&& s);
  inline T&          operator [] (const size_t& i);
  inline             operator bool () const;
  inline             operator T    () const;
  
  const Complex<T>& i() const;
  
  inline T  abs() const;
  inline T  arg() const;
  inline T& real();
  inline T& imag();
  inline const T& real() const;
  inline const T& imag() const;
  T _real;
  T _imag;
};

template <typename T> inline Complex<T>::Complex() {
  ;
}

template <typename T> inline Complex<T>::Complex(const Complex<T>& src) {
  *this = src;
}

template <typename T> inline Complex<T>::Complex(Complex<T>&& src) {
  *this = src;
}

template <typename T> inline Complex<T>::Complex(const T& real, const T& imag) {
  _real = real;
  _imag = imag;
  return;
}

template <typename T> inline Complex<T>::Complex(T&& real) {
  const static T zero(0);
  _real = move(real);
  _imag = zero;
  return;
}

template <typename T> inline Complex<T>::Complex(T&& real, T&& imag) {
  _real = move(real);
  _imag = move(imag);
  return;
}

template <typename T> inline Complex<T>::~Complex() {
  ;
}

template <typename T> inline Complex<T> Complex<T>::operator ~ () const {
  return Complex<T>(  _real, - _imag);
}

template <typename T> inline Complex<T> Complex<T>::operator - () const {
  return Complex<T>(- _real, - _imag);
}

template <typename T> inline Complex<T> Complex<T>::operator + (const Complex<T>& s) const {
  auto result(*this);
  return result += s;
}

template <typename T> inline Complex<T>& Complex<T>::operator += (const Complex<T>& s) {
  _real += s._real;
  _imag += s._imag;
  return *this;
}

template <typename T> inline Complex<T> Complex<T>::operator - (const Complex<T>& s) const {
  auto result(*this);
  return result -= s;
}

template <typename T> inline Complex<T>& Complex<T>::operator -= (const Complex<T>& s) {
  _real -= s._real;
  _imag -= s._imag;
  return *this;
}

template <typename T> inline Complex<T>  Complex<T>::operator * (const T& s) const {
  auto result(*this);
  return result *= s;
}

template <typename T> inline Complex<T>& Complex<T>::operator *= (const T& s) {
  _real *= s;
  _imag *= s;
  return *this;
}
  
template <typename T> inline Complex<T> Complex<T>::operator * (const Complex<T>& s) const {
  return Complex<T>(_real * s._real - _imag * s._imag,
                    _real * s._imag + _imag * s._real);
}
 
template <typename T> inline Complex<T>& Complex<T>::operator *= (const Complex<T>& s) {
  return (*this) = (*this) * s;
}

template <typename T> inline Complex<T> Complex<T>::operator / (const T& s) const {
  auto result(*this);
  return result /= s;
}

template <typename T> inline Complex<T>& Complex<T>::operator /= (const T& s) {
  _real /= s;
  _imag /= s;
  return *this;
}

template <typename T> inline Complex<T> Complex<T>::operator / (const Complex<T>& s) const {
  return (*this * (~ s)) / (s._real * s._real + s._imag * s._imag);
}

template <typename T> inline Complex<T>& Complex<T>::operator /= (const Complex<T>& s) {
  return *this = *this / s;
}

template <typename T> inline bool Complex<T>::operator == (const Complex<T>& s) const {
  return !(*this != s);
}

template <typename T> inline bool Complex<T>::operator != (const Complex<T>& s) const {
  return (_real != s._real) || (_imag != s._imag);
}

template <typename T> inline bool Complex<T>::operator ! () const {
  return !_real && !_imag;
}

template <typename T> inline Complex<T> Complex<T>::operator & (const Complex<T>& s) const {
  auto result(*this);
  return result &= s;
}

template <typename T> inline Complex<T>& Complex<T>::operator &= (const Complex<T>& s) {
  _real &= s._real;
  _imag &= s._imag;
  return *this;
}

template <typename T> inline Complex<T> Complex<T>::operator | (const Complex<T>& s) const {
  auto result(*this);
  return result |= s;
}

template <typename T> inline Complex<T>& Complex<T>::operator |= (const Complex<T>& s) {
  _real |= s._real;
  _imag |= s._imag;
  return *this;
}

template <typename T> inline Complex<T> Complex<T>::operator ^ (const Complex<T>& s) const {
  auto result(*this);
  return result ^= s;
}

template <typename T> inline Complex<T>& Complex<T>::operator ^= (const Complex<T>& s) {
  _real ^= s._real;
  _imag ^= s._imag;
  return *this;
}

template <typename T> inline bool Complex<T>::operator && (const Complex<T>& s) const {
  return *this && s;
}

template <typename T> inline bool Complex<T>::operator || (const Complex<T>& s) const {
  return *this || s;
}

template <typename T> inline Complex<T>& Complex<T>::operator =  (const Complex<T>& s) {
  _real = s._real;
  _imag = s._imag;
  return *this;
}

template <typename T> inline Complex<T>& Complex<T>::operator =  (Complex<T>&& s) {
  _real = move(s._real);
  _imag = move(s._imag);
  return *this;
}

template <typename T> inline T& Complex<T>::operator [] (const size_t& i) {
  assert(0 <= i && i < 2);
  if(i)
    return _imag;
  return _real;
}

template <typename T> inline Complex<T>::operator bool () const {
  return ! (! *this);
}

template <typename T> inline Complex<T>::operator T () const {
  return this->_real;
}

template <typename T> const Complex<T>& Complex<T>::i() const {
  const static auto I(Complex<T>(T(0), T(1)));
  return I;
}

template <typename T> T Complex<T>::abs() const {
  return sqrt(_real * _real + _imag * _imag);
}

template <typename T> T Complex<T>::arg() const {
  return atan2(_imag, _real);
}

template <typename T> inline T& Complex<T>::real() {
  return _real;
}

template <typename T> inline T& Complex<T>::imag() {
  return _imag;
}

template <typename T> inline const T& Complex<T>::real() const {
  return _real;
}

template <typename T> inline const T& Complex<T>::imag() const {
  return _imag;
}

template <typename T> std::ostream& operator << (std::ostream& os, const Complex<T>& v) {
  return os << v.real() << "+i" << v.imag();
}


template <typename T> T abs(const Complex<T>& s) {
  return s.abs();
}

template <typename T> T arg(const Complex<T>& s) {
  return s.arg();
}

template <typename T> const T& real(const Complex<T>& s) {
  return s.real();
}

template <typename T> const T& imag(const Complex<T>& s) {
  return s.imag();
}

template <typename T> Complex<T> exp(const Complex<T>& s) {
  return Complex<T>(exp(s.real())) * Complex<T>(cos(s.imag()), sin(s.imag()));
}

template <typename T> Complex<T> log(const Complex<T>& s) {
  // N.B. main branch
  return Complex<T>(log(abs(s)), arg(s));
}

template <typename T> Complex<T> sqrt(const Complex<T>& s) {
  return exp(log(s) * Complex<T>(T(1) / T(2)));
}

template <typename T> Complex<T> csin(const Complex<T>& s) {
  return (exp(Complex<T>(T(0), s)) - exp(Complex<T>(T(0), - s))) / Complex<T>(T(0), T(2));
}

template <typename T> Complex<T> ccos(const Complex<T>& s) {
  return (exp(Complex<T>(T(0), s)) + exp(Complex<T>(T(0), - s))) / T(2);
}

template <typename T> Complex<T> ctan(const Complex<T>& s) {
  return csin(s) / ccos(s);
}

template <typename T> Complex<T> ccsc(const Complex<T>& s) {
  return Complex<T>(T(1)) / csin(s);
}

template <typename T> Complex<T> csec(const Complex<T>& s) {
  return Complex<T>(T(1)) / ccos(s);
}

template <typename T> T ccot(const T& s) {
  return Complex<T>(T(1)) / ctan(s);
}

#define _INTEGER_FLOAT_
#endif

