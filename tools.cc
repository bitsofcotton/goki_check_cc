#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cctype>
#include <assert.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

//#define int int64_t
#define int int32_t
#define _COMPILE_GOKI_
#include "lieonn.hh"
typedef myfloat num_t;

using std::cout;
using std::cerr;
using std::endl;
using std::atoi;
using std::string;
using std::vector;
using std::sort;

#include <stdlib.h>

static inline num_t myatof(const char* v) {
  return num_t(int(std::atof(v) * 10000)) / num_t(int(10000));
}

#undef int
int main(int argc, const char* argv[]) {
//#define int int64_t
#define int int32_t
  if(argc < 2) goto usage;
  if(strcmp(argv[1], "nop") == 0 ||
     strcmp(argv[1], "limit") == 0 ||
     strcmp(argv[1], "bit") == 0 ||
     strcmp(argv[1], "nbit") == 0 ||
     strcmp(argv[1], "slide") == 0 ||
     strcmp(argv[1], "nslide") == 0 ||
     strcmp(argv[1], "collect") == 0 ||
     strcmp(argv[1], "sharpen") == 0 ||
     strcmp(argv[1], "blur") == 0 ||
     strcmp(argv[1], "bump")    == 0 ||
     strcmp(argv[1], "enlarge") == 0 ||
     strcmp(argv[1], "shrink") == 0 ||
     strcmp(argv[1], "flarge") == 0 ||
     strcmp(argv[1], "blink") == 0 ||
     strcmp(argv[1], "w2b")     == 0 ||
     strcmp(argv[1], "b2w")     == 0 ||
     strcmp(argv[1], "b2wd")    == 0 ||
     strcmp(argv[1], "rgb2xyz") == 0 ||
     strcmp(argv[1], "xyz2rgb") == 0) {
    if(argc < 3) goto usage;
    const int recur(4 < argc ? atoi(argv[4]) : 1);
    const int rot(5 < argc ? atoi(argv[5]) : 0);
    vector<SimpleMatrix<num_t> > data;
    if(!loadp2or3<num_t>(data, argv[2]))
      return - 1;
    if(strcmp(argv[1], "collect") == 0)
      for(int i = 0; i < data.size(); i ++)
        data[i] = filter<num_t>(data[i], COLLECT_BOTH, recur, rot);
    else if(strcmp(argv[1], "enlarge") == 0)
      for(int i = 0; i < data.size(); i ++)
        try {
          data[i] = filter<num_t>(filter<num_t>(data[i], ENLARGE_BOTH, recur, rot), CLIP);
        } catch(const char* e) { cerr << e << endl; }
    else if(strcmp(argv[1], "shrink") == 0)
      for(int i = 0; i < data.size(); i ++)
        try {
          data[i] = filter<num_t>(filter<num_t>(data[i], SHRINK_BOTH, recur, rot), CLIP);
        } catch(const char* e) { cerr << e << endl; }
    else if(strcmp(argv[1], "flarge") == 0)
      for(int i = 0; i < data.size(); i ++)
        data[i] = filter<num_t>(data[i], FLARGE_BOTH, recur, rot);
    else if(strcmp(argv[1], "blink") == 0)
      for(int i = 0; i < data.size(); i ++)
        data[i] = filter<num_t>(data[i], BLINK_BOTH, recur, rot);
    else if(strcmp(argv[1], "sharpen") == 0)
      for(int ii = 0; ii < (recur < 1 ? 1 : recur); ii ++)
        for(int i = 0; i < data.size(); i ++)
          data[i] = filter<num_t>(data[i], SHARPEN_BOTH, 1, rot);
    else if(strcmp(argv[1], "blur") == 0)
      for(int i = 0; i < data.size(); i ++)
        data[i] = filter<num_t>(data[i], BLUR_BOTH, recur, rot);
    else if(strcmp(argv[1], "bump") == 0)
      data[0] = data[1] = data[2] = filter<num_t>(rgb2d<num_t>(data), BUMP_BOTH, recur, rot);
    else if(strcmp(argv[1], "rgb2xyz") == 0)
      data = rgb2xyz<num_t>(data);
    else if(strcmp(argv[1], "xyz2rgb") == 0)
      data = xyz2rgb<num_t>(data);
    else if(strcmp(argv[1], "w2b") == 0) {
      for(int i = 0; i < data[0].rows(); i ++)
        for(int j = 0; j < data[0].cols(); j ++)
          if(data[0](i, j) == num_t(1) &&
             data[1](i, j) == num_t(1) &&
             data[2](i, j) == num_t(1))
            data[0](i, j) = data[1](i, j) = data[2](i, j) = num_t(0);
    } else if(strcmp(argv[1], "b2w") == 0) {
      for(int i = 0; i < data[0].rows(); i ++)
        for(int j = 0; j < data[0].cols(); j ++)
          if(data[0](i, j) == num_t(0) &&
             data[1](i, j) == num_t(0) &&
             data[2](i, j) == num_t(0))
            data[0](i, j) = data[1](i, j) = data[2](i, j) = num_t(1);
    } else if(strcmp(argv[1], "b2wd") == 0) {
      vector<SimpleMatrix<num_t> > ddata;
      if(!loadp2or3<num_t>(ddata, argv[4]))
        return - 1;
      for(int i = 0; i < data[0].rows(); i ++)
        for(int j = 0; j < data[0].cols(); j ++)
          if((data[0](i, j) == num_t(0) &&
              data[1](i, j) == num_t(0) &&
              data[2](i, j) == num_t(0)) ||
             (data[0](i, j) == ddata[0](i, j) &&
              data[1](i, j) == ddata[1](i, j) &&
              data[2](i, j) == ddata[2](i, j)) )
            data[0](i, j) = data[1](i, j) = data[2](i, j) = num_t(1);
    } else if(strcmp(argv[1], "bit") == 0 && 0 < recur) {
      for(int i = 0; i < data.size(); i ++) {
        SimpleMatrix<num_t> work(data[i].rows() * recur, data[i].cols());
        work.O();
        for(int j = 0; j < data[i].rows(); j ++)
          for(int k = 0; k < data[i].cols(); k ++)
            for(int m = 0; m < recur; m ++)
              work(j + m * data[i].rows(), k) =
                num_t(1 & int(data[i](j, k) *
                  pow(num_t(int(2)), num_t(int(m + 1))) ) );
        swap(work, data[i]);
      }
    } else if((strcmp(argv[1], "bit") == 0 && recur < 0) || strcmp(argv[1], "nbit") == 0) {
      for(int i = 0; i < data.size(); i ++) {
        SimpleMatrix<num_t> work(data[i].rows() / abs(recur), data[i].cols());
        work.O();
        for(int j = 0; j < work.rows(); j ++)
          for(int k = 0; k < work.cols(); k ++)
            for(int m = 0; m < abs(recur); m ++)
              work(j, k) += strcmp(argv[1], "nbit") == 0 ?
                (sgn<num_t>(data[i](j + m * data[i].rows() / abs(recur), k) -
                  num_t(int(1)) / num_t(int(2))) + num_t(int(1))) /
                    pow(num_t(int(2)), num_t(int(m + 1))) :
                data[i](j + m * data[i].rows() / abs(recur), k) /
                  pow(num_t(int(2)), num_t(int(m + 1)));
        swap(work, data[i]);
      }
    } else if(strcmp(argv[1], "slide") == 0 && 0 < recur)
      for(int i = 0; i < data.size(); i ++)
        for(int j = 0; j < data[i].rows(); j ++)
          for(int k = 0; k < data[i].cols(); k ++)
            data[i](j, k) = num_t((int(data[i](j, k) * num_t(int(256))) << recur) & 0xff) +
              num_t(int(data[i](j, k) * num_t(int(256))) >> (8 - recur));
    else if((strcmp(argv[1], "slide") == 0 && recur < 0) || strcmp(argv[1], "nslide") == 0)
      for(int i = 0, ns = ! strcmp(argv[1], "nslide"); i < data.size(); i ++)
        for(int j = 0; j < data[i].rows(); j ++)
          for(int k = 0; k < data[i].cols(); k ++)
            data[i](j, k) = num_t(ns ? 0 : (int(data[i](j, k) * num_t(int(256))) << (8 - (- recur))) & 0xff) +
              num_t(int(data[i](j, k) * num_t(int(256))) >> (- recur));
    if(strcmp(argv[1], "nslide") == 0) {
      data = normalize<num_t>(data);
      for(int i = 0; i < data.size(); i ++)
        data[i] *= pow(num_t(int(2)), - num_t(abs(recur)));
    }
    if(!savep2or3<num_t>(argv[3],
        strcmp(argv[1], "b2w") != 0 && strcmp(argv[1], "b2wd") != 0 &&
        strcmp(argv[1], "nop") != 0 && strcmp(argv[1], "limit") != 0 &&
        strcmp(argv[1], "bump") != 0 && strcmp(argv[1], "nslide") != 0 ?
        normalize<num_t>(data) : data,
        strcmp(argv[1], "limit") == 0 ? recur : 65535) )
      return - 1;
  } else if(strcmp(argv[1], "reshape") == 0 ||
            strcmp(argv[1], "recolor") == 0 ||
            strcmp(argv[1], "recolor2") == 0 ||
            strcmp(argv[1], "recolor3") == 0) {
    if(argc < 6) goto usage;
    const int count(atoi(argv[2]));
    vector<SimpleMatrix<num_t> > datac, datas;
    if(!loadp2or3<num_t>(datac, argv[3]))
      return - 1;
    if((strcmp(argv[1], "recolor") == 0 ||
        strcmp(argv[1], "recolor3") == 0 ||
        strcmp(argv[1], "reshape") == 0) &&
       ! loadp2or3<num_t>(datas, argv[4]))
      return - 1;
    if(strcmp(argv[1], "reshape") == 0) {
      const SimpleMatrix<num_t> datav(rgb2d<num_t>(datas));
      for(int i = 0; i < datac.size(); i ++)
        datac[i] = reShape<num_t>(datac[i], datav, count, myatof(argv[6]) );
    } else if(strcmp(argv[1], "recolor2") == 0)
      for(int i = 0; i < datac.size(); i ++)
        datac[i] = reColor<num_t>(datac[i], count, myatof(argv[5]));
    else if(strcmp(argv[1], "recolor3") == 0)
      for(int i = 0; i < datac.size(); i ++)
        datac[i] = reColor3<num_t>(datac[i], datas[i], count);
    else {
      vector<SimpleMatrix<num_t> > xyzc(rgb2xyz<num_t>(datac));
      vector<SimpleMatrix<num_t> > xyzs(rgb2xyz<num_t>(datas));
      for(int i = 0; i < xyzc.size(); i ++)
        xyzc[i] = strcmp(argv[1], "recolor") == 0 ?
          reColor<num_t>(xyzc[i], xyzs[i], count, myatof(argv[6])) :
          reColor<num_t>(xyzc[i], count, myatof(argv[5]));
      datac = xyz2rgb<num_t>(xyzc);
    }
    if(!savep2or3<num_t>(argv[strcmp(argv[1], "recolor2") == 0 ? 4 : 5], normalize<num_t>(datac)))
      return - 1;
  } else if(strcmp(argv[1], "obj" ) == 0) {
    vector<SimpleMatrix<num_t> > data, mask;
    if(argc < 4) goto usage;
    if(!loadp2or3<num_t>(data, argv[2]))
      return - 1;
    SimpleMatrix<num_t> sd(rgb2d<num_t>(data));
    const int rows(sd.rows());
    const int cols(sd.cols());
          vector<SimpleVector<num_t> > points(getTileVec<num_t>(move(sd)));
    if(argv[1][0] == 'O')
      for(int i = 0; i < points.size(); i ++) points[i][2] = - points[i][2];
    saveobj<num_t>(points, num_t(rows), num_t(cols),
                   mesh2<num_t>(points), argv[3]);
    saveMTL<num_t>(argv[3], (string(argv[3]) + string(".mtl")).c_str());
  } else if(strcmp(argv[1], "tilt") == 0 ||
            strcmp(argv[1], "sbox") == 0) {
    if(argc < 8) goto usage;
    const int index(atoi(argv[2]));
    const int Mindex(atoi(argv[3]));
    num_t psi(myatof(argv[4]));
    vector<SimpleMatrix<num_t> > data, bump;
    vector<SimpleVector<num_t> > points;
    vector<SimpleVector<int>   > polys;
    if(!loadp2or3<num_t>(data, argv[5]))
      return - 2;
    const string fn(argv[5]);
    if(!loadp2or3<num_t>(bump, argv[6]))
      return - 2;
    assert(strlen("tilt" ) == strlen("sbox" ));
    bump[0] = rgb2d<num_t>(bump);
    const SimpleMatrix<num_t> tilt0(tilt<num_t>(bump[0] * num_t(int(0)),
      triangles<num_t>(makeRefMatrix<num_t>(data[0], 1), bump[0],
        strncmp(argv[1], "sbox", strlen("sbox")) == 0 ?
          match_t<num_t>() :
          tiltprep<num_t>(data[0], index, Mindex, psi) ),
        strncmp(argv[1], "sbox", strlen("sbox")) == 0 ?
          num_t(index) / num_t(Mindex) *
            num_t(min(data[0].rows(), data[0].cols())) :
          - num_t(1000000)
      ));
    for(int j = 0; j < data.size(); j ++)
      data[j] = pullRefMatrix<num_t>(tilt0, 1, data[j]);
    if(! savep2or3<num_t>(argv[7], data))
      return - 1;
  } else if(strcmp(argv[1], "match") == 0) {
    if(argc < 10) goto usage;
    const int nsub(atoi(argv[2]));
    const int nemph(atoi(argv[3]));
    const int vboxdst(atoi(argv[4]));
    const int vboxsrc(atoi(argv[5]));
    vector<SimpleMatrix<num_t> > in0, in1, bump0, bump1;
    if(!loadp2or3<num_t>(in0, argv[6]))
      return - 2;
    if(!loadp2or3<num_t>(in1, argv[7]))
      return - 2;
    if(!loadp2or3<num_t>(bump0, argv[8]))
      return - 2;
    if(!loadp2or3<num_t>(bump1, argv[9]))
      return - 2;
    const string outbase(argv[10]);
    const vector<SimpleVector<num_t> > shape0(vboxdst < 0
      ? getHesseVec<num_t>(rgb2d<num_t>(bump0), abs(vboxdst))
      : getTileVec<num_t>(rgb2d<num_t>(bump0), abs(vboxdst)));
    const vector<SimpleVector<num_t> > shape1(vboxsrc < 0
      ? getHesseVec<num_t>(rgb2d<num_t>(bump1), abs(vboxsrc))
      : getTileVec<num_t>(rgb2d<num_t>(bump1), abs(vboxsrc)));
    const vector<match_t<num_t> > m(matchPartialR<num_t>(shape0, shape1, nsub));
    vector<SimpleMatrix<num_t> > out;
    out.resize(3);
    const SimpleMatrix<num_t> rin0(makeRefMatrix<num_t>(in0[0], 1));
    const SimpleMatrix<num_t> rin1(makeRefMatrix<num_t>(in1[0], 1 + rin0.rows() * rin0.cols()));
    vector<vector<SimpleVector<int> > > mhull0, mhull1;
    mhull0.reserve(m.size());
    mhull1.reserve(m.size());
    for(int i = 0; i < m.size(); i ++) {
      mhull0.emplace_back(mesh2<num_t>(shape0, m[i].dst));
      mhull1.emplace_back((~ m[i]).hullConv(mhull0[i]));
    }
    out[0] = SimpleMatrix<num_t>(in0[0].rows(), in0[0].cols()).O();
    for(int i = 0; i < m.size(); i ++)
      out[0] += showMatch<num_t>(draw<num_t>(in0[0] * num_t(0),
                     shape0, mhull0[i]), shape0, mhull0[i]);
    out[1] = out[2] = out[0];
    savep2or3<num_t>((outbase + string("-repl0.ppm")).c_str(), normalize<num_t>(out) );
    out[0] = SimpleMatrix<num_t>(in1[0].rows(), in1[0].cols()).O();
    for(int i = 0; i < m.size(); i ++)
      out[0] += showMatch<num_t>(draw<num_t>(in1[0] * num_t(0),
                   shape1, mhull1[i]), shape1, mhull1[i]);
    out[1] = out[2] = out[0];
    savep2or3<num_t>((outbase + string("-repl1.ppm")).c_str(), normalize<num_t>(out) );
    out[0] = SimpleMatrix<num_t>(in0[0].rows(), in0[0].cols()).O();
    for(int i = 0; i < m.size(); i ++)
      out[0] += showMatch<num_t>(draw<num_t>(in0[0] * num_t(0),
                   m[i].transform(shape1), mhull1[i]),
                   m[i].transform(shape1), mhull1[i]);
    out[1] = out[2] = out[0];
    savep2or3<num_t>((outbase + string("-repl2.ppm")).c_str(), normalize<num_t>(out) );
    for(int i = 0; i < nemph; i ++) {
      const num_t iemph(num_t(i) / num_t(nemph));
      SimpleMatrix<num_t> reref(rin1.rows(), rin1.cols());
      reref.O();
      for(int i = 0; i < m.size(); i ++) {
        const SimpleMatrix<num_t> rd(draw<num_t>(rin1, shape1,
          takeShape<num_t>(shape1, shape0, m[i], iemph), mhull1[i]));
        for(int j = 0; j < min(reref.rows(), rd.rows()); j ++)
          for(int k = 0; k < min(reref.cols(), rd.cols()); k ++)
            if(rd(j, k) != num_t(0)) reref(j, k) = rd(j, k);
      }
      for(int idx = 0; idx < out.size(); idx ++)
        out[idx] = pullRefMatrix<num_t>(reref, 1 + rin0.rows() * rin0.cols(), in1[idx]);
      savep2or3<num_t>((outbase + string("-") + to_string(i) + string("-") +
                       to_string(nemph) + string(".ppm")).c_str(), out);
    }
  } else if(strcmp(argv[1], "cat") == 0) {
    if(argc < 3) goto usage;
    vector<vector<SimpleMatrix<num_t> > > in;
    in.resize(argc - 2);
    for(int i = 2; i < argc; i ++) {
      vector<SimpleMatrix<num_t> > ibuf;
      if(!loadp2or3<num_t>(ibuf, argv[i]))
        return - 2;
      in[i - 2] = move(ibuf);
    }
    vector<SimpleMatrix<num_t> > glay;
    glay.reserve(in.size());
    for(int i = 0; i < in.size(); i ++)
      glay.emplace_back(in[i].size() == 1 ? in[i][0] : rgb2d<num_t>(in[i]));
    const vector<vector<int> > cat(catImage<num_t>(glay));
    for(int i = 0; i < cat.size(); i ++) {
      for(int j = 0; j < cat[i].size(); j ++)
        cout << argv[2 + cat[i][j]] << endl;
      cout << endl;
    }
  } else if(strcmp(argv[1], "habit") == 0) {
    vector<SimpleVector<num_t> > pdst, psrc;
    vector<SimpleVector<int>   > poldst, polsrc;
    if(argc < 5 || !loadobj<num_t>(pdst, poldst, argv[2]) ||
                   !loadobj<num_t>(psrc, polsrc, argv[3])) goto usage;
    num_t Mx(0), My(0);
    for(int i = 0; i < pdst.size(); i ++) {
      My = max(num_t(My), abs(pdst[i][0]));
      Mx = max(num_t(Mx), abs(pdst[i][1]));
    }
    for(int i = 0; i < psrc.size(); i ++) {
      My = max(num_t(My), abs(psrc[i][0]));
      Mx = max(num_t(Mx), abs(psrc[i][1]));
    }
    saveobj<num_t>(takeShape<num_t>(pdst, psrc,
      matchPartialR<num_t>(pdst, psrc)[0],
      num_t(1) / num_t(2)), My, Mx, poldst, argv[4]);
#if !defined(_OLDCPP_)
  } else if(strcmp(argv[1], "omake") == 0) {
    vector<vector<num_t> > data;
    string header;
    loaddat<num_t>(argv[3], header, data);
    SimpleMatrix<num_t> buf(atoi(argv[4]), atoi(argv[4]));
    const SimpleMatrix<complex(num_t)> mdft(dft<num_t>(buf.rows()));
    const SimpleMatrix<complex(num_t)> midft(dft<num_t>(- buf.rows()));
    for(int i0 = 1; i0 < data.size(); i0 ++) {
      for(int i = 0; i <= data[i0].size() / buf.rows() / buf.rows(); i ++) {
        for(int k = 0; k < buf.cols(); k ++)
          for(int j = 0; j < buf.rows(); j ++) {
            const int idx(i * buf.rows() * buf.rows() + k * buf.rows() + j);
            buf(j, k) = idx < data[i0].size() ? data[i0][idx] : num_t(0);
          }
        SimpleMatrix<num_t> buf2;
        if(strcmp(argv[2], "diff") == 0)
          buf2 = (midft * (
            filter<num_t>(mdft.template real<num_t>() * buf, COLLECT_BOTH).template cast<complex(num_t) >() +
            filter<num_t>(mdft.template imag<num_t>() * buf, COLLECT_BOTH).template cast<complex(num_t) >() * (complex(num_t))(num_t(0), num_t(1)) ) ).template real<num_t>();
        else if(strcmp(argv[2], "sharpen") == 0)
          buf2 = (midft * (
            filter<num_t>(mdft.template real<num_t>() * buf, SHARPEN_X).template cast<complex(num_t) >() +
            filter<num_t>(mdft.template imag<num_t>() * buf, SHARPEN_X).template cast<complex(num_t) >() * (complex(num_t))(num_t(0), num_t(1)) ) ).template real<num_t>();
        else if(strcmp(argv[2], "bump") == 0)
          buf2 = (midft * (
            filter<num_t>(mdft.template real<num_t>() * buf, BUMP_BOTH).template cast<complex(num_t) >() +
            filter<num_t>(mdft.template imag<num_t>() * buf, BUMP_BOTH).template cast<complex(num_t) >() * (complex(num_t))(num_t(0), num_t(1)) ) ).template real<num_t>();
        for(int k = 0; k < buf2.cols(); k ++)
          for(int j = 0; j < buf2.rows(); j ++) {
            const int idx(i * buf.rows() * buf.rows() + k * buf.rows() + j);
            if(idx < data[i0].size())
              data[i0][idx] = buf2(j, k);
            else
              break;
          }
      }
    }
    num_t M(0);
    for(int i = 1; i < data.size(); i ++) {
      vector<num_t> sdata(data[i]);
      for(int i = 0; i < sdata.size(); i ++)
        sdata[i] = abs(sdata[i]);
      sort(sdata.begin(), sdata.end());
      M = max(M, sdata[sdata.size() * 7 / 8] * num_t(4));
    }
    cout << header;
    for(int i = 0; i < data[0].size(); i ++) {
      cout << data[0][i] << " ";
      for(int j = 1; j < data.size(); j ++)
        cout << (i < data[j].size() ? data[j][i] / M : num_t(0)) << " ";
      cout << endl;
    }
#endif
  } else goto usage;
  return 0;
 usage:
  cout << "Usage:" << endl;
  cout << argv[0] << " (collect|sharpen|blur|bump|enlarge|shrink|flarge|blink|nop|limit|bit|nbit|slide|nslide) <input.ppm> <output.ppm> <recur> <rot>" << endl;
  cout << argv[0] << " cat <input0-4.ppm> ..." << endl;
  cout << argv[0] << " (tilt|sbox) <index> <max_index> <psi> <input.ppm> <input-bump.ppm> <output.ppm>" << endl;
  cout << argv[0] << " obj <input.ppm> <output.obj>" << endl;
  cout << argv[0] << " match <nsub> <nemph> <vbox_dst> <vbox_src> <dst.ppm> <src.ppm> <dst-bump.ppm> <src-bump.ppm> <output-basename>" << endl;
  cout << argv[0] << " reshape <num_shape_per_color> <input_color.ppm> <input_shape.ppm> <output.ppm>" << endl;
  cout << argv[0] << " recolor <num_shape_per_color> <input_color.ppm> <input_shape.ppm> <output.ppm> <intensity>" << endl;
  cout << argv[0] << " recolor2 <num_shape_per_color> <input_color.ppm> <output.ppm> <intensity>" << endl;
  cout << argv[0] << " recolor3 <num_shape_per_color> <input_color.ppm> <input_shape.ppm> <output.ppm>" << endl;
  cout << argv[0] << " habit <in0.obj> <in1.obj> <out.obj>" << endl;
  return - 3;
}

