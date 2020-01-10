#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <assert.h>

#if defined(_WITH_MPFR_)
#include <mpreal.h>
typedef mpfr::mpreal num_t;
using std::sqrt;
using mpfr::pow;
using mpfr::log;
using mpfr::isfinite;
#elif defined(_WITH_NO_FLOAT_)
#include "ifloat.hh"
typedef SimpleFloat<uint16_t, uint32_t, 16, char> num_t;
//typedef SimpleFloat<uint32_t, uint64_t, 32, short> num_t;
#else
typedef double num_t;
#endif

#if defined(_WITHOUT_EIGEN_)
#if defined(_WITH_NO_FLOAT_)
template <typename T> using complex = Complex<T>;
#else
#include <complex>
using std::complex;
#endif
#include <cstring>
#include "simplelin.hh"
#else
#include <Eigen/Core>
#include <Eigen/LU>
using std::complex;
#endif

#if ! defined(_WITH_NO_FLOAT_)
using std::sqrt;
using std::exp;
using std::log;
using std::pow;
using std::sin;
using std::cos;
using std::tan;
using std::atan2;
using std::ceil;
#endif

using std::max;
using std::min;

#if defined(_WITH_GLTF2_)
#include <fx/gltf.h>
#endif

namespace goki {
using std::pair;
using std::make_pair;
using std::vector;
#include "p0.hh"
#include "fileio.hh"
#include "enlarge.hh"
#include "match.hh"
#include "redig.hh"
};
using namespace goki;
using std::cout;
using std::endl;

void usage() {
  cout << "Usage:" << endl;
  cout << "gokicheck enlarge <ratio>  <input.ppm> <output.ppm>" << endl;
  cout << "gokicheck cenl    <ratio>  <inputdst.ppm> <inputsrc.ppm> <output.ppm>" << endl;
  cout << "gokicheck pextend <pixels> <input.ppm> <output.ppm>" << endl;
  cout << "gokicheck collect <input.ppm> <output.ppm>" << endl;
  cout << "gokicheck bump    <input.ppm> <output.ppm>" << endl;
  cout << "gokicheck bump2   <input0.ppm> <input1.ppm> <delta_pixels> <output.ppm>" << endl;
  cout << "gokicheck reshape <num_shape_per_color> <input_color.ppm> <input_shape.ppm> <output.ppm>" << endl;
  cout << "gokicheck obj     <shift_x_pixels> <gather_pixels> <zratio> <input.ppm> <mask.ppm>? <output.obj>" << endl;
  cout << "gokicheck obj     stand <gather_pixels> <thin> <ratio> <zratio> <input.ppm> <mask.ppm>? <output.obj>" << endl;
  cout << "gokicheck tilt    <index> <max_index> <psi> <shift_x_pixels> <input.ppm> <input-bump.(ppm|obj)> <output.ppm>" << endl;
  cout << "gokicheck sbox    <index> <max_index> <input.ppm> <input-bump.(ppm|obj)> <output.ppm>" << endl;
  cout << "gokicheck draw    <input-mask.ppm> <input-obj.(obj|gltf)> <output.ppm>" << endl;
  cout << "gokicheck drawr   <input-mask.ppm> <input-obj.(obj|gltf)> <output.ppm>" << endl;
  cout << "gokicheck drawm   <input-mask.ppm> <input-obj.(obj|gltf)> <output.ppm>" << endl;
  cout << "gokicheck match   <num_of_res_shown> <num_of_hidden_match> <vbox_dst> <vbox_src> <dst.ppm> <src.ppm> <dst-bump.(ppm|obj)> <src-bump.(ppm|obj|gltf)> (<dst-mask.ppm> <src-mask.ppm>)? <output-basename>" << endl;
  cout << "gokicheck matcho  <match> <emph> <vbox_dst> <vbox_src> <dst.ppm> <src.ppm> <dst-bump.(ppm|obj)> <src-bump.(ppm|obj|gltf)> (<dst-mask.ppm> <src-mask.ppm>)? <output-basename>" << endl;
  cout << "gokicheck habit   <in0.obj> <in1.obj> (<index> <max_index> <psi>)? <out.obj>" << endl;
  return;
}

template <typename T> void saveMatches(const std::string& outbase, const match_t<T>& match, const std::vector<typename simpleFile<T>::Vec3>& shape0, const std::vector<typename simpleFile<T>::Vec3>& shape1, const typename simpleFile<T>::Mat in0[3], const typename simpleFile<T>::Mat in1[3], const typename simpleFile<T>::Mat& bump0, const typename simpleFile<T>::Mat& bump1, std::vector<typename simpleFile<T>::Veci3>& ohull0, std::vector<typename simpleFile<T>::Veci3>& ohull1, const T& emph = T(.5)) {
  assert(in0[0].rows() == in1[0].rows() && in0[0].cols() == in1[0].cols());
  reDig<T>      redig;
  simpleFile<T> file;
  typename simpleFile<T>::Mat outs[3];
  const auto rin0(redig.makeRefMatrix(in0[0], 1));
  const auto rin1(redig.makeRefMatrix(in1[0], 1 + rin0.rows() * rin0.cols()));
  const auto mhull0(redig.delaunay2(shape0, match.dstpoints));
  const auto mhull1((~ match).hullConv(mhull0));
  for(int idx = 0; idx < 3; idx ++)
    outs[idx] = redig.showMatch(redig.replace(in0[idx] * T(0), shape0, ohull0),
                                shape0, mhull0);
  file.savep2or3((outbase + std::string("-repl0.ppm")).c_str(), outs, false);
  for(int idx = 0; idx < 3; idx ++)
    outs[idx] = redig.showMatch(redig.replace(in1[idx] * T(0),
                                  match.transform(shape1), ohull1),
                                match.transform(shape1), mhull1);
  file.savep2or3((outbase + std::string("-repl1.ppm")).c_str(), outs, false);
  const auto reref(redig.emphasis(rin0, rin1, bump1, shape0, shape1,
                                  match, mhull0, mhull1, emph));
  for(int idx = 0; idx < 3; idx ++)
    outs[idx] = redig.showMatch(in0[idx], shape0, mhull0);
  file.savep2or3((outbase + std::string("-c0.ppm")).c_str(), outs, false);
  for(int idx = 0; idx < 3; idx ++)
    outs[idx] = redig.showMatch(redig.pullRefMatrix(reref,
                    1 + rin0.rows() * rin0.cols(), in1[idx]),
                    match.transform(shape1), mhull1);
  file.savep2or3((outbase + std::string("-c1.ppm")).c_str(), outs, false);
  for(int idx = 0; idx < 3; idx ++)
    outs[idx] = redig.pullRefMatrix(reref, 1 + rin0.rows() * rin0.cols(), in1[idx]);
  file.savep2or3((outbase + std::string("-c2.ppm")).c_str(), outs, false);
  file.saveobj(redig.takeShape(shape0, shape1,   match, mhull0, mhull1, emph),
               outs[0].rows(), outs[0].cols(),
               ohull0, (outbase + std::string("-emph0.obj")).c_str());
  file.saveobj(redig.takeShape(shape1, shape0, ~ match, mhull1, mhull0, emph),
               outs[0].rows(), outs[0].cols(),
               ohull1, (outbase + std::string("-emph1.obj")).c_str());
/*
  file.saveglTF((outbase + std::string(".gltf")).c_str(),
                shape0, mhull0, ~ match, center, bone);
*/
  return;
}

template <typename T> void resizeDst2(typename simpleFile<T>::Mat mout[3], typename simpleFile<T>::Mat& bump1, typename simpleFile<T>::Mat& mmout1, const typename simpleFile<T>::Mat data[3], const typename simpleFile<T>::Mat& bump, const typename simpleFile<T>::Mat& mask, const int& h, const int& w) {
  mmout1 = bump1 = mout[0] = mout[1] = mout[2] = typename simpleFile<T>::Mat(h, w);
  for(int i = 0; i < min(mout[0].rows(), data[0].rows()); i ++) {
    for(int j = 0; j < min(mout[0].cols(), data[0].cols()); j ++) {
      mout[0](i, j) = data[0](i, j);
      mout[1](i, j) = data[1](i, j);
      mout[2](i, j) = data[2](i, j);
      mmout1(i, j)  = mask(i, j);
      bump1(i, j)   = bump(i, j);
    }
    for(int j = min(mout[0].cols(), data[0].cols()); j < mout[0].cols(); j ++)
      mmout1(i, j) = bump1(i, j) = mout[0](i, j) = mout[1](i, j) = mout[2](i, j) = T(0);
  }
  for(int i = min(mout[0].rows(), data[0].rows()); i < mout[0].rows(); i ++)
    for(int j = 0; j < mout[0].cols(); j ++)
      mmout1(i, j) = bump1(i, j) = mout[0](i, j) = mout[1](i, j) = mout[2](i, j) = T(0);
  return;
}

int main(int argc, const char* argv[]) {
  if(argc < 2) {
    usage();
    return 0;
  }
  simpleFile<num_t> file;
  reDig<num_t>      redig;
  if(strcmp(argv[1], "test") == 0) {
    typename simpleFile<num_t>::Mat data[3];
    if(!file.loadp2or3(data, argv[2]))
      return -1;
    if(!file.savep2or3(argv[3], data, ! true))
      return - 1;
  } else if(strcmp(argv[1], "sharpen") == 0 ||
     strcmp(argv[1], "cenl")    == 0 ||
     strcmp(argv[1], "pextend") == 0 ||
     strcmp(argv[1], "light")   == 0) {
    if(argc < 5 || (strcmp(argv[1], "cenl") == 0 && argc < 6)) {
      usage();
      return 0;
    }
    const auto ratio(std::atoi(argv[2]));
    typename simpleFile<num_t>::Mat data[3];
    if(!file.loadp2or3(data, argv[3]))
      return - 1;
    if(strcmp(argv[1], "sharpen") == 0) {
      Filter<num_t> enlarger;
      for(int i = 0; i < 3; i ++) {
        for(int j = 0; j < ratio; j ++)
          data[i] = enlarger.compute(data[i], enlarger.SHARPEN_BOTH);
        data[i] = enlarger.compute(data[i], enlarger.CLIP);
      }
    } else if(strcmp(argv[1], "pextend") == 0) {
      Filter<num_t> extender;
      extender.plen = ratio;
      for(int i = 0; i < 3; i ++)
        data[i] = extender.compute(extender.compute(data[i], extender.EXTEND_BOTH), extender.CLIP);
    } else if(strcmp(argv[1], "cenl") == 0) {
      typename simpleFile<num_t>::Mat datas[3];
      Filter<num_t> cenl;
      if(!file.loadp2or3(datas, argv[4]))
        return - 1;
      const auto t(num_t(ratio) / num_t(100));
      for(int i = 0; i < 3; i ++)
        data[i] = cenl.compute(data[i] * t + datas[i] * (num_t(1) - t), cenl.CLIP);
      if(!file.savep2or3(argv[5], data, ! true))
        return - 1;
      return 0;
    } else if(strcmp(argv[1], "light") == 0) {
      Filter<num_t> enlarger;
      typename simpleFile<num_t>::Mat yrev[3];
      for(int i = 0; i < 3; i ++) {
        yrev[i].resize(data[i].rows(), data[i].cols());
        for(int j = 0; j < data[i].rows(); j ++)
          yrev[i].row(j) = data[i].row(data[i].rows() - 1 - j);
      }
      for(int i = 0; i < 3; i ++) {
        for(int j = 0; j < ratio; j ++) {
          data[i] = enlarger.compute(data[i], enlarger.SHARPEN_Y);
          yrev[i] = enlarger.compute(yrev[i], enlarger.SHARPEN_Y);
        }
        data[i] = enlarger.compute(data[i], enlarger.CLIP);
        yrev[i] = enlarger.compute(yrev[i], enlarger.CLIP);
      }
      for(int i = 0; i < 3; i ++)
        for(int j = 0; j < data[i].rows(); j ++)
          data[i].row(j) += yrev[i].row(yrev[i].rows() - 1 - j);
      redig.normalize(data, 1.);
    }
    if(!file.savep2or3(argv[4], data, ! true, 65535))
      return - 1;
  } else if(strcmp(argv[1], "collect") == 0 ||
            strcmp(argv[1], "bump")  == 0) {
    if(argc < 4) {
      usage();
      return 0;
    }
    typename simpleFile<num_t>::Mat data[3];
    if(!file.loadp2or3(data, argv[2]))
      return - 1;
    if(strcmp(argv[1], "collect") == 0) {
      Filter<num_t> detect;
      for(int i = 0; i < 3; i ++)
        data[i] = detect.compute(detect.compute(data[i], detect.COLLECT_BOTH), detect.CLIP);
    } else if(strcmp(argv[1], "bump") == 0) {
#if defined(_WITH_MPFR_)
      num_t::set_default_prec(_WITH_MPFR_);
#endif
      Filter<num_t> bump;
      data[0] = data[1] = data[2] = bump.compute(redig.rgb2d(data).template cast<num_t>(), bump.BUMP_BOTH);
    }
    redig.normalize(data, num_t(1));
    if(!file.savep2or3(argv[3], data, ! true))
      return - 1;
  } else if(strcmp(argv[1], "bump2") == 0) {
    if(argc < 5) {
      usage();
      return 0;
    }
    typename simpleFile<num_t>::Mat data0[3], data1[3], out[3];
    if(!file.loadp2or3(data0, argv[2]))
      return - 1;
    if(!file.loadp2or3(data1, argv[3]))
      return - 1;
    const auto pixels(std::atof(argv[4]));
    Filter<num_t> bump;
    out[0] = out[1] = out[2] = redig.autoLevel(- bump.bump2(redig.rgb2d(data0).transpose(), redig.rgb2d(data1).transpose(), pixels).transpose(), data0[0].rows() + data0[0].cols());
    redig.normalize(out, num_t(1));
    if(!file.savep2or3(argv[5], out, ! true))
      return - 1;
  } else if(strcmp(argv[1], "reshape") == 0) {
    if(argc < 6) {
      usage();
      return 0;
    }
    int count(std::atoi(argv[2]));
    typename simpleFile<num_t>::Mat datac[3], datas[3];
    if(!file.loadp2or3(datac, argv[3]))
      return - 1;
    if(!file.loadp2or3(datas, argv[4]))
      return - 1;
    for(int i = 0; i < 3; i ++)
      datac[i] = redig.reShape(datac[i], datas[i], count);
    redig.normalize(datac, num_t(1));
    if(!file.savep2or3(argv[5], datac, ! true))
      return - 1;
  } else if(strcmp(argv[1], "obj") == 0) {
    typename simpleFile<num_t>::Mat data[3], mask[3];
    num_t xoffset(0);
    int    vbox(2);
    num_t addstand(0);
    num_t ratio(1);
    num_t zratio(1);
    int    sidx(0);
    if(strcmp(argv[2], "stand") == 0) {
      if(argc < 9) {
        usage();
        return - 1;
      }
      vbox     = std::atoi(argv[3]);
      addstand = std::atof(argv[4]);
      ratio    = std::atof(argv[5]);
      if(!file.loadp2or3(data, argv[7]))
        return - 1;
      zratio = std::atof(argv[6]);
      if(9 < argc) {
        if(!file.loadp2or3(mask, argv[8]))
          return - 1;
        sidx = 9;
      } else {
        mask[0] = mask[1] = mask[2] = data[0] * num_t(0);
        sidx = 8;
      }
    } else {
      if(argc < 7) {
        usage();
        return - 1;
      }
      xoffset = std::atof(argv[2]);
      vbox    = std::atoi(argv[3]);
      if(!file.loadp2or3(data, argv[5]))
        return - 1;
      zratio  = num_t(std::atof(argv[4])) * sqrt(num_t(data[0].rows() * data[0].cols()));
      if(7 < argc) {
        if(!file.loadp2or3(mask, argv[6]))
          return - 1;
        sidx = 7;
      } else
        sidx = 6;
    }
    std::vector<typename simpleFile<num_t>::Vec3>  points;
    std::vector<typename simpleFile<num_t>::Veci3> facets;
    redig.initialize(vbox);
    redig.getTileVec(data[0], points, facets);
    const auto edges(redig.getEdges(mask[0], points));
    if(edges.size())
      redig.maskVectors(points, facets, mask[0]);
    for(int i = 0; i < points.size(); i ++)
      points[i] *= ratio;
    num_t M(points[0][2]), m(points[0][2]);
    for(int i = 1; i < points.size(); i ++) {
      M = std::max(points[i][2], M);
      m = std::min(points[i][2], m);
    }
    if(M == m) M += num_t(1);
    for(int i = 0; i < points.size(); i ++)
      points[i][2] *= zratio / (M - m);
    file.saveobj(points, ratio * num_t(data[0].rows()), ratio * num_t(data[0].cols()), facets, argv[sidx], edges, addstand, xoffset);
  } else if(strcmp(argv[1], "tilt") == 0 ||
            strcmp(argv[1], "sbox") == 0) {
    if((strcmp(argv[1], "tilt") == 0 && argc < 9) ||
       (strcmp(argv[1], "sbox") == 0 && argc < 7)) {
      usage();
      return - 1;
    }
    int    index(std::atoi(argv[2]));
    int    Mindex(std::atoi(argv[3]));
    num_t psi(0);
    int    offsetx(0);
    int    ipidx(6);
    int    iidx(7);
    int    oidx(8);
    if(strcmp(argv[1], "tilt") == 0) {
      psi = std::atof(argv[4]);
      offsetx = std::atoi(argv[5]);
    } else {
      ipidx = 4;
      iidx  = 5;
      oidx  = 6;
    }
    typename simpleFile<num_t>::Mat data[3], bump[3], out[3];
    std::vector<typename simpleFile<num_t>::Vec3>  points;
    std::vector<typename simpleFile<num_t>::Veci3> polys;
    if(!file.loadp2or3(data, argv[ipidx]))
      return - 2;
    const std::string fn(argv[iidx]);
    bool is_obj(false);
    if(fn[fn.size() - 1] == 'm') {
      if(!file.loadp2or3(bump, argv[iidx]))
        return - 2;
    } else if(fn[fn.size() - 1] == 'j') {
      if(!file.loadobj(points, polys, argv[iidx]))
        return - 2;
      is_obj = true;
    } else
      return - 2;
    typename simpleFile<num_t>::Mat tilt0;
    if(strcmp(argv[1], "sbox") == 0) {
      const match_t<num_t> mtilt;
      if(is_obj)
        tilt0 = redig.tilt(redig.makeRefMatrix(data[0], 1), redig.tiltprep(points, polys, redig.makeRefMatrix(data[0], 1), mtilt), num_t(index) / num_t(Mindex));
      else
        tilt0 = redig.tilt(redig.makeRefMatrix(data[0], 1), bump[0], mtilt, num_t(index) / num_t(Mindex));
    } else {
      auto mtilt(redig.tiltprep(data[0], index, Mindex, psi));
      mtilt.offset[2] += offsetx * data[0].cols();
      if(is_obj)
        tilt0 = redig.tilt(redig.makeRefMatrix(data[0], 1), redig.tiltprep(points, polys, redig.makeRefMatrix(data[0], 1), mtilt));
      else
        tilt0 = redig.tilt(redig.makeRefMatrix(data[0], 1), bump[0], mtilt);
    }
    for(int j = 0; j < 3; j ++)
      out[j] = redig.pullRefMatrix(tilt0, 1, data[j]);
    if(!file.savep2or3(argv[oidx], out, ! true))
      return - 1;
  } else if(strcmp(argv[1], "draw") == 0 ||
            strcmp(argv[1], "drawr") == 0 ||
            strcmp(argv[1], "drawm") == 0) {
    if(argc < 5) {
      usage();
      return - 1;
    }
    typename simpleFile<num_t>::Mat data[3];
    if(!file.loadp2or3(data, argv[2]))
      return - 2;
    std::vector<typename simpleFile<num_t>::Vec3>  datapoly;
    std::vector<typename simpleFile<num_t>::Veci3> polynorms;
    const std::string fn(argv[3]);
    if(fn[fn.size() - 1] == 'j') {
      if(!file.loadobj(datapoly, polynorms, argv[3]))
        return - 2;
    } else if(fn[fn.size() - 1] == 'f') {
      std::vector<typename simpleFile<num_t>::Vec3> center;
      std::vector<std::vector<typename simpleFile<num_t>::Veci4> > bone;
      if(!file.loadglTF(datapoly, polynorms, center, bone, argv[3]))
        return - 2;
    } else {
      usage();
      return - 1;
    }
    assert(datapoly.size());
    num_t My(datapoly[0][0]);
    num_t Mx(datapoly[0][1]);
    num_t my(My);
    num_t mx(Mx);
    for(int i = 0; i < datapoly.size(); i ++) {
      My = max(My, datapoly[i][0]);
      Mx = max(Mx, datapoly[i][1]);
      my = min(my, datapoly[i][0]);
      mx = min(mx, datapoly[i][1]);
    }
    match_t<num_t> m;
    const num_t ratio(sqrt(num_t(data[0].rows() * data[0].cols())) / num_t(max(My - my, Mx - mx)));
    m.ratio     *= ratio;
    m.offset[0] -= my * ratio;
    m.offset[1] -= mx * ratio;
    for(int i = 0; i < datapoly.size(); i ++)
      datapoly[i] = m.transform(datapoly[i]);
    typename simpleFile<num_t>::Mat res[3];
    std::vector<int> idx;
    for(int j = 0; j < datapoly.size(); j ++)
      idx.push_back(j);
    if(strcmp(argv[1], "draw") == 0)
      res[0] = res[1] = res[2] = redig.showMatch(data[0] * num_t(0), datapoly, polynorms, num_t(120));
    else if(strcmp(argv[1], "drawr") == 0) {
      auto mwork(data[0]);
      for(int i = 0; i < mwork.rows(); i ++)
        for(int j = 0; j < mwork.cols(); j ++)
          mwork(i, j) = num_t(1);
      auto prep(redig.tiltprep(datapoly, polynorms, mwork, match_t<num_t>()));
      for(int i = 0; i < prep.size(); i ++)
        prep[i].c = (prep[i].p(2, 0) + prep[i].p(2, 1) + prep[i].p(2, 2)) / num_t(3);
      res[0] = res[1] = res[2] = redig.tilt(data[0] * num_t(0), prep);
    } else {
      auto mwork(data[0]);
      for(int i = 0; i < mwork.rows(); i ++)
        for(int j = 0; j < mwork.cols(); j ++)
          mwork(i, j) = num_t(1);
      res[0] = res[1] = res[2] = mwork - redig.tilt(data[0] * num_t(0), redig.tiltprep(datapoly, polynorms, mwork, match_t<num_t>()));
    }
    redig.normalize(res, num_t(1));
    file.savep2or3(argv[4], res, true);
  } else if(strcmp(argv[1], "matcho") == 0 ||
            strcmp(argv[1], "match") == 0) {
    if(argc < 10) {
      usage();
      return - 1;
    }
    int fnout(10);
    int nshow(0);
    int nhid(0);
    num_t emph(0);
    match_t<num_t> m;
    if(strcmp(argv[1], "match") == 0) {
      nshow = std::atoi(argv[2]);
      nhid  = std::atoi(argv[3]);
    } else {
      std::ifstream input;
      input.open(argv[2]);
      if(input.is_open()) {
        try {
          input >> m;
          cerr << m;
        } catch(...) {
          usage();
          return - 2;
        }
      }
      input.close();
      emph = std::atof(argv[3]);
    }
    int vboxdst(std::atoi(argv[4]));
    int vboxsrc(std::atoi(argv[5]));
    typename simpleFile<num_t>::Mat data[3], data1[3], bdata[3], bdata1[3], mdata[3], mdata1[3], mout[3], mmout1, bump1;
    std::vector<typename simpleFile<num_t>::Veci3> delau0, delau1;
    std::vector<typename simpleFile<num_t>::Vec3>  shape0, shape1;
    if(!file.loadp2or3(data, argv[6]))
      return - 2;
    if(!file.loadp2or3(data1, argv[7]))
      return - 2;
    if(!file.loadp2or3(bdata, argv[8]))
      return - 2;
    const std::string fn(argv[9]);
    if(fn[fn.size() - 1] == 'j') {
      if(!file.loadobj(shape1, delau1, argv[9]))
        return - 2;
      bdata1[0] = bdata1[1] = bdata1[2] = data1[0] * num_t(0);
    } else if(fn[fn.size() - 1] == 'f') {
      std::vector<typename simpleFile<num_t>::Vec3> center;
      std::vector<std::vector<typename simpleFile<num_t>::Veci4> > bone;
      if(!file.loadglTF(shape1, delau1, center, bone, argv[9]))
        return - 2;
      bdata1[0] = bdata1[1] = bdata1[2] = data1[0] * num_t(0);
    } else if(fn[fn.size() - 1] == 'm') {
      if(!file.loadp2or3(bdata1, argv[9]))
        return - 2;
    } else {
      usage();
      return - 1;
    }
    if(11 < argc) {
      if(!file.loadp2or3(mdata, argv[10]))
        return - 2;
      if(!file.loadp2or3(mdata1, argv[11]))
        return - 2;
      fnout = 12;
    } else {
      mdata[0]  = mdata[1]  = mdata[2]  = data[0]  * num_t(0);
      mdata1[0] = mdata1[1] = mdata1[2] = data1[0] * num_t(0);
    }
    resizeDst2<num_t>(mout, bump1, mmout1, data1, bdata1[0], mdata1[0], data[0].rows(), data[0].cols());
    const auto& bump0(bdata[0]);
    redig.initialize(vboxdst);
    redig.getTileVec(bump0, shape0, delau0);
    redig.maskVectors(shape0, delau0, mdata[0]);
    if(fn[fn.size() - 1] == 'm') {
      redig.initialize(vboxsrc);
      redig.getTileVec(bump1, shape1, delau1);
      redig.maskVectors(shape1, delau1, mmout1);
    }
    if(strcmp(argv[1], "matcho") == 0)
      saveMatches<num_t>(std::string(argv[fnout]), m, shape0, shape1, data, mout, bump0, bump1, delau0, delau1, emph);
    else { 
      matchPartial<num_t> statmatch;
      auto matches(statmatch.match(shape0, shape1));
      if(fn[fn.size() - 1] == 'm')
        matches = statmatch.elim(matches, data, mout, bump1, shape1);
      matches.resize(min(int(matches.size()), nhid));
      std::cerr << matches.size() << "pending" << std::endl;
      for(int n = 0; n < min(int(matches.size()), nshow); n ++) {
        std::ofstream output;
        output.open((std::string(argv[fnout]) + std::to_string(n + 1) + std::string(".txt")).c_str());
        if(output.is_open()) {
          try {
            output << matches[n];
          } catch(...) {
            ;
          }
        }
        output.close();
        std::cerr << "Matchingsub: " << n << " / " << matches.size() << std::flush;
        std::vector<int> dstbuf(matches[n].dstpoints);
        std::vector<int> srcbuf(matches[n].srcpoints);
        std::vector<typename simpleFile<num_t>::Vec3> shape0a;
        std::vector<typename simpleFile<num_t>::Vec3> shape1a;
        std::sort(dstbuf.begin(), dstbuf.end());
        std::sort(srcbuf.begin(), srcbuf.end());
        for(int j = 0; j < shape0.size(); j ++)
          if(!binary_search(dstbuf.begin(), dstbuf.end(), j))
            shape0a.push_back(shape0[j]);
        for(int j = 0; j < shape1.size(); j ++)
          if(!binary_search(srcbuf.begin(), srcbuf.end(), j))
            shape1a.push_back(shape1[j]);
        matchPartial<num_t> pstatmatch;
        pstatmatch.ndiv /= 2;
        auto pmatches(pstatmatch.match(shape0a, shape1a));
        if(fn[fn.size() - 1] == 'm')
          pmatches = pstatmatch.elim(pmatches, data, mout, bump1, shape1a);
        pmatches.resize(min(int(pmatches.size()), nhid));
        for(int m = 0; m < min(int(pmatches.size()), nshow); m ++) {
          std::ofstream output;
          output.open((std::string(argv[fnout]) + std::to_string(n + 1) + std::string("-") + std::to_string(m + 1) + std::string(".txt")).c_str());
          if(output.is_open()) {
            try {
              output << pmatches[m];
            } catch(...) {
              ;
            }
          }
          output.close();
        }
      }
    }
  } else if(strcmp(argv[1], "habit") == 0) {
    std::vector<typename simpleFile<num_t>::Vec3>  pdst,   psrc;
    std::vector<typename simpleFile<num_t>::Veci3> poldst, polsrc;
    if(argc < 5 || !file.loadobj(pdst, poldst, argv[2]) ||
                   !file.loadobj(psrc, polsrc, argv[3])) {
      usage();
      return - 2;
    }
    num_t Mx(0), My(0);
    for(int i = 0; i < pdst.size(); i ++) {
      My = max(num_t(My), abs(pdst[i][0]));
      Mx = max(num_t(Mx), abs(pdst[i][1]));
    }
    for(int i = 0; i < psrc.size(); i ++) {
      My = max(num_t(My), abs(psrc[i][0]));
      Mx = max(num_t(Mx), abs(psrc[i][1]));
    }
    if(argc > 7) {
      const auto m(redig.tiltprep(typename simpleFile<num_t>::Mat(int(My), int(Mx)), - std::atoi(argv[4]), std::atoi(argv[5]), std::atof(argv[6])));
      file.saveobj(redig.takeShape(pdst, psrc, m, poldst, polsrc, num_t(1) / num_t(2)),
                   My, Mx, poldst, argv[7]);
    } else {
      matchPartial<num_t> statmatch;
      const auto m(statmatch.match(pdst, psrc)[0]);
      file.saveobj(redig.takeShape(pdst, psrc, m, poldst, polsrc, num_t(1) / num_t(2)),
                   My, Mx, poldst, argv[4]);
    }
  } else if(strcmp(argv[1], "omake") == 0) {
    std::vector<num_t> in;
    std::string line;
    while (std::getline(std::cin, line)) {
      std::stringstream sl(line);
      in.push_back(0);
      sl >> in[in.size() - 1];
    }
    simpleFile<num_t>::Mat buf(std::atoi(argv[2]), std::atoi(argv[2]));
    Filter<num_t> filter;
    num_t M(0);
    std::vector<num_t> rbuf;
    if(strcmp(argv[3], "sharpen") == 0)
      rbuf.resize(in.size() * 2, 0.);
    else
      rbuf.resize(in.size(), 0.);
    const auto dft0(filter.seed(buf.rows(), false));
    const auto idft0(filter.seed(buf.rows(), true));
    for(int i = 0; i <= in.size() / buf.rows() / buf.rows(); i ++) {
      for(int j = 0; j < buf.rows(); j ++)
        for(int k = 0; k < buf.cols(); k ++)
          buf(j, k) = in[i * buf.rows() * buf.rows() + k * buf.rows() + j];
      Filter<num_t>::MatU dft(dft0 * buf.template cast<complex<num_t> >());
      if(strcmp(argv[3], "sharpen") == 0) {
#if defined(_WITHOUT_EIGEN_)
        const auto rp(filter.compute(dft.template real<num_t>(), filter.SHARPEN_X));
        const auto ip(filter.compute(dft.template imag<num_t>(), filter.SHARPEN_X));
        const auto buf2((idft0.template cast<complex<num_t> >() * (rp.template cast<complex<num_t> >() + ip.template cast<complex<num_t> >() * sqrt(complex<num_t>(- num_t(1))))).template real<num_t>());
#else
        const auto rp(filter.compute(dft.real(), filter.SHARPEN_X));
        const auto ip(filter.compute(dft.imag(), filter.SHARPEN_X));
        const auto buf2((idft0.template cast<complex<num_t> >() * (rp.template cast<complex<num_t> >() + ip.template cast<complex<num_t> >() * sqrt(complex<num_t>(- num_t(l))))).real());
#endif
        for(int j = 0; j < buf2.rows(); j ++)
          for(int k = 0; k < buf2.cols(); k ++) {
            M = max(M, abs(buf2(j, k)));
            rbuf[i * buf.rows() * buf.rows() * 2 + k * buf.rows() + j] = buf2(j, k);
          }
        continue;
      } else {
        Filter<num_t>::Mat buf2;
        if(strcmp(argv[3], "diff") == 0){
#if defined(_WITHOUT_EIGEN_)
          const auto rp(filter.compute(dft.template real<num_t>(), filter.DETECT_X));
          const auto ip(filter.compute(dft.template imag<num_t>(), filter.DETECT_X));
          buf2 = (idft0.template cast<complex<num_t> >() * (rp.template cast<complex<num_t> >() + ip.template cast<complex<num_t> >() * sqrt(complex<num_t>(- 1)))).template real<num_t>();
#else
          const auto rp(filter.compute(dft.real(), filter.DETECT_X));
          const auto ip(filter.compute(dft.imag(), filter.DETECT_X));
          buf2 = (idft0 * (rp.template cast<complex<num_t> >() + sqrt(complex<num_t>(- num_t(1))) * ip.template cast<complex<num_t> >())).real();
#endif
        } else if(strcmp(argv[3], "bump") == 0) {
#if defined(_WITHOUT_EIGEN_)
          const auto rp(filter.compute(dft.template real<num_t>(), filter.BUMP_X));
          const auto ip(filter.compute(dft.template imag<num_t>(), filter.BUMP_X));
          buf2 = (idft0.template cast<complex<num_t> >() * (rp.template cast<complex<num_t> >() + ip.template cast<complex<num_t> >() * sqrt(complex<num_t>(- 1)))).template real<num_t>();
#else
          const auto rp(filter.compute(dft.real(), filter.BUMP_X));
          const auto ip(filter.compute(dft.imag(), filter.BUMP_X));
          buf2 = (idft0 * (rp.template cast<complex<num_t> >() + sqrt(complex<num_t>(- num_t(1))) * ip.template cast<complex<num_t> >())).real();
#endif
        }
        for(int j = 0; j < buf2.rows(); j ++)
          for(int k = 0; k < buf2.cols(); k ++) {
            M = max(M, abs(buf2(j, k)));
            rbuf[i * buf.rows() * buf.rows() + k * buf.rows() + j] = buf2(j, k);
          }
      }
    }
    for(int i = 0; i < rbuf.size(); i ++)
      std::cout << rbuf[i] / M << "\r" << std::endl;
  } else {
    usage();
    return - 1;
  }
  return 0;
}

