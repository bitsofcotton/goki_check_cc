#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <assert.h>

#if defined(_WITH_NO_FLOAT_)

#include "ifloat.hh"
template <typename T> using complex = Complex<T>;
//typedef SimpleFloat<uint16_t, uint32_t, 16, char> num_t;
typedef SimpleFloat<uint32_t, uint64_t, 32, short> num_t;

#else

using std::sqrt;
using std::exp;
using std::log;
using std::pow;
using std::sin;
using std::cos;
using std::tan;
using std::atan2;
using std::ceil;

#include <complex>
using std::complex;

#  if defined(_WITH_MPFR_)

#include <mpreal.h>
typedef mpfr::mpreal num_t;
using std::sqrt;
using mpfr::pow;
using mpfr::log;
using mpfr::isfinite;

#  else

typedef double num_t;

#  endif
#endif

#if defined(_WITHOUT_EIGEN_)
#include "simplelin.hh"
#else
#include "simplelin.hh"
#include <Eigen/Core>
#include <Eigen/LU>
#endif

using std::max;
using std::min;
using std::pair;
using std::make_pair;
using std::vector;

namespace goki {
#include "p0.hh"
#include "p1.hh"
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
  cout << "gokicheck (collect|sharpen|bump|illust|enlarge|pextend) <input.ppm> <output.ppm>" << endl;
  cout << "gokicheck pcopy <vbox> <thresh> <zratio> <num_of_emph> <outbase> <inputdst.ppm> <inputdst-bump.ppm> <inputsrc.ppm> <inputsrc-bump.ppm>" << endl;
  cout << "gokicheck ppred <vbox> <thresh> <zratio> <num_of_emph> <outbase> <input0.ppm> <input0-bump.ppm> ..." << endl;
  cout << "gokicheck pred  <output.ppm> <input0.ppm> ..." << endl;
  cout << "gokicheck obj   <gather_pixels> <ratio> <zratio> <thin> <input.ppm> <mask.ppm>? <output.obj>" << endl;
  cout << "gokicheck (tilt|sbox)    <index> <max_index> (<psi>|<zratio>) <input.ppm> <input-bump.(ppm|obj)> <output.ppm>" << endl;
  cout << "gokicheck (match0|match) <num_of_res_shown> <num_of_hidden_match> <vbox_dst> <vbox_src> <zratio> <dst.ppm> <src.ppm> <dst-bump.(ppm|obj)> <src-bump.(ppm|obj)> (<dst-mask.ppm> <src-mask.ppm>)? <output-basename>" << endl;
  cout << "gokicheck matcho  <match> <nemph> <vbox_dst> <vbox_src> <zratio> <dst.ppm> <src.ppm> <dst-bump.(ppm|obj)> <src-bump.(ppm|obj)> (<dst-mask.ppm> <src-mask.ppm>)? <output-basename>" << endl;
  cout << "gokicheck habit   <in0.obj> <in1.obj> (<index> <max_index> <psi>)? <out.obj>" << endl;
  cout << "gokicheck reshape <num_shape_per_color> <input_color.ppm> <input_shape.ppm> <output.ppm>" << endl;
  return;
}

int main(int argc, const char* argv[]) {
  if(argc < 2) {
    usage();
    return 0;
  }
#if defined(_WITH_MPFR_)
  num_t::set_default_prec(_WITH_MPFR_);
#endif
  simpleFile<num_t> file;
  reDig<num_t>      redig;
  Filter<num_t>     filter;
  if(strcmp(argv[1], "collect") == 0 ||
     strcmp(argv[1], "enlarge") == 0 ||
     strcmp(argv[1], "pextend") == 0 ||
     strcmp(argv[1], "sharpen") == 0 ||
     strcmp(argv[1], "bump")    == 0 ||
     strcmp(argv[1], "illust")  == 0 ||
     strcmp(argv[1], "b2w")     == 0 ||
     strcmp(argv[1], "b2wd")    == 0) {
    if(argc < 4) {
      usage();
      return 0;
    }
    if(4 < argc)
      filter = Filter<num_t>(std::atoi(argv[4]));
    typename simpleFile<num_t>::Mat data[3];
    if(!file.loadp2or3(data, argv[2]))
      return - 1;
    if(strcmp(argv[1], "collect") == 0)
      for(int i = 0; i < 3; i ++)
        data[i] = filter.compute(data[i], filter.COLLECT_BOTH);
    else if(strcmp(argv[1], "enlarge") == 0)
      for(int i = 0; i < 3; i ++)
        data[i] = filter.compute(data[i], filter.ENLARGE_BOTH);
    else if(strcmp(argv[1], "pextend") == 0)
      for(int i = 0; i < 3; i ++)
        data[i] = filter.compute(data[i], filter.EXTEND_BOTH);
    else if(strcmp(argv[1], "sharpen") == 0)
      for(int i = 0; i < 3; i ++)
        data[i] = filter.compute(filter.compute(data[i], filter.SHARPEN_BOTH), filter.CLIP);
    else if(strcmp(argv[1], "bump") == 0)
      data[0] = data[1] = data[2] = redig.autoLevel(filter.compute(filter.compute(redig.rgb2l(data), filter.BUMP_BOTH), filter.INTEG_BOTH), 4 * (data[0].rows() + data[0].cols()));
    else if(strcmp(argv[1], "illust") == 0) {
      const auto datav(redig.rgb2l(data));
      const auto bump(redig.autoLevel(filter.compute(filter.compute(datav, filter.BUMP_BOTH), filter.INTEG_BOTH), 4 * (data[0].rows() + data[0].cols())));
      const auto region(redig.reShape(bump, datav, 256));
      typename simpleFile<num_t>::Mat checked(datav.rows(), datav.cols());
      for(int i = 0; i < checked.rows(); i ++)
        for(int j = 0; j < checked.cols(); j ++)
          checked(i, j) = false;
      for(int loop = 0; ; loop ++) {
        std::cerr << "l" << std::flush;
        num_t test(- 1);
        auto  mask(region);
        for(int i = 0; i < mask.rows(); i ++)
          for(int j = 0; j < mask.cols(); j ++)
            if(! checked(i, j)) {
              test = mask(i, j);
              break;
            }
        if(test < num_t(0))
          break;
        for(int i = 0; i < mask.rows(); i ++)
          for(int j = 0; j < mask.cols(); j ++)
            mask(i, j) = num_t(mask(i, j) == test && ! checked(i, j) ? 0 : 1);
        for(int i = 0; i < datav.rows(); i ++)
          for(int j = 0; j < datav.cols(); j ++) 
            if(! checked(i, j) && mask(i, j) < num_t(1) / num_t(2)) {
              std::vector<std::pair<int, int> > store;
              auto wchk(checked);
              redig.floodfill(wchk, store, mask, i, j);
              std::sort(store.begin(), store.end());
              num_t avgr(0);
              num_t avgg(0);
              num_t avgb(0);
              int   cnt(0);
              for(int ii = 0; ii < store.size(); ii ++) {
                if(store[ii].first < 0 || datav.rows() <= store[ii].first)
                  continue;
                for(int jj = store[ii].second;
                        jj < (ii + 1 < store.size() &&
                          store[ii].first == store[ii + 1].first ?
                          store[ii + 1].second : store[ii].second + 1);
                        jj ++) {
                  if(jj < 0 || datav.cols() <= jj)
                    continue;
                  avgr += data[0](store[ii].first, jj); 
                  avgg += data[1](store[ii].first, jj); 
                  avgb += data[2](store[ii].first, jj); 
                  cnt ++;
                }
                if(ii + 1 < store.size() &&
                   store[ii].first == store[ii + 1].first) ii ++;
              }
              if(! cnt) {
                checked(i, j) = true;
                continue;
              }
              avgr /= num_t(cnt);
              avgg /= num_t(cnt);
              avgb /= num_t(cnt);
              for(int ii = 0; ii < store.size(); ii ++) {
                if(store[ii].first < 0 || datav.rows() <= store[ii].first)
                  continue;
                for(int jj = store[ii].second;
                        jj < (ii + 1 < store.size() &&
                          store[ii].first == store[ii + 1].first ?
                          store[ii + 1].second : store[ii].second + 1);
                        jj ++) {
                  if(jj < 0 || datav.cols() <= jj)
                    continue;
                  data[0](store[ii].first, jj) = avgr;
                  data[1](store[ii].first, jj) = avgg;
                  data[2](store[ii].first, jj) = avgb;
                  checked(store[ii].first, jj) = true;
                }
                if(ii + 1 < store.size() &&
                   store[ii].first == store[ii + 1].first) ii ++;
              }
            }
      }
    }
    else if(strcmp(argv[1], "b2w") == 0) {
      for(int i = 0; i < data[0].rows(); i ++)
        for(int j = 0; j < data[0].cols(); j ++)
          if(data[0](i, j) == num_t(0) &&
             data[1](i, j) == num_t(0) &&
             data[1](i, j) == num_t(0))
            data[0](i, j) = data[1](i, j) = data[2](i, j) = num_t(1);
    } else if(strcmp(argv[1], "b2wd") == 0) {
      simpleFile<num_t>::Mat ddata[3];
      if(!file.loadp2or3(ddata, argv[4]))
        return - 1;
      for(int i = 0; i < data[0].rows(); i ++)
        for(int j = 0; j < data[0].cols(); j ++)
          if((data[0](i, j) == num_t(0) &&
              data[1](i, j) == num_t(0) &&
              data[1](i, j) == num_t(0)) ||
             (data[0](i, j) == ddata[0](i, j) &&
              data[1](i, j) == ddata[1](i, j) &&
              data[2](i, j) == ddata[2](i, j)) )
            data[0](i, j) = data[1](i, j) = data[2](i, j) = num_t(1);
    }
    if(strcmp(argv[1], "sharpen") != 0 && strcmp(argv[1], "b2w") != 0 &&
       strcmp(argv[1], "b2wd") != 0)
      redig.normalize(data, num_t(1));
    if(!file.savep2or3(argv[3], data, ! true, strcmp(argv[1], "pextend") == 0 ? 255 : 65535))
      return - 1;
  } else if(strcmp(argv[1], "reshape") == 0) {
    if(argc < 6) {
      usage();
      return 0;
    }
    const auto count(std::atoi(argv[2]));
    typename simpleFile<num_t>::Mat datac[3], datas[3];
    if(!file.loadp2or3(datac, argv[3]))
      return - 1;
    if(!file.loadp2or3(datas, argv[4]))
      return - 1;
    const auto datav(redig.rgb2l(datas));
    for(int i = 0; i < 3; i ++)
      datac[i] = redig.reShape(datac[i], datav, count);
    redig.normalize(datac, num_t(1));
    if(!file.savep2or3(argv[5], datac, ! true))
      return - 1;
  } else if(strcmp(argv[1], "obj") == 0) {
    typename simpleFile<num_t>::Mat data[3], mask[3];
    if(argc < 8) {
      usage();
      return - 1;
    }
    const auto  vbox(std::atoi(argv[2]));
    const num_t ratio(std::atof(argv[3]));
    const num_t zratio(std::atof(argv[4]));
    const num_t thin(std::atof(argv[5]));
    if(!file.loadp2or3(data, argv[6]))
      return - 1;
    int sidx(7);
    if(8 < argc) {
      if(!file.loadp2or3(mask, argv[7]))
        return - 1;
      sidx ++;
    } else
      mask[0] = mask[1] = mask[2] = data[0] * num_t(0);
    std::vector<typename simpleFile<num_t>::Vec3>  points;
    std::vector<typename simpleFile<num_t>::Veci3> facets;
    redig.initialize(vbox, zratio);
    redig.getTileVec(data[0], points, facets);
    const auto edges(redig.getEdges(mask[0], points));
    if(edges.size())
      redig.maskVectors(points, facets, mask[0]);
    for(int i = 0; i < points.size(); i ++)
      points[i] *= ratio;
    file.saveobj(points, ratio * num_t(data[0].rows()),
                         ratio * num_t(data[0].cols()),
                 facets, argv[sidx], edges, thin);
    file.saveMTL(argv[7], (std::string(argv[sidx]) + std::string(".mtl")).c_str());
  } else if(strcmp(argv[1], "tilt") == 0 ||
            strcmp(argv[1], "sbox") == 0) {
    if(argc < 8) {
      usage();
      return - 1;
    }
    const auto index(std::atoi(argv[2]));
    const auto Mindex(std::atoi(argv[3]));
    num_t psi(0);
    num_t zratio(1);
    if(strcmp(argv[1], "tilt") == 0)
      psi = std::atof(argv[4]);
    else
      zratio = std::atof(argv[4]);
    typename simpleFile<num_t>::Mat data[3], bump[3];
    std::vector<typename simpleFile<num_t>::Vec3>  points;
    std::vector<typename simpleFile<num_t>::Veci3> polys;
    if(!file.loadp2or3(data, argv[5]))
      return - 2;
    const std::string fn(argv[6]);
    bool is_obj(false);
    if(fn[fn.size() - 1] == 'm') {
      if(!file.loadp2or3(bump, argv[6]))
        return - 2;
    } else if(fn[fn.size() - 1] == 'j') {
      if(!file.loadobj(points, polys, argv[6]))
        return - 2;
      is_obj = true;
    } else
      return - 2;
    typename simpleFile<num_t>::Mat tilt0;
    const auto mtilt(strcmp(argv[1], "sbox") == 0 ? match_t<num_t>() :
                     redig.tiltprep(data[0], index, Mindex, psi));
    const auto depth(strcmp(argv[1], "sbox") == 0 ?
                     num_t(index) / num_t(Mindex) * zratio *
                       sqrt(num_t(data[0].rows() * data[0].cols())) :
                     - num_t(100000000));
    if(is_obj)
      tilt0 = redig.tilt(redig.makeRefMatrix(data[0], 1), redig.tiltprep(points, polys, redig.makeRefMatrix(data[0], 1), mtilt), depth);
    else
      tilt0 = redig.tilt(redig.makeRefMatrix(data[0], 1), bump[0], mtilt, depth);
    for(int j = 0; j < 3; j ++)
      data[j] = redig.pullRefMatrix(tilt0, 1, data[j]);
    if(!file.savep2or3(argv[7], data, ! true))
      return - 1;
  } else if(strcmp(argv[1], "matcho") == 0 ||
            strcmp(argv[1], "match") == 0 ||
            strcmp(argv[1], "match0") == 0) {
    if(argc < 11) {
      usage();
      return - 1;
    }
    int fnout(11);
    int nshow(0);
    int nhid(0);
    int nemph(0);
    match_t<num_t> m;
    if(strcmp(argv[1], "match") == 0 || strcmp(argv[1], "match0") == 0) {
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
      nemph = std::atoi(argv[3]);
    }
    const auto  vboxdst(std::atoi(argv[4]));
    const auto  vboxsrc(std::atoi(argv[5]));
    const num_t zratio(std::atof(argv[6]));
    typename simpleFile<num_t>::Mat in0[3], in1[3], bump0orig[3], bump1orig[3], mask0orig[3], mask1orig[3];
    std::vector<typename simpleFile<num_t>::Veci3> delau0, delau1;
    std::vector<typename simpleFile<num_t>::Vec3>  shape0, shape1;
    if(!file.loadp2or3(in0, argv[7]))
      return - 2;
    if(!file.loadp2or3(in1, argv[8]))
      return - 2;
    if(!file.loadp2or3(bump0orig, argv[9]))
      return - 2;
    const std::string fn(argv[10]);
    if(fn[fn.size() - 1] == 'j') {
      if(!file.loadobj(shape1, delau1, argv[10]))
        return - 2;
      bump1orig[0] = bump1orig[1] = bump1orig[2] = in1[0] * num_t(0);
    } else if(fn[fn.size() - 1] == 'm') {
      if(!file.loadp2or3(bump1orig, argv[10]))
        return - 2;
    } else {
      usage();
      return - 1;
    }
    if(13 < argc) {
      if(!file.loadp2or3(mask0orig, argv[11]))
        return - 2;
      if(!file.loadp2or3(mask1orig, argv[12]))
        return - 2;
      fnout = 13;
    } else {
      mask0orig[0] = mask0orig[1] = mask0orig[2] = in0[0] * num_t(0);
      mask1orig[0] = mask1orig[1] = mask1orig[2] = in1[0] * num_t(0);
    }
    const auto& bump0(bump0orig[0]);
    const auto& bump1(bump1orig[0]);
    const auto& mask0(mask0orig[0]);
    const auto& mask1(mask1orig[0]);
    redig.initialize(vboxdst, zratio);
    redig.getTileVec(bump0, shape0, delau0);
    redig.maskVectors(shape0, delau0, mask0);
    if(fn[fn.size() - 1] == 'm') {
      redig.initialize(vboxsrc, zratio);
      redig.getTileVec(bump1, shape1, delau1);
      redig.maskVectors(shape1, delau1, mask1);
    }
    if(strcmp(argv[1], "matcho") == 0) {
      typename simpleFile<num_t>::Mat outs[3];
      const auto rin0(redig.makeRefMatrix(in0[0], 1));
      const auto rin1(redig.makeRefMatrix(in1[0], 1 + rin0.rows() * rin0.cols()));
      const auto mhull0(redig.mesh2(shape0, m.dstpoints));
      const auto mhull1((~ m).hullConv(mhull0));
      const std::string outbase(argv[fnout]);
      for(int idx = 0; idx < 3; idx ++)
        outs[idx] = redig.showMatch(redig.draw(in0[idx] * num_t(0),
                                      shape0, delau0), shape0, mhull0);
      file.savep2or3((outbase + std::string("-repl0.ppm")).c_str(), outs, false);
      for(int idx = 0; idx < 3; idx ++)
        outs[idx] = redig.showMatch(redig.draw(in1[idx] * num_t(0),
                                      m.transform(shape1), delau1),
                                    m.transform(shape1), mhull1);
      file.savep2or3((outbase + std::string("-repl1.ppm")).c_str(), outs, false);
      for(int idx = 0; idx < 3; idx ++)
        outs[idx] = redig.showMatch(in0[idx], shape0, mhull0);
      file.savep2or3((outbase + std::string(".ppm")).c_str(), outs, false);
      for(int i = 0; i < nemph; i ++) {
        const auto iemph(num_t(i) / num_t(nemph));
        const auto shape(redig.takeShape(shape1, shape0, ~ m, iemph));
        const auto reref(redig.draw(rin1, shape1, shape, delau1));
        for(int idx = 0; idx < 3; idx ++)
          outs[idx] = redig.pullRefMatrix(reref, 1 + rin0.rows() * rin0.cols(), in1[idx]);
          // outs[idx] = redig.draw(in1[idx] * num_t(0), shape, delau1);
        file.savep2or3((outbase + std::string("-") +
                                  std::to_string(i) +
                                  std::string("-") +
                                  std::to_string(nemph) +
                                  std::string(".ppm")).c_str(), outs, false);
        file.saveobj(redig.takeShape(shape0, shape1, m,   iemph),
                     outs[0].rows(), outs[0].cols(), delau0,
                     (outbase + std::string("-emph0-") +
                                std::to_string(i) +
                                std::string("-") +
                                std::to_string(nemph) +
                                std::string(".obj")).c_str());
        file.saveobj(shape, outs[0].rows(), outs[0].cols(), delau1,
                     (outbase + std::string("-emph1-") +
                                std::to_string(i) +
                                std::string("-") +
                                std::to_string(nemph) +
                                std::string(".obj")).c_str());
      }
    } else { 
      matchPartial<num_t> statmatch;
      auto matches(statmatch.match(shape0, shape1, strcmp(argv[1], "match0") == 0));
      matches.resize(min(int(matches.size()), nhid));
/*
      if(fn[fn.size() - 1] == 'm')
        matches = statmatch.elim(matches, in0, in1, bump1, shape1);
*/
      std::cerr << matches.size() << "pending" << std::endl;
      for(int n = 0; n < min(int(matches.size()), nshow); n ++) {
        std::ofstream output;
        output.open((std::string(argv[fnout]) + std::to_string(n + 1) +
                     std::string(".txt")).c_str());
        if(output.is_open()) {
          try {
            output << matches[n];
          } catch(...) {
            ;
          }
        }
        output.close();
      }
    }
  } else if(strcmp(argv[1], "pred") == 0) {
    if(argc < 5) {
      usage();
      return - 1;
    }
    std::vector<std::vector<typename simpleFile<num_t>::Mat> > in;
    in.resize(argc - 3);
    for(int i = 3; i < argc; i ++) {
      typename simpleFile<num_t>::Mat ibuf[3];
      if(!file.loadp2or3(ibuf, argv[i]))
        return - 2;
      const auto ii(i - 3);
      in[ii].resize(3);
      in[ii][0] = const_cast<typename simpleFile<num_t>::Mat &&>(ibuf[0]);
      in[ii][1] = const_cast<typename simpleFile<num_t>::Mat &&>(ibuf[1]);
      in[ii][2] = const_cast<typename simpleFile<num_t>::Mat &&>(ibuf[2]);
      assert(in[ii][0].rows() == in[0][0].rows());
      assert(in[ii][0].cols() == in[0][0].cols());
    }
    const auto idx(in.size() - 1);
    typename simpleFile<num_t>::Mat out[3];
    for(int i = 0; i < 3; i ++)
      out[i].resize(in[idx][0].rows(), in[idx][0].cols());
    for(int y = 0; y < out[0].rows(); y ++) {
      for(int x = 0; x < out[0].cols(); x ++) {
        P0B<num_t> p0(in.size() / 2);
        auto p1(p0);
        auto p2(p0);
        for(int k = 1; k < in.size(); k ++) {
          out[0](y, x) = p0.next(in[k][0](y, x));
          out[1](y, x) = p1.next(in[k][1](y, x));
          out[2](y, x) = p2.next(in[k][2](y, x));
        }
      }
    }
    redig.normalize(out, 1.);
    file.savep2or3((std::string(argv[2])).c_str(), out, ! true);
  } else if(strcmp(argv[1], "ppred") == 0 ||
            strcmp(argv[1], "pcopy") == 0) {
    if(argc < 11 || ! (argc & 1)) {
      usage();
      return - 1;
    }
    const auto  vbox(std::atoi(argv[2]));
    const auto  thresh(std::atof(argv[3]));
    const num_t zratio(std::atof(argv[4]));
    std::vector<std::vector<typename simpleFile<num_t>::Mat> > in;
    std::vector<typename simpleFile<num_t>::Mat> inb;
    in.resize((argc - 7) / 2);
    inb.resize((argc - 7) / 2);
    int i;
    for(i = 7; i < argc; ) {
      typename simpleFile<num_t>::Mat ibuf[3];
      if(!file.loadp2or3(ibuf, argv[i]))
        exit(- 2);
      const auto ii((i - 7) / 2);
      in[ii].resize(3);
      in[ii][0] = const_cast<typename simpleFile<num_t>::Mat &&>(ibuf[0]);
      in[ii][1] = const_cast<typename simpleFile<num_t>::Mat &&>(ibuf[1]);
      in[ii][2] = const_cast<typename simpleFile<num_t>::Mat &&>(ibuf[2]);
      i ++;
      if(!file.loadp2or3(ibuf, argv[i]))
        exit(- 2);
      inb[ii] = redig.rgb2l(ibuf);
      i ++;
      assert(in[ii][0].rows() == inb[ii].rows());
      assert(in[ii][0].cols() == inb[ii].cols());
      assert(in[ii][0].rows() == in[0][0].rows());
      assert(in[ii][0].cols() == in[0][0].cols());
    }
    std::cerr << "i" << std::flush;
    redig.initialize(vbox, zratio);
    std::vector<std::vector<typename simpleFile<num_t>::Veci3> > delau;
    std::vector<std::vector<typename simpleFile<num_t>::Vec3> >  shape;
    std::vector<std::vector<typename simpleFile<num_t>::Vec3> >  center;
    std::vector<std::vector<std::vector<int> > > attend;
    std::vector<std::vector<num_t> > centerr;
    delau.resize(in.size());
    shape.resize(in.size());
    center.resize(in.size());
    attend.resize(in.size());
    centerr.resize(in.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < in.size(); i ++) {
      auto lredig(redig);
      lredig.getTileVec(inb[i], shape[i], delau[i]);
      lredig.getBone(inb[i], shape[i], center[i], centerr[i], attend[i], thresh);
#if !defined(_OPENMP)
      assert(center[i].size() == centerr[i].size());
      assert(center[i].size() == attend[i].size());
#endif
      std::cerr << center[i].size() << ":" << std::flush;
    }
    const auto idx(strcmp(argv[1], "ppred") == 0 ? center.size() - 1 : 0);
    if(strcmp(argv[1], "ppred") == 0) {
      for(int i = 0; i < in.size() - 1; i ++) {
        center[i] = redig.copyBone(center[idx], centerr[idx], center[i], centerr[i]);
        assert(center[i].size() == center[idx].size());
        std::cerr << "." << std::flush;
      }
    } else
      center[1] = redig.copyBone(center[idx], centerr[idx], center[1], centerr[1]);
    std::vector<std::vector<std::pair<int, int> > > a2xy;
    a2xy.resize(attend[idx].size());
    for(int i = 0; i < attend[idx].size(); i ++) {
      const auto h(in[idx][0].rows() / vbox + 1);
      const auto w(in[idx][0].cols() / vbox + 1);
      a2xy[i].reserve(attend[idx][i].size());
      for(int j = 0; j < attend[idx][i].size(); j ++) {
        const auto& a0ij(attend[idx][i][j]);
              auto& a2xyij(a2xy[i][j]);
        a2xyij.first  = (a0ij / w) * vbox;
        a2xyij.second = (a0ij % w) * vbox;
      }
    }
    std::vector<typename simpleFile<num_t>::Vec3> outcenter;
    typename simpleFile<num_t>::Mat pin[3];
    if(strcmp(argv[1], "ppred") == 0) {
      std::cerr << "p" << std::flush;
      outcenter.resize(center[idx].size());
      for(int i = 0; i < center[idx].size(); i ++) {
        P0B<num_t> p0(center.size() / 2);
        auto p1(p0);
        auto p2(p0);
        outcenter[i] = typename simpleFile<num_t>::Vec3(3);
        for(int j = 1; j < center.size(); j ++) {
          outcenter[i][0] = p0.next(center[j][i][0]);
          outcenter[i][1] = p1.next(center[j][i][1]);
          outcenter[i][2] = p2.next(center[j][i][2]);
        }
      }
      for(int i = 0; i < 3; i ++)
        pin[i].resize(in[idx][0].rows(), in[idx][0].cols());
      for(int i = 0; i < attend[idx].size(); i ++) {
        for(int j = 0; j < attend[idx][i].size(); j ++) {
          P0B<num_t> p0(center.size() / 2);
          auto p1(p0);
          auto p2(p0);
          const auto yy(filter.getImgPt(a2xy[i][j].first,  in[idx][0].rows()));
          const auto xx(filter.getImgPt(a2xy[i][j].second, in[idx][0].cols()));
          for(int k = 1; k < center.size(); k ++) {
            const auto yf(filter.getImgPt(a2xy[i][j].first  + int(center[k][i][0] - center[idx][i][0]), in[idx][0].rows()));
            const auto xf(filter.getImgPt(a2xy[i][j].second + int(center[k][i][1] - center[idx][i][1]), in[idx][0].cols()));
            pin[0](yy, xx) = p0.next(in[k][0](yf, xf));
            pin[1](yy, xx) = p1.next(in[k][1](yf, xf));
            pin[2](yy, xx) = p2.next(in[k][2](yf, xf));
          }
          const auto n0(pin[0](yy, xx));
          const auto n1(pin[1](yy, xx));
          const auto n2(pin[2](yy, xx));
          for(int y0 = yy; y0 < min(yy + vbox, int(in[idx][0].rows()) - 1); y0 ++)
            for(int x0 = xx; x0 < min(xx + vbox, int(in[idx][0].cols()) - 1); x0 ++) {
              pin[0](y0, x0) = n0;
              pin[1](y0, x0) = n1;
              pin[2](y0, x0) = n2;
            }
        }
      }
    } else {
      outcenter = center[1];
      pin[0]    = in[idx][0];
      pin[1]    = in[idx][1];
      pin[2]    = in[idx][2];
    }
    typename simpleFile<num_t>::Mat out[3];
    const auto rin0(redig.makeRefMatrix(in[idx][0], 1));
    for(int j = 0; j < std::atoi(argv[5]); j ++) {
      const auto iemph(num_t(j + 1) / num_t(std::atoi(argv[5])));
      const auto newshape(redig.takeShape(shape[idx], center[idx], outcenter, attend[idx], iemph));
      const auto reref(redig.draw(rin0, shape[idx], newshape, delau[idx]));
      for(int ii = 0; ii < 3; ii ++)
        out[ii] = filter.compute(redig.pullRefMatrix(reref, 1, (in[idx][ii] * (num_t(1) - iemph) + pin[ii] * iemph)), filter.CLIP);
      file.savep2or3((std::string(argv[6]) + std::to_string(j) + std::string(".ppm")).c_str(), out, ! true);
      file.saveobj(newshape, out[0].rows(), out[0].cols(), delau[0],
                   (std::string(argv[6]) + std::to_string(j) + std::string(".obj")).c_str());
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
      file.saveobj(redig.takeShape(pdst, psrc, m, num_t(1) / num_t(2)),
                   My, Mx, poldst, argv[7]);
    } else {
      matchPartial<num_t> statmatch;
      const auto m(statmatch.match(pdst, psrc)[0]);
      file.saveobj(redig.takeShape(pdst, psrc, m, num_t(1) / num_t(2)),
                   My, Mx, poldst, argv[4]);
    }
  } else if(strcmp(argv[1], "omake") == 0) {
    std::vector<std::vector<num_t> > data;
    std::string header;
    file.loaddat(argv[3], header, data);
    simpleFile<num_t>::Mat buf(std::atoi(argv[4]), std::atoi(argv[4]));
    const auto dft(filter.seed(buf.rows(), false));
    const auto idft(filter.seed(buf.rows(), true));
    for(int i0 = 1; i0 < data.size(); i0 ++) {
      for(int i = 0; i <= data[i0].size() / buf.rows() / buf.rows(); i ++) {
        for(int k = 0; k < buf.cols(); k ++)
          for(int j = 0; j < buf.rows(); j ++) {
            const auto idx(i * buf.rows() * buf.rows() + k * buf.rows() + j);
            buf(j, k) = idx < data[i0].size() ? data[i0][idx] : num_t(0);
          }
        Filter<num_t>::Mat buf2;
        if(strcmp(argv[2], "diff") == 0){
#if defined(_WITHOUT_EIGEN_)
          buf2 = (idft * (
            filter.compute(dft.template real<num_t>() * buf, filter.DETECT_X).template cast<complex<num_t> >() +
            filter.compute(dft.template imag<num_t>() * buf, filter.DETECT_X).template cast<complex<num_t> >() * complex<num_t>(num_t(0), num_t(1)) ) ).template real<num_t>();
#else
          buf2 = (idft * (
            filter.compute(dft.real() * buf, filter.DETECT_X).template cast<complex<num_t> >() +
            filter.compute(dft.imag() * buf, filter.DETECT_X).template cast<complex<num_t> >() * complex<num_t>(num_t(0), num_t(1)) ) ).real();
#endif
        } else if(strcmp(argv[2], "sharpen") == 0) {
#if defined(_WITHOUT_EIGEN_)
          buf2 = (idft * (
            filter.compute(dft.template real<num_t>() * buf, filter.SHARPEN_X).template cast<complex<num_t> >() +
            filter.compute(dft.template imag<num_t>() * buf, filter.SHARPEN_X).template cast<complex<num_t> >() * complex<num_t>(num_t(0), num_t(1)) ) ).template real<num_t>();
#else
          buf2 = (idft * (
            filter.compute(dft.real() * buf, filter.SHARPEN_X).template cast<complex<num_t> >() +
            filter.compute(dft.imag() * buf, filter.SHARPEN_X).template cast<complex<num_t> >() * complex<num_t>(num_t(0), num_t(1)) ) ).real();
#endif
        } else if(strcmp(argv[2], "bump") == 0) {
#if defined(_WITHOUT_EIGEN_)
          buf2 = (idft * (
            filter.compute(dft.template real<num_t>() * buf, filter.BUMP_X).template cast<complex<num_t> >() +
            filter.compute(dft.template imag<num_t>() * buf, filter.BUMP_X).template cast<complex<num_t> >() * complex<num_t>(num_t(0), num_t(1)) ) ).template real<num_t>();
#else
          buf2 = (idft * (
            filter.compute(dft.real() * buf, filter.BUMP_X).template cast<complex<num_t> >() +
            filter.compute(dft.imag() * buf, filter.BUMP_X).template cast<complex<num_t> >() * complex<num_t>(num_t(0), num_t(1)) ) ).real();
#endif
        }
        for(int k = 0; k < buf2.cols(); k ++)
          for(int j = 0; j < buf2.rows(); j ++) {
            const auto idx(i * buf.rows() * buf.rows() + k * buf.rows() + j);
            if(idx < data[i0].size())
              data[i0][idx] = buf2(j, k);
            else
              break;
          }
      }
    }
    num_t M(0);
    for(int i = 1; i < data.size(); i ++) {
      auto sdata(data[i]);
      for(int i = 0; i < sdata.size(); i ++)
        sdata[i] = abs(sdata[i]);
      std::sort(sdata.begin(), sdata.end());
      M = max(M, sdata[sdata.size() * 7 / 8] * num_t(4));
    }
    std::cout << header;
    for(int i = 0; i < data[0].size(); i ++) {
      std::cout << data[0][i] << " ";
      for(int j = 1; j < data.size(); j ++)
        std::cout << (i < data[j].size() ? data[j][i] / M : num_t(0)) << " ";
      std::cout << std::endl;
    }
  } else {
    usage();
    return - 1;
  }
  return 0;
}

