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
//typedef SimpleFloat<uint32_t, uint64_t, 32, int16_t> num_t;
typedef SimpleFloat<uint32_t, uint64_t, 32, Signed<DUInt<uint64_t, 64>, 128> > num_t;

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

typedef long double num_t;

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

#include <omp.h>
namespace goki {
#include "p0.hh"
#include "p1.hh"
#include "decompose.hh"
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
  cout << "gokicheck (collect|sharpen|bump|enlarge|pextend) <input.ppm> <output.ppm>" << endl;
  cout << "gokicheck ppred <vbox> <thresh> <zratio> <num_of_emph> <outbase> <input0.ppm> <input0-bump.ppm> ..." << endl;
  cout << "gokicheck pred  <output.ppm> <input0.ppm> ..." << endl;
  cout << "gokicheck obj   <gather_pixels> <ratio> <zratio> <thin> <input.ppm> <mask.ppm>? <output.obj>" << endl;
  cout << "gokicheck (tilt|sbox)    <index> <max_index> (<psi>|<zratio>) <input.ppm> <input-bump.(ppm|obj)> <output.ppm>" << endl;
  cout << "gokicheck (match0|match) <num_of_res_shown> <num_of_hidden_match> <vbox_dst> <vbox_src> <zratio> <dst.ppm> <src.ppm> <dst-bump.(ppm|obj)> <src-bump.(ppm|obj)> (<dst-mask.ppm> <src-mask.ppm>)? <output-basename>" << endl;
  cout << "gokicheck matcho  <match> <nemph> <vbox_dst> <vbox_src> <zratio> <dst.ppm> <src.ppm> <dst-bump.(ppm|obj)> <src-bump.(ppm|obj)> (<dst-mask.ppm> <src-mask.ppm>)? <output-basename>" << endl;
  cout << "gokicheck (rmatch0|rmatch) <num_of_res_shown> <num_of_sub_match> <vbox_dst> <vbox_src> <zratio> <dst.ppm> <src.ppm> <dst-bump.(ppm|obj)> <src-bump.(ppm|obj)> (<dst-mask.ppm> <src-mask.ppm>)? <output-basename>" << endl;
  cout << "gokicheck habit   <in0.obj> <in1.obj> (<index> <max_index> <psi>)? <out.obj>" << endl;
  cout << "gokicheck reshape <num_shape_per_color> <input_color.ppm> <input_shape.ppm> <output.ppm>" << endl;
  cout << "gokicheck recolor <num_shape_per_color> <input_color.ppm> <input_shape.ppm> <output.ppm>" << endl;
  cout << "gokicheck recolor2 <num_shape_per_color> <input_color.ppm> <intensity> <output.ppm>" << endl;
  cout << "gokicheck retrace <num_shape_per_point> <inputdst.ppm> <inputsrc.ppm> <output.ppm> <intensity>" << endl;
  cout << "gokicheck retrace2 <num_shape_per_point> <inputdst.ppm> <output.ppm> <intensity>" << endl;
  cout << "gokicheck reimage <num_shape_per_point> <inputdst.ppm> <inputsrc.ppm> <output.ppm> <intensity>" << endl;
  cout << "gokicheck reimage2 <num_shape_per_point> <inputdst.ppm> <output.ppm> <intensity>" << endl;
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
  P0<num_t>         p;
  if(strcmp(argv[1], "collect") == 0 ||
     strcmp(argv[1], "enlarge") == 0 ||
     strcmp(argv[1], "pextend") == 0 ||
     strcmp(argv[1], "sharpen") == 0 ||
     strcmp(argv[1], "bump")    == 0 ||
     strcmp(argv[1], "w2b")     == 0 ||
     strcmp(argv[1], "b2w")     == 0 ||
     strcmp(argv[1], "b2wd")    == 0) {
    if(argc < 4) {
      usage();
      return 0;
    }
    if(4 < argc && strcmp(argv[1], "b2wd"))
      filter = Filter<num_t>(std::atoi(argv[4]));
    const auto rot(5 < argc ? std::atoi(argv[5]) : 1);
    typename simpleFile<num_t>::Mat data[3];
    if(!file.loadp2or3(data, argv[2]))
      return - 1;
    if(strcmp(argv[1], "collect") == 0)
      for(int i = 0; i < 3; i ++)
        data[i] = filter.rotcompute(data[i], filter.COLLECT_BOTH, rot);
    else if(strcmp(argv[1], "enlarge") == 0)
      for(int i = 0; i < 3; i ++)
        data[i] = filter.rotcompute(data[i], filter.ENLARGE_BOTH, rot);
    else if(strcmp(argv[1], "pextend") == 0)
      for(int i = 0; i < 3; i ++)
        data[i] = filter.compute(filter.rotcompute(data[i], filter.EXTEND_BOTH, rot), filter.CLIP);
    else if(strcmp(argv[1], "sharpen") == 0)
      for(int i = 0; i < 3; i ++)
        data[i] = filter.compute(filter.rotcompute(data[i], filter.SHARPEN_BOTH, rot), filter.CLIP);
    else if(strcmp(argv[1], "bump") == 0)
      data[0] = data[1] = data[2] = redig.autoLevel(filter.rotcompute(redig.rgb2d(data), filter.BUMP_BOTH, rot), 4 * (data[0].rows() + data[0].cols()));
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
      simpleFile<num_t>::Mat ddata[3];
      if(!file.loadp2or3(ddata, argv[4]))
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
    }
    if(strcmp(argv[1], "sharpen") != 0 && strcmp(argv[1], "b2w") != 0 &&
       strcmp(argv[1], "b2wd") != 0)
      redig.normalize(data, num_t(1));
    if(!file.savep2or3(argv[3], data, ! true, strcmp(argv[1], "pextend") == 0 ? 255 : 65535))
      return - 1;
  } else if(strcmp(argv[1], "reshape") == 0 ||
            strcmp(argv[1], "recolor") == 0 ||
            strcmp(argv[1], "recolor2") == 0) {
    if(argc < 6) {
      usage();
      return 0;
    }
    const auto count(std::atoi(argv[2]));
    typename simpleFile<num_t>::Mat datac[3], datas[3];
    if(!file.loadp2or3(datac, argv[3]))
      return - 1;
    if((strcmp(argv[1], "recolor") == 0 ||
        strcmp(argv[1], "reshape") == 0) &&
       ! file.loadp2or3(datas, argv[4]))
      return - 1;
    if(strcmp(argv[1], "reshape") == 0) {
      const auto datav(redig.rgb2d(datas));
      for(int i = 0; i < 3; i ++)
        datac[i] = redig.reShape(datac[i], datav, count);
    } else {
      typename simpleFile<num_t>::Mat xyzc[3], xyzs[3];
      redig.rgb2xyz(xyzc, datac);
      redig.rgb2xyz(xyzs, datas);
      for(int i = 0; i < 3; i ++)
        xyzc[i] = strcmp(argv[1], "recolor") == 0 ?
          redig.reColor(xyzc[i], xyzs[i], count) :
          redig.reColor(xyzc[i], count, num_t(std::atof(argv[4])));
      redig.xyz2rgb(datac, xyzc);
    }
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
    file.saveMTL(argv[sidx], (std::string(argv[sidx]) + std::string(".mtl")).c_str());
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
                     - num_t(1000000));
    if(is_obj)
      tilt0 = redig.tilt(redig.makeRefMatrix(data[0], 1), redig.tiltprep(points, polys, redig.makeRefMatrix(data[0], 1), mtilt), depth);
    else
      tilt0 = redig.tilt(redig.makeRefMatrix(data[0], 1), bump[0], mtilt, depth);
    for(int j = 0; j < 3; j ++)
      data[j] = redig.pullRefMatrix(tilt0, 1, data[j]);
    if(!file.savep2or3(argv[7], data, ! true))
      return - 1;
  } else if(strcmp(argv[1], "matcho") == 0 ||
            strcmp(argv[1], "match")  == 0 ||
            strcmp(argv[1], "match0") == 0 ||
            strcmp(argv[1], "rmatch") == 0 ||
            strcmp(argv[1], "rmatch0") == 0) {
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
      nhid = std::atoi(argv[3]);
    } else if(strcmp(argv[1], "rmatch") == 0 ||
              strcmp(argv[1], "rmatch0") == 0) {
      nemph = std::atoi(argv[2]);
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
    if(strcmp(argv[1], "rmatch") == 0 || strcmp(argv[1], "rmatch0") == 0 ||
       strcmp(argv[1], "matcho") == 0) {
      std::vector<match_t<num_t> > mm;
      std::vector<std::vector<typename simpleFile<num_t>::Veci3> > sdelau0, sdelau1;
      std::vector<std::vector<typename simpleFile<num_t>::Vec3>  > sshape0, sshape1;
      if(strcmp(argv[1], "rmatch")  == 0 || strcmp(argv[1], "rmatch0") == 0) {
        matchPartial<num_t> statmatch;
        sdelau0.emplace_back(delau0);
        sdelau1.emplace_back(delau1);
        sshape0.emplace_back(shape0);
        sshape1.emplace_back(shape1);
        for(int i = 1; i < nhid; i ++) {
          const auto mmm(statmatch.match(sshape0[i - 1], sshape1[i - 1], strcmp(argv[1], "rmatch0") == 0));
          // const auto mmm(statmatch.elim(statmatch.match(sshape0[i - 1], sshape1[i - 1], strcmp(argv[1], "rmatch0") == 0), in0, in1, bump1, sshape1[i - 1]));
/*
          auto mmm(statmatch.match(sshape0[i - 1], sshape1[i - 1], strcmp(argv[1], "rmatch0") == 0));
          if(nhid < mmm.size()) mmm.resize(nhid);
          mmm = statmatch.elim(mmm, in0, in1, bump1, sshape1[i - 1]);
*/
          if(! mmm.size()) break;
          mm.emplace_back(mmm[0]);
          auto dstsort(mm[i - 1].dstpoints);
          auto srcsort(mm[i - 1].srcpoints);
          std::sort(dstsort.begin(), dstsort.end());
          std::sort(srcsort.begin(), srcsort.end());
          sdelau0.emplace_back(std::vector<typename simpleFile<num_t>::Veci3>());
          sdelau1.emplace_back(std::vector<typename simpleFile<num_t>::Veci3>());
          sshape0.emplace_back(std::vector<typename simpleFile<num_t>::Vec3>());
          sshape1.emplace_back(std::vector<typename simpleFile<num_t>::Vec3>());
          std::vector<int> revdst;
          std::vector<int> revsrc;
          for(int j = 0; j < sshape0[i - 1].size(); j ++) {
            if(std::binary_search(dstsort.begin(), dstsort.end(), j)) {
              revdst.emplace_back(- 1);
              continue;
            }
            sshape0[i].emplace_back(sshape0[i - 1][j]);
            revdst.emplace_back(sshape0[i].size() - 1);
          }
          for(int j = 0; j < sshape1[i - 1].size(); j ++) {
            if(std::binary_search(srcsort.begin(), srcsort.end(), j)) {
              revsrc.emplace_back(- 1);
              continue;
            }
            sshape1[i].emplace_back(sshape1[i - 1][j]);
            revsrc.emplace_back(sshape1[i].size() - 1);
          }
          for(int j = 0; j < sdelau0[i - 1].size(); j ++) {
            auto sd0(sdelau0[i - 1][j]);
            for(int k = 0; k < sd0.size(); k ++) {
              sd0[k] = revdst[sd0[k]];
              if(sd0[k] < 0) goto nsd0;
            }
            sdelau0[i].emplace_back(sd0);
           nsd0:
            ;
          }
          for(int j = 0; j < sdelau1[i - 1].size(); j ++) {
            auto sd1(sdelau1[i - 1][j]);
            for(int k = 0; k < sd1.size(); k ++) {
              sd1[k] = revsrc[sd1[k]];
              if(sd1[k] < 0) goto nsd1;
            }
            sdelau1[i].emplace_back(sd1);
           nsd1:
            ;
          }
          if(! sdelau0[i].size() || ! sdelau1[i].size() || ! sshape0[i].size() || ! sshape1[i].size())
            break;
        }
      } else {
        mm.emplace_back(m);
        sshape0.emplace_back(shape0);
        sshape1.emplace_back(shape1);
        sdelau0.emplace_back(delau0);
        sdelau1.emplace_back(delau1);
      }
      typename simpleFile<num_t>::Mat outs[3];
      const auto rin0(redig.makeRefMatrix(in0[0], 1));
      const auto rin1(redig.makeRefMatrix(in1[0], 1 + rin0.rows() * rin0.cols()));
      std::vector<std::vector<typename simpleFile<num_t>::Veci3> > mhull0;
      std::vector<std::vector<typename simpleFile<num_t>::Veci3> > mhull1;
      for(int i = 0; i < mm.size(); i ++) {
        mhull0.emplace_back(redig.mesh2(sshape0[i], mm[i].dstpoints));
        mhull1.emplace_back((~ mm[i]).hullConv(mhull0[i]));
      }
      const std::string outbase(argv[fnout]);
      for(int idx = 0; idx < 3; idx ++) {
        for(int i = 0; i < mm.size(); i ++)
          if(i)
            outs[idx] += redig.showMatch(redig.draw(in0[idx] * num_t(0),
                                         sshape0[i], sdelau0[i]),
                           sshape0[i], sdelau0[i]);
          else
            outs[idx]  = redig.showMatch(redig.draw(in0[idx] * num_t(0),
                                         sshape0[i], sdelau0[i]),
                           sshape0[i], sdelau0[i]);
      }
      redig.normalize(outs, 1.);
      file.savep2or3((outbase + std::string("-repl0.ppm")).c_str(), outs, false);
      for(int idx = 0; idx < 3; idx ++) {
        for(int i = 0; i < mm.size(); i ++)
          if(i)
            outs[idx] += redig.showMatch(redig.draw(in1[idx] * num_t(0),
                                       mm[i].transform(sshape1[i]), sdelau1[i]),
                                       mm[i].transform(sshape1[i]), mhull1[i]);
          else
            outs[idx]  = redig.showMatch(redig.draw(in1[idx] * num_t(0),
                                       mm[i].transform(sshape1[i]), sdelau1[i]),
                                       mm[i].transform(sshape1[i]), mhull1[i]);
      }
      redig.normalize(outs, 1.);
      file.savep2or3((outbase + std::string("-repl1.ppm")).c_str(), outs, false);
      for(int idx = 0; idx < 3; idx ++) {
        for(int i = 0; i < mm.size(); i ++)
          if(i)
            outs[idx] += redig.showMatch(in0[idx], sshape0[i], sdelau0[i]);
          else
            outs[idx]  = redig.showMatch(in0[idx], sshape0[i], sdelau0[i]);
      }
      redig.normalize(outs, 1.);
      file.savep2or3((outbase + std::string(".ppm")).c_str(), outs, false);
      for(int i = 0; i < nemph; i ++) {
        const auto iemph(num_t(i) / num_t(nemph));
        simpleFile<num_t>::Mat reref;
        for(int i = 0; i < mm.size(); i ++) {
          const auto rd(redig.draw(rin1, sshape1[i],
                          redig.takeShape(sshape1[i], sshape0[i],
                            ~ mm[i], iemph), sdelau1[i]));
          if(i)
            for(int j = 0; j < reref.rows(); j ++)
              for(int k = 0; k < reref.cols(); k ++) {
                if(rd(j, k) != num_t(0)) reref(j, k) = rd(j, k);
              }
          else
            reref = rd;
        }
        for(int idx = 0; idx < 3; idx ++)
          outs[idx] = redig.pullRefMatrix(reref, 1 + rin0.rows() * rin0.cols(), in1[idx]);
          // outs[idx] = redig.draw(in1[idx] * num_t(0), shape, delau1);
        file.savep2or3((outbase + std::string("-") +
                                  std::to_string(i) +
                                  std::string("-") +
                                  std::to_string(nemph) +
                                  std::string(".ppm")).c_str(), outs, false);
/*
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
*/
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
    const auto& comp(p.next(in.size()));
    for(int y = 0; y < out[0].rows(); y ++) {
      for(int x = 0; x < out[0].cols(); x ++) {
        out[0](y, x) = out[1](y, x) = out[2](y, x) = num_t(0);
        for(int k = 0; k < in.size(); k ++)
          for(int kk = 0; kk < 3; kk ++)
            out[kk](y, x) += in[k][kk](y, x) * comp[k];
      }
    }
    redig.normalize(out, 1.);
    file.savep2or3((std::string(argv[2])).c_str(), out, ! true);
  } else if(strcmp(argv[1], "ppred") == 0 ||
            strcmp(argv[1], "ppredr") == 0) {
    if(argc < 11 || (argc & 1)) {
      usage();
      return - 1;
    }
    const auto  vbox(std::atoi(argv[2]));
    const auto  thresh(std::atof(argv[3]));
    const num_t zratio(std::atof(argv[4]));
    std::vector<std::vector<typename simpleFile<num_t>::Mat> > in;
    std::vector<typename simpleFile<num_t>::Mat> inb;
    in.resize((argc - 8) / 2);
    inb.resize((argc - 8) / 2);
    int i;
    for(i = 8; i < argc; ) {
      typename simpleFile<num_t>::Mat ibuf[3];
      if(!file.loadp2or3(ibuf, argv[i]))
        exit(- 2);
      const auto ii((i - 8) / 2);
      in[ii].resize(3);
      in[ii][0] = const_cast<typename simpleFile<num_t>::Mat &&>(ibuf[0]);
      in[ii][1] = const_cast<typename simpleFile<num_t>::Mat &&>(ibuf[1]);
      in[ii][2] = const_cast<typename simpleFile<num_t>::Mat &&>(ibuf[2]);
      i ++;
      if(!file.loadp2or3(ibuf, argv[i]))
        exit(- 2);
      inb[ii] = redig.rgb2d(ibuf);
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
    std::vector<typename simpleFile<num_t>::Vec3> outcenter;
    std::vector<typename simpleFile<num_t>::Mat>  pout;
    typename simpleFile<num_t>::Mat out[3];
    const auto rin0(redig.makeRefMatrix(in[0][0], 1));
    {
            auto center1(center);
      const auto idx(center.size() - 1);
      std::vector<std::vector<typename simpleFile<num_t>::Vec3> > lcenter;
      std::vector<std::vector<typename simpleFile<num_t>::Mat> > lin;
      lcenter.reserve(std::atoi(argv[6]) * 2);
      lin.reserve(std::atoi(argv[6]) * 2);
      for(int i = max(int(in.size()) - std::atoi(argv[6]), 0), ii = 0;
              i < in.size(); i ++, ii ++) {
        lcenter.emplace_back(i == idx ? center[i] : redig.copyBone(center[idx], centerr[idx], center[i], centerr[i]));
        lin.emplace_back(in[i]);
        assert(lcenter[ii].size() == center[idx].size());
        std::cerr << "." << std::flush;
      }
      redig.complement(pout, outcenter, lin, lcenter, attend[attend.size() - 1],
                       redig.getReverseLookup(attend[attend.size() - 1],
                         in[in.size() - 1][0]),
                       num_t(lcenter.size()));
      in.emplace_back(pout);
      center.emplace_back(outcenter);
      centerr.emplace_back(centerr[centerr.size() - 1]);
    }
    const auto center0(center);
    for(int idx = 0; idx < center.size() - 1; idx ++) {
      std::vector<std::vector<typename simpleFile<num_t>::Vec3> > lcenter;
      std::vector<std::vector<typename simpleFile<num_t>::Mat> > lin;
      lcenter.reserve(std::atoi(argv[6]) * 2);
      lin.reserve(std::atoi(argv[6]) * 2);
      for(int i = max(idx - std::atoi(argv[6]), 0), ii = 0;
              i < min(idx + std::atoi(argv[6]) + 1, int(in.size()));
              i ++, ii ++) {
        lcenter.emplace_back(i == idx ? center0[i] : redig.copyBone(center0[idx], centerr[idx], center0[i], centerr[i]));
        lin.emplace_back(in[i]);
        assert(lcenter[ii].size() == center0[idx].size());
        std::cerr << "." << std::flush;
      }
      const auto a2xy(redig.getReverseLookup(attend[idx], in[idx][0]));
      std::vector<typename simpleFile<num_t>::Vec3> leftc;
      std::vector<typename simpleFile<num_t>::Vec3> rightc;
      redig.complement(pout, leftc, lin, lcenter, attend[idx], a2xy,
        num_t(idx - max(0, idx - std::atoi(argv[6]))) +
          num_t(1) / num_t(2 * std::atoi(argv[5])) );
      redig.complement(pout, rightc, lin, lcenter, attend[idx], a2xy,
        num_t(idx - max(0, idx - std::atoi(argv[6]))) +
          num_t(std::atoi(argv[5]) * 2 - 1) / num_t(2 * std::atoi(argv[5])) );
      for(int i = - std::atoi(argv[5]); i < std::atoi(argv[5]); i ++) {
        redig.complement(pout, outcenter, lin, lcenter, attend[idx], a2xy,
          num_t(idx - max(0, idx - std::atoi(argv[6]))) +
            (num_t(i) + num_t(1) / num_t(2)) / num_t(std::atoi(argv[5])) );
        const auto newshape(redig.takeShape(shape[idx], leftc, rightc, attend[idx], num_t(i) / num_t(2 * std::atoi(argv[5]))));
        const auto reref(redig.draw(rin0, shape[idx], newshape, delau[idx]));
        for(int ii = 0; ii < 3; ii ++)
          out[ii] = filter.compute(redig.pullRefMatrix(reref, 1, pout[ii]), filter.CLIP);
        file.savep2or3((std::string(argv[7]) + std::to_string((2 * idx + 1) * std::atoi(argv[5]) + i) + std::string(".ppm")).c_str(), out, ! true);
        file.saveobj(newshape, out[0].rows(), out[0].cols(), delau[0],
                     (std::string(argv[7]) + std::to_string((2 * idx + 1) * std::atoi(argv[5]) + i) + std::string(".obj")).c_str());
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
      file.saveobj(redig.takeShape(pdst, psrc, m, num_t(1) / num_t(2)),
                   My, Mx, poldst, argv[7]);
    } else {
      matchPartial<num_t> statmatch;
      const auto m(statmatch.match(pdst, psrc)[0]);
      file.saveobj(redig.takeShape(pdst, psrc, m, num_t(1) / num_t(2)),
                   My, Mx, poldst, argv[4]);
    }
  } else if(strcmp(argv[1], "retrace") == 0 ||
            strcmp(argv[1], "reimage") == 0) {
    if(argc < 7) {
      usage();
      return - 1;
    }
    redig.initialize(1);
    typename simpleFile<num_t>::Mat dst[3], src[3], out[3];
    if(!file.loadp2or3(dst, argv[3]))
      exit(- 2);
    if(!file.loadp2or3(src, argv[4]))
      exit(- 2);
    if(strcmp(argv[1], "retrace") == 0)
      out[0] = out[1] = out[2] =
        redig.reTrace(redig.normalize(redig.rgb2d(dst), num_t(1)),
          redig.normalize(redig.rgb2d(src), num_t(1)),
          num_t(std::atof(argv[6])), std::atoi(argv[2]));
    else
      for(int i = 0; i < 3; i ++)
        out[i] = redig.reImage(dst[i], src[i],
          num_t(std::atof(argv[6])), std::atoi(argv[2]));
    redig.normalize(out, 1.);
    if(!file.savep2or3(argv[5], out, ! true, 255))
      return - 3;
  } else if(strcmp(argv[1], "retrace2") == 0 ||
            strcmp(argv[1], "reimage2") == 0) {
    if(argc < 6) {
      usage();
      return - 1;
    }
    redig.initialize(1);
    typename simpleFile<num_t>::Mat dst[3], out[3];
    if(!file.loadp2or3(dst, argv[3]))
      exit(- 2);
    if(strcmp(argv[1], "retrace2") == 0)
      out[0] = out[1] = out[2] =
        redig.reTrace(redig.normalize(redig.rgb2d(dst), num_t(1)),
          num_t(std::atof(argv[5])), std::atoi(argv[2]));
    else {
      for(int i = 0; i < 3; i ++)
        out[i] = redig.reImage(dst[i],
          num_t(std::atof(argv[5])), std::atoi(argv[2]));
      redig.autoLevel(out, 4 * (out[0].rows() + out[0].cols()));
    }
    redig.normalize(out, 1.);
    if(!file.savep2or3(argv[4], out, ! true, 255))
      return - 3;
  } else if(strcmp(argv[1], "omake") == 0) {
    std::vector<std::vector<num_t> > data;
    std::string header;
    file.loaddat(argv[3], header, data);
    simpleFile<num_t>::Mat buf(std::atoi(argv[4]), std::atoi(argv[4]));
    const auto& dft(p.seed(buf.rows()));
    const auto& idft(p.seed(- buf.rows()));
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

