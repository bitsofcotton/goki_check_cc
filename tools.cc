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

#if defined(_WITH_NO_FLOAT)

#include "ifloat.hh"
template <typename T> using complex = Complex<T>;
typedef SimpleFloat<uint16_t, uint32_t, 16, char> num_t;
//typedef SimpleFloat<uint32_t, uint64_t, 32, short> num_t;

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
  cout << "gokicheck collect <input.ppm> <output.ppm>" << endl;
  cout << "gokicheck light   <n_recursive> <input.ppm> <output.ppm>" << endl;
  cout << "gokicheck bump    <n_integrate> <delta_pixels> <input.ppm> <output.ppm>" << endl;
  cout << "gokicheck pextend <pixels> <input.ppm> <output.ppm>" << endl;
  cout << "gokicheck reshape <num_shape_per_color> <input_color.ppm> <input_shape.ppm> <output.ppm>" << endl;
  cout << "gokicheck div     <dummy> <input0.ppm> <input1.ppm> <output.ppm>" << endl;
  cout << "gokicheck obj     <gather_pixels> <ratio> <zratio> <thin> <input.ppm> <mask.ppm>? <output.obj>" << endl;
  cout << "gokicheck tilt    <index> <max_index> <psi> <input.ppm> <input-bump.(ppm|obj)> <output.ppm>" << endl;
  cout << "gokicheck sbox    <index> <max_index> <input.ppm> <input-bump.(ppm|obj)> <output.ppm>" << endl;
  cout << "gokicheck match0  <num_of_res_shown> <num_of_hidden_match> <vbox_dst> <vbox_src> <dst.ppm> <src.ppm> <dst-bump.(ppm|obj)> <src-bump.(ppm|obj)> (<dst-mask.ppm> <src-mask.ppm>)? <output-basename>" << endl;
  cout << "gokicheck match   <num_of_res_shown> <num_of_hidden_match> <vbox_dst> <vbox_src> <dst.ppm> <src.ppm> <dst-bump.(ppm|obj)> <src-bump.(ppm|obj)> (<dst-mask.ppm> <src-mask.ppm>)? <output-basename>" << endl;
  cout << "gokicheck matcho  <match> <nemph> <vbox_dst> <vbox_src> <dst.ppm> <src.ppm> <dst-bump.(ppm|obj)> <src-bump.(ppm|obj)> (<dst-mask.ppm> <src-mask.ppm>)? <output-basename>" << endl;
  cout << "gokicheck pose    <vbox> <thresh> <posex.txt> <posey.txt> <input.ppm> <input-bump.ppm>" << endl;
  cout << "gokicheck poso    <vbox> <thresh> <posex.txt> <posey.txt> <input.ppm> <input-bump.ppm> <num_of_res_shown> <output-base>" << endl;
  cout << "gokicheck habit   <in0.obj> <in1.obj> (<index> <max_index> <psi>)? <out.obj>" << endl;
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
     strcmp(argv[1], "light")   == 0 ||
     strcmp(argv[1], "bump")    == 0 ||
     strcmp(argv[1], "bumpi")   == 0) {
    const auto f_col( ! strcmp(argv[1], "collect") ||
                      ! strcmp(argv[1], "enlarge") );
    if((f_col && argc < 4) || ((! f_col) && argc < 5)) {
      usage();
      return 0;
    }
    const auto ratio(f_col ? 0. : std::atof(argv[2]));
    const auto inidx(f_col ? 2  : 3);
    typename simpleFile<num_t>::Mat data[3];
    if(!file.loadp2or3(data, argv[inidx]))
      return - 1;
    if(strcmp(argv[1], "collect") == 0) {
      for(int i = 0; i < 3; i ++)
        data[i] = filter.compute(data[i], filter.COLLECT_BOTH);
    } else if(strcmp(argv[1], "enlarge") == 0) {
      for(int i = 0; i < 3; i ++)
        data[i] = filter.compute(data[i], filter.ENLARGE_BOTH);
    } else if(strcmp(argv[1], "pextend") == 0) {
      filter.plen = ratio;
      typename simpleFile<num_t>::Mat xyz[3];
      redig.rgb2xyz(xyz, data);
      for(int i = 0; i < 3; i ++)
        xyz[i] = filter.compute(xyz[i], filter.EXTEND_BOTH);
      redig.xyz2rgb(data, xyz);
      for(int i = 0; i < 3; i ++)
        data[i] = filter.compute(data[i], filter.CLIP);
    } else if(strcmp(argv[1], "light") == 0) {
      filter.lrecur = ratio;
      for(int i = 0; i < 3; i ++)
        data[i] = filter.compute(data[i], filter.SHARPEN_BOTH);
    } else if(strcmp(argv[1], "bump") == 0) {
      filter.dist  = 1;
      filter.bumpi = ratio;
      data[0] = data[1] = data[2] = filter.compute(redig.rgb2l(data), filter.BUMP_BOTH);
    } else if(strcmp(argv[1], "bumpi") == 0) {
      filter.dist  = 60;
      filter.bumpi = ratio;
      data[0] = data[1] = data[2] = filter.compute(redig.rgb2l(data), filter.BUMP_BOTH);
    }
    if(strcmp(argv[1], "pextend") != 0) {
      if(strcmp(argv[1], "bump") != 0 && strcmp(argv[1], "bumpi") != 0)
        redig.autoLevel(data, 3 * 2 * (data[0].rows() + data[0].cols()));
      redig.normalize(data, num_t(1));
    }
    if(!file.savep2or3(argv[inidx + 1], data, ! true, 65535))
      return - 1;
  } else if(strcmp(argv[1], "reshape") == 0 ||
            strcmp(argv[1], "div")   == 0) {
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
    if(strcmp(argv[1], "reshape") == 0) {
      const auto datav(redig.rgb2l(datas));
      for(int i = 0; i < 3; i ++)
        datac[i] = redig.reShape(datac[i], datav, count);
    } else if(strcmp(argv[1], "div") == 0) {
      for(int i = 0; i < 3; i ++) {
        assert(datas[i].rows() == datac[i].rows() &&
               datas[i].cols() == datac[i].cols());
        datas[i] = filter.compute(datas[i], filter.BCLIP);
        datac[i] = filter.compute(datac[i], filter.BCLIP);
      }
      for(int i = 0; i < datas[0].rows(); i ++)
        for(int j = 0; j < datas[0].cols(); j ++)
          for(int k = 0; k < 3; k ++)
            datac[k](i, j) *= datas[k](i, j) / datac[k](i, j);
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
    redig.initialize(vbox);
    redig.getTileVec(data[0], points, facets);
    const auto edges(redig.getEdges(mask[0], points));
    if(edges.size())
      redig.maskVectors(points, facets, mask[0]);
    for(int i = 0; i < points.size(); i ++) {
      points[i]    *= ratio;
      points[i][2] *= zratio * sqrt(data[0].rows() * data[0].cols());
    }
    file.saveobj(points, ratio * num_t(data[0].rows()),
                         ratio * num_t(data[0].cols()),
                 facets, argv[sidx], edges, thin);
    std::cout << std::string(argv[sidx]) + std::string(".mtl") << std::endl;
    file.saveMTL(argv[7], (std::string(argv[sidx]) + std::string(".mtl")).c_str());
  } else if(strcmp(argv[1], "tilt") == 0 ||
            strcmp(argv[1], "sbox") == 0) {
    if((strcmp(argv[1], "tilt") == 0 && argc < 8) ||
       (strcmp(argv[1], "sbox") == 0 && argc < 7)) {
      usage();
      return - 1;
    }
    const auto index(std::atoi(argv[2]));
    const auto Mindex(std::atoi(argv[3]));
    num_t psi(0);
    int   ipidx(5);
    if(strcmp(argv[1], "tilt") == 0)
      psi = std::atof(argv[4]);
    else
      ipidx --;
    typename simpleFile<num_t>::Mat data[3], bump[3];
    std::vector<typename simpleFile<num_t>::Vec3>  points;
    std::vector<typename simpleFile<num_t>::Veci3> polys;
    if(!file.loadp2or3(data, argv[ipidx]))
      return - 2;
    const std::string fn(argv[ipidx + 1]);
    bool is_obj(false);
    if(fn[fn.size() - 1] == 'm') {
      if(!file.loadp2or3(bump, argv[ipidx + 1]))
        return - 2;
    } else if(fn[fn.size() - 1] == 'j') {
      if(!file.loadobj(points, polys, argv[ipidx + 1]))
        return - 2;
      is_obj = true;
    } else
      return - 2;
    typename simpleFile<num_t>::Mat tilt0;
    const auto mtilt(strcmp(argv[1], "sbox") == 0 ? match_t<num_t>() :
                     redig.tiltprep(data[0], index, Mindex, psi));
    const auto depth(strcmp(argv[1], "sbox") == 0 ?
                     num_t(index) / num_t(Mindex) :
                     - num_t(1e8));
    if(is_obj)
      tilt0 = redig.tilt(redig.makeRefMatrix(data[0], 1), redig.tiltprep(points, polys, redig.makeRefMatrix(data[0], 1), mtilt), depth);
    else
      tilt0 = redig.tilt(redig.makeRefMatrix(data[0], 1), bump[0], mtilt, depth);
    for(int j = 0; j < 3; j ++)
      data[j] = redig.pullRefMatrix(tilt0, 1, data[j]);
    if(!file.savep2or3(argv[ipidx + 2], data, ! true))
      return - 1;
  } else if(strcmp(argv[1], "matcho") == 0 ||
            strcmp(argv[1], "match") == 0 ||
            strcmp(argv[1], "match0") == 0) {
    if(argc < 10) {
      usage();
      return - 1;
    }
    int fnout(10);
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
    const auto vboxdst(std::atoi(argv[4]));
    const auto vboxsrc(std::atoi(argv[5]));
    typename simpleFile<num_t>::Mat in0[3], in1[3], bump0orig[3], bump1orig[3], mask0orig[3], mask1orig[3];
    std::vector<typename simpleFile<num_t>::Veci3> delau0, delau1;
    std::vector<typename simpleFile<num_t>::Vec3>  shape0, shape1;
    if(!file.loadp2or3(in0, argv[6]))
      return - 2;
    if(!file.loadp2or3(in1, argv[7]))
      return - 2;
    if(!file.loadp2or3(bump0orig, argv[8]))
      return - 2;
    const std::string fn(argv[9]);
    if(fn[fn.size() - 1] == 'j') {
      if(!file.loadobj(shape1, delau1, argv[9]))
        return - 2;
      bump1orig[0] = bump1orig[1] = bump1orig[2] = in1[0] * num_t(0);
    } else if(fn[fn.size() - 1] == 'm') {
      if(!file.loadp2or3(bump1orig, argv[9]))
        return - 2;
    } else {
      usage();
      return - 1;
    }
    if(11 < argc) {
      if(!file.loadp2or3(mask0orig, argv[10]))
        return - 2;
      if(!file.loadp2or3(mask1orig, argv[11]))
        return - 2;
      fnout = 12;
    } else {
      mask0orig[0] = mask0orig[1] = mask0orig[2] = in0[0] * num_t(0);
      mask1orig[0] = mask1orig[1] = mask1orig[2] = in1[0] * num_t(0);
    }
    const auto& bump0(bump0orig[0]);
    const auto& bump1(bump1orig[0]);
    const auto& mask0(mask0orig[0]);
    const auto& mask1(mask1orig[0]);
    redig.initialize(vboxdst);
    redig.getTileVec(bump0, shape0, delau0);
    redig.maskVectors(shape0, delau0, mask0);
    if(fn[fn.size() - 1] == 'm') {
      redig.initialize(vboxsrc);
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
        outs[idx] = redig.showMatch(redig.replace(in0[idx] * num_t(0),
                                      shape0, delau0), shape0, mhull0);
      file.savep2or3((outbase + std::string("-repl0.ppm")).c_str(), outs, false);
      for(int idx = 0; idx < 3; idx ++)
        outs[idx] = redig.showMatch(redig.replace(in1[idx] * num_t(0),
                                      m.transform(shape1), delau1),
                                    m.transform(shape1), mhull1);
      file.savep2or3((outbase + std::string("-repl1.ppm")).c_str(), outs, false);
      for(int idx = 0; idx < 3; idx ++)
        outs[idx] = redig.showMatch(in0[idx], shape0, mhull0);
      file.savep2or3((outbase + std::string(".ppm")).c_str(), outs, false);
      for(int i = 0; i < nemph; i ++) {
        const auto iemph(num_t(i) / num_t(nemph));
        const auto reref(redig.replace(rin1, rin0, shape1, shape0,
                                       ~ m, delau0, iemph) );
        for(int idx = 0; idx < 3; idx ++)
          outs[idx] = redig.pullRefMatrix(reref, 1, in0[idx]) *
                        (num_t(1) - iemph) +
                      redig.pullRefMatrix(reref, 1 + rin0.rows() * rin0.cols(),
                                          in1[idx]) * iemph;
        file.savep2or3((outbase + std::string("-") +
                                  std::to_string(iemph) +
                                  std::string(".ppm")).c_str(), outs, false);
        file.saveobj(redig.takeShape(shape0, shape1, m,   iemph),
                     outs[0].rows(), outs[0].cols(),
                     delau0, (outbase + std::string("-emph0-") +
                                        std::to_string(iemph) +
                                        std::string(".obj")).c_str());
        file.saveobj(redig.takeShape(shape1, shape0, ~ m, iemph),
                     outs[0].rows(), outs[0].cols(),
                     delau1, (outbase + std::string("-emph1-") +
                                        std::to_string(iemph) +
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
  } else if(strcmp(argv[1], "pmerge") == 0) {
    if(argc < 8) {
      usage();
      return - 1;
    }
    const auto vbox(std::atoi(argv[2]));
    const auto thresh(std::atof(argv[3]));
    typename simpleFile<num_t>::Mat in0[3], in1[3];
    if(!file.loadp2or3(in0, argv[4]))
      return - 2;
    if(!file.loadp2or3(in1, argv[5]))
      return - 2;
    redig.initialize(vbox);
    std::vector<typename simpleFile<num_t>::Veci3> delaux0, delauy0;
    std::vector<typename simpleFile<num_t>::Veci3> delaux1, delauy1;
    std::vector<typename simpleFile<num_t>::Vec3>  shapex0, shapey0;
    std::vector<typename simpleFile<num_t>::Vec3>  shapex1, shapey1;
    std::vector<typename simpleFile<num_t>::Vec3>  centerx0, centery0;
    std::vector<typename simpleFile<num_t>::Vec3>  centerx1, centery1;
    std::vector<std::vector<int> > attendx0, attendy0;
    std::vector<std::vector<int> > attendx1, attendy1;
    const auto bump0(redig.rgb2l(in0));
    const auto bump1(redig.rgb2l(in1));
    redig.getTileVec(bump0, shapex0, delaux0);
    redig.getBones1d(bump0, shapex0, centerx0, attendx0, thresh);
    redig.getTileVec(bump1, shapex1, delaux1);
    redig.getBones1d(bump1, shapex1, centerx1, attendx1, thresh);
    const auto bump0t(bump0.transpose());
    const auto bump1t(bump1.transpose());
    redig.getTileVec(bump0t, shapey0, delauy0);
    redig.getBones1d(bump0t, shapey0, centery0, attendy0, thresh);
    redig.getTileVec(bump1t, shapey1, delauy1);
    redig.getBones1d(bump1t, shapey1, centery1, attendy1, thresh);
    const auto centerx(redig.copyBone(centerx0, attendx0, centerx1, attendx1));
    const auto centery(redig.copyBone(centery0, attendy0, centery1, attendy1));
    std::ofstream outputx, outputy;
    outputx.open(argv[6]);
    outputy.open(argv[7]);
    if(outputx.is_open() && outputy.is_open()) {
      try {
        for(int i = 0; i < centerx.size(); i ++)
          outputx << centerx[i][0] << " " << centerx[i][1] << " " << centerx[i][2] << std::endl;
        for(int i = 0; i < centery.size(); i ++)
          outputy << centery[i][0] << " " << centery[i][1] << " " << centery[i][2] << std::endl;
      } catch(...) {
        ;
      }
    }
    outputx.close();
    outputy.close();
  } else if(strcmp(argv[1], "pose") == 0 ||
            strcmp(argv[1], "poso") == 0) {
    if(argc < (strcmp(argv[1], "poso") == 0 ? 9 : 8)) {
      usage();
      return - 1;
    }
    const auto vbox(std::atoi(argv[2]));
    const auto thresh(std::atof(argv[3]));
    std::vector<typename simpleFile<num_t>::Vec3> outcenterx, outcentery;
    if(strcmp(argv[1], "poso") == 0) {
      std::ifstream inputx, inputy;
      inputx.open(argv[4]);
      inputy.open(argv[5]);
      try {
        std::string buf;
        while(std::getline(inputx, buf)) {
          std::stringstream sbuf(buf);
          typename simpleFile<num_t>::Vec3 work(3);
          sbuf >> work[0];
          sbuf >> work[1];
          sbuf >> work[2];
          outcenterx.push_back(work);
        }
        while(std::getline(inputy, buf)) {
          std::stringstream sbuf(buf);
          typename simpleFile<num_t>::Vec3 work(3);
          sbuf >> work[0];
          sbuf >> work[1];
          sbuf >> work[2];
          outcentery.push_back(work);
        }
        std::cerr << outcenterx.size() << "x" << outcentery.size() << "bone points" << std::endl;
      } catch(...) {
        usage();
        return - 2;
      }
      inputx.close();
      inputy.close();
    }
    typename simpleFile<num_t>::Mat in[3], bump[3];
    std::vector<typename simpleFile<num_t>::Veci3> delaux, delauy;
    std::vector<typename simpleFile<num_t>::Vec3>  shapex, shapey;
    std::vector<typename simpleFile<num_t>::Vec3>  centerx, centery;
    std::vector<std::vector<int> > attendx, attendy;
    if(!file.loadp2or3(in, argv[6]))
      return - 2;
    if(!file.loadp2or3(bump, argv[7]))
      return - 2;
    redig.initialize(vbox);
    auto bump0(redig.rgb2l(bump));
    redig.getTileVec(bump0, shapex, delaux);
    redig.getBones1d(bump0, shapex, centerx, attendx, thresh);
    bump0 = bump0.transpose();
    redig.getTileVec(bump0, shapey, delauy);
    redig.getBones1d(bump0, shapey, centery, attendy, thresh);
    if(strcmp(argv[1], "pose") == 0) {
      std::ofstream outputx, outputy;
      outputx.open(argv[4]);
      outputy.open(argv[5]);
      if(outputx.is_open() && outputy.is_open()) {
        try {
          for(int i = 0; i < centerx.size(); i ++)
            outputx << centerx[i][0] << " " << centerx[i][1] << " " << centerx[i][2] << std::endl;
          for(int i = 0; i < centery.size(); i ++)
            outputy << centery[i][0] << " " << centery[i][1] << " " << centery[i][2] << std::endl;
        } catch(...) {
          ;
        }
      }
      outputx.close();
      outputy.close();
    } else {
      assert(centery.size() == outcentery.size());
      for(int i = 0; i < shapey.size(); i ++) {
        const auto buf(shapey[i]);
        shapey[i][0] = buf[1];
        shapey[i][1] = buf[0];
      }
      for(int i = 0; i < centery.size(); i ++) {
        const auto buf(centery[i]);
        const auto obuf(outcentery[i]);
        centery[i][0] = buf[1];
        centery[i][1] = buf[0];
        outcentery[i][0] = obuf[1];
        outcentery[i][1] = obuf[0];
      }
      // N.B. bump0 is transposed, attendx for shapey.
      const auto Mx(bump0.rows() / vbox + 1);
      const auto My(bump0.cols() / vbox + 1);
      for(int i = 0; i < attendx.size(); i ++)
        for(int j = 0; j < attendx[i].size(); j ++)
          attendx[i][j] = (attendx[i][j] % Mx) * My + (attendx[i][j] / Mx);
      const auto rin0(redig.makeRefMatrix(in[0], 1));
      for(int j = 0; j < std::atoi(argv[8]); j ++) {
        typename simpleFile<num_t>::Mat out[3];
        const auto iemph(num_t(j + 1) / num_t(std::atoi(argv[8])));
        const auto reref(redig.draw(rin0, shapey,
            redig.takeShape(
              redig.takeShape(shapey, centery, outcentery, attendy, iemph),
              centerx, outcenterx, attendx, iemph), delauy));
        for(int idx = 0; idx < 3; idx ++)
          out[idx] = redig.pullRefMatrix(reref, 1, in[idx]);
        if(!file.savep2or3((std::string(argv[9]) + std::to_string(j) + std::string(".ppm")).c_str(), out, ! true))
          return - 1;
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
        } else if(strcmp(argv[2], "light") == 0) {
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

