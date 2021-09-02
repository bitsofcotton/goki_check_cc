#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cctype>
#include <random>
#include <assert.h>

//#define int int64_t
#define int int32_t
#include "lieonn.hh"
typedef myfloat num_t;
#include "p0.hh"
#include "p1.hh"
#include "decompose.hh"
#include "catg.hh"
#include "fileio.hh"
#include "enlarge.hh"
#include "match.hh"
#include "redig.hh"

using std::cout;
using std::cerr;
using std::endl;
using std::atoi;
using std::string;
using std::to_string;
using std::vector;
using std::sort;
using std::binary_search;
using std::make_pair;

#include <stdlib.h>

void usage() {
  cout << "Usage:" << endl;
  cout << "gokicheck (collect|integ|sharpen|bump|enlarge|flarge|pextend|blink|lpf|represent) <input.ppm> <output.ppm> <recur> <rot>" << endl;
  cout << "gokicheck bumpc <psi> <gather_pixels> <zratio> <color.ppm> <bump0.ppm> <output.ppm>" << endl;
  cout << "gokicheck (pred|lenl) <output.ppm> <input0.ppm> ..." << endl;
  cout << "gokicheck (cat|composite) <output.ppm> <input0.ppm> <input0-represent.ppm> ..." << endl;
  cout << "gokicheck obj   <gather_pixels> <ratio> <zratio> <thin> <input.ppm> <mask.ppm>? <output.obj>" << endl;
  cout << "gokicheck (tilt|sbox)    <index> <max_index> (<psi>|<zratio>) <input.ppm> <input-bump.(ppm|obj)> <output.ppm>" << endl;
  cout << "gokicheck match <nsub> <nemph> <vbox_dst> <vbox_src> <zratio> <dst.ppm> <src.ppm> <dst-bump.(ppm|obj)> <src-bump.(ppm|obj)> (<dst-mask.ppm> <src-mask.ppm>)? <output-basename>" << endl;
  cout << "gokicheck habit   <in0.obj> <in1.obj> (<index> <max_index> <psi>)? <out.obj>" << endl;
  cout << "gokicheck reshape <num_shape_per_color> <input_color.ppm> <input_shape.ppm> <output.ppm>" << endl;
  cout << "gokicheck recolor <num_shape_per_color> <input_color.ppm> <input_shape.ppm> <output.ppm> <intensity>" << endl;
  cout << "gokicheck recolor2 <num_shape_per_color> <input_color.ppm> <output.ppm> <intensity>" << endl;
  cout << "gokicheck recolor3 <num_shape_per_color> <input_color.ppm> <input_shape.ppm> <output.ppm>" << endl;
  cout << "gokicheck retrace <num_shape_per_point> <inputdst.ppm> <inputsrc.ppm> <output.ppm> <intensity>" << endl;
  cout << "gokicheck retrace2 <num_shape_per_point> <inputdst.ppm> <output.ppm> <intensity>" << endl;
  cout << "gokicheck newtrace <num_shape_per_point> <size> <output.ppm>" << endl;
  cout << "gokicheck reimage <num_shape_per_point> <inputdst.ppm> <inputsrc.ppm> <output.ppm> <intensity>" << endl;
  cout << "gokicheck reimage2 <num_shape_per_point> <inputdst.ppm> <output.ppm> <intensity>" << endl;
  cout << "gokicheck penl <opt.ppm> <in.ppm> <output.ppm>" << endl;
  return;
}

#undef int
int main(int argc, const char* argv[]) {
//#define int int64_t
#define int int32_t
  if(argc < 2) {
    usage();
    return 0;
  }
  simpleFile<num_t> file;
  reDig<num_t>      redig;
  if(strcmp(argv[1], "collect") == 0 ||
     strcmp(argv[1], "integ") == 0 ||
     strcmp(argv[1], "rot") == 0 ||
     strcmp(argv[1], "enlarge") == 0 ||
     strcmp(argv[1], "flarge") == 0 ||
     strcmp(argv[1], "pextend") == 0 ||
     strcmp(argv[1], "blink") == 0 ||
     strcmp(argv[1], "lpf") == 0 ||
     strcmp(argv[1], "sharpen") == 0 ||
     strcmp(argv[1], "bump")    == 0 ||
     strcmp(argv[1], "represent") == 0 ||
     strcmp(argv[1], "w2b")     == 0 ||
     strcmp(argv[1], "b2w")     == 0 ||
     strcmp(argv[1], "b2wd")    == 0) {
    if(argc < 3) {
      usage();
      return 0;
    }
    const auto recur(4 < argc ? atoi(argv[4]) : 1);
    const auto n(5 < argc ? atoi(argv[5]) : 1);
    typename simpleFile<num_t>::Mat data[3];
    if(!file.loadp2or3(data, argv[2]))
      return - 1;
    if(strcmp(argv[1], "collect") == 0)
      for(int i = 0; i < 3; i ++)
        data[i] = filter<num_t>(data[i], COLLECT_BOTH, n, recur);
    else if(strcmp(argv[1], "integ") == 0)
      for(int i = 0; i < 3; i ++)
        data[i] = filter<num_t>(data[i], INTEG_BOTH, n, recur);
    else if(strcmp(argv[1], "rot") == 0)
      for(int i = 0; i < 3; i ++)
        data[i] = center<num_t>(rotate<num_t>(rotate<num_t>(data[i], atan(num_t(1)) / num_t(n)), - atan(num_t(1)) / num_t(n)), data[i]);
    else if(strcmp(argv[1], "enlarge") == 0)
      for(int i = 0; i < 3; i ++)
        data[i] = filter<num_t>(filter<num_t>(data[i], ENLARGE_BOTH, n, recur), CLIP);
    else if(strcmp(argv[1], "flarge") == 0)
      for(int i = 0; i < 3; i ++)
        data[i] = filter<num_t>(data[i], FLARGE_BOTH, n, recur);
    else if(strcmp(argv[1], "pextend") == 0)
      for(int i = 0; i < 3; i ++)
        data[i] = filter<num_t>(data[i], EXTEND_BOTH, 1, recur);
    else if(strcmp(argv[1], "blink") == 0)
      for(int i = 0; i < 3; i ++)
        data[i] = filter<num_t>(data[i], BLINK_BOTH, n, recur);
    else if(strcmp(argv[1], "lpf") == 0)
      for(int i = 0; i < 3; i ++)
        data[i] = filter<num_t>(data[i], LPF_BOTH, n, recur);
    else if(strcmp(argv[1], "sharpen") == 0)
      for(int i = 0; i < 3; i ++)
        data[i] = filter<num_t>(data[i], SHARPEN_BOTH, n, recur);
    else if(strcmp(argv[1], "bump") == 0)
      data[0] = data[1] = data[2] = filter<num_t>(redig.rgb2d(data), BUMP_BOTH, n, recur);
    else if(strcmp(argv[1], "represent") == 0)
      data[0] = data[1] = data[2] = filter<num_t>(redig.rgb2d(data), REPRESENT, n, recur);
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
    if(strcmp(argv[1], "b2w") != 0 &&
       strcmp(argv[1], "b2wd") != 0)
      redig.normalize(data, num_t(1));
    if(!file.savep2or3(argv[3], data, ! true, strcmp(argv[1], "pextend") == 0 ? 255 : 65535))
      return - 1;
  } else if(strcmp(argv[1], "penl") == 0) {
    if(argc < 5) {
      usage();
      return 0;
    }
    typename simpleFile<num_t>::Mat opt[3];
    typename simpleFile<num_t>::Mat data[3];
    if(!file.loadp2or3(opt, argv[2]))
      return - 1;
    if(!file.loadp2or3(data, argv[3]))
      return - 1;
    typename simpleFile<num_t>::Mat out[3], cnt;
    out[0] = out[1] = out[2] = cnt = typename simpleFile<num_t>::Mat(data[0].rows() * 5 / 4 + 1, data[0].cols() * 5 / 4 + 1).O();
    typename simpleFile<num_t>::Mat cr(5, 5);
    for(int i = 0; i < cr.rows(); i ++)
      for(int j = 0; j < cr.cols(); j ++)
        cr(i, j) = num_t(1);
    for(int i = 0; i < 3; i ++)
      for(int j = 0; j < data[i].rows() - 5; j += 2)
        for(int k = 0; k < data[i].cols() - 5; k += 2) {
          num_t m(data[i](j, k));
          num_t M(data[i](j, k));
          for(int ii = 0; ii < 4; ii ++)
            for(int jj = 0; jj < 4; jj ++) {
              m = min(m, data[i](ii, jj));
              M = max(M, data[i](ii, jj));
            }
          out[i].setMatrix(j * 5 / 4, k * 5 / 4,
            out[i].subMatrix(j * 5 / 4, k * 5 / 4, 5, 5) +
             redig.compImage(redig.normalize(
              data[i].subMatrix(j, k, 4, 4), num_t(1) / num_t(2)), opt[0]) *
             (M == m ? num_t(1) : M - m)
          );
          cnt.setMatrix(j * 5 / 4, k * 5 / 4,
            cnt.subMatrix(j * 5 / 4, k * 5 / 4, 5, 5) + cr);
        }
    for(int i = 0; i < 3; i ++) {
      for(int j = 0; j < data[i].rows(); j ++)
        for(int k = 0; k < data[i].cols(); k ++)
          out[i](j, k) /= max(num_t(1), cnt(j, k));
      out[i] = filter<num_t>(out[i], CLIP, 1, 1);
    }
    if(!file.savep2or3(argv[4], out, ! true, 255))
      return - 1;
  } else if(strcmp(argv[1], "bumpc") == 0) {
    if(argc < 8) {
      usage();
      return 0;
    }
    typename simpleFile<num_t>::Mat datac[3], bump[3];
    if(!file.loadp2or3(datac, argv[5]))
      return -1;
    if(!file.loadp2or3(bump, argv[6]))
      return -1;
    redig.initialize(atoi(argv[3]), std::atof(argv[4]));
    bump[0] = bump[1] = bump[2] = redig.bump(redig.rgb2d(datac), redig.rgb2d(bump), std::atof(argv[2]));
    redig.normalize(bump, num_t(1));
    if(!file.savep2or3(argv[7], bump, ! true) )
      return - 2;
  } else if(strcmp(argv[1], "reshape") == 0 ||
            strcmp(argv[1], "recolor") == 0 ||
            strcmp(argv[1], "recolor2") == 0 ||
            strcmp(argv[1], "recolor3") == 0) {
    if(argc < 6) {
      usage();
      return 0;
    }
    const auto count(atoi(argv[2]));
    typename simpleFile<num_t>::Mat datac[3], datas[3];
    if(!file.loadp2or3(datac, argv[3]))
      return - 1;
    if((strcmp(argv[1], "recolor") == 0 ||
        strcmp(argv[1], "recolor3") == 0 ||
        strcmp(argv[1], "reshape") == 0) &&
       ! file.loadp2or3(datas, argv[4]))
      return - 1;
    if(strcmp(argv[1], "reshape") == 0) {
      const auto datav(redig.rgb2d(datas));
      for(int i = 0; i < 3; i ++)
        datac[i] = redig.reShape(datac[i], datav, count, std::atof(argv[6]));
    } else if(strcmp(argv[1], "recolor3") == 0)
      for(int i = 0; i < 3; i ++)
        datac[i] = redig.reColor3(datac[i], datas[i], count);
    else {
      typename simpleFile<num_t>::Mat xyzc[3], xyzs[3];
      redig.rgb2xyz(xyzc, datac);
      redig.rgb2xyz(xyzs, datas);
      for(int i = 0; i < 3; i ++)
        xyzc[i] = strcmp(argv[1], "recolor") == 0 ?
          redig.reColor(xyzc[i], xyzs[i], count, num_t(std::atof(argv[6]))) :
          redig.reColor(xyzc[i], count, num_t(std::atof(argv[5])));
      redig.xyz2rgb(datac, xyzc);
    }
    redig.normalize(datac, num_t(1));
    if(!file.savep2or3(argv[strcmp(argv[1], "recolor2") == 0 ? 4 : 5], datac, ! true))
      return - 1;
  } else if(strcmp(argv[1], "obj") == 0) {
    typename simpleFile<num_t>::Mat data[3], mask[3];
    if(argc < 8) {
      usage();
      return - 1;
    }
    const auto  vbox(atoi(argv[2]));
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
    vector<typename simpleFile<num_t>::Vec>  points;
    vector<typename simpleFile<num_t>::Veci> facets;
    redig.initialize(vbox, zratio);
    redig.getTileVec(data[0], points, facets);
    for(int i = 0; i < points.size(); i ++)
      points[i] *= ratio;
    const auto ff(redig.floodfill(redig.rgb2d(mask), points)[0]);
    file.saveobj(points, ratio * num_t(data[0].rows()),
                         ratio * num_t(data[0].cols()),
                 8 < argc ? redig.mesh2(points, ff) : facets, argv[sidx],
                 redig.edge(points, ff), thin);
    file.saveMTL(argv[sidx], (string(argv[sidx]) + string(".mtl")).c_str());
  } else if(strcmp(argv[1], "tilt") == 0 ||
            strcmp(argv[1], "sbox") == 0) {
    if(argc < 8) {
      usage();
      return - 1;
    }
    const auto index(atoi(argv[2]));
    const auto Mindex(atoi(argv[3]));
    num_t psi(0);
    num_t zratio(1);
    if(strcmp(argv[1], "tilt") == 0)
      psi = std::atof(argv[4]);
    else
      zratio = std::atof(argv[4]);
    typename simpleFile<num_t>::Mat data[3], bump[3];
    vector<typename simpleFile<num_t>::Vec>  points;
    vector<typename simpleFile<num_t>::Veci> polys;
    if(!file.loadp2or3(data, argv[5]))
      return - 2;
    const string fn(argv[6]);
    if(!file.loadp2or3(bump, argv[6]))
      return - 2;
    typename simpleFile<num_t>::Mat tilt0;
    const auto mtilt(strcmp(argv[1], "sbox") == 0 ? match_t<num_t>() :
                     redig.tiltprep(data[0], index, Mindex, psi));
    const auto depth(strcmp(argv[1], "sbox") == 0 ?
                     num_t(index) / num_t(Mindex) * zratio *
                       sqrt(num_t(data[0].rows() * data[0].cols())) :
                     - num_t(1000000));
    tilt0 = redig.tilt(redig.makeRefMatrix(data[0], 1), bump[0], mtilt, depth);
    for(int j = 0; j < 3; j ++)
      data[j] = redig.pullRefMatrix(tilt0, 1, data[j]);
    if(!file.savep2or3(argv[7], data, ! true))
      return - 1;
  } else if(strcmp(argv[1], "match") == 0) {
    if(argc < 11) {
      usage();
      return - 1;
    }
          int  fnout(11);
    const auto nsub(atoi(argv[2]));
    const auto nemph(atoi(argv[3]));
    const auto vboxdst(atoi(argv[4]));
    const auto vboxsrc(atoi(argv[5]));
    const num_t zratio(std::atof(argv[6]));
    typename simpleFile<num_t>::Mat in0[3], in1[3], bump0orig[3], bump1orig[3], mask0orig[3], mask1orig[3];
    vector<typename simpleFile<num_t>::Veci> delau0, delau1;
    vector<typename simpleFile<num_t>::Vec>  shape0, shape1;
    if(!file.loadp2or3(in0, argv[7]))
      return - 2;
    if(!file.loadp2or3(in1, argv[8]))
      return - 2;
    if(!file.loadp2or3(bump0orig, argv[9]))
      return - 2;
    const string fn(argv[10]);
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
    redig.initialize(vboxsrc, zratio);
    redig.getTileVec(bump1, shape1, delau1);
    redig.maskVectors(shape1, delau1, mask1);
    vector<vector<typename simpleFile<num_t>::Veci> > sdelau0, sdelau1;
    vector<vector<typename simpleFile<num_t>::Vec>  > sshape0, sshape1;
    sdelau0.emplace_back(delau0);
    sdelau1.emplace_back(delau1);
    sshape0.emplace_back(shape0);
    sshape1.emplace_back(shape1);
    vector<match_t<num_t> > mm;
    for(int i = 1; i <= nsub; i ++) {
      const auto m(matchPartial<num_t>(sshape0[i - 1], sshape1[i - 1]));
      mm.emplace_back(m);
      sdelau0.emplace_back(vector<typename simpleFile<num_t>::Veci>());
      sdelau1.emplace_back(vector<typename simpleFile<num_t>::Veci>());
      sshape0.emplace_back(vector<typename simpleFile<num_t>::Vec>());
      sshape1.emplace_back(vector<typename simpleFile<num_t>::Vec>());
      auto dstsort(m.dstpoints);
      auto srcsort(m.srcpoints);
      std::sort(dstsort.begin(), dstsort.end());
      std::sort(srcsort.begin(), srcsort.end());
      vector<int> ds0;
      vector<int> ds1;
      ds0.resize(sshape0[i - 1].size(), - 1);
      ds1.resize(sshape1[i - 1].size(), - 1);
      for(int j = 0, jj = 0; j < sshape0[i - 1].size(); j ++) {
        if(binary_search(dstsort.begin(), dstsort.end(), j))
          continue;
        sshape0[i].emplace_back(sshape0[i - 1][j]);
        ds0[j] = jj ++;
      }
      for(int j = 0, jj = 0; j < sshape1[i - 1].size(); j ++) {
        if(binary_search(srcsort.begin(), srcsort.end(), j))
          continue;
        sshape1[i].emplace_back(sshape1[i - 1][j]);
        ds1[j] = jj ++;
      }
      sdelau0[i].reserve(sdelau0[i - 1].size());
      sdelau1[i].reserve(sdelau1[i - 1].size());
      for(int j = 0; j < sdelau0[i - 1].size(); j ++) {
        auto sd0(sdelau0[i - 1][j]);
        for(int k = 0; k < sd0.size(); k ++)
          if((sd0[k] = ds0[sd0[k]]) < 0) goto nxtsd0;
        sdelau0[i].emplace_back(move(sd0));
       nxtsd0:
        ;
      }
      for(int j = 0; j < sdelau1[i - 1].size(); j ++) {
        auto sd1(sdelau1[i - 1][j]);
        for(int k = 0; k < sd1.size(); k ++)
          if((sd1[k] = ds1[sd1[k]]) < 0) goto nxtsd1;
        sdelau1[i].emplace_back(move(sd1));
       nxtsd1:
        ;
      }
      if(! sdelau0[i].size() || ! sdelau1[i].size() || ! sshape0[i].size() || ! sshape1[i].size())
        break;
    }
    typename simpleFile<num_t>::Mat outs[3];
    const auto rin0(redig.makeRefMatrix(in0[0], 1));
    const auto rin1(redig.makeRefMatrix(in1[0], 1 + rin0.rows() * rin0.cols()));
    vector<vector<typename simpleFile<num_t>::Veci> > mhull0;
    vector<vector<typename simpleFile<num_t>::Veci> > mhull1;
    for(int i = 0; i < mm.size(); i ++) {
      mhull0.emplace_back(redig.mesh2(sshape0[i], mm[i].dstpoints));
      mhull1.emplace_back((~ mm[i]).hullConv(mhull0[i]));
    }
    const string outbase(argv[fnout]);
    for(int idx = 0; idx < 3; idx ++)
      for(int i = 0; i < mm.size(); i ++) {
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
    file.savep2or3((outbase + string("-repl0.ppm")).c_str(), outs, false);
    for(int idx = 0; idx < 3; idx ++)
      for(int i = 0; i < mm.size(); i ++) {
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
    file.savep2or3((outbase + string("-repl1.ppm")).c_str(), outs, false);
    for(int idx = 0; idx < 3; idx ++) {
      for(int i = 0; i < mm.size(); i ++)
        if(i)
          outs[idx] += redig.showMatch(in0[idx], sshape0[i], sdelau0[i]);
        else
          outs[idx]  = redig.showMatch(in0[idx], sshape0[i], sdelau0[i]);
    }
    redig.normalize(outs, 1.);
    file.savep2or3((outbase + string(".ppm")).c_str(), outs, false);
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
      file.savep2or3((outbase + string("-") +
                                to_string(i) +
                                string("-") +
                                to_string(nemph) +
                                string(".ppm")).c_str(), outs, false);
/*
      file.saveobj(redig.takeShape(shape0, shape1, m,   iemph),
                   outs[0].rows(), outs[0].cols(), delau0,
                   (outbase + string("-emph0-") +
                              to_string(i) +
                              string("-") +
                              to_string(nemph) +
                              string(".obj")).c_str());
      file.saveobj(shape, outs[0].rows(), outs[0].cols(), delau1,
                   (outbase + string("-emph1-") +
                              to_string(i) +
                              string("-") +
                              to_string(nemph) +
                              string(".obj")).c_str());
*/
    }
  } else if(strcmp(argv[1], "pred") == 0 ||
            strcmp(argv[1], "lenl") == 0 ||
            strcmp(argv[1], "cat") == 0 ||
            strcmp(argv[1], "composite") == 0) {
    if(argc < 4) {
      usage();
      return - 1;
    }
    vector<vector<typename simpleFile<num_t>::Mat> > in;
    in.resize(argc - 3);
    for(int i = 3; i < argc; i ++) {
      typename simpleFile<num_t>::Mat ibuf[3];
      if(!file.loadp2or3(ibuf, argv[i]))
        return - 2;
      const auto ii(i - 3);
      in[ii].resize(3);
      for(int j = 0; j < 3; j ++)
        in[ii][j] = std::move(ibuf[j]);
      if(strcmp(argv[1], "cat") == 0 && ! (ii & 1)) {
        assert(in[ii][0].rows() == in[0][0].rows());
        assert(in[ii][0].cols() == in[0][0].cols());
      }
    }
    const auto idx(in.size() - 1);
    typename simpleFile<num_t>::Mat out[3];
    if(strcmp(argv[1], "pred") == 0) {
      for(int i = 0; i < 3; i ++)
        out[i].resize(in[idx][0].rows(), in[idx][0].cols());
      const auto comp(nextP0<num_t>(in.size() * 2 + 1));
      for(int y = 0; y < out[0].rows(); y ++) {
        for(int x = 0; x < out[0].cols(); x ++) {
          out[0](y, x) = out[1](y, x) = out[2](y, x) = num_t(0);
          for(int k = 0; k < in.size(); k ++)
            for(int m = 0; m < 3; m ++) {
              out[m](y, x) += atan(in[k][m](y, x)) * comp[2 * k];
              if(k < in.size() - 1)
                out[m](y, x) += (atan(in[k][m](y, x)) +
                  atan(in[k + 1][m](y, x))) / num_t(2) * comp[2 * k + 1];
            }
          for(int m = 0; m < 3; m ++)
            out[m](y, x) = tan(out[m](y, x));
        }
      }
      for(int i = 0; i < 3; i ++)
        out[i] = filter<num_t>(out[i], CLIP);
      file.savep2or3(argv[2], out, ! true);
    } else if(strcmp(argv[1], "lenl") == 0) {
      std::vector<std::pair<typename simpleFile<num_t>::Mat, typename simpleFile<num_t>::Mat> > pair;
      const auto d5(dft<num_t>(5).subMatrix(0, 0, 4, 5));
      const auto d4(dft<num_t>(- 4));
      const auto taylc(d4 * d5);
      const auto tayl(taylc.template real<num_t>());
      for(int i = 0; i < in.size(); i ++)
        for(int j = 0; j < 3; j ++)
          for(int k = 0; k < in[i][j].rows() - 5; k += 2)
            for(int kk = 0; kk < in[i][j].cols() - 5; kk += 2) {
              const auto work(redig.normalize(in[i][j].subMatrix(k, kk, 5, 5), num_t(1) / num_t(2)));
              pair.emplace_back(std::make_pair(work, redig.normalize(tayl * work * tayl.transpose(), num_t(1) / num_t(2))));
            }
      out[0] = out[1] = out[2] = redig.optImage(pair);
      file.savep2or3(argv[2], out, ! true, 65535);
    } else if(strcmp(argv[1], "cat") == 0) {
      vector<typename simpleFile<num_t>::Mat> rep;
      vector<typename simpleFile<num_t>::Mat> glay;
      glay.reserve(in.size());
      for(int i = 0; i < in.size(); i ++) {
        typename simpleFile<num_t>::Mat inn[3];
        for(int j = 0; j < 3; j ++)
          inn[j] = const_cast<typename simpleFile<num_t>::Mat &&>(in[i][j]);
        (i & 1 ? rep : glay).emplace_back(redig.rgb2d(inn));
      }
      const auto cat(redig.catImage(rep, glay));
      for(int i = 0; i < cat.size(); i ++) {
        out[0] = out[1] = out[2] = cat[i];
        redig.normalize(out, 1.);
        file.savep2or3((string(argv[2]) + string("-") + to_string(i) + string(".ppm")).c_str(), out, ! true);
      }
    } else if(strcmp(argv[1], "composite") == 0) {
      vector<typename simpleFile<num_t>::Mat> glay;
      glay.reserve(in.size());
      for(int i = 0; i < in.size(); i ++) {
        typename simpleFile<num_t>::Mat inn[3];
        for(int j = 0; j < 3; j ++)
          inn[j] = std::move(in[i][j]);
        glay.emplace_back(redig.rgb2d(inn));
      }
      const auto composite(redig.compositeImage(glay));
      for(int i = 0; i < composite.size(); i ++) {
        out[0] = out[1] = out[2] = composite[i];
        redig.normalize(out, 1.);
        file.savep2or3((string(argv[2]) + string("-") + to_string(i) + string(".ppm")).c_str(), out, ! true);
      }
    }
  } else if(strcmp(argv[1], "habit") == 0) {
    vector<typename simpleFile<num_t>::Vec>  pdst,   psrc;
    vector<typename simpleFile<num_t>::Veci> poldst, polsrc;
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
      const auto m(redig.tiltprep(typename simpleFile<num_t>::Mat(int(My), int(Mx)), - atoi(argv[4]), atoi(argv[5]), std::atof(argv[6])));
      file.saveobj(redig.takeShape(pdst, psrc, m, num_t(1) / num_t(2)),
                   My, Mx, poldst, argv[7]);
    } else {
      const auto m(matchPartial<num_t>(pdst, psrc));
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
          num_t(std::atof(argv[6])), atoi(argv[2]));
    else
      for(int i = 0; i < 3; i ++)
        out[i] = redig.reImage(dst[i], src[i],
          num_t(std::atof(argv[6])), atoi(argv[2]));
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
          num_t(std::atof(argv[5])), atoi(argv[2]));
    else {
      for(int i = 0; i < 3; i ++)
        out[i] = redig.reImage(dst[i],
          num_t(std::atof(argv[5])), atoi(argv[2]));
      redig.autoLevel(out, 4 * (out[0].rows() + out[0].cols()));
    }
    redig.normalize(out, 1.);
    if(!file.savep2or3(argv[4], out, ! true, 255))
      return - 3;
  } else if(strcmp(argv[1], "newtrace") == 0) {
    if(argc < 5) {
      usage();
      return -1;
    }
    SimpleVector<num_t> m(atoi(argv[2]));
    Decompose<num_t> dec(m.size());
    auto n(m);
    auto f(m);
    m[0] = n[0] = num_t(atoi(argv[2]) * 2);
    // XXX (f[0] := 1 causes flat result.):
    f[0] = num_t(0);
    for(int i = 1; i < m.size(); i ++) {
      m[i] = m[i - 1] + (num_t(arc4random_uniform(3000)) - num_t(1500)) / num_t(1500);
      n[i] = n[i - 1] + (num_t(arc4random_uniform(3000)) - num_t(1500)) / num_t(1500);
      f[i] = num_t(1);
    }
    f /= sqrt(f.dot(f));
    const auto pp(make_pair(make_pair(atoi(argv[3]),
      atoi(argv[3])), make_pair(0, 0)));
    typename simpleFile<num_t>::Mat M[3];
    M[0] = M[1] = M[2] = redig.normalize(redig.applyTrace(
      make_pair(dec.synth(m, f), dec.synth(n, f)),
      make_pair(pp, pp)), num_t(1));
    file.savep2or3(argv[4], M, true, 255);
  } else if(strcmp(argv[1], "omake") == 0) {
    vector<vector<num_t> > data;
    string header;
    file.loaddat(argv[3], header, data);
    simpleFile<num_t>::Mat buf(atoi(argv[4]), atoi(argv[4]));
    const auto mdft(dft<num_t>(buf.rows()));
    const auto midft(dft<num_t>(- buf.rows()));
    for(int i0 = 1; i0 < data.size(); i0 ++) {
      for(int i = 0; i <= data[i0].size() / buf.rows() / buf.rows(); i ++) {
        for(int k = 0; k < buf.cols(); k ++)
          for(int j = 0; j < buf.rows(); j ++) {
            const auto idx(i * buf.rows() * buf.rows() + k * buf.rows() + j);
            buf(j, k) = idx < data[i0].size() ? data[i0][idx] : num_t(0);
          }
        typename simpleFile<num_t>::Mat buf2;
        if(strcmp(argv[2], "diff") == 0)
          buf2 = (midft * (
            filter<num_t>(mdft.template real<num_t>() * buf, DETECT_X).template cast<complex<num_t> >() +
            filter<num_t>(mdft.template imag<num_t>() * buf, DETECT_X).template cast<complex<num_t> >() * complex<num_t>(num_t(0), num_t(1)) ) ).template real<num_t>();
        else if(strcmp(argv[2], "sharpen") == 0)
          buf2 = (midft * (
            filter<num_t>(mdft.template real<num_t>() * buf, SHARPEN_X).template cast<complex<num_t> >() +
            filter<num_t>(mdft.template imag<num_t>() * buf, SHARPEN_X).template cast<complex<num_t> >() * complex<num_t>(num_t(0), num_t(1)) ) ).template real<num_t>();
        else if(strcmp(argv[2], "bump") == 0)
          buf2 = (midft * (
            filter<num_t>(mdft.template real<num_t>() * buf, BUMP_BOTH).template cast<complex<num_t> >() +
            filter<num_t>(mdft.template imag<num_t>() * buf, BUMP_BOTH).template cast<complex<num_t> >() * complex<num_t>(num_t(0), num_t(1)) ) ).template real<num_t>();
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
  } else {
    usage();
    return - 1;
  }
  return 0;
}

