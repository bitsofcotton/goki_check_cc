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
  cout << "gokicheck (pred|lenl|composite) <output.ppm> <input0.ppm> ..." << endl;
  cout << "gokicheck (cat|catr) <input0.ppm> ..." << endl;
  cout << "gokicheck obj <gather_pixels> <ratio> <zratio> <input.ppm> <output.obj>" << endl;
  cout << "gokicheck (tilt|sbox) <index> <max_index> (<psi>|<zratio>) <input.ppm> <input-bump.ppm> <output.ppm>" << endl;
  cout << "gokicheck match <nsub> <nemph> <vbox_dst> <vbox_src> <zratio> <dst.ppm> <src.ppm> <dst-bump.ppm> <src-bump.ppm> <output-basename>" << endl;
  cout << "gokicheck habit   <in0.obj> <in1.obj> <out.obj>" << endl;
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
     strcmp(argv[1], "b2wd")    == 0 ||
     strcmp(argv[1], "flatten") == 0) {
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
    } else if(strcmp(argv[1], "flatten") == 0) {
            auto row(data[0].row(0));
            auto col(data[0].col(0));
      for(int k = 0; k < 3; k ++)
        for(int i = 1; i < data[k].rows(); i ++)
          row += data[k].row(i);
      for(int k = 0; k < 3; k ++)
        for(int i = 1; i < data[k].cols(); i ++)
          col += data[k].col(i);
      row /= num_t(3 * data[0].rows());
      col /= num_t(3 * data[0].cols());
      // Sum(f(x) + t * x) == 0. <=> Sum(f(x)) + t * x(x + 1) / 2 == 0
      // t = - 2 * Sum(f(x)) / x / (x + 1)
      num_t rt(0);
      num_t ct(0);
      for(int i = 0; i < row.size(); i ++)
        rt += row[i];
      for(int i = 0; i < col.size(); i ++)
        ct += col[i];
      rt = num_t(2) * rt / num_t(row.size() * (row.size() - 1));
      ct = num_t(2) * ct / num_t(col.size() * (col.size() - 1));
      for(int i = 0; i < data[0].rows(); i ++)
        for(int j = 0; j < data[0].cols(); j ++)
          for(int k = 0; k < 3; k ++)
            data[k](i, j) += ct * num_t(i) + rt * num_t(j);
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
    if(argc < 7) {
      usage();
      return - 1;
    }
    const auto  vbox(atoi(argv[2]));
    const num_t ratio(std::atof(argv[3]));
    const num_t zratio(std::atof(argv[4]));
    if(!file.loadp2or3(data, argv[5]))
      return - 1;
    redig.initialize(abs(vbox), zratio);
          auto points(vbox < 0 ? redig.getHesseVec(redig.rgb2d(data))
                               : redig.getTileVec(redig.rgb2d(data)));
    for(int i = 0; i < points.size(); i ++)
      points[i] *= ratio;
    file.saveobj(points, ratio * num_t(data[0].rows()),
                         ratio * num_t(data[0].cols()),
                 redig.mesh2(points), argv[6]);
    file.saveMTL(argv[6], (string(argv[6]) + string(".mtl")).c_str());
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
    tilt0 = redig.tilt(redig.makeRefMatrix(data[0], 1), bump[0],
      strcmp(argv[1], "sbox") == 0 ? match_t<num_t>() :
        redig.tiltprep(data[0], index, Mindex, psi),
      strcmp(argv[1], "sbox") == 0 ?
        num_t(index) / num_t(Mindex) * zratio *
          sqrt(num_t(data[0].rows() * data[0].cols())) :
        - num_t(1000000)
      );
    for(int j = 0; j < 3; j ++)
      data[j] = redig.pullRefMatrix(tilt0, 1, data[j]);
    if(!file.savep2or3(argv[7], data, ! true))
      return - 1;
  } else if(strcmp(argv[1], "match") == 0) {
    if(argc < 11) {
      usage();
      return - 1;
    }
    const auto nsub(atoi(argv[2]));
    const auto nemph(atoi(argv[3]));
    const auto vboxdst(atoi(argv[4]));
    const auto vboxsrc(atoi(argv[5]));
    const num_t zratio(std::atof(argv[6]));
    typename simpleFile<num_t>::Mat in0[3], in1[3], bump0[3], bump1[3];
    if(!file.loadp2or3(in0, argv[7]))
      return - 2;
    if(!file.loadp2or3(in1, argv[8]))
      return - 2;
    if(!file.loadp2or3(bump0, argv[9]))
      return - 2;
    if(!file.loadp2or3(bump1, argv[10]))
      return - 2;
    const string outbase(argv[11]);
    redig.initialize(abs(vboxdst), zratio);
    const auto shape0(vboxdst < 0
      ? redig.getHesseVec(redig.rgb2d(bump0))
      : redig.getTileVec(redig.rgb2d(bump0)));
    redig.initialize(abs(vboxsrc), zratio);
    const auto shape1(vboxsrc < 0
      ? redig.getHesseVec(redig.rgb2d(bump1))
      : redig.getTileVec(redig.rgb2d(bump1)));
    auto m(shape0.size() < shape1.size()
      ? matchPartial< num_t>(shape0, shape1, nsub)
      : matchPartialR<num_t>(shape0, shape1, nsub));
    typename simpleFile<num_t>::Mat outs[3];
    const auto rin0(redig.makeRefMatrix(in0[0], 1));
    const auto rin1(redig.makeRefMatrix(in1[0], 1 + rin0.rows() * rin0.cols()));
    vector<vector<typename simpleFile<num_t>::Veci> > mhull0;
    vector<vector<typename simpleFile<num_t>::Veci> > mhull1;
    for(int i = 0; i < m.size(); i ++) {
      mhull0.emplace_back(redig.mesh2(shape0, m[i].dst));
      mhull1.emplace_back((~ m[i]).hullConv(mhull0[i]));
    }
    outs[0] = SimpleMatrix<num_t>(in0[0].rows(), in0[0].cols()).O();
    for(int i = 0; i < m.size(); i ++)
      outs[0] += redig.showMatch(redig.draw(in0[0] * num_t(0),
                     shape0, mhull0[i]), shape0, mhull0[i]);
    outs[1] = outs[2] = outs[0];
    redig.normalize(outs, 1.);
    file.savep2or3((outbase + string("-repl0.ppm")).c_str(), outs, false);
    outs[0] = SimpleMatrix<num_t>(in1[0].rows(), in1[0].cols()).O();
    for(int i = 0; i < m.size(); i ++)
      outs[0] += redig.showMatch(redig.draw(in1[0] * num_t(0),
                   shape1, mhull1[i]), shape1, mhull1[i]);
    outs[1] = outs[2] = outs[0];
    redig.normalize(outs, 1.);
    file.savep2or3((outbase + string("-repl1.ppm")).c_str(), outs, false);
    outs[0] = SimpleMatrix<num_t>(in0[0].rows(), in0[0].cols()).O();
    for(int i = 0; i < m.size(); i ++)
      outs[0] += redig.showMatch(redig.draw(in0[0] * num_t(0),
                   m[i].transform(shape1), mhull1[i]),
                   m[i].transform(shape1), mhull1[i]);
    outs[1] = outs[2] = outs[0];
    redig.normalize(outs, 1.);
    file.savep2or3((outbase + string("-repl2.ppm")).c_str(), outs, false);
    for(int i = 0; i < nemph; i ++) {
      const auto iemph(num_t(i) / num_t(nemph));
      simpleFile<num_t>::Mat reref(rin1.rows(), rin1.cols());
      reref.O();
      for(int i = 0; i < m.size(); i ++) {
        const auto rd(redig.draw(rin1, shape1,
          redig.takeShape(shape1, shape0, ~ m[i], iemph), mhull1[i]));
        for(int j = 0; j < min(reref.rows(), rd.rows()); j ++)
          for(int k = 0; k < min(reref.cols(), rd.cols()); k ++)
            if(rd(j, k) != num_t(0)) reref(j, k) = rd(j, k);
      }
      for(int idx = 0; idx < 3; idx ++)
        outs[idx] = redig.pullRefMatrix(reref, 1 + rin0.rows() * rin0.cols(), in1[idx]);
      file.savep2or3((outbase + string("-") + to_string(i) + string("-") +
                      to_string(nemph) + string(".ppm")).c_str(), outs, false);
    }
  } else if(strcmp(argv[1], "pred") == 0 ||
            strcmp(argv[1], "lenl") == 0 ||
            strcmp(argv[1], "cat") == 0 ||
            strcmp(argv[1], "catr") == 0 ||
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
    }
    const auto idx(in.size() - 1);
    typename simpleFile<num_t>::Mat out[3];
    if(strcmp(argv[1], "pred") == 0) {
      for(int i = 0; i < 3; i ++) {
        out[i].resize(in[idx][0].rows(), in[idx][0].cols());
        out[i].O();
      }
      const auto& comp(pnext<num_t>(in.size() * 2 + 1));
      for(int y = 0; y < out[0].rows(); y ++)
        for(int x = 0; x < out[0].cols(); x ++)
          for(int k = 0; k < in.size(); k ++)
            for(int m = 0; m < 3; m ++) {
              out[m](y, x) += in[k][m](y, x) * comp[2 * k];
              if(k < in.size() - 1)
                out[m](y, x) += (in[k][m](y, x) +
                  in[k + 1][m](y, x)) / num_t(2) * comp[2 * k + 1];
            }
      redig.normalize(out, num_t(1));
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
              pair.emplace_back(std::make_pair(redig.normalize(tayl * work * tayl.transpose(), num_t(1) / num_t(2)), work));
            }
      out[0] = out[1] = out[2] = redig.optImage(pair);
      file.savep2or3(argv[2], out, ! true, 65535);
    } else if(strcmp(argv[1], "cat") == 0 ||
              strcmp(argv[1], "catr") == 0) {
      vector<typename simpleFile<num_t>::Mat> glay;
      glay.reserve(strcmp(argv[1], "catr") == 0 ? in.size()
                    : in.size() * min(int(in[0][0].rows()), 1 + 5 + 1));
      const auto in00sz(in[0][0].rows());
      for(int i = 0; i < in.size(); i ++) {
        typename simpleFile<num_t>::Mat inn[3];
        for(int j = 0; j < 3; j ++)
          inn[j] = std::move(in[i][j]);
        if(strcmp(argv[1], "catr") == 0)
          glay.emplace_back(redig.rgb2d(inn));
        else {
          auto work(redig.rgb2d(inn));
          assert(work.rows() == in00sz);
          for(int j = 0; j < min(int(work.rows()), 1 + 5 + 1); j ++) {
            SimpleMatrix<num_t> rr(1, work.cols());
            rr.row(0) = std::move(work.row(j));
            glay.emplace_back(rr);
          }
        }
      }
      const auto cat(redig.catImage(glay));
      for(int i = 0; i < cat.size(); i ++) {
        for(int j = 0; j < cat[i].size(); j ++)
          std::cout << argv[3 + (strcmp(argv[1], "catr") == 0 ? cat[i][j]
                         : cat[i][j] / min(in00sz, 1 + 5 + 1)) ] << std::endl;
        std::cout << std::endl;
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
    file.saveobj(redig.takeShape(pdst, psrc,
      matchPartial<num_t>(pdst, psrc)[0],
      num_t(1) / num_t(2)), My, Mx, poldst, argv[4]);
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
            filter<num_t>(mdft.template real<num_t>() * buf, DETECT_BOTH).template cast<complex<num_t> >() +
            filter<num_t>(mdft.template imag<num_t>() * buf, DETECT_BOTH).template cast<complex<num_t> >() * complex<num_t>(num_t(0), num_t(1)) ) ).template real<num_t>();
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

