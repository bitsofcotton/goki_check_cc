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

#if defined(_OPENMP)
#include <omp.h>
#endif

//#define int int64_t
#define int int32_t
#include "lieonn.hh"
typedef myfloat num_t;
#include "p0.hh"
#include "p1.hh"
#include "decompose.hh"
#include "catg.hh"
#include "goki.hh"

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
  cout << "gokicheck (collect|sharpen|bump|enlarge|flarge|pextend|blink|represent) <input.ppm> <output.ppm> <recur>" << endl;
  cout << "gokicheck (pred|lenl|composite) <output.ppm> <input0.ppm> ..." << endl;
  cout << "gokicheck (cat|catr) <input0.ppm> ..." << endl;
  cout << "gokicheck obj <gather_pixels> <ratio> <input.ppm> <output.obj>" << endl;
  cout << "gokicheck (tilt|sbox) <index> <max_index> <psi> <input.ppm> <input-bump.ppm> <output.ppm>" << endl;
  cout << "gokicheck match <nsub> <nemph> <vbox_dst> <vbox_src> <dst.ppm> <src.ppm> <dst-bump.ppm> <src-bump.ppm> <output-basename>" << endl;
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
  if(strcmp(argv[1], "nop") == 0 ||
     strcmp(argv[1], "collect") == 0 ||
     strcmp(argv[1], "integ") == 0 ||
     strcmp(argv[1], "diffraw") == 0 ||
     strcmp(argv[1], "integraw") == 0 ||
     strcmp(argv[1], "enlarge") == 0 ||
     strcmp(argv[1], "flarge") == 0 ||
     strcmp(argv[1], "pextend") == 0 ||
     strcmp(argv[1], "blink") == 0 ||
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
    const auto rot(5 < argc ? atoi(argv[5]) : 0);
    vector<SimpleMatrix<num_t> > data;
    if(!loadp2or3<num_t>(data, argv[2]))
      return - 1;
    if(strcmp(argv[1], "collect") == 0)
      for(int i = 0; i < data.size(); i ++)
        data[i] = filter<num_t>(data[i], COLLECT_BOTH, recur, rot);
    else if(strcmp(argv[1], "integ") == 0)
      for(int i = 0; i < data.size(); i ++)
        data[i] = filter<num_t>(data[i], INTEG_BOTH, recur, rot);
    else if(strcmp(argv[1], "diffraw") == 0) {
      auto avgdiffL(dft<num_t>(data[0].rows()));
      auto avgdiffR(dft<num_t>(data[0].cols()));
      for(int i = 1; i < avgdiffL.rows(); i ++)
        avgdiffL.row(i)  *= complex<num_t>(num_t(int(0)), - num_t(int(2))) * atan(num_t(int(1))) * num_t(int(4)) * num_t(i) / num_t(int(avgdiffL.rows()));
      for(int i = 1; i < avgdiffR.rows(); i ++)
        avgdiffR.row(i)  *= complex<num_t>(num_t(int(0)), - num_t(int(2))) * atan(num_t(int(1))) * num_t(int(4)) * num_t(i) / num_t(int(avgdiffR.rows()));
      const auto dL((dft<num_t>(- data[0].rows()) * avgdiffL).template real<num_t>());
      const auto dR((dft<num_t>(- data[0].cols()) * avgdiffR).template real<num_t>());
      for(int i = 0; i < data.size(); i ++)
        data[i] = dL * data[i] * dR;
    } else if(strcmp(argv[1], "integraw") == 0)
      for(int i = 0; i < data.size(); i ++)
        data[i] = diff<num_t>(- data[i].rows()) * data[i] * diff<num_t>(- data[i].cols()).transpose();
    else if(strcmp(argv[1], "enlarge") == 0)
      for(int i = 0; i < data.size(); i ++)
        data[i] = filter<num_t>(filter<num_t>(data[i], ENLARGE_BOTH, recur, rot), CLIP);
    else if(strcmp(argv[1], "flarge") == 0)
      for(int i = 0; i < data.size(); i ++)
        data[i] = filter<num_t>(data[i], FLARGE_BOTH, recur, rot);
    else if(strcmp(argv[1], "pextend") == 0)
      for(int i = 0; i < data.size(); i ++)
        data[i] = filter<num_t>(data[i], EXTEND_BOTH, recur);
    else if(strcmp(argv[1], "blink") == 0)
      for(int i = 0; i < data.size(); i ++)
        data[i] = filter<num_t>(data[i], BLINK_BOTH, recur, rot);
    else if(strcmp(argv[1], "sharpen") == 0)
      for(int i = 0; i < data.size(); i ++)
        data[i] = filter<num_t>(data[i], SHARPEN_BOTH);
    else if(strcmp(argv[1], "represent") == 0)
      data[0] = data[1] = data[2] = filter<num_t>(rgb2d<num_t>(data), REPRESENT, recur);
    else if(strcmp(argv[1], "bump") == 0)
      data[0] = data[1] = data[2] = autoLevel<num_t>(filter<num_t>(rgb2d<num_t>(data), BUMP_BOTH, recur, rot), data[0].rows() + data[0].cols());
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
    }
    if(!savep2or3<num_t>(argv[3],
        strcmp(argv[1], "b2w") != 0 && strcmp(argv[1], "b2wd") != 0 &&
        strcmp(argv[1], "nop") != 0
        ? normalize<num_t>(data) : data,
        ! true, 65535))
      return - 1;
  } else if(strcmp(argv[1], "penl") == 0) {
    if(argc < 5) {
      usage();
      return 0;
    }
    vector<SimpleMatrix<num_t> > opt, data;
    if(!loadp2or3<num_t>(opt, argv[2]))
      return - 1;
    if(!loadp2or3<num_t>(data, argv[3]))
      return - 1;
    vector<SimpleMatrix<num_t> > out;
    SimpleMatrix<num_t> cnt;
    out.resize(3);
    out[0] = out[1] = out[2] = cnt = SimpleMatrix<num_t>(data[0].rows() * 5 / 4 + 1, data[0].cols() * 5 / 4 + 1).O();
    SimpleMatrix<num_t> cr(5, 5);
    for(int i = 0; i < cr.rows(); i ++)
      for(int j = 0; j < cr.cols(); j ++)
        cr(i, j) = num_t(1);
    for(int i = 0; i < data.size(); i ++)
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
             compImage<num_t>(normalize<num_t>(
              data[i].subMatrix(j, k, 4, 4), num_t(1) / num_t(2)), opt[0]) *
             (M == m ? num_t(1) : M - m)
          );
          cnt.setMatrix(j * 5 / 4, k * 5 / 4,
            cnt.subMatrix(j * 5 / 4, k * 5 / 4, 5, 5) + cr);
        }
    for(int i = 0; i < data.size(); i ++) {
      for(int j = 0; j < data[i].rows(); j ++)
        for(int k = 0; k < data[i].cols(); k ++)
          out[i](j, k) /= max(num_t(1), cnt(j, k));
      out[i] = filter<num_t>(out[i], CLIP);
    }
    if(!savep2or3<num_t>(argv[4], out, ! true, 65535))
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
    vector<SimpleMatrix<num_t> > datac, datas;
    if(!loadp2or3<num_t>(datac, argv[3]))
      return - 1;
    if((strcmp(argv[1], "recolor") == 0 ||
        strcmp(argv[1], "recolor3") == 0 ||
        strcmp(argv[1], "reshape") == 0) &&
       ! loadp2or3<num_t>(datas, argv[4]))
      return - 1;
    if(strcmp(argv[1], "reshape") == 0) {
      const auto datav(rgb2d<num_t>(datas));
      for(int i = 0; i < datac.size(); i ++)
        datac[i] = reShape<num_t>(datac[i], datav, count, std::atof(argv[6]));
    } else if(strcmp(argv[1], "recolor3") == 0)
      for(int i = 0; i < datac.size(); i ++)
        datac[i] = reColor3<num_t>(datac[i], datas[i], count);
    else {
      auto xyzc(rgb2xyz<num_t>(datac));
      auto xyzs(rgb2xyz<num_t>(datas));
      for(int i = 0; i < xyzc.size(); i ++)
        xyzc[i] = strcmp(argv[1], "recolor") == 0 ?
          reColor<num_t>(xyzc[i], xyzs[i], count, num_t(std::atof(argv[6]))) :
          reColor<num_t>(xyzc[i], count, num_t(std::atof(argv[5])));
      datac = xyz2rgb<num_t>(xyzc);
    }
    if(!savep2or3<num_t>(argv[strcmp(argv[1], "recolor2") == 0 ? 4 : 5], normalize<num_t>(datac), ! true))
      return - 1;
  } else if(strcmp(argv[1], "obj") == 0) {
    vector<SimpleMatrix<num_t> > data, mask;
    if(argc < 6) {
      usage();
      return - 1;
    }
    const auto  vbox(atoi(argv[2]));
    const num_t ratio(std::atof(argv[3]));
    if(!loadp2or3<num_t>(data, argv[4]))
      return - 1;
    auto points(vbox < 0 ? getHesseVec<num_t>(rgb2d<num_t>(data), abs(vbox))
                         : getTileVec<num_t>(rgb2d<num_t>(data), abs(vbox)));
    for(int i = 0; i < points.size(); i ++)
      points[i] *= ratio;
    saveobj<num_t>(points, ratio * num_t(data[0].rows()),
                           ratio * num_t(data[0].cols()),
                   mesh2<num_t>(points), argv[5]);
    saveMTL<num_t>(argv[5], (string(argv[5]) + string(".mtl")).c_str());
  } else if(strcmp(argv[1], "tilt") == 0 ||
            strcmp(argv[1], "sbox") == 0) {
    if(argc < 8) {
      usage();
      return - 1;
    }
    const auto index(atoi(argv[2]));
    const auto Mindex(atoi(argv[3]));
    num_t psi(0);
    psi = std::atof(argv[4]);
    vector<SimpleMatrix<num_t> > data, bump;
    vector<SimpleVector<num_t> > points;
    vector<SimpleVector<int>   > polys;
    if(!loadp2or3<num_t>(data, argv[5]))
      return - 2;
    const string fn(argv[6]);
    if(!loadp2or3<num_t>(bump, argv[6]))
      return - 2;
    const auto tilt0(tilt<num_t>(makeRefMatrix<num_t>(data[0], 1), bump[0],
      strcmp(argv[1], "sbox") == 0 ? match_t<num_t>() :
        tiltprep<num_t>(data[0], index, Mindex, psi),
      strcmp(argv[1], "sbox") == 0 ?
        num_t(index) / num_t(Mindex) *
          num_t(min(data[0].rows(), data[0].cols())) :
        - num_t(1000000)
      ));
    for(int j = 0; j < data.size(); j ++)
      data[j] = pullRefMatrix<num_t>(tilt0, 1, data[j]);
    if(!savep2or3<num_t>(argv[7], data, ! true))
      return - 1;
  } else if(strcmp(argv[1], "match") == 0) {
    if(argc < 10) {
      usage();
      return - 1;
    }
    const auto nsub(atoi(argv[2]));
    const auto nemph(atoi(argv[3]));
    const auto vboxdst(atoi(argv[4]));
    const auto vboxsrc(atoi(argv[5]));
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
    const auto shape0(vboxdst < 0
      ? getHesseVec<num_t>(rgb2d<num_t>(bump0), abs(vboxdst))
      : getTileVec<num_t>(rgb2d<num_t>(bump0), abs(vboxdst)));
    const auto shape1(vboxsrc < 0
      ? getHesseVec<num_t>(rgb2d<num_t>(bump1), abs(vboxsrc))
      : getTileVec<num_t>(rgb2d<num_t>(bump1), abs(vboxsrc)));
    auto m(shape0.size() < shape1.size()
      ? matchPartial< num_t>(shape0, shape1, nsub)
      : matchPartialR<num_t>(shape0, shape1, nsub));
    vector<SimpleMatrix<num_t> > out;
    out.resize(3);
    const auto rin0(makeRefMatrix<num_t>(in0[0], 1));
    const auto rin1(makeRefMatrix<num_t>(in1[0], 1 + rin0.rows() * rin0.cols()));
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
    savep2or3<num_t>((outbase + string("-repl0.ppm")).c_str(), normalize<num_t>(out), false);
    out[0] = SimpleMatrix<num_t>(in1[0].rows(), in1[0].cols()).O();
    for(int i = 0; i < m.size(); i ++)
      out[0] += showMatch<num_t>(draw<num_t>(in1[0] * num_t(0),
                   shape1, mhull1[i]), shape1, mhull1[i]);
    out[1] = out[2] = out[0];
    savep2or3<num_t>((outbase + string("-repl1.ppm")).c_str(), normalize<num_t>(out), false);
    out[0] = SimpleMatrix<num_t>(in0[0].rows(), in0[0].cols()).O();
    for(int i = 0; i < m.size(); i ++)
      out[0] += showMatch<num_t>(draw<num_t>(in0[0] * num_t(0),
                   m[i].transform(shape1), mhull1[i]),
                   m[i].transform(shape1), mhull1[i]);
    out[1] = out[2] = out[0];
    savep2or3<num_t>((outbase + string("-repl2.ppm")).c_str(), normalize<num_t>(out), false);
    for(int i = 0; i < nemph; i ++) {
      const auto iemph(num_t(i) / num_t(nemph));
      SimpleMatrix<num_t> reref(rin1.rows(), rin1.cols());
      reref.O();
      for(int i = 0; i < m.size(); i ++) {
        const auto rd(draw<num_t>(rin1, shape1,
          takeShape<num_t>(shape1, shape0, ~ m[i], iemph), mhull1[i]));
        for(int j = 0; j < min(reref.rows(), rd.rows()); j ++)
          for(int k = 0; k < min(reref.cols(), rd.cols()); k ++)
            if(rd(j, k) != num_t(0)) reref(j, k) = rd(j, k);
      }
      for(int idx = 0; idx < out.size(); idx ++)
        out[idx] = pullRefMatrix<num_t>(reref, 1 + rin0.rows() * rin0.cols(), in1[idx]);
      savep2or3<num_t>((outbase + string("-") + to_string(i) + string("-") +
                       to_string(nemph) + string(".ppm")).c_str(), out, false);
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
    vector<vector<SimpleMatrix<num_t> > > in;
    in.resize(argc - (strcmp(argv[1], "pred") == 0 ? 4 : 3));
    for(int i = strcmp(argv[1], "pred") == 0 ? 4 : 3; i < argc; i ++) {
      vector<SimpleMatrix<num_t> > ibuf;
      if(!loadp2or3<num_t>(ibuf, argv[i]))
        return - 2;
      in[i - (strcmp(argv[1], "pred") == 0 ? 4 : 3)] = std::move(ibuf);
    }
    const auto idx(in.size() - 1);
    vector<SimpleMatrix<num_t> > out;
    out.resize(3);
    if(strcmp(argv[1], "pred") == 0) {
      for(int i = 0; i < 3; i ++)
        out[i].resize(in[idx][0].rows(), in[idx][0].cols());
      std::vector<std::vector<SimpleMatrix<num_t> > > mout;
      mout.resize(in.size() / 12, out);
      auto nout(mout);
      for(int i = 0; i < out.size(); i ++)
        for(int j = 0; j < out[i].rows(); j ++) {
          SimpleMatrix<num_t> m(in.size(), out[i].cols());
          for(int k = 0; k < in.size(); k ++)
            m.row(k) = std::move(in[k][i].row(j));
          auto ext(filter<num_t>(m, EXTEND_Y, std::atoi(argv[3])));
          for(int k = 0; k < std::atoi(argv[3]); k ++)
            ext += filter<num_t>(m, EXTEND_Y, std::atoi(argv[3]));
          for(int k = 0; k < mout.size(); k ++) {
            mout[k][i].row(j) = std::move(ext.row(k));
            nout[k][i].row(j) = std::move(ext.row(ext.rows() - 1 - k));
          }
        }
      for(int k = 0; k < mout.size(); k ++) {
        savep2or3<num_t>((std::string(argv[2]) + std::string("-p-") + std::to_string(k)).c_str(), normalize<num_t>(mout[k]), ! true, 65535);
        savep2or3<num_t>((std::string(argv[2]) + std::string("-n-") + std::to_string(k)).c_str(), normalize<num_t>(nout[k]), ! true, 65535);
      }
    } else if(strcmp(argv[1], "lenl") == 0) {
      vector<std::pair<SimpleMatrix<num_t>, SimpleMatrix<num_t> > > pair;
      const auto d5(dft<num_t>(5).subMatrix(0, 0, 4, 5));
      const auto d4(dft<num_t>(- 4));
      const auto taylc(d4 * d5);
      const auto tayl(taylc.template real<num_t>());
      for(int i = 0; i < in.size(); i ++)
        for(int j = 0; j < in[i].size(); j ++)
          for(int k = 0; k < in[i][j].rows() - 5; k += 2)
            for(int kk = 0; kk < in[i][j].cols() - 5; kk += 2) {
              const auto work(normalize<num_t>(in[i][j].subMatrix(k, kk, 5, 5), num_t(1) / num_t(2)));
              pair.emplace_back(std::make_pair(normalize<num_t>(tayl * work * tayl.transpose(), num_t(1) / num_t(2)), work));
            }
      out[0] = out[1] = out[2] = optImage<num_t>(pair);
      savep2or3<num_t>(argv[2], out, ! true, 65535);
    } else if(strcmp(argv[1], "cat") == 0 ||
              strcmp(argv[1], "catr") == 0) {
      vector<SimpleMatrix<num_t> > glay;
      glay.reserve(strcmp(argv[1], "catr") == 0 ? in.size()
                    : in.size() * min(int(in[0][0].rows()), 1 + 5 + 1));
      const auto& in00sz(in[0][0].rows());
      for(int i = 0; i < in.size(); i ++) {
        if(strcmp(argv[1], "catr") == 0)
          glay.emplace_back(rgb2d<num_t>(in[i]));
        else {
          auto work(rgb2d<num_t>(in[i]));
          assert(work.rows() == in00sz);
          for(int j = 0; j < min(int(work.rows()), 1 + 5 + 1); j ++) {
            SimpleMatrix<num_t> rr(1, work.cols());
            rr.row(0) = std::move(work.row(j));
            glay.emplace_back(rr);
          }
        }
      }
      const auto cat(catImage<num_t>(glay));
      for(int i = 0; i < cat.size(); i ++) {
        for(int j = 0; j < cat[i].size(); j ++)
          std::cout << argv[3 + (strcmp(argv[1], "catr") == 0 ? cat[i][j]
                         : cat[i][j] / min(in00sz, 1 + 5 + 1)) ] << std::endl;
        std::cout << std::endl;
      }
    } else if(strcmp(argv[1], "composite") == 0) {
      vector<SimpleMatrix<num_t> > glay;
      glay.reserve(in.size());
      for(int i = 0; i < in.size(); i ++)
        glay.emplace_back(rgb2d<num_t>(in[i]));
      const auto composite(compositeImage<num_t>(glay));
      for(int i = 0; i < composite.size(); i ++) {
        out[0] = out[1] = out[2] = composite[i];
        savep2or3<num_t>((string(argv[2]) + string("-") + to_string(i) + string(".ppm")).c_str(), normalize<num_t>(out), ! true);
      }
    }
  } else if(strcmp(argv[1], "habit") == 0) {
    vector<SimpleVector<num_t> > pdst, psrc;
    vector<SimpleVector<int>   > poldst, polsrc;
    if(argc < 5 || !loadobj<num_t>(pdst, poldst, argv[2]) ||
                   !loadobj<num_t>(psrc, polsrc, argv[3])) {
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
    saveobj<num_t>(takeShape<num_t>(pdst, psrc,
      matchPartial<num_t>(pdst, psrc)[0],
      num_t(1) / num_t(2)), My, Mx, poldst, argv[4]);
  } else if(strcmp(argv[1], "retrace") == 0 ||
            strcmp(argv[1], "reimage") == 0) {
    if(argc < 7) {
      usage();
      return - 1;
    }
    vector<SimpleMatrix<num_t> > dst, src, out;
    out.resize(3);
    if(!loadp2or3<num_t>(dst, argv[3]))
      exit(- 2);
    if(!loadp2or3<num_t>(src, argv[4]))
      exit(- 2);
    if(strcmp(argv[1], "retrace") == 0)
      out[0] = out[1] = out[2] =
        reTrace<num_t>(normalize<num_t>(rgb2d<num_t>(dst)),
          normalize<num_t>(rgb2d<num_t>(src)),
          num_t(std::atof(argv[6])), atoi(argv[2]));
    else
      for(int i = 0; i < out.size(); i ++)
        out[i] = reImage<num_t>(dst[i], src[i],
          num_t(std::atof(argv[6])), atoi(argv[2]));
    if(!savep2or3<num_t>(argv[5], normalize<num_t>(out), ! true, 255))
      return - 3;
  } else if(strcmp(argv[1], "retrace2") == 0 ||
            strcmp(argv[1], "reimage2") == 0) {
    if(argc < 6) {
      usage();
      return - 1;
    }
    vector<SimpleMatrix<num_t> > dst, out;
    out.resize(3);
    if(!loadp2or3<num_t>(dst, argv[3]))
      exit(- 2);
    if(strcmp(argv[1], "retrace2") == 0)
      out[0] = out[1] = out[2] =
        reTrace<num_t>(normalize<num_t>(rgb2d<num_t>(dst)),
          num_t(std::atof(argv[5])), atoi(argv[2]));
    else {
      for(int i = 0; i < out.size(); i ++)
        out[i] = reImage<num_t>(dst[i],
          num_t(std::atof(argv[5])), atoi(argv[2]));
      autoLevel<num_t>(out, 4 * (out[0].rows() + out[0].cols()));
    }
    if(!savep2or3<num_t>(argv[4], normalize<num_t>(out), ! true, 255))
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
    // N.B. (f[0] := 1 causes flat result.):
    f[0] = num_t(0);
    for(int i = 1; i < m.size(); i ++) {
      m[i] = m[i - 1] + (num_t(arc4random_uniform(3000)) - num_t(1500)) / num_t(1500);
      n[i] = n[i - 1] + (num_t(arc4random_uniform(3000)) - num_t(1500)) / num_t(1500);
      f[i] = num_t(1);
    }
    f /= sqrt(f.dot(f));
    const auto pp(make_pair(make_pair(atoi(argv[3]),
      atoi(argv[3])), make_pair(0, 0)));
    vector<SimpleMatrix<num_t> > M;
    M.resize(3);
    M[0] = M[1] = M[2] = normalize<num_t>(applyTrace<num_t>(
      make_pair(dec.synth(m, f), dec.synth(n, f)),
      make_pair(pp, pp)) );
    savep2or3<num_t>(argv[4], M, true, 255);
  } else if(strcmp(argv[1], "omake") == 0) {
    vector<vector<num_t> > data;
    string header;
    loaddat<num_t>(argv[3], header, data);
    SimpleMatrix<num_t> buf(atoi(argv[4]), atoi(argv[4]));
    const auto mdft(dft<num_t>(buf.rows()));
    const auto midft(dft<num_t>(- buf.rows()));
    for(int i0 = 1; i0 < data.size(); i0 ++) {
      for(int i = 0; i <= data[i0].size() / buf.rows() / buf.rows(); i ++) {
        for(int k = 0; k < buf.cols(); k ++)
          for(int j = 0; j < buf.rows(); j ++) {
            const auto idx(i * buf.rows() * buf.rows() + k * buf.rows() + j);
            buf(j, k) = idx < data[i0].size() ? data[i0][idx] : num_t(0);
          }
        SimpleMatrix<num_t> buf2;
        if(strcmp(argv[2], "diff") == 0)
          buf2 = (midft * (
            filter<num_t>(mdft.template real<num_t>() * buf, COLLECT_BOTH).template cast<complex<num_t> >() +
            filter<num_t>(mdft.template imag<num_t>() * buf, COLLECT_BOTH).template cast<complex<num_t> >() * complex<num_t>(num_t(0), num_t(1)) ) ).template real<num_t>();
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

