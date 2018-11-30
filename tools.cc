#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <assert.h>

#if defined(_WITHOUT_EIGEN_)
#include <complex>
#include "simplelin.hh"
#else
#include <Eigen/Core>
#include <Eigen/LU>
#endif

#if defined(_WITH_GLTF2_)
#include <fx/gltf.h>
#endif

namespace goki {
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
  cout << "gokicheck pextend <pixels> <input.ppm> <output.ppm>" << endl;
  cout << "gokicheck collect <input.ppm> <output.ppm>" << endl;
  cout << "gokicheck idetect <input.ppm> <output.ppm>" << endl;
  cout << "gokicheck bump    <input.ppm> <output.ppm>" << endl;
  cout << "gokicheck bumpe   <input.ppm> <output.ppm>" << endl;
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
  for(int idx = 0; idx < 3; idx ++)
    outs[idx] = redig.replace(in0[idx] * T(0), shape0, ohull0) +
                redig.replace(in1[idx] * T(0), match.transform(shape1), ohull1);
  file.savep2or3((outbase + std::string("-match.ppm")).c_str(), outs, false);
  const auto reref(redig.emphasis(rin0, rin1, bump1, shape0, shape1,
                                  match, ohull0, ohull1, emph));
  const auto mhull0(redig.delaunay2(shape0, match.dstpoints));
  const auto mhull1(match.hull(match.srcpoints, match.reverseHull(match.dstpoints, mhull0)));
  for(int idx = 0; idx < 3; idx ++)
    outs[idx] = redig.showMatch(in0[idx], shape0, mhull0) * (T(1) - emph) +
      redig.showMatch(redig.pullRefMatrix(reref, 1 + rin0.rows() * rin0.cols(),
                        in1[idx]), shape1, mhull1) * emph;
  file.savep2or3((outbase + std::string("-cover.ppm")).c_str(), outs, false);
  file.saveobj(redig.takeShape(shape0, shape1,   match, mhull0, mhull1, emph),
               ohull0, (outbase + std::string("-emph.obj")).c_str());
  file.saveobj(redig.takeShape(shape1, shape0, ~ match, mhull1, mhull0, emph),
               ohull1, (outbase + std::string("-emphr.obj")).c_str());
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
  simpleFile<double> file;
  reDig<double>      redig;
  if(strcmp(argv[1], "enlarge") == 0 ||
     strcmp(argv[1], "pextend") == 0) {
    if(argc < 5) {
      usage();
      return 0;
    }
    const auto ratio(std::atoi(argv[2]));
    typename simpleFile<double>::Mat data[3];
    if(!file.loadp2or3(data, argv[3]))
      return - 1;
    if(strcmp(argv[1], "enlarge") == 0) {
      enlarger2ex<double> enlarger, denlarger;
      for(int j = 0; j < ratio; j ++)
        for(int i = 0; i < 3; i ++) {
          const auto xye(enlarger.compute(data[i], enlarger.ENLARGE_BOTH));
          data[i] = xye + redig.tilt45(denlarger.compute(redig.tilt45(data[i], false), denlarger.ENLARGE_BOTH), true, xye) * sqrt(xye.rows() * xye.cols()) / double(xye.rows() + xye.cols());
        }
    } else if(strcmp(argv[1], "pextend") == 0) {
      enlarger2ex<double> extender;
      for(int i = 0; i < 3; i ++)
        for(int j = 0; j < ratio; j ++)
          data[i] = extender.compute(data[i], extender.EXTEND_Y);
/*
      typename simpleFile<double>::Mat result[3];
      for(int j = 0; j < 3; j ++) {
        result[j] = typename simpleFile<double>::Mat(data[j].rows() + ratio * 2, data[j].cols());
        for(int k = 0; k < result[j].rows(); k ++)
          for(int l = 0; l < result[j].cols(); l ++)
            result[j](k, l) = 0.;
        for(int k = 0; k < data[j].rows(); k ++)
          result[j].row(k + ratio) = data[j].row(k);
      }
      for(int j = 0; j < ratio; j ++) {
        for(int i = 0; i < 3; i ++) {
          typename simpleFile<double>::Mat work(data[i].rows() / (j + 1), data[i].cols());
          for(int k = 0; k < data[i].rows() / (j + 1); k ++)
            work.row(k) = data[i].row(k * (j + 1));
          result[i].row(ratio - j - 1) = extender.compute(extender.compute(work, extender.REVERSE_Y), extender.EXTEND_Y0).row(work.rows());
          work = typename simpleFile<double>::Mat(data[i].rows() / (j + 1), data[i].cols());
          for(int k = 0; k < data[i].rows() / (j + 1); k ++)
            work.row(k) = data[i].row(k * (j + 1) - (data[i].rows() / (j + 1) - 1) * (j + 1) + data[i].rows() - 1);
          result[i].row(result[i].rows() - ratio + j) = extender.compute(work, extender.EXTEND_Y0).row(work.rows());
        }
      }
      for(int i = 0; i < 3; i ++)
        data[i] = result[i];
*/
    }
    redig.normalize(data, 1.);
    if(!file.savep2or3(argv[4], data, ! true))
      return - 1;
  } else if(strcmp(argv[1], "collect") == 0 ||
            strcmp(argv[1], "idetect") == 0 ||
            strcmp(argv[1], "bump")  == 0 ||
            strcmp(argv[1], "bumpe") == 0) {
    if(argc < 4) {
      usage();
      return 0;
    }
    typename simpleFile<double>::Mat data[3];
    if(!file.loadp2or3(data, argv[2]))
      return - 1;
    if(strcmp(argv[1], "collect") == 0) {
      enlarger2ex<double> detect, ddetect;
      for(int i = 0; i < 3; i ++) {
        const auto xye(detect.compute(data[i], detect.COLLECT_BOTH));
        data[i] = xye + redig.tilt45(ddetect.compute(redig.tilt45(data[i], false), ddetect.COLLECT_BOTH), true, xye);
      }
    } else if(strcmp(argv[1], "idetect") == 0) {
      enlarger2ex<double> idetect, didetect;
      for(int i = 0; i < 3; i ++) {
        const auto xye(idetect.compute(data[i], idetect.IDETECT_BOTH));
        data[i] = xye + redig.tilt45(didetect.compute(redig.tilt45(data[i], false), didetect.IDETECT_BOTH), true, xye);
      }
    } else if(strcmp(argv[1], "bump") == 0) {
      enlarger2ex<double> bump;
      auto xye(bump.compute(redig.rgb2d(data), bump.BUMP_BOTH));
      data[0] = data[1] = data[2] = redig.autoLevel(xye + redig.tilt45(bump.compute(redig.tilt45(redig.rgb2d(data), false), bump.BUMP_BOTH), true, xye), (xye.rows() + xye.cols()) * 2);
    } else if(strcmp(argv[1], "bumpe") == 0) {
      enlarger2ex<double> bump;
      const auto data0(bump.compute(redig.normalize(redig.rgb2d(data), 1.), bump.DEDGE));
      auto xye(bump.compute(data0, bump.BUMP_BOTH));
      data[0] = data[1] = data[2] = redig.autoLevel(xye + redig.tilt45(bump.compute(redig.tilt45(data0, false), bump.BUMP_BOTH), true, xye), (xye.rows() + xye.cols()) * 2);
    }
    redig.normalize(data, 1.);
    if(!file.savep2or3(argv[3], data, ! true))
      return - 1;
  } else if(strcmp(argv[1], "obj") == 0) {
    typename simpleFile<double>::Mat data[3], mask[3];
    double xoffset(0);
    int    vbox(2);
    double addstand(0);
    double ratio(1);
    double zratio(1);
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
        mask[0] = mask[1] = mask[2] = data[0] * 0.;
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
      zratio  = std::atof(argv[4]) * sqrt(double(data[0].rows() * data[0].cols()));
      if(7 < argc) {
        if(!file.loadp2or3(mask, argv[6]))
          return - 1;
        sidx = 7;
      } else
        sidx = 6;
    }
    std::vector<typename simpleFile<double>::Vec3>  points;
    std::vector<typename simpleFile<double>::Veci3> facets;
    redig.initialize(vbox);
    redig.getTileVec(data[0], points, facets);
    const auto edges(redig.getEdges(mask[0], points));
    if(edges.size())
      redig.maskVectors(points, facets, mask[0]);
    for(int i = 0; i < points.size(); i ++)
      points[i] *= ratio;
    double M(points[0][2]), m(points[0][2]);
    for(int i = 1; i < points.size(); i ++) {
      M = std::max(points[i][2], M);
      m = std::min(points[i][2], m);
    }
    if(M == m) M += 1.;
    for(int i = 0; i < points.size(); i ++)
      points[i][2] *= zratio / (M - m);
    file.saveobj(points, facets, argv[sidx], edges, addstand, xoffset);
  } else if(strcmp(argv[1], "tilt") == 0 ||
            strcmp(argv[1], "sbox") == 0) {
    if((strcmp(argv[1], "tilt") == 0 && argc < 9) ||
       (strcmp(argv[1], "sbox") == 0 && argc < 7)) {
      usage();
      return - 1;
    }
    int    index(std::atoi(argv[2]));
    int    Mindex(std::atoi(argv[3]));
    double psi(0);
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
    typename simpleFile<double>::Mat data[3], bump[3], out[3];
    std::vector<typename simpleFile<double>::Vec3>  points;
    std::vector<typename simpleFile<double>::Veci3> polys;
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
    typename simpleFile<double>::Mat tilt0;
    if(strcmp(argv[1], "sbox") == 0) {
      const match_t<double> mtilt;
      if(is_obj)
        tilt0 = redig.tilt(redig.makeRefMatrix(data[0], 1), redig.tiltprep(points, polys, redig.makeRefMatrix(data[0], 1), mtilt), index / (double)Mindex);
      else
        tilt0 = redig.tilt(redig.makeRefMatrix(data[0], 1), bump[0], mtilt, index / (double)Mindex);
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
    typename simpleFile<double>::Mat data[3];
    if(!file.loadp2or3(data, argv[2]))
      return - 2;
    std::vector<typename simpleFile<double>::Vec3>  datapoly;
    std::vector<typename simpleFile<double>::Veci3> polynorms;
    const std::string fn(argv[3]);
    if(fn[fn.size() - 1] == 'j') {
      if(!file.loadobj(datapoly, polynorms, argv[3]))
        return - 2;
    } else if(fn[fn.size() - 1] == 'f') {
      std::vector<typename simpleFile<double>::Vec3> center;
      std::vector<std::vector<typename simpleFile<double>::Veci4> > bone;
      if(!file.loadglTF(datapoly, polynorms, center, bone, argv[3]))
        return - 2;
    } else {
      usage();
      return - 1;
    }
    assert(datapoly.size());
    double My(datapoly[0][0]);
    double Mx(datapoly[0][1]);
    double my(My);
    double mx(Mx);
    for(int i = 0; i < datapoly.size(); i ++) {
      My = max(My, datapoly[i][0]);
      Mx = max(Mx, datapoly[i][1]);
      my = min(my, datapoly[i][0]);
      mx = min(mx, datapoly[i][1]);
    }
    match_t<double> m;
    const double ratio(sqrt(double(data[0].rows() * data[0].cols())) / max(My - my, Mx - mx));
    m.ratio     *= ratio;
    m.offset[0] -= my * ratio;
    m.offset[1] -= mx * ratio;
    for(int i = 0; i < datapoly.size(); i ++)
      datapoly[i] = m.transform(datapoly[i]);
    typename simpleFile<double>::Mat res[3];
    std::vector<int> idx;
    for(int j = 0; j < datapoly.size(); j ++)
      idx.push_back(j);
    if(strcmp(argv[1], "draw") == 0)
      res[0] = res[1] = res[2] = redig.showMatch(data[0] * 0., datapoly, polynorms, double(120));
    else if(strcmp(argv[1], "drawr") == 0) {
      auto mwork(data[0]);
      for(int i = 0; i < mwork.rows(); i ++)
        for(int j = 0; j < mwork.cols(); j ++)
          mwork(i, j) = 1.;
      auto prep(redig.tiltprep(datapoly, polynorms, mwork, match_t<double>()));
      for(int i = 0; i < prep.size(); i ++)
        prep[i].c = (prep[i].p(2, 0) + prep[i].p(2, 1) + prep[i].p(2, 2)) / 3.;
      res[0] = res[1] = res[2] = redig.tilt(data[0] * 0., prep);
    } else {
      auto mwork(data[0]);
      for(int i = 0; i < mwork.rows(); i ++)
        for(int j = 0; j < mwork.cols(); j ++)
          mwork(i, j) = 1.;
      res[0] = res[1] = res[2] = mwork - redig.tilt(data[0] * 0., redig.tiltprep(datapoly, polynorms, mwork, match_t<double>()));
    }
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
    float emph(0);
    match_t<double> m;
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
    typename simpleFile<double>::Mat data[3], data1[3], bdata[3], bdata1[3], mdata[3], mdata1[3], mout[3], mmout1, bump1;
    std::vector<typename simpleFile<double>::Veci3> delau0, delau1;
    std::vector<typename simpleFile<double>::Vec3>  shape0, shape1;
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
      bdata1[0] = bdata1[1] = bdata1[2] = data1[0] * 0.;
    } else if(fn[fn.size() - 1] == 'f') {
      std::vector<typename simpleFile<double>::Vec3> center;
      std::vector<std::vector<typename simpleFile<double>::Veci4> > bone;
      if(!file.loadglTF(shape1, delau1, center, bone, argv[9]))
        return - 2;
      bdata1[0] = bdata1[1] = bdata1[2] = data1[0] * 0.;
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
      mdata[0]  = mdata[1]  = mdata[2]  = data[0]  * 0.;
      mdata1[0] = mdata1[1] = mdata1[2] = data1[0] * 0.;
    }
    resizeDst2<double>(mout, bump1, mmout1, data1, bdata1[0], mdata1[0], data[0].rows(), data[0].cols());
    auto& bump0(bdata[0]);
    redig.initialize(vboxdst);
    redig.getTileVec(bump0, shape0, delau0);
    redig.maskVectors(shape0, delau0, mdata[0]);
    redig.initialize(vboxsrc);
    if(fn[fn.size() - 1] == 'm') {
      redig.getTileVec(bump1, shape1, delau1);
      redig.maskVectors(shape1, delau1, mmout1);
    }
    if(strcmp(argv[1], "matcho") == 0)
      saveMatches<double>(std::string(argv[fnout]), m, shape0, shape1, data, mout, bump0, bump1, delau0, delau1, emph);
    else { 
      matchPartial<double> statmatch;
      auto matches(statmatch.match(shape0, shape1));
      // matches = statmatch.elim(matches, data, mout, bump1, shape1);
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
        std::vector<typename simpleFile<double>::Vec3> shape0a;
        std::vector<typename simpleFile<double>::Vec3> shape1a;
        std::sort(dstbuf.begin(), dstbuf.end());
        std::sort(srcbuf.begin(), srcbuf.end());
        for(int j = 0; j < shape0.size(); j ++)
          if(!binary_search(dstbuf.begin(), dstbuf.end(), j))
            shape0a.push_back(shape0[j]);
        for(int j = 0; j < shape1.size(); j ++)
          if(!binary_search(srcbuf.begin(), srcbuf.end(), j))
            shape1a.push_back(shape1[j]);
        matchPartial<double> pstatmatch;
        pstatmatch.ndiv /= 2;
        auto pmatches(pstatmatch.match(shape0a, shape1a));
        // pmatches = pstatmatch.elim(pmatches, data, mout, bump1, shape1a);
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
    std::vector<typename simpleFile<double>::Vec3>  pdst,   psrc;
    std::vector<typename simpleFile<double>::Veci3> poldst, polsrc;
    if(argc < 5 || !file.loadobj(pdst, poldst, argv[2]) ||
                   !file.loadobj(psrc, polsrc, argv[3])) {
      usage();
      return - 2;
    }
    double Mx(0), My(0);
    for(int i = 0; i < pdst.size(); i ++) {
      My = max(My, std::abs(pdst[i][0]));
      Mx = max(Mx, std::abs(pdst[i][1]));
    }
    for(int i = 0; i < psrc.size(); i ++) {
      My = max(My, std::abs(psrc[i][0]));
      Mx = max(Mx, std::abs(psrc[i][1]));
    }
    if(argc > 7) {
      const auto m(redig.tiltprep(typename simpleFile<double>::Mat(int(My), int(Mx)), - std::atoi(argv[4]), std::atoi(argv[5]), std::atof(argv[6])));
      file.saveobj(redig.takeShape(pdst, psrc, m, poldst, polsrc, double(.5)),
                   poldst, argv[7]);
    } else {
      matchPartial<double> statmatch;
      const auto m(statmatch.match(pdst, psrc)[0]);
      file.saveobj(redig.takeShape(pdst, psrc, m, poldst, polsrc, double(.5)),
                   poldst, argv[4]);
    }
  } else if(strcmp(argv[1], "omake") == 0) {
    std::vector<double> in;
    std::string line;
    while (std::getline(std::cin, line)) {
      std::stringstream sl(line);
      in.push_back(0);
      sl >> in[in.size() - 1];
    }
    simpleFile<double>::Mat buf(std::atoi(argv[2]), std::atoi(argv[2]));
    enlarger2ex<double> filter;
    double M(0);
    std::vector<double> rbuf;
    if(strcmp(argv[3], "enlarge") == 0)
      rbuf.resize(in.size() * 2, 0.);
    else
      rbuf.resize(in.size(), 0.);
    const auto dft0(filter.seed(buf.rows(), false));
    const auto idft0(filter.seed(buf.rows(), true));
    for(int i = 0; i <= in.size() / buf.rows() / buf.rows(); i ++) {
      for(int j = 0; j < buf.rows(); j ++)
        for(int k = 0; k < buf.cols(); k ++)
          buf(j, k) = in[i * buf.rows() * buf.rows() + k * buf.rows() + j];
      enlarger2ex<double>::MatU dft(dft0 * buf.template cast<complex<double> >());
      if(strcmp(argv[3], "enlarge") == 0) {
#if defined(_WITHOUT_EIGEN_)
        const auto rp(filter.compute(dft.template real<double>(), filter.ENLARGE_X));
        const auto ip(filter.compute(dft.template imag<double>(), filter.ENLARGE_X));
        const auto buf2((idft0.template cast<complex<double> >() * (rp.template cast<complex<double> >() + ip.template cast<complex<double> >() * sqrt(complex<double>(- 1.)))).template real<double>());
#else
        const auto rp(filter.compute(dft.real(), filter.ENLARGE_X));
        const auto ip(filter.compute(dft.imag(), filter.ENLARGE_X));
        const auto buf2((idft0.template cast<complex<double> >() * (rp.template cast<complex<double> >() + ip.template cast<complex<double> >() * sqrt(complex<double>(- 1.)))).real());
#endif
        for(int j = 0; j < buf2.rows(); j ++)
          for(int k = 0; k < buf2.cols(); k ++) {
            M = max(M, abs(buf2(j, k)));
            rbuf[i * buf.rows() * buf.rows() * 2 + k * buf.rows() + j] = buf2(j, k);
          }
        continue;
      } else {
        enlarger2ex<double>::Mat buf2;
        if(strcmp(argv[3], "diff") == 0){
#if defined(_WITHOUT_EIGEN_)
          const auto rp(filter.compute(dft.template real<double>(), filter.DETECT_X));
          const auto ip(filter.compute(dft.template imag<double>(), filter.DETECT_X));
          buf2 = (idft0.template cast<complex<double> >() * (rp.template cast<complex<double> >() + ip.template cast<complex<double> >() * sqrt(complex<double>(- 1)))).template real<double>();
#else
          const auto rp(filter.compute(dft.real(), filter.DETECT_X));
          const auto ip(filter.compute(dft.imag(), filter.DETECT_X));
          buf2 = (idft0 * (rp.template cast<complex<double> >() + sqrt(complex<double>(- 1.)) * ip.template cast<complex<double> >())).real();
#endif
        } else if(strcmp(argv[3], "integ") == 0){
#if defined(_WITHOUT_EIGEN_)
          const auto rp(filter.compute(dft.template real<double>(), filter.IDETECT_X));
          const auto ip(filter.compute(dft.template imag<double>(), filter.IDETECT_X));
          buf2 = (idft0.template cast<complex<double> >() * (rp.template cast<complex<double> >() + ip.template cast<complex<double> >() * sqrt(complex<double>(- 1)))).template real<double>();
#else
          const auto rp(filter.compute(dft.real(), filter.IDETECT_X));
          const auto ip(filter.compute(dft.imag(), filter.IDETECT_X));
          buf2 = (idft0 * (rp.template cast<complex<double> >() + sqrt(complex<double>(- 1.)) * ip.template cast<complex<double> >())).real();
#endif
        } else if(strcmp(argv[3], "bump") == 0) {
#if defined(_WITHOUT_EIGEN_)
          const auto rp(filter.compute(dft.template real<double>(), filter.BUMP_X));
          const auto ip(filter.compute(dft.template imag<double>(), filter.BUMP_X));
          buf2 = (idft0.template cast<complex<double> >() * (rp.template cast<complex<double> >() + ip.template cast<complex<double> >() * sqrt(complex<double>(- 1)))).template real<double>();
#else
          const auto rp(filter.compute(dft.real(), filter.BUMP_X));
          const auto ip(filter.compute(dft.imag(), filter.BUMP_X));
          buf2 = (idft0 * (rp.template cast<complex<double> >() + sqrt(complex<double>(- 1.)) * ip.template cast<complex<double> >())).real();
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

