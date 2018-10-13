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
#include "simplelin.hh"
#include <complex>
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
  cout << "gokicheck obj     <shift_x_pixels> <gather_pixels> <zratio> <input.ppm> <mask.ppm>? <output.obj>" << endl;
  cout << "gokicheck obj     stand <gather_pixels> <thin> <ratio> <zratio> <input.ppm> <mask.ppm>? <output.obj>" << endl;
  cout << "gokicheck tilt    <index> <max_index> <psi> <shift_x_pixels> <input.ppm> <input-bump.(ppm|obj)> <output.ppm>" << endl;
  cout << "gokicheck sbox    <index> <max_index> <input.ppm> <input-bump.(ppm|obj)> <output.ppm>" << endl;
  cout << "gokicheck draw    <input-mask.ppm> <input-obj.(obj|gltf)> <output.ppm>" << endl;
  cout << "gokicheck match   <num_of_res_shown> <num_of_hidden_match> <num_emph> <vbox_dst> <vbox_src> <dst.ppm> <src.ppm> <dst-bump.(ppm|obj)> <src-bump.(ppm|obj|gltf)> (<dst-mask.ppm> <src-mask.ppm>)? <output-basename>" << endl;
  cout << "gokicheck habit   <in0.obj> <in1.obj> (<index> <max_index> <psi>)? <out.obj>" << endl;
  return;
}

template <typename T> void saveMatches(const std::string& outbase, const match_t<T>& match, const std::vector<typename simpleFile<T>::Vec3>& shape0, const std::vector<typename simpleFile<T>::Vec3>& shape1, const typename simpleFile<T>::Mat in0[3], const typename simpleFile<T>::Mat in1[3], const typename simpleFile<T>::Mat& bump0, const typename simpleFile<T>::Mat& bump1, const std::vector<T>& emph) {
  reDig<T> redig;
  simpleFile<T> file;
  std::cerr << "(" << match.rdepth << ", " << match.ratio << ")" << std::endl;
  
  typename simpleFile<T>::Mat outs[3], outs2[3];
  const auto mhull0(redig.delaunay2(shape0, match.dstpoints));
  const auto mhull1(match.hull(match.srcpoints, match.reverseHull(match.dstpoints, mhull0)));
  
  typename simpleFile<T>::Mat sin1[3];
  for(int idx = 0; idx < 3; idx ++)
    sin1[idx] = redig.showMatch(in1[idx], shape1, mhull1);
  redig.normalize(sin1, 1.);
  const auto tilt0(redig.tilt(redig.makeRefMatrix(in0[0], 1), bump1, match));
  for(int idx = 0; idx < 3; idx ++)
    outs[idx] = redig.pullRefMatrix(tilt0, 1, sin1[idx]);
  std::string outfile(outbase + std::string("-src.ppm"));
  file.savep2or3(outfile.c_str(), outs, false);
  
  for(int idx = 0; idx < 3; idx ++)
    outs[idx] = redig.showMatch(in0[idx], shape0, mhull0);
  redig.normalize(outs, 1.);
  outfile = outbase + std::string("-dst.ppm");
  file.savep2or3(outfile.c_str(), outs, false);
  
  for(int idx = 0; idx < 3; idx ++)
    outs[idx] = redig.replace(in0[idx], match.transform(shape1), match, mhull1);
  outfile = outbase + std::string("-repl.ppm");
  file.savep2or3(outfile.c_str(), outs, false);
  
  const auto rin0(redig.makeRefMatrix(in0[0], 1));
  const auto rin1(redig.makeRefMatrix(in1[0], 1 + rin0.rows() * rin0.cols()));
  for(int kk = 0; kk < emph.size(); kk ++) {
    const auto reref(redig.emphasis(rin0, rin1, bump1, shape0, shape1, match, mhull0, mhull1, emph[kk]));
    for(int idx = 0; idx < 3; idx ++) {
      outs[idx] = (in0[idx] * (emph.size() - 1 - kk) + redig.pullRefMatrix(reref, 1 + rin0.rows() * rin0.cols(), sin1[idx]) * kk) / double(emph.size() - 1);
      outs2[idx] = redig.pullRefMatrix(reref, 1 + rin0.rows() * rin0.cols(), in1[idx]);
    }
    outfile = outbase + std::string("-emph-") + std::to_string(kk) + std::string(".ppm");
    file.savep2or3(outfile.c_str(), outs, false);
    outfile = outbase + std::string("-emph2-") + std::to_string(kk) + std::string(".ppm");
    file.savep2or3(outfile.c_str(), outs2, false);
  }
  
  file.saveobj(redig.takeShape(shape0, shape1, match, mhull0, mhull1,
                               emph[emph.size() - 1]), mhull0,
               (outbase + std::string("-emph.obj")).c_str());
  file.saveobj(redig.takeShape(shape1, shape0, ~ match, mhull1, mhull0,
                               emph[emph.size() - 1]), mhull1,
               (outbase + std::string("-emphr.obj")).c_str());
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
      const int count(sqrt(sqrt(double(data[0].rows() * data[0].cols()))));
      for(int j = 0; j < count; j ++)
        for(int i = 0; i < 3; i ++)
          data[i] = extender.compute(data[i], extender.EXTEND_BOTH);
    }
    redig.normalize(data, 1.);
    if(!file.savep2or3(argv[4], data, ! true))
      return - 1;
  } else if(strcmp(argv[1], "collect") == 0 ||
            strcmp(argv[1], "idetect") == 0 ||
            strcmp(argv[1], "bump") == 0) {
    if(argc < 4) {
      usage();
      return 0;
    }
    typename simpleFile<double>::Mat data[3];
    if(!file.loadp2or3(data, argv[2]))
      return - 1;
    if(strcmp(argv[1], "collect") == 0) {
      // collect.
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
      const auto xye(bump.compute(redig.rgb2d(data), bump.BUMP_BOTH));
      data[0] = data[1] = data[2] = redig.autoLevel(xye + redig.tilt45(bump.compute(redig.tilt45(redig.rgb2d(data), false), bump.BUMP_BOTH), true, xye), 4 * (xye.rows() + xye.cols()));
    }
    redig.normalize(data, 1.);
    if(!file.savep2or3(argv[3], data, ! true))
      return - 1;
  } else if(strcmp(argv[1], "obj") == 0) {
    typename simpleFile<double>::Mat data[3], mask[3];
    int    xoffset(0);
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
      zratio = std::atof(argv[6]) * sqrt(double(data[0].rows() * data[0].cols()));
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
      xoffset = std::atoi(argv[2]);
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
    auto mtilt(redig.tiltprep(data[0], index, Mindex, psi));
    mtilt.offset[2] += offsetx * data[0].cols();
    typename simpleFile<double>::Mat tilt0;
    if(is_obj)
      tilt0 = redig.tilt(redig.makeRefMatrix(data[0], 1), redig.tiltprep(points, polys, redig.makeRefMatrix(data[0], 1), mtilt));
    else
      tilt0 = redig.tilt(redig.makeRefMatrix(data[0], 1), bump[0], mtilt);
    for(int j = 0; j < 3; j ++)
      out[j] = redig.pullRefMatrix(tilt0, 1, data[j]);
    if(!file.savep2or3(argv[oidx], out, ! true))
      return - 1;
  } else if(strcmp(argv[1], "draw") == 0) {
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
    res[0] = res[1] = res[2] = redig.showMatch(data[0] * 0., datapoly, polynorms, double(120));
    file.savep2or3(argv[4], res, true);
  } else if(strcmp(argv[1], "match") == 0) {
    if(argc < 11) {
      usage();
      return - 1;
    }
    int fnout(11);
    int nshow(std::atoi(argv[2]));
    int nhid(std::atoi(argv[3]));
    int nemph(std::atoi(argv[4]));
    int vboxdst(std::atoi(argv[5]));
    int vboxsrc(std::atoi(argv[6]));
    typename simpleFile<double>::Mat data[3], data1[3], bdata[3], bdata1[3], mdata[3], mdata1[3], mout[3], mmout1, bump1;
    std::vector<typename simpleFile<double>::Veci3> delau0, delau1;
    std::vector<typename simpleFile<double>::Vec3>  shape0, shape1;
    if(!file.loadp2or3(data, argv[7]))
      return - 2;
    if(!file.loadp2or3(data1, argv[8]))
      return - 2;
    if(!file.loadp2or3(bdata, argv[9]))
      return - 2;
    const std::string fn(argv[10]);
    if(fn[fn.size() - 1] == 'j') {
      if(!file.loadobj(shape1, delau1, argv[10]))
        return - 2;
      bdata1[0] = bdata1[1] = bdata1[2] = data1[0] * 0.;
    } else if(fn[fn.size() - 1] == 'f') {
      std::vector<typename simpleFile<double>::Vec3> center;
      std::vector<std::vector<typename simpleFile<double>::Veci4> > bone;
      if(!file.loadglTF(shape1, delau1, center, bone, argv[10]))
        return - 2;
      bdata1[0] = bdata1[1] = bdata1[2] = data1[0] * 0.;
    } else if(fn[fn.size() - 1] == 'm') {
      if(!file.loadp2or3(bdata1, argv[10]))
        return - 2;
    } else {
      usage();
      return - 1;
    }
    if(12 < argc) {
      if(!file.loadp2or3(mdata, argv[11]))
        return - 2;
      if(!file.loadp2or3(mdata1, argv[12]))
        return - 2;
      fnout = 13;
    } else {
      mdata[0]  = mdata[1]  = mdata[2]  = data[0]  * 0.;
      mdata1[0] = mdata1[1] = mdata1[2] = data1[0] * 0.;
    }
    resizeDst2<double>(mout, bump1, mmout1, data1, bdata1[0], mdata1[0], data[0].rows(), data[0].cols());
    std::vector<double> emph;
    for(int i = 0; i <= nemph; i ++)
      emph.push_back(double(i) / nemph);
    auto& bump0(bdata[0]);
    redig.initialize(vboxdst);
    redig.getTileVec(bump0, shape0, delau0);
    redig.maskVectors(shape0, delau0, mdata[0]);
    redig.initialize(vboxsrc);
    if(fn[fn.size() - 1] == 'm')
      redig.getTileVec(bump1, shape1, delau1);
    redig.maskVectors(shape1, delau1, mmout1);
    matchPartial<double> statmatch;
    auto matches(statmatch.match(shape0, shape1));
    // matches = statmatch.elim(matches, data, mout, bump1, shape1);
    matches.resize(min(int(matches.size()), nhid));
    for(int n = 0; n < min(int(matches.size()), nshow); n ++) {
      std::cerr << "Matchingsub: " << n << " / " << matches.size();
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
        std::cerr << "Writing " << m << " / " << pmatches.size() << " - " << n << " / " << matches.size();
        saveMatches<double>(std::string(argv[fnout]) + std::to_string(n + 1) + std::string("-") + std::to_string(m + 1), pmatches[m], shape0a, shape1a, data, mout, bump0, bump1, emph);
      }
      std::cerr << "Writing " << n << " / " << matches.size();
      saveMatches<double>(std::string(argv[fnout]) + std::to_string(n + 1), matches[n], shape0, shape1, data, mout, bump0, bump1, emph);
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
    matchPartial<double> statmatch;
    if(argc > 7) {
      const auto m0(redig.tiltprep(typename simpleFile<double>::Mat(int(My), int(Mx)), std::atoi(argv[4]), std::atoi(argv[5]), std::atof(argv[6])));
            auto m1(redig.tiltprep(typename simpleFile<double>::Mat(int(My), int(Mx)), std::atoi(argv[4]), std::atoi(argv[5]), std::atof(argv[6]) * 2));
      statmatch.complementMatch(m1, pdst, psrc);
      file.saveobj(m0.transform(redig.takeShape(pdst, psrc, m1, poldst, polsrc, double(.5))),
                   poldst, argv[7]);
    } else {
      const auto m(statmatch.match(pdst, psrc)[0]);
      file.saveobj(redig.takeShape(pdst, psrc, m, poldst, polsrc, double(.5)),
                   poldst, argv[4]);
    }
  } else {
    usage();
    return - 1;
  }
  return 0;
}

