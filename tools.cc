#include <cstdio>
#include <cmath>
#include <iostream>

#if defined(_WITHOUT_EIGEN_)
#include "simplelin.hh"
#else
#include <Eigen/Core>
#endif

#include "fileio.hh"
#include "enlarge.hh"
#include "match.hh"
#include "redig.hh"

using namespace std;

// XXX: configure me:
const int    nshow(4);
const int    nshowh(16);
const int    nmatchhid(20);
const int    nemph(8);
const double Memph(.25);
const int    vbox(16);
const double aroffset(.01);
const double tiltrots(.001);
const double tiltrote(.1);
const double psi(.025);
const int    Mpoly(8000);

void usage() {
  cout << "Usage:" << endl;
  cout << "./tools enlarge <ratio>  <input.ppm> <output.ppm>" << endl;
  cout << "./tools pextend <pixels> <input.ppm> <output.ppm>" << endl;
  cout << "./tools collect <input.ppm> <output.ppm>" << endl;
  cout << "./tools idetect <input.ppm> <output.ppm>" << endl;
  cout << "./tools bump    <input.ppm> <output.ppm>" << endl;
  cout << "./tools obj     <shift_x_pixels> <gather_pixels> <zratio> <input.ppm> <mask.ppm>? <output.obj>" << endl;
  cout << "./tools obj     stand <gather_pixels> <thin> <ratio> <zratio> <input.ppm> <mask.ppm>? <output.obj>" << endl;
  cout << "./tools tilt    <index> <max_index> <psi> <shift_x_pixels> <input.ppm> <input-bump.(ppm|obj)> <output.ppm>" << endl;
  cout << "./tools draw    <input-mask.ppm> <input-obj.(obj|gltf)> <output.ppm>" << endl;
  cout << "./tools match   <inp..." << endl;
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
    bool   addstand(false);
    double thin(1);
    double ratio(1);
    double zratio(1);
    int    sidx(0);
    if(strcmp(argv[2], "stand") == 0) {
      if(argc < 9) {
        usage();
        return - 1;
      }
      addstand = true;
      vbox   = std::atoi(argv[3]);
      thin   = std::atof(argv[4]);
      ratio  = std::atof(argv[5]);
      if(!file.loadp2or3(data, argv[7]))
        return - 1;
      zratio = std::atof(argv[6]) * sqrt(double(data[0].rows() * data[0].cols()));
      if(9 < argc) {
        if(!file.loadp2or3(mask, argv[8]))
          return - 1;
        sidx = 9;
      } else
        sidx = 8;
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
    auto edges(redig.getEdges(mask[0], points));
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
    if(0 < mask[0].rows() && 0 < mask[0].cols())
      file.saveobj(points, facets, argv[sidx], addstand, false, edges, 2., xoffset);
    else
      file.saveobj(points, facets, argv[sidx], addstand, false, edges, 2., xoffset);
  } else if(strcmp(argv[1], "tilt") == 0) {
    if(argc < 9) {
      usage();
      return - 1;
    }
    int    index(std::atoi(argv[2]));
    int    Mindex(std::atoi(argv[3]));
    double psi(std::atof(argv[4]));
    int    offsetx(std::atoi(argv[5]));
    typename simpleFile<double>::Mat data[3], bump[3], out[3];
    std::vector<typename simpleFile<double>::Vec3>  points;
    std::vector<typename simpleFile<double>::Veci3> polys;
    if(!file.loadp2or3(data, argv[6]))
      return - 2;
    const std::string fn(argv[7]);
    bool is_obj(false);
    if(fn[fn.size() - 1] == 'm') {
      if(!file.loadp2or3(bump, argv[7]))
        return - 2;
    } else if(fn[fn.size() - 1] == 'j') {
      if(!file.loadobj(points, polys, argv[7]))
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
    if(!file.savep2or3(argv[8], out, ! true))
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
    match_t<double> m;
    m.offset[0] += min(data[0].rows(), data[0].cols()) / 2.;
    m.offset[1] += min(data[0].rows(), data[0].cols()) / 2.;
    m.ratio     *= min(data[0].rows(), data[0].cols()) / 2.;
    for(int i = 0; i < datapoly.size(); i ++)
      datapoly[i] = m.transform(datapoly[i]);
    typename simpleFile<double>::Mat res[3];
    std::vector<int> idx;
    for(int j = 0; j < datapoly.size(); j ++)
      idx.push_back(j);
    res[0] = res[1] = res[2] = redig.showMatch(data[0] * 0., datapoly, polynorms, double(120));
    file.savep2or3(argv[4], res, true);
  } else if(strcmp(argv[1], "match") == 0) {
    typename simpleFile<double>::Mat data[3], data1[3], bdata[3], bdata1[3], mdata[3], mdata1[3], mout[3], mmout1, bump1;
    if(!file.loadp2or3(data, argv[3]))
      return - 2;
    if(!file.loadp2or3(data1, argv[4]))
      return - 2;
    if(!file.loadp2or3(bdata, argv[5]))
      return - 2;
    if(!file.loadp2or3(bdata1, argv[6]))
      return - 2;
    if(!file.loadp2or3(mdata, argv[7]))
      return - 2;
    if(!file.loadp2or3(mdata1, argv[8]))
      return - 2;
    std::vector<double> emph;
    for(int i = 0; i <= nemph; i ++)
      emph.push_back(double(i) / nemph * Memph);
    resizeDst2<double>(mout, bump1, mmout1, data1, bdata1[0], mdata1[0], data[0].rows(), data[0].cols());
    std::vector<typename simpleFile<double>::Veci3> delau0, delau1;
    std::vector<typename simpleFile<double>::Vec3>  shape0, shape1;
    auto& bump0(bdata[0]);
    redig.getTileVec(bump0, shape0, delau0);
    redig.getTileVec(bump1, shape1, delau1);
    redig.maskVectors(shape0, delau0, mdata[0]);
    redig.maskVectors(shape1, delau1, mmout1);
    matchPartial<double> statmatch;
    //if(mode != 9)
      statmatch.ndiv = statmatch.ndiv * 1.2;
    auto matches(statmatch.match(shape0, shape1));
    matches.resize(min(int(matches.size()), nmatchhid));
    matches = statmatch.elim(matches, data, mout, bump1, shape1);
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
      pmatches.resize(min(int(pmatches.size()), nmatchhid));
      pmatches = pstatmatch.elim(pmatches, data, mout, bump1, shape1a);
      for(int m = 0; m < min(int(pmatches.size()), nshow); m ++) {
        std::cerr << "Writing " << m << " / " << pmatches.size() << " - " << n << " / " << matches.size();
        saveMatches<double>(std::string(argv[3]) + std::to_string(n + 1) + std::string("-") + std::to_string(m + 1), pmatches[m], shape0a, shape1a, data, mout, bump0, bump1, emph);
      }
      std::cerr << "Writing " << n << " / " << matches.size();
      saveMatches<double>(std::string(argv[3]) + std::to_string(n + 1), matches[n], shape0, shape1, data, mout, bump0, bump1, emph);
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
    if(argc > 9) {
      auto m(redig.tiltprep(typename simpleFile<double>::Mat(int(My), int(Mx)), std::atoi(argv[6]), std::atoi(argv[7]), std::atof(argv[8])));
      statmatch.complementMatch(m, pdst, psrc);
      file.saveobj(redig.takeShape(pdst, psrc, m, poldst, polsrc, double(.5)),
                   poldst, argv[9]);
    } else {
      const auto m(statmatch.match(pdst, psrc)[0]);
      file.saveobj(redig.takeShape(pdst, psrc, m, poldst, polsrc, double(.5)),
                   poldst, argv[9]);
    }
  } else {
    usage();
    return - 1;
  }
  return 0;
/*
      // 3d - 2d match.
      std::vector<typename simpleFile<double>::Vec3> datapoly;
      std::vector<typename simpleFile<double>::Veci3> polynorms;
      typename simpleFile<double>::Mat bump0[3], mask0[3];
      if(!file.loadobj(datapoly, polynorms, argv[4]))
        return - 2;
      if(!file.loadp2or3(bump0, argv[5]))
        return - 2;
      if(!file.loadp2or3(mask0, argv[6]))
        return - 2;
      if(datapoly.size() > Mpoly) {
    matchPartial<double> statmatch;
    const auto match(statmatch.match(pdst, psrc));
    for(int i = 0; i < nshow; i ++) {
      const auto& m(match[i]);
      file.saveobj(redig.takeShape(pdst, psrc, m, poldst, polsrc, double(.5)),
                   poldst, (argv[3] + std::string("-emph-") + to_string(i) +
                                      std::string(".obj")).c_str());
    }
      auto m(redig.tiltprep(in, std::atoi(argv[6]), std::atoi(argv[7]), - psi2));
      statmatch.complementMatch(m, pdst, psrc);
      const double psi2(std::atof(argv[8]));
      double Mx(0), My(0);
      for(int i = 0; i < pdst.size(); i ++) {
        My = max(My, std::abs(pdst[i][0]));
        Mx = max(Mx, std::abs(pdst[i][1]));
      }
      for(int i = 0; i < psrc.size(); i ++) {
        My = max(My, std::abs(psrc[i][0]));
        Mx = max(Mx, std::abs(psrc[i][1]));
      }
      typename simpleFile<double>::Mat in(int(My + 1.), int(Mx + 1.));
      matchPartial<double> statmatch;
      auto m(redig.tiltprep(in, std::atoi(argv[6]), std::atoi(argv[7]), - psi2));
      statmatch.complementMatch(m, pdst, psrc);
      file.saveobj(redig.takeShape(pdst, psrc, m, poldst, polsrc, double(.5)),
                   poldst, (argv[3] + std::string("-emph") +
                                      std::string(".obj")).c_str());
      // 3d - 2d match.
      std::vector<typename simpleFile<double>::Vec3> datapoly;
      std::vector<typename simpleFile<double>::Veci3> polynorms;
      typename simpleFile<double>::Mat bump0[3], mask0[3];
      if(!file.loadobj(datapoly, polynorms, argv[4]))
        return - 2;
      if(!file.loadp2or3(bump0, argv[5]))
        return - 2;
      if(!file.loadp2or3(mask0, argv[6]))
        return - 2;
      if(datapoly.size() > Mpoly) {
        std::cerr << "Too many vertices." << std::endl;
        return - 2;
      }
      std::vector<double> emph;
      for(int i = 0; i <= nemph; i ++)
        emph.push_back(double(i) / nemph * Memph);
      std::vector<typename simpleFile<double>::Vec3> shape;
      std::vector<typename simpleFile<double>::Veci3> delau;
      auto& bump(bump0[0]);
      redig.getTileVec(bump, shape, delau);
      redig.maskVectors(shape, delau, mask0[0]);
      typename simpleFile<double>::Mat zero[3];
      for(int i = 0; i < 3; i ++)
        zero[i] = bump * double(0);
      matchPartial<double> statmatch;
      const auto matches(statmatch.match(shape, datapoly));
      for(int n = 0; n < min(int(matches.size()), nshow); n ++) {
        std::cerr << "Writing " << n << " / " << matches.size();
        auto match(matches[n]);
        match_t<double> m0;
        match.rot    = m0.rot;
        match.offset = m0.offset;
        saveMatches<double>(std::string(argv[3]) + std::to_string(n + 1), match, shape, matches[n].transform(datapoly), data, zero, bump, zero[0], emph);
      }
    }
    return 0;
  case 11:
    {
      // 2d - 2d match with hidden 3d model.
      typename simpleFile<double>::Mat data1[3], bdata[3], bdata1[3], mdata[3], mdata1[3], mout[3], mmout1, bump1;
      std::vector<typename simpleFile<double>::Vec3> datapoly;
      std::vector<typename simpleFile<double>::Veci3> polynorms;
      if(!file.loadp2or3(data1, argv[4]))
        return - 2;
      if(!file.loadp2or3(bdata, argv[5]))
        return - 2;
      if(!file.loadp2or3(bdata1, argv[6]))
        return - 2;
      if(!file.loadp2or3(mdata, argv[7]))
        return - 2;
      if(!file.loadp2or3(mdata1, argv[8]))
        return - 2;
      if(!file.loadobj(datapoly, polynorms, argv[9]))
        return - 2;
      if(datapoly.size() > Mpoly) {
        std::cerr << "Too many vertices." << std::endl;
        return - 2;
      }
      std::vector<double> emph;
      for(int i = 0; i <= nemph; i ++)
        emph.push_back(double(i) / nemph * Memph);
      resizeDst2<double>(mout, bump1, mmout1, data1, bdata1[0], mdata1[0], data[0].rows(), data[0].cols());
      std::vector<typename simpleFile<double>::Veci3> delau0, delau1;
      std::vector<typename simpleFile<double>::Vec3> shape0, shape1;
      auto& bump0(bdata[0]);
      redig.getTileVec(bump0, shape0, delau0);
      redig.getTileVec(bump1, shape1, delau1);
      redig.maskVectors(shape0, delau0, mdata[0]);
      redig.maskVectors(shape1, delau1, mmout1);
      matchPartial<double> statmatch;
      const auto match0(statmatch.match(shape0, datapoly));
      const auto match1(statmatch.match(shape1, datapoly));
      std::vector<match_t<double> > matches;
      for(int n = 0; n < min(int(match0.size()), nshowh); n ++)
        for(int m = 0; m < min(int(match1.size()), nshowh); m ++)
          matches.push_back(match0[n] / match1[m]);
      matches = statmatch.elim(matches, data, mout, bump1, shape1);
      matches = statmatch.elim(matches, data, mout, bump1, shape1);
      for(int n = 0; n < min(int(matches.size()), nshow); n ++) {
        const auto& relmatch(matches[n]);
        std::cerr << "Writing " << n << " / " << matches.size();
        saveMatches<double>(std::string(argv[3]) + std::string("-") + std::to_string(n), relmatch, shape0, shape1, data, mout, bump0, bump1, emph);
      }
    }
    return 0;
  case 15:
    // match 3d with bone.
    {
      std::vector<std::vector<typename simpleFile<double>::Vec3> >  datapoly;
      std::vector<std::vector<typename simpleFile<double>::Veci3> > polynorms;
      std::vector<std::vector<typename simpleFile<double>::Veci4> > bone;
      std::vector<typename simpleFile<double>::Vec3>  center;
      typename simpleFile<double>::Mat bump0[3], mask0[3];
      if(!file.loadglTF(datapoly, polynorms, center, bone, argv[4]))
        return - 2;
      if(!file.loadp2or3(bump0, argv[5]))
        return - 2;
      if(!file.loadp2or3(mask0, argv[6]))
        return - 2;
      std::vector<double> emph;
      for(int i = 0; i <= nemph; i ++)
        emph.push_back(double(i) / nemph * Memph);
      std::vector<typename simpleFile<double>::Vec3>  shape;
      std::vector<typename simpleFile<double>::Veci3> delau;
      auto& bump(bump0[0]);
      redig.getTileVec(bump, shape, delau);
      redig.maskVectors(shape, delau, mask0[0]);
      typename simpleFile<double>::Mat zero[3];
      for(int i = 0; i < 3; i ++)
        zero[i] = bump * double(0);
      matchWhole<double> wholematch;
      const auto matches(wholematch.match(shape, datapoly, center, bone));
      for(int n = 0; n < min(int(matches.size()), nshow); n ++) {
        std::cerr << "Writing " << n << " / " << matches.size();
        for(int m = 0; m < matches[n].size(); m ++)
          if(matches[n][m].dstpoints.size()) {
            auto match(matches[n][m]);
            match_t<double> m0;
            match.rot    = m0.rot;
            match.offset = m0.offset;
            saveMatches<double>(std::string(argv[3]) + std::to_string(n + 1) + std::string("-") + std::to_string(m), match, shape, matches[n][m].transform(datapoly[m]), data, zero, bump, zero[0], emph);
          }
      }
    }
    return 0;
  case 28:
    // match3d3d
    {
      std::vector<typename simpleFile<double>::Vec3>  datapoly0;
      std::vector<typename simpleFile<double>::Veci3> polynorms0;
      std::vector<std::vector<typename simpleFile<double>::Vec3> >  datapoly;
      std::vector<std::vector<typename simpleFile<double>::Veci3> > polynorms;
      std::vector<std::vector<typename simpleFile<double>::Veci4> > bone;
      std::vector<typename simpleFile<double>::Vec3>  center;
      typename simpleFile<double>::Mat bump0[3], mask0[3];
      if(!file.loadglTF(datapoly, polynorms, center, bone, argv[4]))
        return - 2;
      if(!file.loadobj(datapoly0, polynorms0, argv[5]))
        return - 2;
      match_t<double> m;
      m.offset[0] += min(data[0].rows(), data[0].cols()) / 2.;
      m.offset[1] += min(data[0].rows(), data[0].cols()) / 2.;
      m.ratio     *= min(data[0].rows(), data[0].cols()) / 2.;
      for(int i = 0; i < datapoly0.size(); i ++)
        datapoly0[i] = m.transform(datapoly0[i]);
      for(int i = 0; i < datapoly.size(); i ++)
        for(int j = 0; j < datapoly[i].size(); j ++)
          datapoly[i][j] = m.transform(datapoly[i][j]);
      typename simpleFile<double>::Mat zero[3];
      for(int i = 0; i < 3; i ++)
        zero[i] = data[0] * double(0);
      std::vector<double> emph;
      emph.push_back(0.); emph.push_back(1.);
      matchWhole<double> wholematch;
      const auto matches(wholematch.match(datapoly0, datapoly, center, bone));
      for(int n = 0; n < min(int(matches.size()), nshow); n ++) {
        std::cerr << "Writing " << n << " / " << matches.size();
        for(int m = 0; m < matches[n].size(); m ++)
          if(matches[n][m].dstpoints.size())
            saveMatches<double>(std::string(argv[3]) + std::to_string(n + 1) + std::string("-") + std::to_string(m), matches[n][m], datapoly0, datapoly[m], zero, zero, zero[0], zero[0], emph);
      }
    }
*/
}

