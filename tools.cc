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
const int    vbox0(2);
const int    vbox(16);
const double rz(1. / 6.);
const double offsetx(.1);
const double aroffset(.01);
const double arrot(.025);
const int    M_TILT(32);
const int    M_TILTROT(8);
const double tiltrots(.001);
const double tiltrote(.1);
const double psi(.025);
const double psi2(.1);
const int    Mpoly(2000);

void usage() {
  cout << "Usage: tools (enlarge|enlarge4|pextend|collect|idetect|bump|obj|arobj|bump2|rbump2|tilt|tilt2|tilt3|tilt4|tilt5|tiltp|match|match3d|match3dbone|match2dh3d|match2dh3dbone|pmatch|pmatch3d|match3d3d|maskobj|habit|habit2|drawobj|drawgltf) <input filename>.p[gp]m <output filename>.p[gp]m <args>?" << endl;
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
    outs[idx] = redig.replace(in0[idx], shape1, match, mhull1);
  outfile = outbase + std::string("-repl.ppm");
  file.savep2or3(outfile.c_str(), outs, false);
  
  const auto rin0(redig.makeRefMatrix(in0[0], 1));
  const auto rin1(redig.makeRefMatrix(in1[0], 1 + rin0.rows() * rin0.cols()));
  for(int kk = 0; kk < emph.size(); kk ++) {
    const auto reref(redig.emphasis(rin0, rin1, bump1, shape0, shape1, match, mhull0, mhull1, emph[kk]));
    for(int idx = 0; idx < 3; idx ++) {
      outs[idx] = (in0[idx] * (emph.size() - 1 - kk) + redig.pullRefMatrix(reref, 1 + rin0.rows() * rin0.cols(), sin1[idx]) * kk) / double(emph.size() - 1);
      outs2[idx] = (in0[idx] * (emph.size() - 1 - kk) + redig.pullRefMatrix(reref, 1 + rin0.rows() * rin0.cols(), in1[idx]) * kk) / double(emph.size() - 1);
/*
      outs2[idx] = in0[idx];
      const auto rework(redig.pullRefMatrix(reref, 1 + rin0.rows() * rin0.cols(), in1[idx]));
      for(int i = 0; i < reref.rows(); i ++)
        for(int j = 0; j < reref.cols(); j ++)
          if(rin0.rows() * rin0.cols() + 1 < reref(i, j))
            outs2[idx](i, j) = rework(i, j);
*/
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
  if(argc < 4) {
    usage();
    return 0;
  }
  int mode = - 1;
  if(strcmp(argv[1], "enlarge") == 0)
    mode = 0;
  else if(strcmp(argv[1], "enlarge4") == 0)
    mode = 18;
  else if(strcmp(argv[1], "pextend") == 0)
    mode = 23;
  else if(strcmp(argv[1], "bump") == 0)
    mode = 2;
  else if(strcmp(argv[1], "obj") == 0)
    mode = 7;
  else if(strcmp(argv[1], "bump2") == 0)
    mode = 3;
  else if(strcmp(argv[1], "rbump2") == 0)
    mode = 5;
  else if(strcmp(argv[1], "collect") == 0)
    mode = 4;
  else if(strcmp(argv[1], "idetect") == 0)
    mode = 1;
  else if(strcmp(argv[1], "tilt") == 0)
    mode = 6;
  else if(strcmp(argv[1], "match") == 0)
    mode = 9;
  else if(strcmp(argv[1], "match3d") == 0)
    mode = 10;
  else if(strcmp(argv[1], "match3dbone") == 0)
    mode = 15;
  else if(strcmp(argv[1], "match2dh3d") == 0)
    mode = 11;
  else if(strcmp(argv[1], "match2dh3dbone") == 0)
    mode = 16;
  else if(strcmp(argv[1], "pmatch") == 0)
    mode = 25;
  else if(strcmp(argv[1], "pmatch3d") == 0)
    mode = 26;
  else if(strcmp(argv[1], "match3d3d") == 0)
    mode = 28;
  else if(strcmp(argv[1], "maskobj") == 0)
    mode = 12;
  else if(strcmp(argv[1], "maskobj2") == 0)
    mode = 27;
  else if(strcmp(argv[1], "arobj") == 0)
    mode = 17;
  else if(strcmp(argv[1], "habit") == 0)
    mode = 13;
  else if(strcmp(argv[1], "tilt2") == 0)
    mode = 14;
  else if(strcmp(argv[1], "tilt3") == 0)
    mode = 8;
  else if(strcmp(argv[1], "tilt4") == 0)
    mode = 24;
  else if(strcmp(argv[1], "tilt5") == 0)
    mode = 20;
  else if(strcmp(argv[1], "tiltp") == 0)
    mode = 22;
  else if(strcmp(argv[1], "habit2") == 0)
    mode = 19;
  else if(strcmp(argv[1], "drawobj") == 0)
    mode = 29;
  else if(strcmp(argv[1], "drawgltf") == 0)
    mode = 21;
  if(mode < 0) {
    usage();
    return - 1;
  }
  reDig<double> redig;
  simpleFile<double> file;
  typename simpleFile<double>::Mat data[3];
  if(!file.loadp2or3(data, argv[2]))
    return - 2;
  redig.initialize(vbox, rz);
  redig.normalize(data, 1.);
  switch(mode) {
  case 0:
    {
      // enlarge.
      enlarger2ex<double> enlarger, denlarger;
      for(int i = 0; i < 3; i ++) {
        const auto xye(enlarger.compute(data[i], enlarger.ENLARGE_BOTH));
        data[i] = xye + redig.tilt45(denlarger.compute(redig.tilt45(data[i], false), denlarger.ENLARGE_BOTH), true, xye) * sqrt(xye.rows() * xye.cols()) / double(xye.rows() + xye.cols());
      }
    }
    break;
  case 18:
    {
      // enlarge4.
      enlarger2ex<double> enlarger, denlarger;
      for(int j = 0; j < 4; j ++)
        for(int i = 0; i < 3; i ++) {
          const auto xye(enlarger.compute(data[i], enlarger.ENLARGE_BOTH));
          data[i]  = xye + redig.tilt45(denlarger.compute(redig.tilt45(data[i], false), denlarger.ENLARGE_BOTH), true, xye);
          data[i] /= 2.;
        }
    }
    break;
  case 23:
    {
      // extend.
      enlarger2ex<double> extender, extender2;
      const int count(sqrt(sqrt(double(data[0].rows() * data[0].cols()))));
      for(int j = 0; j < count; j ++)
        for(int i = 0; i < 3; i ++)
          data[i] = redig.reversey(redig.reversey(extender2.compute(redig.reversey(redig.reversey(extender.compute(data[i], extender.EXTEND_BOTH)).transpose()), extender2.EXTEND_BOTH)).transpose());
    }
    break;
  case 4:
    {
      // collect.
      enlarger2ex<double> detect, ddetect;
      for(int i = 0; i < 3; i ++) {
        const auto xye(detect.compute(data[i], detect.COLLECT_BOTH));
        data[i] = xye + redig.tilt45(ddetect.compute(redig.tilt45(data[i], false), ddetect.COLLECT_BOTH), true, xye);
      }
    }
    break;
  case 1:
    {
      // idetect.
      enlarger2ex<double> idetect, didetect;
      for(int i = 0; i < 3; i ++) {
        const auto xye(idetect.compute(data[i], idetect.IDETECT_BOTH));
        data[i] = xye + redig.tilt45(didetect.compute(redig.tilt45(data[i], false), didetect.IDETECT_BOTH), true, xye);
      }
    }
    break;
  case 2:
    {
      // bump.
      enlarger2ex<double> bump;
      const auto xye(bump.compute(redig.rgb2d(data), bump.BUMP_BOTH));
      data[1] = data[2] = data[0] = xye;
      // data[1] = data[2] = data[0] = xye + redig.tilt45(bump.compute(redig.tilt45(redig.rgb2d(data), false), bump.BUMP_BOTH), true, xye);
    }
    break;
  case 7:
  case 17:
    {
      // obj.
      std::vector<typename simpleFile<double>::Vec3> points;
      std::vector<typename simpleFile<double>::Veci3> facets;
      redig.initialize(vbox0, rz);
      redig.getTileVec(data[0], points, facets);
      if(mode == 7)
        file.saveobj(points, facets, argv[3], false);
      else {
        file.saveobj(points, facets, (std::string(argv[3]) + std::string("-L.obj")).c_str(), true, false, std::vector<std::vector<int> >(), 2., - aroffset);
        file.saveobj(points, facets, (std::string(argv[3]) + std::string("-R.obj")).c_str(), true, false, std::vector<std::vector<int> >(), 2.,   aroffset);
      }
    }
    return 0;
  case 3:
    {
      // bump2 for training output.
      enlarger2ex<double> bump;
      auto X(data[0] * .49000/.17697 + data[1] * .31000/.17697  + data[2] * .20000/.17697);
      auto Z(                          data[1] * .010000/.17697 + data[2] * .99000/.17697);
      data[0] = X / 8.;
      data[1] = bump.compute(redig.rgb2l(data), bump.BUMP_BOTH) / 8.;
      data[2] = Z / 8.;
    }
    // with no auto-level.
    if(!file.savep2or3(argv[3], data, ! true))
      return - 3;
    return 0;
  case 5:
    {
      // reverse bump2.
      auto R(  data[0] * .41847    - data[1] * .15866   - data[2] * .082835);
      auto G(- data[0] * .091169   + data[1] * .25243   + data[2] * .015708);
      auto B(  data[0] * .00092090 - data[1] * .0025498 + data[2] * .17860);
      data[0] = R;
      data[1] = G;
      data[2] = B;
    }
    break;
  case 6:
    {
      // tilt.
      typename simpleFile<double>::Mat bump[3], out[3];
      std::vector<typename simpleFile<double>::Vec3>  points;
      std::vector<typename simpleFile<double>::Veci3> polys;
      const std::string fn(argv[4]);
      bool is_obj(false);
      if(fn[fn.size() - 1] == 'm') {
        if(!file.loadp2or3(bump, argv[4]))
          return - 2;
      } else if(fn[fn.size() - 1] == 'j') {
        if(!file.loadobj(points, polys, argv[4]))
          return - 2;
        is_obj = true;
      } else
        return - 2;
      for(int i = 0; i < M_TILT; i ++) {
        const auto mtilt(redig.tiltprep(data[0], i, M_TILT, psi));
        typename simpleFile<double>::Mat tilt0;
        if(is_obj)
          tilt0 = redig.tilt(redig.makeRefMatrix(data[0], 1), redig.tiltprep(points, polys, redig.makeRefMatrix(data[0], 1), mtilt));
        else
          tilt0 = redig.tilt(redig.makeRefMatrix(data[0], 1), bump[0], mtilt);
        for(int j = 0; j < 3; j ++)
          out[j] = redig.pullRefMatrix(tilt0, 1, data[j]);
        std::string outfile(argv[3]);
        outfile += std::string("-") + std::to_string(i) + std::string(".ppm");
        file.savep2or3(outfile.c_str(), out, false);
      }
      return 0;
    }
    break;
  case 14:
    // tilt2
    {
      typename simpleFile<double>::Mat bump[3], data2[2][3], out[3];
      std::vector<typename simpleFile<double>::Vec3>  points[2];
      std::vector<typename simpleFile<double>::Veci3> polys;
      const std::string fn(argv[4]);
      bool is_obj(false);
      if(fn[fn.size() - 1] == 'm') {
        if(!file.loadp2or3(bump, argv[4]))
          return - 2;
        data2[0][0] = data2[0][1] = data2[0][2] = data2[1][0] = data2[1][1] = data2[1][2] = bump[0] * 0.;
        for(int i = 0; i < bump[0].rows(); i ++)
          for(int j = 0; j < bump[0].cols(); j ++) {
            const int jj(- offsetx * bump[0].cols() + j);
            if(0 <= jj && jj < bump[0].cols()) {
              bump[0](i, j) = bump[2](i, jj);
              for(int k = 0; k < 3; k ++)
                data2[0][k](i, j) = data[k](i, j);
            } else {
              bump[0](i, j) = 0.;
              for(int k = 0; k < 3; k ++)
                data2[0][k](i, j) = 0.;
            }
            const int jj2(offsetx * bump[0].cols() + j);
            if(0 <= jj2 && jj2 < bump[0].cols()) {
              bump[1](i, j) = bump[2](i, jj2);
              for(int k = 0; k < 3; k ++)
                data2[1][k](i, j) = data[k](i, j);
            } else {
              bump[1](i, j) = 0.;
              for(int k = 0; k < 3; k ++)
                data2[1][k](i, j) = 0.;
            }
          }
      } else if(fn[fn.size() - 1] == 'j') {
        if(!file.loadobj(points[0], polys, argv[4]))
          return - 2;
        points[1] = points[0];
        for(int i = 0; i < points[0].size(); i ++) {
          points[0][i][1] -= offsetx * data[0].cols();
          points[1][i][1] += offsetx * data[0].cols();
        }
        is_obj = true;
      } else
        return - 2;
      for(int i = 0; i < 2; i ++) {
        const auto mtilt(redig.tiltprep(data[0], i, 2, psi));
        typename simpleFile<double>::Mat tilt0;
        if(is_obj)
          tilt0 = redig.tilt(redig.makeRefMatrix(data[0], 1), redig.tiltprep(points[i % 2], polys, redig.makeRefMatrix(data[0], 1), mtilt));
        else
          tilt0 = redig.tilt(redig.makeRefMatrix(data[0], 1), bump[i % 2], mtilt);
        for(int j = 0; j < 3; j ++)
          if(is_obj)
            out[j] = redig.pullRefMatrix(tilt0, 1, data[j]);
          else
            out[j] = redig.pullRefMatrix(tilt0, 1, data2[i % 2][j]);
        std::string outfile(argv[3]);
        const char* names[2] = {"-L.ppm", "-R.ppm"};
        outfile += std::string(names[i % 2]);
        file.savep2or3(outfile.c_str(), out, false);
      }
      return 0;
    }
    break;
  case 8:
    {
      // tilt3.
      typename simpleFile<double>::Mat bump[3], out[3];
      std::vector<typename simpleFile<double>::Vec3>  points;
      std::vector<typename simpleFile<double>::Veci3> polys;
      const std::string fn(argv[4]);
      bool is_obj(false);
      if(fn[fn.size() - 1] == 'm') {
        if(!file.loadp2or3(bump, argv[4]))
          return - 2;
      } else if(fn[fn.size() - 1] == 'j') {
        if(!file.loadobj(points, polys, argv[4]))
          return - 2;
        is_obj = true;
      } else
        return - 2;
      for(int i = 0; i < 4; i ++) {
        const auto mtilt(redig.tiltprep(data[0], i, 4, psi2));
        typename simpleFile<double>::Mat tilt0;
        if(is_obj)
          tilt0 = redig.tilt(redig.makeRefMatrix(data[0], 1), redig.tiltprep(points, polys, redig.makeRefMatrix(data[0], 1), mtilt));
        else
          tilt0 = redig.tilt(redig.makeRefMatrix(data[0], 1), bump[0], mtilt);
        for(int j = 0; j < 3; j ++)
          out[j] = redig.pullRefMatrix(tilt0, 1, data[j]);
        std::string outfile(argv[3]);
        outfile += std::string("-") + std::to_string(i) + std::string(".ppm");
        file.savep2or3(outfile.c_str(), out, false);
      }
      return 0;
    }
    break;
  case 24:
  case 20:
    {
      // tilt4.
      // tilt5.
      typename simpleFile<double>::Mat bump[3], out[3];
      std::vector<typename simpleFile<double>::Vec3>  points;
      std::vector<typename simpleFile<double>::Veci3> polys;
      const std::string fn(argv[4]);
      bool is_obj(false);
      if(fn[fn.size() - 1] == 'm') {
        if(!file.loadp2or3(bump, argv[4]))
          return - 2;
      } else if(fn[fn.size() - 1] == 'j') {
        if(!file.loadobj(points, polys, argv[4]))
          return - 2;
        is_obj = true;
      } else
        return - 2;
      for(int i = 0; i < M_TILT; i ++) {
        const auto mtilt(redig.tiltprep(data[0], 0, 2, double(mode == 20 ? 2 : 1) * psi2 * ((M_TILT - 1) / 2. - i) / ((M_TILT - 1) / 2.)));
        typename simpleFile<double>::Mat tilt0;
        if(is_obj)
          tilt0 = redig.tilt(redig.makeRefMatrix(data[0], 1), redig.tiltprep(points, polys, redig.makeRefMatrix(data[0], 1), mtilt));
        else
          tilt0 = redig.tilt(redig.makeRefMatrix(data[0], 1), bump[0], mtilt);
        for(int j = 0; j < 3; j ++)
          out[j] = redig.pullRefMatrix(tilt0, 1, data[j]);
        std::string outfile(argv[3]);
        outfile += std::string("-") + std::to_string(i) + std::string(".ppm");
        file.savep2or3(outfile.c_str(), out, false);
      }
      return 0;
    }
    break;
  case 9:
  case 25:
    {
      // 2d - 2d match with hidden calculated 3d.
      // and pmatch.
      typename simpleFile<double>::Mat data1[3], bdata[3], bdata1[3], mdata[3], mdata1[3], mout[3], mmout1, bump1;
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
      matchPartialPartial<double> statmatch;
      if(mode != 9)
        statmatch.ndiv = statmatch.ndiv * 1.2;
      auto matches(statmatch.match(shape0, shape1));
      matches.resize(min(int(matches.size()), nmatchhid));
      matches = statmatch.elim(matches, data, mout, bump1, shape1);
      if(mode == 9)
        for(int n = 0; n < min(int(matches.size()), nshow); n ++) {
          std::cerr << "Writing " << n << " / " << matches.size();
          saveMatches<double>(std::string(argv[3]) + std::to_string(n + 1), matches[n], shape0, shape1, data, mout, bump0, bump1, emph);
        }
      else
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
          matchPartialPartial<double> pstatmatch;
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
    }
    return 0;
  case 10:
    {
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
      matchPartialPartial<double> statmatch;
      const auto matches(statmatch.match(shape, datapoly));
      for(int n = 0; n < min(int(matches.size()), nshow); n ++) {
        std::cerr << "Writing " << n << " / " << matches.size();
        std::vector<typename simpleFile<double>::Vec3> mdatapoly;
        mdatapoly.reserve(datapoly.size());
        auto match(matches[n]);
        for(int k = 0; k < datapoly.size(); k ++)
          mdatapoly.push_back(match.transform(datapoly[k]) / match.ratio);
        match.rot    *= match.rot.transpose();
        match.offset *= double(0);
        saveMatches<double>(std::string(argv[3]) + std::to_string(n + 1), match, shape, mdatapoly, data, zero, bump, zero[0], emph);
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
      matchPartialPartial<double> statmatch;
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
  case 12:
  case 27:
    {
      // maskobj.
      // maskobj2.
      std::vector<typename simpleFile<double>::Vec3> points;
      std::vector<typename simpleFile<double>::Veci3> polys;
      if((mode == 12 && argc < 6) || (mode == 27 && argc < 5) ||
         !file.loadobj(points, polys, argv[3])) {
        usage();
        return - 2;
      }
      redig.initialize(vbox0, rz);
      redig.maskVectors(points, polys, data[0]);
      auto edges(redig.getEdges(data[0], points));
      if(mode == 27) {
        file.saveobj(points, polys, argv[4], false);
        return - 3;
      }
      double ratio;
      std::stringstream stream(argv[5]);
      stream >> ratio;
      for(int i = 0; i < points.size(); i ++)
        points[i] *= ratio;
      std::stringstream stream2(argv[6]);
      stream2 >> ratio;
      double M(points[0][2]), m(points[0][2]);
      for(int i = 1; i < points.size(); i ++) {
        M = std::max(points[i][2], M);
        m = std::min(points[i][2], m);
      }
      if(M == m) M += 1.;
      for(int i = 0; i < points.size(); i ++)
        points[i][2] *= ratio / (M - m);
      file.saveobj(points, polys, argv[4], false, true, edges);
    }
    return 0;
  case 13:
    // habit.
    {
      std::vector<typename simpleFile<double>::Vec3> pdst, psrc;
      std::vector<typename simpleFile<double>::Veci3> poldst, polsrc;
      if(argc < 5 || !file.loadobj(pdst, poldst, argv[4]) ||
                     !file.loadobj(psrc, polsrc, argv[4])) {
        usage();
        return - 2;
      }
      matchPartialPartial<double> statmatch;
      const auto match(statmatch.match(pdst, psrc));
      for(int i = 0; i < nshow; i ++) {
        const auto mhull0(redig.delaunay2(pdst, match[i].dstpoints));
        const auto mhull1(match[i].hull(match[i].srcpoints, match[i].reverseHull(match[i].dstpoints, mhull0)));
        file.saveobj(redig.takeShape(pdst, psrc, match[i], mhull0, mhull1, Memph),
                     poldst, (argv[3] + std::string("-emph-") + to_string(i) +
                                        std::string(".obj")).c_str());
      }
    }
    return 0;
  case 19:
    // habit2.
    {
      std::vector<typename simpleFile<double>::Vec3>  pdst,   psrc;
      std::vector<typename simpleFile<double>::Veci3> poldst, polsrc;
      if(argc < 8 || !file.loadobj(pdst, poldst, argv[4]) ||
                     !file.loadobj(psrc, polsrc, argv[5])) {
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
      typename simpleFile<double>::Mat in(int(My + 1.), int(Mx + 1.));
      matchPartialPartial<double> statmatch;
      auto m(redig.tiltprep(in, std::atoi(argv[6]), std::atoi(argv[7]), psi2 * 2.));
      statmatch.complementMatch(m, pdst, psrc, statmatch.makeG(pdst), statmatch.makeG(psrc));
      const auto mhull0(redig.delaunay2(pdst, m.dstpoints));
      const auto mhull1(m.hull(m.srcpoints, m.reverseHull(m.dstpoints, mhull0)));
      file.saveobj(redig.takeShape(pdst, psrc, m, mhull0, mhull1, double(.5)),
                   poldst, (argv[3] + std::string("-emph") +
                                      std::string(".obj")).c_str());
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
      matchWholePartial<double> wholematch;
      const auto matches(wholematch.match(shape, datapoly, center, bone));
      for(int n = 0; n < min(int(matches.size()), nshow); n ++) {
        std::cerr << "Writing " << n << " / " << matches.size();
        for(int m = 0; m < matches[n].size(); m ++)
          if(matches[n][m].dstpoints.size())
            saveMatches<double>(std::string(argv[3]) + std::to_string(n + 1) + std::string("-") + std::to_string(m), matches[n][m], shape, datapoly[m], data, zero, bump, zero[0], emph);
      }
    }
    break;
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
      matchWholePartial<double> wholematch;
      const auto matches(wholematch.match(datapoly0, datapoly, center, bone));
      for(int n = 0; n < min(int(matches.size()), nshow); n ++) {
        std::cerr << "Writing " << n << " / " << matches.size();
        for(int m = 0; m < matches[n].size(); m ++)
          if(matches[n][m].dstpoints.size())
            saveMatches<double>(std::string(argv[3]) + std::to_string(n + 1) + std::string("-") + std::to_string(m), matches[n][m], datapoly0, datapoly[m], zero, zero, zero[0], zero[0], emph);
      }
    }
    return 0;
  case 29:
    // draw obj.
    {
      std::vector<typename simpleFile<double>::Vec3>  datapoly;
      std::vector<typename simpleFile<double>::Veci3> polynorms;
      if(!file.loadobj(datapoly, polynorms, argv[4]))
        return - 2;
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
      file.savep2or3((std::string(argv[3]) + std::string("-obj.ppm")).c_str(), res, true);
    }
    return 0;
  case 21:
    // draw gltf.
    {
      std::vector<std::vector<typename simpleFile<double>::Vec3> >  datapoly;
      std::vector<std::vector<typename simpleFile<double>::Veci3> > polynorms;
      std::vector<std::vector<typename simpleFile<double>::Veci4> > bone;
      std::vector<typename simpleFile<double>::Vec3>  center;
      if(!file.loadglTF(datapoly, polynorms, center, bone, argv[4]))
        return - 2;
      match_t<double> m;
      m.offset[0] += min(data[0].rows(), data[0].cols()) / 2.;
      m.offset[1] += min(data[0].rows(), data[0].cols()) / 2.;
      m.ratio     *= min(data[0].rows(), data[0].cols()) / 2.;
      for(int i = 0; i < datapoly.size(); i ++)
        for(int j = 0; j < datapoly[i].size(); j ++)
          datapoly[i][j] = m.transform(datapoly[i][j]);
      for(int i = 0; i < datapoly.size(); i ++) {
        if(!datapoly[i].size()) continue;
        cerr << i << ": " << endl;
        for(int j = 0; j < datapoly[i].size(); j ++) {
          for(int k = 0; k < datapoly[i][j].size(); k ++)
            cerr << datapoly[i][j][k] << ", ";
          cerr << endl;
        }
        typename simpleFile<double>::Mat res[3];
        std::vector<int> idx;
        for(int j = 0; j < datapoly[i].size(); j ++)
          idx.push_back(j);
        res[0] = res[1] = res[2] = redig.showMatch(data[0] * 0., datapoly[i], polynorms[i], double(120));
        file.savep2or3((std::string(argv[3]) + std::string("-") + std::to_string(i) + std::string(".ppm")).c_str(), res, true);
      }
    }
    return 0;
  case 22:
    // tiltp
    {
      typename simpleFile<double>::Mat bump[3], data2[2][3], out[3];
      std::vector<typename simpleFile<double>::Vec3>  points[2];
      std::vector<typename simpleFile<double>::Veci3> polys;
      const std::string fn(argv[4]);
      bool is_obj(false);
      if(fn[fn.size() - 1] == 'm') {
        if(!file.loadp2or3(bump, argv[4]))
          return - 2;
        data2[0][0] = data2[0][1] = data2[0][2] = data2[1][0] = data2[1][1] = data2[1][2] = bump[0] * 0.;
        for(int i = 0; i < bump[0].rows(); i ++)
          for(int j = 0; j < bump[0].cols(); j ++) {
            const int jj(- offsetx * bump[0].cols() + j);
            if(0 <= jj && jj < bump[0].cols()) {
              bump[0](i, j) = bump[2](i, jj);
              for(int k = 0; k < 3; k ++)
                data2[0][k](i, j) = data[k](i, j);
            } else {
              bump[0](i, j) = 0.;
              for(int k = 0; k < 3; k ++)
                data2[0][k](i, j) = 0.;
            }
            const int jj2(offsetx * bump[0].cols() + j);
            if(0 <= jj2 && jj2 < bump[0].cols()) {
              bump[1](i, j) = bump[2](i, jj2);
              for(int k = 0; k < 3; k ++)
                data2[1][k](i, j) = data[k](i, j);
            } else {
              bump[1](i, j) = 0.;
              for(int k = 0; k < 3; k ++)
                data2[1][k](i, j) = 0.;
            }
          }
      } else if(fn[fn.size() - 1] == 'j') {
        if(!file.loadobj(points[0], polys, argv[4]))
          return - 2;
        points[1] = points[0];
        for(int i = 0; i < points[0].size(); i ++) {
          points[0][i][1] -= offsetx * data[0].cols();
          points[1][i][1] += offsetx * data[0].cols();
        }
        is_obj = true;
      } else
        return - 2;
      for(int k = 0; k < M_TILTROT; k ++) {
        for(int i = 0; i < 2; i ++) {
          const auto mtilt(redig.tiltprep(data[0], i, 2, (tiltrote * k + tiltrots * (M_TILTROT - k - 1)) / M_TILTROT));
          typename simpleFile<double>::Mat tilt0;
          if(is_obj)
            tilt0 = redig.tilt(redig.makeRefMatrix(data[0], 1), redig.tiltprep(points[i % 2], polys, redig.makeRefMatrix(data[0], 1), mtilt));
          else
            tilt0 = redig.tilt(redig.makeRefMatrix(data[0], 1), bump[i % 2], mtilt);
          for(int j = 0; j < 3; j ++)
            if(is_obj)
              out[j] = redig.pullRefMatrix(tilt0, 1, data[j]);
            else
              out[j] = redig.pullRefMatrix(tilt0, 1, data2[i % 2][j]);
          std::string outfile(argv[3]);
          outfile += std::string("-") + std::to_string(k);
          const char* names[2] = {"-L.ppm", "-R.ppm"};
          outfile += std::string(names[i % 2]);
          file.savep2or3(outfile.c_str(), out, false);
        }
      }
      return 0;
    }
    break;
  default:
    usage();
    return - 1;
  }
  redig.normalize(data, 1.);
  if(!file.savep2or3(argv[3], data, ! true))
    return - 3;
  return 0;
}

