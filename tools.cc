#include <cstdio>
#include <cmath>
#include <iostream>
#include <Eigen/Core>
#include "fileio.hh"
#include "enlarge.hh"
#include "match.hh"
#include "redig.hh"

using namespace std;

// XXX: configure me:
const int    nshow(4);
const int    nshowh(16);
const int    nemph(2);
const double Memph(.25);
const int    vbox0(2);
const int    vbox(16);
const double rz(1 / 3.);
const double offsetx(.1);
const int    M_TILT(32);
const double psi(.025);
const int    Mpoly(2000);

void usage() {
  cout << "Usage: tools (enlarge|collect|idetect|bump|obj|bump2|rbump2|tilt|tilt2|match|match3d|match3dbone|match2dh3d|match2dh3dbone|maskobj|habit) <input filename>.p[gp]m <output filename>.p[gp]m <args>?" << endl;
  return;
}

template <typename T> void saveMatches(const std::string& outbase, const match_t<T>& match, const std::vector<Eigen::Matrix<T, 3, 1> >& shape0, const std::vector<Eigen::Matrix<T, 3, 1> >& shape1, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> in0[3], const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> in1[3], const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& bump0, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& bump1, const std::vector<T>& emph) {
  // generated bump map is inverted.
  reDig<T> redig;
  std::cerr << "(" << match.rdepth << ", " << match.ratio << ")" << std::endl;
  
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> outs[3];
  const auto mhull0(redig.delaunay2(shape0, match.dstpoints));
  const auto mhull1(match.hull(match.srcpoints, match.reverseHull(match.dstpoints, mhull0)));
  
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> sin1[3];
  for(int idx = 0; idx < 3; idx ++)
    sin1[idx] = redig.showMatch(in1[idx], shape1, mhull1);
  redig.normalize(sin1, 1.);
  for(int idx = 0; idx < 3; idx ++)
    outs[idx] = redig.tilt(sin1[idx], bump1, match);
  std::string outfile(outbase + std::string("-src.ppm"));
  savep2or3<T>(outfile.c_str(), outs, false);
  
  for(int idx = 0; idx < 3; idx ++)
    outs[idx] = redig.showMatch(in0[idx], shape0, mhull0);
  redig.normalize(outs, 1.);
  outfile = outbase + std::string("-dst.ppm");
  savep2or3<T>(outfile.c_str(), outs, false);
  
  for(int idx = 0; idx < 3; idx ++)
    outs[idx] = redig.replace(in0[idx], shape1, match, mhull1);
  outfile = outbase + std::string("-repl.ppm");
  savep2or3<T>(outfile.c_str(), outs, false);
  
  const auto rin0(redig.makeRefMatrix(in0[0], 1));
  const auto rin1(redig.makeRefMatrix(in1[0], 1 + rin0.rows() * rin0.cols()));
  for(int kk = 0; kk < emph.size(); kk ++) {
    const auto reref(redig.emphasis(rin0, rin1, bump1, shape0, shape1, match, mhull0, mhull1, emph[kk]));
    for(int idx = 0; idx < 3; idx ++) 
      outs[idx] = (in0[idx] + redig.pullRefMatrix(reref, 1 + rin0.rows() * rin0.cols(), sin1[idx])) / T(2);
    outfile = outbase + std::string("-emph-") + std::to_string(kk) + std::string(".ppm");
    savep2or3<T>(outfile.c_str(), outs, false);
  }
  
  saveobj(redig.takeShape(shape0, shape1, match, mhull0, mhull1,
                          emph[emph.size() - 1]), mhull0,
          (outbase + std::string("-emph.obj")).c_str());
  saveobj(redig.takeShape(shape1, shape0, ~ match, mhull1, mhull0,
                          emph[emph.size() - 1]), mhull1,
          (outbase + std::string("-emphr.obj")).c_str());
  return;
}

template <typename T> void resizeDst2(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mout[3], Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& bump1, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mmout1, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> data[3], const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& bump, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mask, const int& h, const int& w) {
  mmout1 = bump1 = mout[0] = mout[1] = mout[2] = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(h, w);
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
  else if(strcmp(argv[1], "maskobj") == 0)
    mode = 12;
  else if(strcmp(argv[1], "habit") == 0)
    mode = 13;
  else if(strcmp(argv[1], "tilt2") == 0)
    mode = 14;
  if(mode < 0) {
    usage();
    return - 1;
  }
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> data[3];
  if(!loadp2or3<double>(data, argv[2]))
    return - 2;
  reDig<double> redig;
  redig.initialize(vbox, rz);
  switch(mode) {
  case 0:
    {
      // enlarge.
      enlarger2ex<double> enlarger, denlarger;
      for(int i = 0; i < 3; i ++) {
        const auto xye(enlarger.compute(data[i], enlarger.ENLARGE_BOTH));
        data[i] = xye + redig.tilt45(denlarger.compute(redig.tilt45(data[i], false), denlarger.ENLARGE_BOTH), true, xye);
      }
    }
    break;
  case 4:
    {
      // collect.
      enlarger2ex<double> detect, ddetect;
      for(int i = 0; i < 3; i ++) {
        const auto xye(detect.compute(data[i], detect.COLLECT_BOTH));
        data[i] = xye * (data[i].rows() + data[i].cols()) + redig.tilt45(ddetect.compute(redig.tilt45(data[i], false), ddetect.COLLECT_BOTH), true, xye) * sqrt(double(data[i].rows() * data[i].cols()));
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
      const auto xye(bump.compute(redig.rgb2l(data), bump.BUMP_BOTH));
      // data[0] = xye;
      data[0] = redig.autoLevel(xye, data[0].rows() + data[0].cols());
      // data[0] = xye + redig.tilt45(bump.compute(redig.tilt45(redig.rgb2l(data), false), bump.BUMP_BOTH), true, xye);
      data[1] = data[2] = data[0];
    }
    break;
  case 7:
    {
      // obj.
      std::vector<Eigen::Matrix<double, 3, 1> > points;
      std::vector<Eigen::Matrix<int,    3, 1> > facets;
      redig.initialize(vbox0, rz);
      redig.getTileVec(data[0], points, facets);
      saveobj(points, facets, argv[3]);
    }
    return 0;
  case 3:
    {
      // bump2 for training output.
      enlarger2ex<double> bump;
      auto X(.49000/.17697 * data[0] + .31000/.17697  * data[1] + .20000/.17697 * data[2]);
      auto Z(                          .010000/.17697 * data[1] + .99000/.17697 * data[2]);
      data[0] = X     / 8.;
      data[1] = bump.compute(redig.rgb2l(data), bump.BUMP_BOTH) / 8.;
      data[2] = Z     / 8.;
    }
    // with no auto-level.
    if(!savep2or3<double>(argv[3], data, ! true))
      return - 3;
    return 0;
  case 5:
    {
      // reverse bump2.
      auto R(  .41847    * data[0] - .15866   * data[1] - .082835 * data[2]);
      auto G(- .091169   * data[0] + .25243   * data[1] + .015708 * data[2]);
      auto B(  .00092090 * data[0] - .0025498 * data[1] + .17860  * data[2]);
      data[0] = R;
      data[1] = G;
      data[2] = B;
    }
    break;
  case 6:
    {
      // tilt.
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> bump[3], out[3];
      if(!loadp2or3<double>(bump, argv[4]))
        return - 2;
      auto zero(bump[0] * double(0));
      for(int i = 0; i < M_TILT; i ++) {
        for(int j = 0; j < 3; j ++)
          out[j] = redig.tilt(data[j], bump[0], i, M_TILT, psi);
        std::string outfile(argv[3]);
        outfile += std::string("-") + std::to_string(i) + std::string(".ppm");
        savep2or3<double>(outfile.c_str(), out, false);
      }
      return 0;
    }
    break;
  case 14:
    // tilt2
    {
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> bump[3], out[3];
      if(!loadp2or3<double>(bump, argv[4]))
        return - 2;
      for(int i = 0; i < bump[0].rows(); i ++)
        for(int j = 0; j < bump[0].cols(); j ++) {
          const int jj(- offsetx * bump[0].cols() + j);
          if(0 <= jj && jj < bump[0].cols())
            bump[0](i, j) = bump[2](i, jj);
          else
            bump[0](i, j) = 0.;
          const int jj2(offsetx * bump[0].cols() + j);
          if(0 <= jj2 && jj2 < bump[0].cols())
            bump[1](i, j) = bump[2](i, jj2);
          else
            bump[1](i, j) = 0.;
        }
      auto zero(bump[0] * double(0));
      for(int i = 0; i < 2; i ++) {
        for(int j = 0; j < 3; j ++)
          out[j] = redig.tilt(data[j], bump[i % 2], i, 2, psi);
        std::string outfile(argv[3]);
        const char* names[2] = {"-L.ppm", "-R.ppm"};
        outfile += std::string(names[i % 2]);
        savep2or3<double>(outfile.c_str(), out, false);
      }
      return 0;
    }
    break;
  case 9:
    {
      // 2d - 2d match with hidden calculated 3d.
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> data1[3], bdata[3], bdata1[3], mdata[3], mdata1[3], mout[3], mmout1, bump1;
      if(!loadp2or3<double>(data1, argv[4]))
        return - 2;
      if(!loadp2or3<double>(bdata, argv[5]))
        return - 2;
      if(!loadp2or3<double>(bdata1, argv[6]))
        return - 2;
      if(!loadp2or3<double>(mdata, argv[7]))
        return - 2;
      if(!loadp2or3<double>(mdata1, argv[8]))
        return - 2;
      std::vector<double> emph;
      for(int i = 0; i <= nemph; i ++)
        emph.push_back(double(i) / nemph * Memph);
      resizeDst2(mout, bump1, mmout1, data1, bdata1[0], mdata1[0], data[0].rows(), data[0].cols());
      std::vector<Eigen::Matrix<int,    3, 1> > sute;
      std::vector<Eigen::Matrix<double, 3, 1> > shape0, shape1;
      auto& bump0(bdata[0]);
      redig.getTileVec(bump0, shape0, sute);
      redig.getTileVec(bump1, shape1, sute);
      std::vector<int> id0, id1;
      for(int i = 0; i < shape0.size(); i ++)
        id0.push_back(i);
      for(int i = 0; i < shape1.size(); i ++)
        id1.push_back(i);
      redig.maskVectors(shape0, redig.delaunay2(shape0, id0), mdata[0]);
      redig.maskVectors(shape1, redig.delaunay2(shape1, id1), mmout1);
      matchPartialPartial<double> statmatch;
      const auto matches(statmatch.elim(statmatch.match(shape0, shape1), data, mout, bump1, shape1));
      for(int n = 0; n < min(int(matches.size()), nshow); n ++) {
        std::cerr << "Writing " << n << " / " << matches.size();
        saveMatches<double>(std::string(argv[3]) + std::to_string(n + 1), matches[n], shape0, shape1, data, mout, bump0, bump1, emph);
      }
    }
    return 0;
  case 10:
    {
      // 3d - 2d match.
      std::vector<Eigen::Matrix<double, 3, 1> > datapoly;
      std::vector<Eigen::Matrix<int, 3, 1> >    polynorms;
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> bump0[3], mask0[3];
      if(!loadobj<double>(datapoly, polynorms, argv[4]))
        return - 2;
      if(!loadp2or3<double>(bump0, argv[5]))
        return - 2;
      if(!loadp2or3<double>(mask0, argv[6]))
        return - 2;
      if(datapoly.size() > Mpoly) {
        std::cerr << "Too many vertices." << std::endl;
        return - 2;
      }
      std::vector<double> emph;
      for(int i = 0; i <= nemph; i ++)
        emph.push_back(double(i) / nemph * Memph);
      std::vector<Eigen::Matrix<double, 3, 1> > shape;
      std::vector<Eigen::Matrix<int,    3, 1> > sute;
      auto& bump(bump0[0]);
      redig.getTileVec(bump, shape, sute);
      std::vector<int> id;
      for(int i = 0; i < shape.size(); i ++)
        id.push_back(i);
      auto delau(redig.delaunay2(shape, id));
      redig.maskVectors(shape, delau, mask0[0]);
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> zero[3];
      for(int i = 0; i < 3; i ++)
        zero[i] = bump * double(0);
      matchPartialPartial<double> statmatch;
      const auto matches(statmatch.match(shape, datapoly));
      for(int n = 0; n < min(int(matches.size()), nshow); n ++) {
        std::cerr << "Writing " << n << " / " << matches.size();
        std::vector<Eigen::Matrix<double, 3, 1> > mdatapoly;
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
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> data1[3], bdata[3], bdata1[3], mdata[3], mdata1[3], mout[3], mmout1, bump1;
      std::vector<Eigen::Matrix<double, 3, 1> > datapoly;
      std::vector<Eigen::Matrix<int, 3, 1> >    polynorms;
      if(!loadp2or3<double>(data1, argv[4]))
        return - 2;
      if(!loadp2or3<double>(bdata, argv[5]))
        return - 2;
      if(!loadp2or3<double>(bdata1, argv[6]))
        return - 2;
      if(!loadp2or3<double>(mdata, argv[7]))
        return - 2;
      if(!loadp2or3<double>(mdata1, argv[8]))
        return - 2;
      if(!loadobj<double>(datapoly, polynorms, argv[9]))
        return - 2;
      if(datapoly.size() > Mpoly) {
        std::cerr << "Too many vertices." << std::endl;
        return - 2;
      }
      std::vector<double> emph;
      for(int i = 0; i <= nemph; i ++)
        emph.push_back(double(i) / nemph * Memph);
      resizeDst2(mout, bump1, mmout1, data1, bdata1[0], mdata1[0], data[0].rows(), data[0].cols());
      std::vector<Eigen::Matrix<int,    3, 1> > sute;
      std::vector<Eigen::Matrix<double, 3, 1> > shape0, shape1;
      auto& bump0(bdata[0]);
      redig.getTileVec(bump0, shape0, sute);
      redig.getTileVec(bump1, shape1, sute);
      std::vector<int> id0, id1;
      for(int i = 0; i < shape0.size(); i ++)
        id0.push_back(i);
      for(int i = 0; i < shape1.size(); i ++)
        id1.push_back(i);
      redig.maskVectors(shape0, redig.delaunay2(shape0, id0), mdata[0]);
      redig.maskVectors(shape1, redig.delaunay2(shape1, id1), mmout1);
      matchPartialPartial<double> statmatch;
      const auto match0(statmatch.match(shape0, datapoly));
      const auto match1(statmatch.match(shape1, datapoly));
      std::vector<match_t<double> > matches;
      for(int n = 0; n < min(int(match0.size()), nshowh); n ++)
        for(int m = 0; m < min(int(match1.size()), nshowh); m ++)
          matches.push_back(match0[n] / match1[m]);
      matches = statmatch.elim(matches, data, mout, bump1, shape1);
      for(int n = 0; n < min(int(matches.size()), nshow); n ++) {
        const auto& relmatch(matches[n]);
        std::cerr << "Writing " << n << " / " << matches.size();
        saveMatches<double>(std::string(argv[3]) + std::string("-") + std::to_string(n), relmatch, shape0, shape1, data, mout, bump0, bump1, emph);
      }
    }
    return 0;
  case 12:
    {
      // maskobj.
      std::vector<Eigen::Matrix<double, 3, 1> > points;
      std::vector<Eigen::Matrix<int, 3, 1> >    polys;
      if(argc < 6 || !loadobj<double>(points, polys, argv[3])) {
        usage();
        return - 2;
      }
      redig.initialize(vbox0, rz);
      redig.maskVectors(points, polys, data[0]);
      auto edges(redig.getEdges(data[0], points));
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
      saveobj(points, polys, argv[4], true, edges);
    }
    return 0;
  case 13:
    // habit.
    {
      std::vector<Eigen::Matrix<double, 3, 1> > pdst,   psrc;
      std::vector<Eigen::Matrix<int,    3, 1> > poldst, polsrc;
      if(argc < 5 || !loadobj<double>(pdst, poldst, argv[4]) ||
                     !loadobj<double>(psrc, polsrc, argv[4])) {
        usage();
        return - 2;
      }
      matchPartialPartial<double> statmatch;
      const auto match(statmatch.match(pdst, psrc));
      for(int i = 0; i < nshow; i ++) {
        const auto mhull0(redig.delaunay2(pdst, match[i].dstpoints));
        const auto mhull1(match[i].hull(match[i].srcpoints, match[i].reverseHull(match[i].dstpoints, mhull0)));
        saveobj(redig.takeShape(pdst, psrc, match[i], mhull0, mhull1, Memph),
                mhull0, (argv[3] + std::string("-emph-") + to_string(i) +
                                   std::string(".obj")).c_str());
      }
    }
    return 0;
  default:
    usage();
    return - 1;
  }
  redig.normalize(data, 1.);
  if(!savep2or3<double>(argv[3], data, ! true))
    return - 3;
  return 0;
}

