#include <cstdio>
#include <cmath>
#include <iostream>
#include <Eigen/Core>
#include "ppm2eigen.hh"
#include "obj2vector.hh"
#include "enlarge.hh"
#include "fisheye.hh"
#include "tilt.hh"
#include "match.hh"

using namespace std;

// XXX: configure me:
const int    nshow(4);
const int    nemph(4);
const int    vbox(4);
const int    M_TILT(32);
const double psi(.95);
const int    Mpoly(2000);

void usage() {
  cout << "Usage: tools (enlarge|collect|idetect|bump|obj|bump2|rbump2|tilt|match|match3d|match2dh3d|maskobj) <input filename>.p[gp]m <output filename>.p[gp]m <args>?" << endl;
  return;
}

template <typename T> void saveMatches(const std::string& outbase, const match_t<T>& match, const std::vector<Eigen::Matrix<T, 3, 1> >& shape0, const std::vector<Eigen::Matrix<T, 3, 1> >& shape1, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> in0[3], const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> in1[3], const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& bump0, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& bump1, const std::vector<T>& emph) {
  reDig<T>  redig;
  tilter<T> tilt;
  // generated bump map is inverted.
  tilt.initialize(- sqrt(double(bump1.rows() * bump1.cols())) / double(6));
  
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> outs[3];
  
  Eigen::Matrix<double, 3, 3> I3;
  Eigen::Matrix<double, 3, 1> zero3;
  for(int i = 0; i < 3; i ++)
    for(int j = 0; j < 3; j ++)
      I3(i, j) = (i == j ? double(1) : double(0));
  zero3[0] = zero3[1] = zero3[2] = double(0);
  
  const auto mhull0(redig.delaunay2(shape0, match.dstpoints));
  const auto mhull1(match.hull(match.srcpoints, match.reverseHull(match.dstpoints, mhull0)));
  
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> sin1[3];
  for(int idx = 0; idx < 3; idx ++) {
    sin1[idx] = redig.showMatch(in1[idx], shape1, mhull1);
    outs[idx] = tilt.tilt(sin1[idx], bump1, match.rot, I3, match.offset, match.ratio, zero3);
  }
  normalize<T>(outs, 1.);
  std::string outfile(outbase + std::string("-src.ppm"));
  savep2or3<T>(outfile.c_str(), outs, false);
  
  for(int idx = 0; idx < 3; idx ++)
    outs[idx] = redig.showMatch(in0[idx], shape0, mhull0);
  normalize<T>(outs, 1.);
  outfile = outbase + std::string("-dst.ppm");
  savep2or3<T>(outfile.c_str(), outs, false);
  
  for(int idx = 0; idx < 3; idx ++)
    outs[idx] = redig.replace(in0[idx], shape1, match, mhull1);
  outfile = outbase + std::string("-repl.ppm");
  savep2or3<T>(outfile.c_str(), outs, false);
  
  const auto rin0(redig.makeRefMatrix(in0[0], 1));
  const auto rin1(redig.makeRefMatrix(in1[0], 1 + rin0.rows() * rin0.cols()));
  for(int kk = 0; kk < emph.size(); kk ++) {
    const auto reref(redig.emphasis(rin0, rin1, bump1, shape0, shape1, match, mhull0, mhull1, emph[kk], tilt));
    for(int idx = 0; idx < 3; idx ++) 
      outs[idx] = (in0[idx] + redig.pullRefMatrix(reref, 1 + rin0.rows() * rin0.cols(), sin1[idx])) / T(2);
    outfile = outbase + std::string("-emph-") + std::to_string(kk) + std::string(".ppm");
    savep2or3<T>(outfile.c_str(), outs, false);
  }
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
  else if(strcmp(argv[1], "match2dh3d") == 0)
    mode = 11;
  else if(strcmp(argv[1], "maskobj") == 0)
    mode = 12;
  if(mode < 0) {
    usage();
    return - 1;
  }
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> data[3];
  if(!loadp2or3<double>(data, argv[2]))
    return - 2;
  switch(mode) {
  case 0:
    {
      // enlarge.
      enlarger2ex<double> enlarger, denlarger;
      std::vector<double> stat;
      for(int i = 0; i < 3; i ++) {
        const auto xye(enlarger.compute(data[i], enlarger.ENLARGE_BOTH));
        data[i] = xye + tilt45(denlarger.compute(tilt45(data[i], false), enlarger.ENLARGE_BOTH), true, xye);
        for(int j = 0; j < data[i].rows(); j ++)
          for(int k = 0; k < data[i].cols(); k ++)
            stat.push_back(data[i](j, k));
      }
      std::sort(stat.begin(), stat.end());
      const int sauto((data[0].rows() + data[0].cols()) * 4 * 4);
      for(int i = 0; i < 3; i ++)
        for(int j = 0; j < data[i].rows(); j ++)
          for(int k = 0; k < data[i].cols(); k ++)
            data[i](j, k) = std::max(std::min(stat[stat.size() - sauto], data[i](j, k)), stat[sauto]);
    }
    break;
  case 4:
    {
      // collect.
      enlarger2ex<double> detect;
      for(int i = 0; i < 3; i ++)
        data[i] = detect.compute(data[i], detect.COLLECT_BOTH);
    }
    break;
  case 1:
    {
      // idetect.
      enlarger2ex<double> idetect;
      for(int i = 0; i < 3; i ++)
        data[i] = idetect.compute(data[i], idetect.IDETECT_BOTH);
    }
    break;
  case 2:
    {
      // bump.
      PseudoBump<double> bump;
      data[0] = bump.getPseudoBump(rgb2l(data));
      data[1] = data[2] = data[0];
    }
    break;
  case 7:
    {
      // obj.
      PseudoBump<double> bump;
      std::vector<Eigen::Matrix<double, 3, 1> > points;
      std::vector<Eigen::Matrix<int,    3, 1> > facets;
      bump.getPseudoVec(data[0], points, facets, vbox);
      saveobj(points, facets, argv[3]);
    }
    return 0;
  case 3:
    {
      // bump2 for training output.
      PseudoBump<double> bump;
      auto X(.49000/.17697 * data[0] + .31000/.17697  * data[1] + .20000/.17697 * data[2]);
      auto Z(                          .010000/.17697 * data[1] + .99000/.17697 * data[2]);
      data[0] = X     / 8.;
      data[1] = bump.getPseudoBump(rgb2l(data)) / 8.;
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
      tilter<double> tilt;
      tilt.initialize(sqrt(double(bump[0].rows() * bump[0].cols())) / double(6));
      for(int i = 0; i < M_TILT; i ++) {
        for(int j = 0; j < 3; j ++)
          out[j] = tilt.tilt(tilt.tilt(data[j], bump[0], i, M_TILT, psi), zero, - i, M_TILT, psi);
        std::string outfile(argv[3]);
        outfile += std::string("-") + std::to_string(i) + std::string(".ppm");
        savep2or3<double>(outfile.c_str(), out, false);
      }
      return 0;
    }
    break;
  case 9:
    {
      // 2d - 2d match with hidden calculated 3d.
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> data1[3], bdata[3], bdata1[3], mout[3], bump1;
      if(!loadp2or3<double>(data1, argv[4]))
        return - 2;
      if(!loadp2or3<double>(bdata, argv[5]))
        return - 2;
      if(!loadp2or3<double>(bdata1, argv[6]))
        return - 2;
      std::vector<double> emph;
      for(int i = 0; i <= nemph; i ++)
        emph.push_back(double(i) / nemph);
      bump1 = mout[0] = mout[1] = mout[2] = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(data[0].rows(), data[0].cols());
      for(int i = 0; i < min(mout[0].rows(), data1[0].rows()); i ++) {
        for(int j = 0; j < min(mout[0].cols(), data1[0].cols()); j ++) {
          mout[0](i, j) = data1[0](i, j);
          mout[1](i, j) = data1[1](i, j);
          mout[2](i, j) = data1[2](i, j);
          bump1(i, j)   = bdata1[0](i, j);
        }
        for(int j = min(mout[0].cols(), data1[0].cols()); j < mout[0].cols(); j ++)
          bump1(i, j) = mout[0](i, j) = mout[1](i, j) = mout[2](i, j) = double(0);
      }
      for(int i = min(mout[0].rows(), data1[0].rows()); i < mout[0].rows(); i ++)
        for(int j = 0; j < mout[0].cols(); j ++)
          bump1(i, j) = mout[0](i, j) = mout[1](i, j) = mout[2](i, j) = double(0);
      PseudoBump<double> bump;
      std::vector<Eigen::Matrix<int,    3, 1> > sute;
      std::vector<Eigen::Matrix<double, 3, 1> > shape0, shape1;
      auto& bump0(bdata[0]);
      bump.getPseudoVec(bump0, shape0, sute, vbox);
      bump.getPseudoVec(bump1, shape1, sute, vbox);
      matchPartialPartial<double> statmatch;
      const auto matches(statmatch.match(shape0, shape1));
      for(int n = 0; n < min(int(matches.size()), nshow); n ++) {
        std::cerr << n << " / " << matches.size() << "(" << matches[n].rdepth << ", " << matches[n].ratio << ")" << endl;
        saveMatches<double>(std::string(argv[3]) + std::to_string(n + 1), matches[n], shape0, shape1, data, mout, bump0, bump1, emph);
      }
    }
    return 0;
  case 10:
    {
      // 3d - 2d match.
      std::vector<Eigen::Matrix<double, 3, 1> > datapoly;
      std::vector<Eigen::Matrix<int, 3, 1> >    polynorms;
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> bump0[3];
      if(!loadobj<double>(datapoly, polynorms, argv[4]))
        return - 2;
      if(!loadp2or3<double>(bump0, argv[5]))
        return - 2;
      if(datapoly.size() > Mpoly) {
        std::cerr << "Too many vertices." << std::endl;
        return - 2;
      }
      std::vector<double> emph;
      for(int i = 0; i <= nemph; i ++)
        emph.push_back(double(i) / nemph);
      PseudoBump<double> bumper;
      std::vector<Eigen::Matrix<double, 3, 1> > shape;
      std::vector<Eigen::Matrix<int,    3, 1> > sute;
      auto& bump(bump0[0]);
      bumper.getPseudoVec(bump, shape, sute, vbox);
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> zero[3];
      for(int i = 0; i < 3; i ++)
        zero[i] = bump * double(0);
      matchPartialPartial<double> statmatch;
      const auto matches(statmatch.match(shape, datapoly));
      for(int n = 0; n < min(int(matches.size()), nshow); n ++) {
        std::cerr << "Writing " << n << " / " << matches.size() << "(" << matches[n].rdepth << ", " << matches[n].ratio << ")" << endl;
        std::vector<Eigen::Matrix<double, 3, 1> > mdatapoly;
        auto match(matches[n]);
        for(int k = 0; k < datapoly.size(); k ++)
          mdatapoly.push_back(match.transform(datapoly[k]) / match.ratio);
        match.offset *= double(0);
        match.rot    *= match.rot.transpose();
        saveMatches<double>(std::string(argv[3]) + std::to_string(n + 1), match, shape, mdatapoly, data, zero, bump, zero[0], emph);
        
        reDig<double> redig;
        const auto hsrc(redig.delaunay2(shape, match.srcpoints));
        const auto hdst(match.hull(match.dstpoints, match.reverseHull(match.srcpoints, hsrc)));
        const auto dp2(redig.takeShape(mdatapoly, shape, ~ match, hsrc, hdst, double(1)));
        saveobj(dp2, polynorms, (std::string(argv[3]) + std::to_string(n + 1) + std::string(".obj")).c_str());;
      }
    }
    return 0;
  case 11:
    {
      // 2d - 2d match with hidden 3d model.
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> data1[3], bdata[3], bdata1[3], mout[3], bump1;
      std::vector<Eigen::Matrix<double, 3, 1> > datapoly;
      std::vector<Eigen::Matrix<int, 3, 1> >    polynorms;
      if(!loadp2or3<double>(data1, argv[4]))
        return - 2;
      if(!loadp2or3<double>(bdata, argv[5]))
        return - 2;
      if(!loadp2or3<double>(bdata1, argv[6]))
        return - 2;
      if(!loadobj<double>(datapoly, polynorms, argv[7]))
        return - 2;
      if(datapoly.size() > Mpoly) {
        std::cerr << "Too many vertices." << std::endl;
        return - 2;
      }
      std::vector<double> emph;
      for(int i = 0; i <= nemph; i ++)
        emph.push_back(double(i) / nemph);
      bump1 = mout[0] = mout[1] = mout[2] = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(data[0].rows(), data[0].cols());
      for(int i = 0; i < min(mout[0].rows(), data1[0].rows()); i ++) {
        for(int j = 0; j < min(mout[0].cols(), data1[0].cols()); j ++) {
          mout[0](i, j) = data1[0](i, j);
          mout[1](i, j) = data1[1](i, j);
          mout[2](i, j) = data1[2](i, j);
          bump1(i, j)   = bdata1[0](i, j);
        }
        for(int j = min(mout[0].cols(), data1[0].cols()); j < mout[0].cols(); j ++)
          bump1(i, j) = mout[0](i, j) = mout[1](i, j) = mout[2](i, j) = double(0);
      }
      for(int i = min(mout[0].rows(), data1[0].rows()); i < mout[0].rows(); i ++)
        for(int j = 0; j < mout[0].cols(); j ++)
          bump1(i, j) = mout[0](i, j) = mout[1](i, j) = mout[2](i, j) = double(0);
      PseudoBump<double> bump;
      std::vector<Eigen::Matrix<int,    3, 1> > sute;
      std::vector<Eigen::Matrix<double, 3, 1> > shape0, shape1;
      auto& bump0(bdata[0]);
      bump.getPseudoVec(bump0, shape0, sute, vbox);
      bump.getPseudoVec(bump1, shape1, sute, vbox);
      matchPartialPartial<double> statmatch;
      const auto match0(statmatch.match(shape0, datapoly));
      const auto match1(statmatch.match(shape1, datapoly));
      for(int n = 0; n < min(int(match0.size()), nshow); n ++)
        for(int m = 0; m < min(int(match1.size()), nshow); m ++) {
          match_t<double> relmatch;
          relmatch.rot    = match0[n].rot    * match1[m].rot.transpose();
          relmatch.ratio  = match0[n].ratio  / match1[m].ratio;
          relmatch.offset = match0[n].offset - relmatch.rot * relmatch.ratio * match1[m].offset;
          relmatch.rdepth = match0[n].rdepth + match1[m].rdepth;
          relmatch.dstpoints = vector<int>();
          relmatch.srcpoints = vector<int>();
          for(int i = 0; i < match0[n].srcpoints.size(); i ++)
            for(int j = 0; j < match1[m].srcpoints.size(); j ++)
              if(match0[n].srcpoints[i] == match1[m].srcpoints[j]) {
                relmatch.dstpoints.push_back(match0[n].dstpoints[i]);
                relmatch.srcpoints.push_back(match1[m].dstpoints[j]);
              }
          std::cerr << "(" << n << ", " << m << ") / (" << match0.size() << ", " << match1.size() << ") : (" << relmatch.rdepth << ", " << relmatch.ratio << ")" << endl;
          saveMatches<double>(std::string(argv[3]) + std::to_string(n + 1) + std::string("-") + std::to_string(m + 1), relmatch, shape0, shape1, data, mout, bump0, bump1, emph);
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
      maskVectors<double>(points, polys, data[0]);
      auto edges(getEdges<double>(data[0], points, vbox));
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
  default:
    usage();
    break;
  }
  normalize<double>(data, 1.);
  if(!savep2or3<double>(argv[3], data, ! true))
    return - 3;
  return 0;
}

