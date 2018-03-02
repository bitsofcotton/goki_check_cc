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
const int    vbox(8);
const int    objvbox(2);
const int    M_TILT(32);
const double psi(.95);
const int    Mpoly(2000);

void usage() {
  cout << "Usage: tools (enlarge|collect|idetect|bump|obj|bump2|rbump2|tilt|match|match3d|match2dh3d|maskobj) <input filename>.p[gp]m <output filename>.p[gp]m <args>?" << endl;
  return;
}

int main(int argc, const char* argv[]) {
  Eigen::Matrix<double, 3, 3> I3;
  Eigen::Matrix<double, 3, 1> zero3;
  for(int i = 0; i < 3; i ++)
    for(int j = 0; j < 3; j ++)
      I3(i, j) = (i == j ? double(1) : double(0));
  zero3[0] = zero3[1] = zero3[2] = double(0);
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
      enlarger2ex<double> enlarger;
      for(int i = 0; i < 3; i ++)
        data[i] = enlarger.compute(data[i], enlarger.ENLARGE_BOTH);
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
      PseudoBump<double>  bump;
      enlarger2ex<double> collector;
      data[0] = bump.getPseudoBump(rgb2l(data));
      data[1] = collector.compute(rgb2l(data), collector.COLLECT_BOTH);
      for(int i = 0; i < data[0].rows(); i ++)
        for(int j = 0; j < data[0].cols(); j ++)
          data[0](i, j) /= double(1) + data[1](i, j);
      data[1] = data[2] = data[0];
    }
    break;
  case 7:
    {
      // obj.
      PseudoBump<double> bump;
      std::vector<Eigen::Matrix<double, 3, 1> > points;
      std::vector<Eigen::Matrix<int,    3, 1> > facets;
      bump.getPseudoVec(data[0], points, facets, objvbox);
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
      tilt.initialize(sqrt(double(bump[0].rows() * bump[0].cols())));
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
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> data1[3];
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> bdata[3], bdata1[3];
      if(!loadp2or3<double>(data1, argv[4]))
        return - 2;
      if(!loadp2or3<double>(bdata, argv[5]))
        return - 2;
      if(!loadp2or3<double>(bdata1, argv[6]))
        return - 2;
      std::vector<double> emph;
      for(int i = 0; i <= nemph; i ++)
        emph.push_back(double(i) / nemph);
      PseudoBump<double> bump;
      std::vector<Eigen::Matrix<int,    3, 1> > sute;
      std::vector<Eigen::Matrix<double, 3, 1> > shape0, shape1;
      auto& bump0(bdata[0]);
      auto& bump1(bdata1[0]);
      bump.getPseudoVec(bump0, shape0, sute, vbox);
      bump.getPseudoVec(bump1, shape1, sute, vbox);
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mout[3], mbump;
      mbump = mout[0] = mout[1] = mout[2] = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(data[0].rows(), data[0].cols());
      for(int i = 0; i < min(mout[0].rows(), data1[0].rows()); i ++) {
        for(int j = 0; j < min(mout[0].cols(), data1[0].cols()); j ++) {
          mout[0](i, j) = data1[0](i, j);
          mout[1](i, j) = data1[1](i, j);
          mout[2](i, j) = data1[2](i, j);
          mbump(i, j)   = bump1(i, j);
        }
        for(int j = min(mout[0].cols(), data1[0].cols()); j < mout[0].cols(); j ++)
          mout[0](i, j) = mout[1](i, j) = mout[2](i, j) = mbump(i, j) = double(0);
      }
      for(int i = min(mout[0].rows(), data1[0].rows()); i < mout[0].rows(); i ++)
        for(int j = 0; j < mout[0].cols(); j ++)
          mout[0](i, j) = mout[1](i, j) = mout[2](i, j) = mbump(i, j) = double(0);
      matchPartialPartial<double> statmatch;
      auto matches(statmatch.match(shape0, shape1));
      for(int n = 0; n < min(int(matches.size()), nshow); n ++) {
        cerr << n << " / " << matches.size() << "(" << matches[n].rdepth << ", " << matches[n].ratio << ")" << endl;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> outs[3], outs2[3], outs3[3], outs4[3], outs5[3];
        auto hull(loadBumpSimpleMesh<double>(shape0, matches[n].dstpoints));
        std::vector<Eigen::Matrix<int, 3, 1> > mhull0, mhull1;
        for(int idx = 0; idx < hull.size(); idx ++) {
          Eigen::Matrix<int, 3, 1> buf;
          for(int idx2 = 0; idx2 < hull[idx].size(); idx2 ++)
            buf[idx2] = matches[n].dstpoints[hull[idx][idx2]];
          mhull0.push_back(buf);
          for(int idx2 = 0; idx2 < hull[idx].size(); idx2 ++)
            buf[idx2] = matches[n].srcpoints[hull[idx][idx2]];
          mhull1.push_back(buf);
        }
        tilter<double> tilt;
        for(int idx = 0; idx < 3; idx ++)
          outs[idx] = tilt.tilt(showMatch<double>(mout[idx], shape1, mhull1), mbump, matches[n].rot, I3, matches[n].offset, matches[n].ratio, zero3);
        normalize<double>(outs, 1.);
        std::string outfile;
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-src.ppm");
        savep2or3<double>(outfile.c_str(), outs, false);
        
        for(int idx = 0; idx < 3; idx ++)
          outs2[idx] = showMatch<double>(data[idx], shape0, mhull0);
        normalize<double>(outs2, 1.);
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-dst.ppm");
        savep2or3<double>(outfile.c_str(), outs2, false);
        
        for(int idx = 0; idx < 3; idx ++)
          outs3[idx] = outs[idx] + outs2[idx];
        normalize<double>(outs3, 1.);
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-match.ppm");
        savep2or3<double>(outfile.c_str(), outs3, false);
        
        match_t<double> rmatchn(~ matches[n]);
        reDig<double>   redig;
        for(int kk = 0; kk < emph.size(); kk ++) {
          for(int idx = 0; idx < 3; idx ++) {
            outs4[idx] = redig.emphasis(mout[idx], mbump, data[idx], shape1, shape0, rmatchn,    hull, emph[kk]);
            outs5[idx] = redig.emphasis(data[idx], bump0, mout[idx], shape0, shape1, matches[n], hull, emph[kk]);
          }
          outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-emphasis-a-") + std::to_string(kk) + std::string(".ppm");
          savep2or3<double>(outfile.c_str(), outs4, false);
          outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-emphasis-b-") + std::to_string(kk) + std::string(".ppm");
          savep2or3<double>(outfile.c_str(), outs5, false);
        }
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
      auto zero(data[0] * double(0));
      matchPartialPartial<double>   statmatch;
      auto matches(statmatch.match(shape, datapoly));
      for(int n = 0; n < min(int(matches.size()), nshow); n ++) {
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> outs[3], outs2[3], outs3[3], outs4[3];
        cerr << "Writing " << n << " / " << matches.size() << "(" << matches[n].rdepth << ", " << matches[n].ratio << ")" << endl;
        
        auto ch(loadBumpSimpleMesh<double>(shape, matches[n].dstpoints));

        std::vector<Eigen::Matrix<double, 3, 1> > shape3d;
        for(int idx = 0; idx < datapoly.size(); idx ++)
          shape3d.push_back(matches[n].rot * matches[n].ratio * datapoly[idx] + matches[n].offset);
        std::vector<Eigen::Matrix<int, 3, 1> > dch;
        for(int idx = 0; idx < ch.size(); idx ++) {
          Eigen::Matrix<int, 3, 1> buf;
          for(int idx2 = 0; idx2 < ch[idx].size(); idx2 ++)
            buf[idx2] = matches[n].srcpoints[ch[idx][idx2]];
          dch.push_back(buf);
        }
        for(int idx = 0; idx < 3; idx ++)
          outs[idx] = showMatch<double>(zero, shape3d, dch);
        normalize<double>(outs, 1.);
        std::string outfile;
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-src.ppm");
        savep2or3<double>(outfile.c_str(), outs, false);
        
        std::vector<Eigen::Matrix<int, 3, 1> > mch;
        for(int idx = 0; idx < ch.size(); idx ++) {
          Eigen::Matrix<int, 3, 1> buf;
          for(int idx2 = 0; idx2 < ch[idx].size(); idx2 ++)
            buf[idx2] = matches[n].dstpoints[ch[idx][idx2]];
          mch.push_back(buf);
        }
        for(int idx = 0; idx < 3; idx ++)
          outs2[idx] = showMatch<double>(data[idx], shape, mch);
        normalize<double>(outs2, 1.);
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-dst.ppm");
        savep2or3<double>(outfile.c_str(), outs2, false);
        
        for(int idx = 0; idx < 3; idx ++)
          outs3[idx] = outs[idx] + outs2[idx];
        normalize<double>(outs3, 1.);
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-match.ppm");
        savep2or3<double>(outfile.c_str(), outs3, false);
        
        match_t<double> rmatchn(~ matches[n]);
        reDig<double>   redig;
        for(int kk = 0; kk < emph.size(); kk ++) {
          for(int idx = 0; idx < 3; idx ++) {
            outs4[idx] = redig.emphasis(data[idx], bump, zero, shape, datapoly, matches[n], ch, emph[kk]);
          }
          outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-emphasis-") + std::to_string(kk) + std::string(".ppm");
          savep2or3<double>(outfile.c_str(), outs4, false);
        }
      }
    }
    return 0;
  case 11:
    {
      // 2d - 2d match with hidden 3d model.
    }
    break;
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
      auto edges(getEdges<double>(data[0], points, objvbox));
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

