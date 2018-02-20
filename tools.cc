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
const double zrs(1.2);
const double zre( .8);
const int    zrl(3);
const int    nshow(8);
const int    nemph(4);
const int    vbox(8);
const double threshr(.75);
const int    M_TILT(16);
const int    Mpoly(2000);

void usage() {
  cout << "Usage: tools (enlarge|bump|obj|bump2|rbump2|collect|tilt|lpoly|match|match3d|match2dh3d|maskobj) <input filename>.p[gp]m <output filename>.p[gp]m <args>?" << endl;
  return;
}

int main(int argc, const char* argv[]) {
  std::ios::sync_with_stdio(false);
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
  else if(strcmp(argv[1], "tilt") == 0)
    mode = 6;
  else if(strcmp(argv[1], "lpoly") == 0)
    mode = 8;
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
        data[i] = enlarger.compute(data[i], enlarger.ENLARGE_3BOTH);
    }
    break;
  case 4:
    {
      // collect.
      enlarger2ex<double> detect;
      for(int i = 0; i < 3; i ++)
        data[i] = detect.compute(data[i], enlarger2ex<double>::COLLECT_BOTH);
    }
    break;
  case 2:
    {
      // bump.
      PseudoBump<double> bump;
      data[0] = data[1] = data[2] = bump.getPseudoBump(rgb2l(data));
    }
    break;
  case 7:
    {
      // obj.
      PseudoBump<double> bump;
      std::vector<Eigen::Matrix<double, 3, 1> > points;
      std::vector<Eigen::Matrix<int,    3, 1> > facets;
      bump.getPseudoVec(data[0], points, facets);
      if(4 < argc) {
        double ratio;
        std::stringstream stream(argv[4]);
        stream >> ratio;
        for(int i = 0; i < points.size(); i ++)
          points[i] *= ratio;
      }
      double lz(0);
      for(int i = 0; i < points.size(); i ++)
        lz = max(lz, abs(points[i][2]));
      saveobj(points, facets, argv[3], true, lz + double(2));
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
      for(int i = 0; i < M_TILT; i ++) {
        for(int j = 0; j < 3; j ++)
          out[j] = tilt.tilt(tilt.tilt(data[j], bump[0], i, M_TILT, .95), zero, - i, M_TILT, .99);
        std::string outfile(argv[3]);
        outfile += std::string("-") + std::to_string(i) + std::string(".ppm");
        savep2or3<double>(outfile.c_str(), out, false);
      }
      return 0;
    }
    break;
  case 8:
    {
      // lowpoly.
      lowFreq<double> lf;
      auto points(lf.getLowFreq(rgb2l(data)));
      std::vector<int> dstpoints;
      for(int i = 0; i < points.size(); i ++)
        dstpoints.push_back(i);
      normalize<double>(data, 1.);
      auto pptr(loadBumpSimpleMesh<double>(points, dstpoints));
      for(int i = 0; i < 3; i ++)
        data[i] = showMatch<double>(data[i], points, pptr);
      std::cout << "Handled points:" << std::endl;
      for(int i = 0; i < points.size(); i ++)
        std::cout << points[i].transpose() << std::endl;
      saveobj(points, pptr, argv[4]);
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
        emph.push_back(double(i) / nemph / double(8));
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
      statmatch.threshr = threshr;
      std::vector<match_t<double> > matches;
      for(double zr = zrs;
          (zrs / zre < double(1) && zr < zre) || 
          (zre / zrs < double(1) && zre < zr);
          zr *= pow(zre / zrs, double(1) / double(zrl)))
        for(double zr2 = zrs;
            (zrs / zre < double(1) && zr2 < zre) ||
            (zre / zrs < double(1) && zre < zr2);
            zr2 *= pow(zre / zrs, double(1) / double(zrl))) {
          std::vector<Eigen::Matrix<double, 3, 1> > sshape0(shape0), sshape1(shape1);
          for(int i = 0; i < sshape0.size(); i ++)
            sshape0[i][2] *= zr;
          for(int i = 0; i < sshape1.size(); i ++)
            sshape1[i][2] *= zr2;
          statmatch.match(sshape0, sshape1, matches);
      }
      double zr(zrs);
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
            outs4[idx] = redig.emphasis(mout[idx], mbump, shape1, shape0, rmatchn,    hull, emph[kk]);
            outs5[idx] = redig.emphasis(data[idx], bump0, shape0, shape1, matches[n], hull, emph[kk]);
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
        emph.push_back(double(i) / nemph / double(8));
      PseudoBump<double> bumper;
      std::vector<Eigen::Matrix<double, 3, 1> > shape;
      std::vector<Eigen::Matrix<int,    3, 1> > sute;
      auto& bump(bump0[0]);
      bumper.getPseudoVec(bump, shape, sute, vbox);
      auto zero(data[0] * double(0));
      matchPartialPartial<double> statmatch;
      std::vector<match_t<double> > matches;
      for(double zr = zrs;
          (zrs / zre < double(1) && zr < zre) || 
          (zre / zrs < double(1) && zre < zr);
          zr *= pow(zre / zrs, double(1) / double(zrl))) {
        std::vector<Eigen::Matrix<double, 3, 1> > sshape(shape);
        for(int i = 0; i < shape.size(); i ++)
          sshape[i][2] *= zr;
        statmatch.match(sshape, datapoly, matches);
      }
      double zr(zrs);
      for(int n = 0; n < min(int(matches.size()), nshow); n ++) {
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> outs[3], outs2[3], outs3[3], outs4[3];
        cerr << "Writing " << n << " / " << matches.size() << "(" << matches[n].rdepth << ", " << matches[n].ratio << ")" << endl;
        
        auto ch(loadBumpSimpleMesh<double>(shape, matches[n].dstpoints));

        vector<Eigen::Matrix<double, 3, 1> > shape3d;
        for(int idx = 0; idx < datapoly.size(); idx ++)
          shape3d.push_back(matches[n].rot * matches[n].ratio * datapoly[idx] + matches[n].offset);
        for(int idx = 0; idx < 3; idx ++)
          outs[idx] = showMatch<double>(zero, shape3d, ch);
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
            outs4[idx] = redig.emphasis(data[idx], bump, shape, datapoly, matches[n], ch, emph[kk]);
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
      std::vector<Eigen::Matrix<double, 3, 1> > datapoly;
      std::vector<Eigen::Matrix<int, 3, 1> >    polynorms;
      if(argc < 5 || !loadobj<double>(datapoly, polynorms, argv[3]))
        return - 2;
      std::vector<int> elim, elimp, after;
      for(int i = 0, ii = 0; i < datapoly.size(); i ++) {
        const int y(std::max(std::min(int(datapoly[i][0]), int(data[0].rows() - 1)), 0));
        const int x(std::max(std::min(int(datapoly[i][1]), int(data[0].cols() - 1)), 0));
        if(data[0](y, x) > .5) {
          elim.push_back(i);
          after.push_back(- 1);
        } else
          after.push_back(ii ++);
      }
      for(int i = 0; i < polynorms.size(); i ++)
        if(std::binary_search(elim.begin(), elim.end(), polynorms[i][0]) ||
           std::binary_search(elim.begin(), elim.end(), polynorms[i][1]) ||
           std::binary_search(elim.begin(), elim.end(), polynorms[i][2]))
          elimp.push_back(i);
      for(int i = 0, j = 0; i < elim.size(); i ++) {
        datapoly.erase(datapoly.begin() + (elim[i] - j), datapoly.begin() + (elim[i] - j + 1));
        j ++;
      }
      for(int i = 0; i < polynorms.size(); i ++)
        for(int j = 0; j < polynorms[i].size(); j ++)
          polynorms[i][j] = after[polynorms[i][j]];
      for(int i = 0, j = 0; i < elimp.size(); i ++) {
        polynorms.erase(polynorms.begin() + (elimp[i] - j), polynorms.begin() + (elimp[i] - j + 1));
        j ++;
      }
      if(5 < argc) {
        double ratio;
        std::stringstream stream(argv[5]);
        stream >> ratio;
        for(int i = 0; i < datapoly.size(); i ++)
          datapoly[i] *= ratio;
      }
      double lz(0);
      for(int i = 0; i < datapoly.size(); i ++)
        lz = max(lz, abs(datapoly[i][2]));
      saveobj(datapoly, polynorms, argv[4], true, lz + double(2));
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

