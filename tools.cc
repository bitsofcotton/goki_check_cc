#include <cstdio>
#include <cmath>
#include <iostream>
#include <Eigen/Core>
#include "ppm2eigen.hh"
#include "obj2vector.hh"
#include "enlarge.hh"
#include "fisheye.hh"
#include "tilt.hh"
#include "scancontext.hh"

using namespace std;

void usage() {
  cout << "Usage: tools (enlarge|bump|collect|idetect|tilt|lpoly|3poly|match|match3d) <input filename>.p[gp]m <output filename>.p[gp]m <args>?" << endl;
  return;
}

int main(int argc, const char* argv[]) {
  std::ios::sync_with_stdio(false);
  if(argc < 4) {
    usage();
    return 0;
  }
  int mode = - 1;
  if(strcmp(argv[1], "enlarge") == 0)
    mode = 0;
  else if(strcmp(argv[1], "bump") == 0)
    mode = 2;
  else if(strcmp(argv[1], "collect") == 0)
    mode = 4;
  else if(strcmp(argv[1], "idetect") == 0)
    mode = 5;
  else if(strcmp(argv[1], "tilt") == 0)
    mode = 6;
  else if(strcmp(argv[1], "lpoly") == 0)
    mode = 8;
  else if(strcmp(argv[1], "3poly") == 0)
    mode = 11;
  else if(strcmp(argv[1], "match") == 0)
    mode = 9;
  else if(strcmp(argv[1], "match3d") == 0)
    mode = 10;
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
      enlarger2ex<double> enlarger;
      for(int i = 0; i < 3; i ++)
        data[i] = enlarger.compute(data[i], enlarger.ENLARGE_3BOTH);
    }
    break;
  case 4:
    {
      enlarger2ex<double> detect;
      for(int i = 0; i < 3; i ++)
        data[i] = detect.compute(data[i], enlarger2ex<double>::COLLECT_BOTH);
    }
    break;
  case 5:
    {
      enlarger2ex<double> idetector;
      for(int i = 0; i < 3; i ++)
        data[i] = idetector.compute(data[i], enlarger2ex<double>::IDETECT_BOTH);
    }
    break;
  case 2:
    {
      PseudoBump<double> bump;
      std::vector<Eigen::Matrix<double, 3, 1> > points;
      std::vector<Eigen::Matrix<int,    3, 1> > delaunay;
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> bumps;
      const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b_mesh(bump.getPseudoBumpVec(rgb2l(data), points, delaunay, bumps));
      data[0] = data[1] = data[2] = bumps;
      std::cout << "Handled points:" << std::endl;
      for(int i = 0; i < points.size(); i ++)
        std::cout << points[i].transpose() << std::endl;
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> outs[3];
      outs[0] = outs[1] = outs[2] = b_mesh;
      normalize<double>(outs, 1.);
      savep2or3<double>(argv[4], outs, false);
      saveobj(points, delaunay, argv[5], double(4));
    }
    break;
  case 6:
    {
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> bump[3], out[3], zero;
      if(!loadp2or3<double>(bump, argv[4]))
        return - 2;
      zero = bump[0] * double(0);
      tilter<double> tilt;
      const int M_TILT(16);
      tilt.initialize(.3125);
      for(int i = 0; i < M_TILT; i ++) {
        for(int j = 0; j < 3; j ++)
          out[j] = tilt.tilt(tilt.tilt(data[j], bump[0], i, M_TILT, .98), zero, - i, M_TILT, .98);
        std::string outfile(argv[3]);
        outfile += std::string("-") + std::to_string(i) + std::string(".ppm");
        savep2or3<double>(outfile.c_str(), out, false);
      }
      return 0;
    }
    break;
  case 8:
    {
      lowFreq<double> lf;
      auto points(lf.getLowFreq(rgb2l(data)));
      std::vector<int> dstpoints;
      for(int i = 0; i < points.size(); i ++)
        dstpoints.push_back(i);
      normalize<double>(data, 1.);
      std::vector<Eigen::Matrix<int, 3, 1> > pptr(loadBumpSimpleMesh<double>(points, dstpoints));
      for(int i = 0; i < 3; i ++)
        data[i] = showMatch<double>(data[i], points, pptr);
      std::cout << "Handled points:" << std::endl;
      for(int i = 0; i < points.size(); i ++)
        std::cout << points[i].transpose() << std::endl;
      saveobj(points, pptr, argv[4]);
    }
    break;
  case 11:
    {
      std::vector<Eigen::Matrix<double, 3, 1> > datapoly;
      std::vector<Eigen::Matrix<int, 3, 1> >    polynorms;
      if(!loadobj<double>(datapoly, polynorms, argv[4]))
        return - 1;
      for(int i = 0; i < datapoly.size(); i ++) {
        datapoly[i][0] += data[0].rows() / 2.;
        datapoly[i][1] += data[0].cols() / 2.;
      }
      for(int i = 0; i < polynorms.size(); i ++)
        std::cout << polynorms[i].transpose() << std::endl;
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> zero(data[0]);
      zero *= double(0);
      for(int idx = 0; idx < 3; idx ++)
        data[idx] = showMatch<>(zero, datapoly, polynorms);
    }
    break;
  case 9:
    {
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> data1[3];
      if(!loadp2or3<double>(data1, argv[4]))
        return - 2;
      // XXX: configure me:
      double zrs(1.5);
      double zre( .5);
      int    zrl(3);
      int    nshow(8);
      double emph(.95);
      PseudoBump<double> bump;
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> sutei;
      std::vector<Eigen::Matrix<int,    3, 1> > sute;
      std::vector<Eigen::Matrix<double, 3, 1> > shape0, shape1;
      auto bump0(bump.getPseudoBumpVec(rgb2l(data),  shape0, sute, sutei));
      auto bump1(bump.getPseudoBumpVec(rgb2l(data1), shape1, sute, sutei));
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mout[3];
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mbump;
      Eigen::Matrix<double, 3, 3> I3;
      mbump = mout[0] = mout[1] = mout[2] = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(data[0].rows(), data[0].cols());
      for(int i = 0; i < 3; i ++)
        for(int j = 0; j < 3; j ++)
          I3(i, j) = (i == j ? double(1) : double(0));
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
      Eigen::Matrix<double, 3, 1> zero3;
      zero3[0] = zero3[1] = zero3[2] = double(0);
      matchPartialPartial<double> statmatch;
      statmatch.threshr = double(.75);
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
        tilter<double> tilt;
        reDig<double>  redig;
        tilt.initialize(zr / (data[0].rows() * data[0].cols()));
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
        
        for(int idx = 0; idx < 3; idx ++)
          outs4[idx] = redig.emphasis(data[idx], bump0, shape0, shape1, matches[n], hull, double(1.) - emph);
        normalize<double>(outs4, 1.);
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-emphasis-0.ppm");
        savep2or3<double>(outfile.c_str(), outs4, false);
        for(int idx = 0; idx < 3; idx ++)
          outs5[idx] = redig.emphasis(data[idx], bump0, shape0, shape1, matches[n], hull, double(1.) + emph);
        normalize<double>(outs5, 1.);
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-emphasis-2.ppm");
        savep2or3<double>(outfile.c_str(), outs5, false);
      }
    }
    return 0;
  case 10:
    {
      std::vector<Eigen::Matrix<double, 3, 1> > datapoly;
      std::vector<Eigen::Matrix<int, 3, 1> >    polynorms;
      if(!loadobj<double>(datapoly, polynorms, argv[4]))
        return - 2;
      // XXX: configure me:
      if(datapoly.size() > 2000) {
        std::cerr << "Too many vertices." << std::endl;
        return - 2;
      }
      double zrs(1.5);
      double zre( .5);
      int    zrl(3);
      int    nshow(8);
      double emph(.95);
      PseudoBump<double> bumper;
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> sutei;
      std::vector<Eigen::Matrix<double, 3, 1> > shape;
      std::vector<Eigen::Matrix<int,    3, 1> > sute;
      auto bump(bumper.getPseudoBumpVec(rgb2l(data), shape, sute, sutei));
      Eigen::Matrix<double, 3, 1> zero3;
      zero3[0] = zero3[1] = zero3[2] = double(0);
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> zero(data[0].rows(), data[0].cols());
      for(int i = 0; i < zero.rows(); i ++)
        for(int j = 0; j < zero.cols(); j ++)
          zero(i, j) = double(0);
      Eigen::Matrix<double, 3, 3> I3;
      for(int i = 0; i < 3; i ++)
        for(int j = 0; j < 3; j ++)
          I3(i, j) = (i == j ? double(1) : double(0));
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
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> outs[3], outs2[3], outs3[3], outs4[3], outs5[3];
        tilter<double> tilt;
        reDig<double>  redig;
        tilt.initialize(zr / sqrt((data[0].rows() * data[0].cols()) / (data[0].rows() * data[0].cols())));
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
        
        for(int idx = 0; idx < 3; idx ++)
          outs4[idx] = redig.emphasis(data[idx], bump, shape, datapoly, matches[n], ch, double(1) - emph);
        normalize<double>(outs4, 1.);
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-emphasis-0.ppm");
        savep2or3<double>(outfile.c_str(), outs4, false);
        for(int idx = 0; idx < 3; idx ++)
          outs5[idx] = redig.emphasis(data[idx], bump, shape, datapoly, matches[n], ch, double(1) + emph);
        normalize<double>(outs5, 1.);
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-emphasis-2.ppm");
        savep2or3<double>(outfile.c_str(), outs5, false);
      }
    }
    return 0;
  default:
    break;
  }
  normalize<double>(data, 1.);
  if(!savep2or3<double>(argv[3], data, ! true))
    return - 3;
  return 0;
}

