#include <cstdio>
#include <cmath>
#include <iostream>
#include <Eigen/Core>
#include "ppm2eigen.hh"
#include "obj2vector.hh"
#include "enlarge.hh"
#include "fisheye.hh"
#include "edgedetect.hh"
#include "tilt.hh"
#include "scancontext.hh"

using namespace std;

void usage() {
  cout << "Usage: tools (enlarge|enlargeds|bump|collect|tilt|lpoly|3poly|match|match3d) <input filename>.p[gp]m <output filename>.p[gp]m <args>?" << endl;
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
  else if(strcmp(argv[1], "enlargeds") == 0)
    mode = 1;
  else if(strcmp(argv[1], "bump") == 0)
    mode = 2;
  else if(strcmp(argv[1], "collect") == 0)
    mode = 4;
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
  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> data[3];
  if(!loadp2or3<float>(data, argv[2]))
    return - 2;
  switch(mode) {
  case 0:
    {
      enlarger2ex<float> enlarger;
      for(int kk = 0; kk < 1; kk ++) {
        for(int i = 0; i < 3; i ++)
          data[i] = enlarger.enlarge2(data[i], enlarger2ex<float>::ENLARGE_QUAD);
        autoLevel<float>(data, 3);
      }
    }
    break;
  case 1:
    {
      enlarger2exds<float> enlarger;
      for(int kk = 0; kk < 1; kk ++) {
        for(int i = 0; i < 3; i ++)
          data[i] = enlarger.enlarge2ds(data[i], enlarger2exds<float>::ENLARGE_QUAD);
        autoLevel<float>(data, 3);
      }
    }
    break;
  case 4:
    {
      edgedetect<float> detect;
      for(int i = 0; i < 3; i ++)
        data[i] = detect.detect(data[i], edgedetect<float>::COLLECT_BOTH);
    }
    break;
  case 2:
    {
      PseudoBump<double> bump;
      std::vector<Eigen::Matrix<double, 3, 1> > points;
      std::vector<Eigen::Matrix<int,    3, 1> > delaunay;
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> bumps;
      const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> b_mesh(bump.getPseudoBumpVec(rgb2l(data).template cast<double>(), points, delaunay, bumps, true).template cast<float>());
      data[0] = data[1] = data[2] = bumps.template cast<float>();
      for(int i = 0; i < points.size(); i ++)
        points[i][2] *= double(20);
      std::cout << "Handled points:" << std::endl;
      for(int i = 0; i < points.size(); i ++)
        std::cout << points[i].transpose() << std::endl;
      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> outs[3];
      outs[0] = outs[1] = outs[2] = b_mesh.template cast<float>();
      normalize<float>(outs, 1.);
      savep2or3<float>(argv[4], outs, false);
      saveobj(points, delaunay, argv[5]);
    }
    break;
  case 6:
    {
      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> bump[3];
      if(!loadp2or3<float>(bump, argv[4]))
        return - 2;
      tilter<float> tilt;
      const int M_TILT = 32;
      for(int i = 0; i < M_TILT; i ++) {
        Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> out[3];
        for(int j = 0; j < 3; j ++)
          out[j] = tilt.tilt(data[j], bump[0], i, M_TILT, .99);
        std::string outfile(argv[3]);
        outfile += std::string("-") + std::to_string(i + 1) + std::string(".ppm");
        savep2or3<float>(outfile.c_str(), out, false);
      }
      return 0;
    }
    break;
  case 8:
    {
      lowFreq<double>    lf;
      std::vector<Eigen::Matrix<double, 3, 1> > points(lf.getLowFreq(rgb2l(data).template cast<double>(), 600));
      std::vector<int> dstpoints;
      for(int i = 0; i < points.size(); i ++)
        dstpoints.push_back(i);
      normalize<float>(data, 1.);
      std::vector<Eigen::Matrix<int, 3, 1> > pptr(loadBumpSimpleMesh<double>(points, dstpoints));
      for(int i = 0; i < 3; i ++)
        data[i] = showMatch<double>(data[i].template cast<double>(), points, pptr).template cast<float>();
      std::cout << "Handled points:" << std::endl;
      for(int i = 0; i < points.size(); i ++)
        std::cout << points[i].transpose() << std::endl;
      saveobj(points, pptr, argv[4]);
    }
    break;
  case 11:
    {
      std::vector<Eigen::Matrix<float, 3, 1> > datapoly;
      std::vector<Eigen::Matrix<int, 3, 1> > polynorms;
      if(!loadobj<float>(datapoly, polynorms, argv[4]))
        return - 1;
      for(int i = 0; i < datapoly.size(); i ++) {
        datapoly[i][0] += data[0].rows() / 2.;
        datapoly[i][1] += data[0].cols() / 2.;
      }
      for(int i = 0; i < polynorms.size(); i ++)
        std::cout << polynorms[i].transpose() << std::endl;
      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> zero(data[0]);
      zero *= float(0);
      for(int idx = 0; idx < 3; idx ++)
        data[idx] = showMatch<float>(zero, datapoly, polynorms);
    }
    break;
  case 9:
    {
      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> data1[3], bump0, bump1;
      if(!loadp2or3<float>(data1, argv[4]))
        return - 2;
      // XXX: configure me.
      float zrs(1.);
      float zre(.25);
      int   zrl(4);
      int   nshow(6);
      float emph(.75);
      PseudoBump<double> bump;
      bump.vmax = 150;
      std::vector<Eigen::Matrix<double, 3, 1> > shape0, shape1;
      std::vector<Eigen::Matrix<int,    3, 1> > delaunay0, delaunay1;
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> sute;
      bump0 = bump.getPseudoBumpVec(rgb2l(data).template cast<double>(), shape0, delaunay0, sute, true).template cast<float>();
      bump.vmax = 60;
      bump1 = bump.getPseudoBumpVec(rgb2l(data1).template cast<double>(), shape1, delaunay1, sute, true).template cast<float>();
      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> mout[3];
      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> mbump;
      Eigen::Matrix<double, 3, 3> I3;
      mbump = mout[0] = mout[1] = mout[2] = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>(data[0].rows(), data[0].cols());
      for(int i = 0; i < 3; i ++)
        for(int j = 0; j < 3; j ++)
          I3(i, j) = (i == j ? float(1) : float(0));
      for(int i = 0; i < min(mout[0].rows(), data1[0].rows()); i ++) {
        for(int j = 0; j < min(mout[0].cols(), data1[0].cols()); j ++) {
          mout[0](i, j) = data1[0](i, j);
          mout[1](i, j) = data1[1](i, j);
          mout[2](i, j) = data1[2](i, j);
          mbump(i, j)   = bump1(i, j);
        }
        for(int j = min(mout[0].cols(), data1[0].cols()); j < mout[0].cols(); j ++)
          mout[0](i, j) = mout[1](i, j) = mout[2](i, j) = mbump(i, j) = float(0);
      }
      
      for(int i = min(mout[0].rows(), data1[0].rows()); i < mout[0].rows(); i ++)
        for(int j = 0; j < mout[0].cols(); j ++)
          mout[0](i, j) = mout[1](i, j) = mout[2](i, j) = mbump(i, j) = float(0);
      Eigen::Matrix<double, 3, 1> zero3;
      zero3[0] = zero3[1] = zero3[2] = float(0);
      matchPartialPartial<double> statmatch;
      std::vector<match_t<double> > matches;
      for(float zr = zrs;
          (zrs / zre < float(1) && zr < zre) || 
          (zre / zrs < float(1) && zre < zr);
          zr *= pow(zre / zrs, float(1) / float(zrl)))
        for(float zr2 = zrs;
            (zrs / zre < float(1) && zr2 < zre) ||
            (zre / zrs < float(1) && zre < zr2);
            zr2 *= pow(zre / zrs, float(1) / float(zrl))) {
          std::vector<Eigen::Matrix<double, 3, 1> > sshape0(shape0), sshape1(shape1);
          for(int i = 0; i < sshape0.size(); i ++)
            sshape0[i][2] *= zr;
          for(int i = 0; i < sshape1.size(); i ++)
            sshape1[i][2] *= zr2;
          std::vector<match_t<double> > lmatches(statmatch.match(sshape0, sshape1));
          for(int i = 0; i < lmatches.size(); i ++)
            matches.push_back(lmatches[i]);
      }
      std::sort(matches.begin(), matches.end(), cmpwrap<double>);
      float zr(zrs);
      for(int n = 0; n < min(int(matches.size()), nshow); n ++) {
        Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> outs[3], outs2[3], outs3[3], outs4[3], outs5[3];
        tilter<double> tilt;
        reDig<double>  redig;
        tilt.initialize(zr * sqrt((data1[0].rows() * data1[0].cols()) / (data[0].rows() * data[0].cols())));
        std::vector<Eigen::Matrix<int, 3, 1> > hull1(loadBumpSimpleMesh<double>(shape1, matches[n].srcpoints));
        std::vector<Eigen::Matrix<int, 3, 1> > mhull0, mhull1;
        for(int idx = 0; idx < hull1.size(); idx ++) {
          Eigen::Matrix<int, 3, 1> buf;
          for(int idx2 = 0; idx2 < hull1[idx].size(); idx2 ++)
            buf[idx2] = matches[n].dstpoints[hull1[idx][idx2]];
          mhull0.push_back(buf);
          for(int idx2 = 0; idx2 < hull1[idx].size(); idx2 ++)
            buf[idx2] = matches[n].srcpoints[hull1[idx][idx2]];
          mhull1.push_back(buf);
        }
        cerr << "Writing " << n << " / " << matches.size() << "(" << float(1) / matches[n].rdepth << ", " << matches[n].rpoints << ", " << matches[n].ratio << ")" << endl;
        for(int idx = 0; idx < 3; idx ++)
          outs[idx] = tilt.tilt(showMatch<double>(mout[idx].template cast<double>(), shape1, mhull1), mbump.template cast<double>(), matches[n].rot, I3, matches[n].offset, matches[n].ratio, zero3).template cast<float>();
        normalize<float>(outs, 1.);
        std::string outfile;
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-src.ppm");
        savep2or3<float>(outfile.c_str(), outs, false);
        
        for(int idx = 0; idx < 3; idx ++)
          outs2[idx] = showMatch<double>(data[idx].template cast<double>(), shape0, mhull0).template cast<float>();
        normalize<float>(outs2, 1.);
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-dst.ppm");
        savep2or3<float>(outfile.c_str(), outs2, false);
        
        for(int idx = 0; idx < 3; idx ++)
          outs3[idx] = outs[idx] + outs2[idx];
        normalize<float>(outs3, 1.);
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-match.ppm");
        savep2or3<float>(outfile.c_str(), outs3, false);
        
        for(int idx = 0; idx < 3; idx ++)
          outs4[idx] = redig.emphasis(data[idx].template cast<double>(), bump0.template cast<double>(), shape0, shape1, matches[n], hull1, float(1.) - emph).template cast<float>();
        normalize<float>(outs4, 1.);
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-emphasis-0.ppm");
        savep2or3<float>(outfile.c_str(), outs4, false);
        for(int idx = 0; idx < 3; idx ++)
          outs5[idx] = redig.emphasis(data[idx].template cast<double>(), bump1.template cast<double>(), shape0, shape1, matches[n], hull1, float(1.) + emph).template cast<float>();
        normalize<float>(outs5, 1.);
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-emphasis-2.ppm");
        savep2or3<float>(outfile.c_str(), outs5, false);
      }
    }
    return 0;
  case 10:
    {
      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> bump;
      std::vector<Eigen::Matrix<double, 3, 1> > datapoly;
      std::vector<Eigen::Matrix<int, 3, 1> >    polynorms;
      if(!loadobj<double>(datapoly, polynorms, argv[4]))
        return - 2;
      // XXX configure me:
      if(datapoly.size() > 2000) {
        std::cerr << "Too many vertices." << std::endl;
        return - 2;
      }
      float zrs(1.);
      float zre(.25);
      int   zrl(4);
      int   nshow(6);
      float emph(.75);
      PseudoBump<double> bumper;
      bumper.vmax = 60;
      std::vector<Eigen::Matrix<double, 3, 1> > shape;
      std::vector<Eigen::Matrix<int, 3, 1> > poly;
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> sute;
      bump = bumper.getPseudoBumpVec(rgb2l(data).template cast<double>(), shape, poly, sute, true).template cast<float>();
      Eigen::Matrix<double, 3, 1> zero3;
      zero3[0] = zero3[1] = zero3[2] = float(0);
      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> zero(data[0].rows(), data[0].cols());
      for(int i = 0; i < zero.rows(); i ++)
        for(int j = 0; j < zero.cols(); j ++)
          zero(i, j) = float(0);
      Eigen::Matrix<double, 3, 3> I3;
      for(int i = 0; i < 3; i ++)
        for(int j = 0; j < 3; j ++)
          I3(i, j) = (i == j ? float(1) : float(0));
      matchPartialPartial<double> statmatch;
      std::vector<match_t<double> > matches;
      for(float zr = zrs;
          (zrs / zre < float(1) && zr < zre) || 
          (zre / zrs < float(1) && zre < zr);
          zr *= pow(zre / zrs, float(1) / float(zrl))) {
        std::vector<Eigen::Matrix<double, 3, 1> > sshape(shape);
        for(int i = 0; i < shape.size(); i ++)
          sshape[i][2] *= zr;
        std::vector<match_t<double> > lmatches(statmatch.match(sshape, datapoly));
        std::copy(lmatches.begin(), lmatches.end(), std::back_inserter(matches));
      }
      float zr(zrs);
      std::sort(matches.begin(), matches.end(), cmpwrap<double>);
      for(int n = 0; n < min(int(matches.size()), nshow); n ++) {
        Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> outs[3], outs2[3], outs3[3], outs4[3], outs5[3];
        tilter<double> tilt;
        reDig<double>  redig;
        tilt.initialize(zr * sqrt((data[0].rows() * data[0].cols()) / (data[0].rows() * data[0].cols())));
        cerr << "Writing " << n << " / " << matches.size() << "(" << float(1) / matches[n].rdepth << ", " << matches[n].rpoints << ", " << matches[n].ratio << ")" << endl;
        
        vector<Eigen::Matrix<double, 3, 1> > shape3d;
        for(int idx = 0; idx < datapoly.size(); idx ++)
          shape3d.push_back(matches[n].rot * matches[n].ratio * datapoly[idx] + matches[n].offset);
        for(int idx = 0; idx < 3; idx ++)
          outs[idx] = showMatch<double>(zero.template cast<double>(), shape3d, polynorms).template cast<float>();
        normalize<float>(outs, 1.);
        std::string outfile;
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-src.ppm");
        savep2or3<float>(outfile.c_str(), outs, false);
        
        std::vector<Eigen::Matrix<int, 3, 1> > ch(loadBumpSimpleMesh<double>(datapoly, matches[n].dstpoints));
        std::vector<Eigen::Matrix<int, 3, 1> > mch;
        for(int idx = 0; idx < ch.size(); idx ++) {
          Eigen::Matrix<int, 3, 1> buf;
          for(int idx2 = 0; idx2 < ch[idx].size(); idx2 ++)
            buf[idx2] = matches[n].dstpoints[ch[idx][idx2]];
          mch.push_back(buf);
        }
        for(int idx = 0; idx < 3; idx ++)
          outs2[idx] = showMatch<double>(data[idx].template cast<double>(), shape, mch).template cast<float>();
        normalize<float>(outs2, 1.);
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-dst.ppm");
        savep2or3<float>(outfile.c_str(), outs2, false);
        
        for(int idx = 0; idx < 3; idx ++)
          outs3[idx] = outs[idx] + outs2[idx];
        normalize<float>(outs3, 1.);
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-match.ppm");
        savep2or3<float>(outfile.c_str(), outs3, false);
        
        for(int idx = 0; idx < 3; idx ++)
          outs4[idx] = redig.emphasis(data[idx].template cast<double>(), bump.template cast<double>(), shape, datapoly, matches[n], ch, float(1) - emph).template cast<float>();
        normalize<float>(outs4, 1.);
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-emphasis-0.ppm");
        savep2or3<float>(outfile.c_str(), outs4, false);
        for(int idx = 0; idx < 3; idx ++)
          outs5[idx] = redig.emphasis(data[idx].template cast<double>(), bump.template cast<double>(), shape, datapoly, matches[n], ch, float(1) + emph).template cast<float>();
        normalize<float>(outs5, 1.);
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-emphasis-2.ppm");
        savep2or3<float>(outfile.c_str(), outs5, false);
      }
    }
    return 0;
  default:
    break;
  }
  normalize<float>(data, 1.);
  if(!savep2or3<float>(argv[3], data, ! true))
    return - 3;
  return 0;
}

