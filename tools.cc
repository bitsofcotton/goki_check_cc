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
      for(int i = 0; i < 3; i ++) {
        data[i] = enlarger.enlarge2(data[i], enlarger2ex<float>::ENLARGE_BOTH);
      }
    }
    break;
  case 1:
    {
      enlarger2exds<float> enlarger;
      for(int i = 0; i < 3; i ++)
        data[i] = enlarger.enlarge2ds(data[i], enlarger2exds<float>::ENLARGE_BOTH);
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
      PseudoBump<float> bump;
      data[0] = bump.getPseudoBump(bump.rgb2l(data), false);
      data[1] = data[0];
      data[2] = data[0];
    }
    break;
  case 6:
    {
      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> bump[3];
      if(!loadp2or3<float>(bump, argv[4]))
        return - 2;
      tilter<float> tilt;
      tilt.initialize(8.);
      const int M_TILT = 32;
      for(int i = 0; i < M_TILT; i ++) {
        Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> out[3];
        for(int j = 0; j < 3; j ++)
          out[j] = tilt.tilt(data[j], bump[0], i, M_TILT, .999);
        std::string outfile(argv[3]);
        outfile += std::string("-") + std::to_string(i + 1) + std::string(".ppm");
        savep2or3<float>(outfile.c_str(), out, false);
      }
      return 0;
    }
    break;
  case 8:
    {
      lowFreq<float>    lf;
      PseudoBump<float> bump;
      edgedetect<float> detect;
      std::vector<Eigen::Matrix<float, 3, 1> > points(lf.getLowFreq(bump.rgb2l(data), 30));
      std::vector<int> dstpoints;
      for(int i = 0; i < points.size(); i ++)
        dstpoints.push_back(i);
      normalize<float>(data, 1.);
      std::vector<Eigen::Matrix<int, 3, 1> > pptr(loadBumpSimpleMesh<float>(points, dstpoints));
      for(int i = 0; i < 3; i ++)
        data[i] = showMatch<float>(data[i], points, pptr);
      std::cout << "Handled points:" << std::endl;
      for(int i = 0; i < points.size(); i ++)
        std::cout << points[i].transpose() << std::endl;
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
      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> data1[3], data2[3], data3[3];
      if(!loadp2or3<float>(data1, argv[4]))
        return - 2;
      if(!loadp2or3<float>(data2, argv[5]))
        return - 2;
      if(!loadp2or3<float>(data3, argv[6]))
        return - 2;
      // XXX: configure me.
      float thresh_para(.9);
      float thresh_len(.2);
      float thresh_points(.125);
      float thresh_r(.125);
      float zrs(1.);
      float zre(.25);
      int   zrl(4);
      // bump to bump match.
      float r_max_theta(.01);
      int   div(20);
      int   nshow(6);
      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> mout[3];
      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> mbump;
      Eigen::Matrix<float, 3, 3> I3;
      float emph(.1);
      mbump = mout[0] = mout[1] = mout[2] = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>(data3[0].rows(), data3[0].cols());
      for(int i = 0; i < 3; i ++)
        for(int j = 0; j < 3; j ++)
          I3(i, j) = (i == j ? float(1) : float(0));
      for(int i = 0; i < min(mout[0].rows(), data2[0].rows()); i ++) {
        for(int j = 0; j < min(mout[0].cols(), data2[0].cols()); j ++) {
          mout[0](i, j) = data2[0](i, j);
          mout[1](i, j) = data2[1](i, j);
          mout[2](i, j) = data2[2](i, j);
          mbump(i, j)   = data1[0](i, j);
        }
        for(int j = min(mout[0].cols(), data2[0].cols()); j < mout[0].cols(); j ++)
          mout[0](i, j) = mout[1](i, j) = mout[2](i, j) = mbump(i, j) = float(0);
      }
      
      for(int i = min(mout[0].rows(), data2[0].rows()); i < mout[0].rows(); i ++)
        for(int j = 0; j < mout[0].cols(); j ++)
          mout[0](i, j) = mout[1](i, j) = mout[2](i, j) = mbump(i, j) = float(0);
      Eigen::Matrix<float, 3, 1> zero3;
      zero3[0] = zero3[1] = zero3[2] = float(0);
      matchPartialPartial<float> statmatch;
      PseudoBump<float> bump;
      lowFreq<float> lf;
      std::vector<Eigen::Matrix<float, 3, 1> > shape0(lf.getLowFreq(bump.rgb2l(data)));
      std::vector<Eigen::Matrix<float, 3, 1> > shape1(lf.getLowFreq(bump.rgb2l(data1)));
      std::vector<match_t<float> > matches;
      for(float zr = zrs;
          (zrs / zre < float(1) && zr < zre) || 
          (zre / zrs < float(1) && zre < zr);
          zr *= pow(zre / zrs, float(1) / float(zrl)))
        for(float zr2 = zrs;
            (zrs / zre < float(1) && zr2 < zre) ||
            (zre / zrs < float(1) && zre < zr2);
            zr2 *= pow(zre / zrs, float(1) / float(zrl))) {
          std::vector<Eigen::Matrix<float, 3, 1> > sshape0, sshape1;
          for(int i = 0; i < shape0.size(); i ++) {
            sshape0.push_back(shape0[i]);
            sshape1.push_back(shape1[i]);
            sshape0[i][2] *= zr;
            sshape1[i][2] *= zr2;
          }
          statmatch.init(sshape0, thresh_para, thresh_len, thresh_points, thresh_r);
          std::vector<match_t<float> > lmatches(statmatch.match(sshape1, div, r_max_theta));
          std::copy(lmatches.begin(), lmatches.end(), std::back_inserter(matches));
      }
      std::sort(matches.begin(), matches.end(), cmpwrap<float>);
      float zr(zrs);
      for(int n = 0; n < min(int(matches.size()), nshow); n ++) {
        Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> outs[3], outs2[3], outs3[3], outs4[3], outs5[3];
        tilter<float> tilt;
        reDig<float>  redig;
        tilt.initialize(zr * sqrt((data1[0].rows() * data1[0].cols()) / (data[0].rows() * data[0].cols())));
        std::vector<Eigen::Matrix<int, 3, 1> > hull1(loadBumpSimpleMesh<float>(shape1, matches[n].srcpoints));
        std::vector<Eigen::Matrix<int, 3, 1> > mhull0, mhull1;
        for(int idx = 0; idx < hull1.size(); idx ++) {
          Eigen::Matrix<int, 3, 1> buf;
          for(int idx2 = 0; idx2 < hull1[idx].size(); idx2 ++)
            buf[idx2] = matches[n].dstpoints[hull1[idx][idx2]];
          mhull0.push_back(buf);
        }
        for(int idx = 0; idx < hull1.size(); idx ++) {
          Eigen::Matrix<int, 3, 1> buf;
          for(int idx2 = 0; idx2 < hull1[idx].size(); idx2 ++)
            buf[idx2] = matches[n].srcpoints[hull1[idx][idx2]];
          mhull1.push_back(buf);
        }
        cerr << "Writing " << n << " / " << matches.size() << "(" << float(1) / matches[n].rdepth << ", " << matches[n].rpoints << ", " << matches[n].ratio << ")" << endl;
        for(int idx = 0; idx < 3; idx ++)
          outs[idx] = tilt.tilt(showMatch<float>(mout[idx], shape1, mhull1), mbump, matches[n].rot, I3, matches[n].offset, matches[n].ratio, zero3);
        normalize<float>(outs, 1.);
        std::string outfile;
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-src.ppm");
        savep2or3<float>(outfile.c_str(), outs, false);
        
        for(int idx = 0; idx < 3; idx ++)
          outs2[idx] = showMatch<float>(data3[idx], shape0, mhull0);
        normalize<float>(outs2, 1.);
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-dst.ppm");
        savep2or3<float>(outfile.c_str(), outs2, false);
        
        for(int idx = 0; idx < 3; idx ++)
          outs3[idx] = outs[idx] + outs2[idx];
        normalize<float>(outs3, 1.);
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-match.ppm");
        savep2or3<float>(outfile.c_str(), outs3, false);
        
        for(int idx = 0; idx < 3; idx ++)
          outs4[idx] = redig.emphasis(data3[idx], data[idx], shape0, shape1, matches[n], hull1, float(1.) - emph);
        normalize<float>(outs4, 1.);
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-emphasis-0.ppm");
        savep2or3<float>(outfile.c_str(), outs4, false);
        for(int idx = 0; idx < 3; idx ++)
          outs5[idx] = redig.emphasis(data3[idx], data[idx], shape0, shape1, matches[n], hull1, float(1.) + emph);
        normalize<float>(outs5, 1.);
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-emphasis-2.ppm");
        savep2or3<float>(outfile.c_str(), outs5, false);
      }
    }
    return 0;
  case 10:
    {
      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> data1[3];
      std::vector<Eigen::Matrix<float, 3, 1> > datapoly;
      std::vector<Eigen::Matrix<int, 3, 1> > polynorms;
      if(!loadp2or3<float>(data1, argv[4]))
        return - 2;
      if(!loadobj<float>(datapoly, polynorms, argv[5]))
        return - 2;
      Eigen::Matrix<float, 3, 1> zero3;
      zero3[0] = zero3[1] = zero3[2] = float(0);
      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> zero(data1[0].rows(), data1[0].cols());
      for(int i = 0; i < zero.rows(); i ++)
        for(int j = 0; j < zero.cols(); j ++)
          zero(i, j) = float(0);
      Eigen::Matrix<float, 3, 3> I3;
      for(int i = 0; i < 3; i ++)
        for(int j = 0; j < 3; j ++)
          I3(i, j) = (i == j ? float(1) : float(0));
      
      // XXX: configure me.
      float thresh_para(.9);
      float thresh_len(.2);
      float thresh_points(.25);
      float thresh_r(.125);
      float zrs(1.);
      float zre(.25);
      int   zrl(4);
      // bump to bump match.
      float r_max_theta(.01);
      int   div(20);
      int   nshow(6);
      float emph(.1);
      matchPartialPartial<float> statmatch;
      PseudoBump<float> bump;
      lowFreq<float> lf;
      std::vector<Eigen::Matrix<float, 3, 1> > shape(lf.getLowFreq(bump.rgb2l(data)));
      std::vector<match_t<float> > matches;
      for(float zr = zrs;
          (zrs / zre < float(1) && zr < zre) || 
          (zre / zrs < float(1) && zre < zr);
          zr *= pow(zre / zrs, float(1) / float(zrl))) {
        std::vector<Eigen::Matrix<float, 3, 1> > sshape(shape);
        for(int i = 0; i < shape.size(); i ++)
          sshape[i][2] *= zr;
        statmatch.init(sshape, thresh_para, thresh_len, thresh_points, thresh_r);
        std::vector<match_t<float> > lmatches(statmatch.match(datapoly, div, r_max_theta));
        std::copy(lmatches.begin(), lmatches.end(), std::back_inserter(matches));
      }
      float zr(zrs);
      std::sort(matches.begin(), matches.end(), cmpwrap<float>);
      for(int n = 0; n < min(int(matches.size()), nshow); n ++) {
        Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> outs[3], outs2[3], outs3[3], outs4[3], outs5[3];
        tilter<float> tilt;
        reDig<float>  redig;
        tilt.initialize(zr * sqrt((data1[0].rows() * data1[0].cols()) / (data[0].rows() * data[0].cols())));
        cerr << "Writing " << n << " / " << matches.size() << "(" << float(1) / matches[n].rdepth << ", " << matches[n].rpoints << ", " << matches[n].ratio << ")" << endl;
        
        vector<Eigen::Matrix<float, 3, 1> > shape3d;
        for(int idx = 0; idx < datapoly.size(); idx ++)
          shape3d.push_back(matches[n].rot * matches[n].ratio * datapoly[idx] + matches[n].offset);
        for(int idx = 0; idx < 3; idx ++)
          outs[idx] = showMatch<float>(zero, shape3d, polynorms);
        normalize<float>(outs, 1.);
        std::string outfile;
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-src.ppm");
        savep2or3<float>(outfile.c_str(), outs, false);
        
        std::vector<Eigen::Matrix<int, 3, 1> > ch(loadBumpSimpleMesh<float>(datapoly, matches[n].srcpoints));
        std::vector<Eigen::Matrix<int, 3, 1> > mch;
        for(int idx = 0; idx < ch.size(); idx ++) {
          Eigen::Matrix<int, 3, 1> buf;
          for(int idx2 = 0; idx2 < ch[idx].size(); idx2 ++)
            buf[idx2] = matches[n].dstpoints[ch[idx][idx2]];
          mch.push_back(buf);
        }
        for(int idx = 0; idx < 3; idx ++)
          outs2[idx] = showMatch<float>(data1[idx], shape, mch);
        normalize<float>(outs2, 1.);
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-dst.ppm");
        savep2or3<float>(outfile.c_str(), outs2, false);
        
        for(int idx = 0; idx < 3; idx ++)
          outs3[idx] = outs[idx] + outs2[idx];
        normalize<float>(outs3, 1.);
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-match.ppm");
        savep2or3<float>(outfile.c_str(), outs3, false);
        
        for(int idx = 0; idx < 3; idx ++)
          outs4[idx] = redig.emphasis(data1[idx], data[idx], shape, datapoly, matches[n], ch, float(1) - emph);
        normalize<float>(outs4, 1.);
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-emphasis-0.ppm");
        savep2or3<float>(outfile.c_str(), outs4, false);
        for(int idx = 0; idx < 3; idx ++)
          outs5[idx] = redig.emphasis(data1[idx], data[idx], shape, datapoly, matches[n], ch, float(1) + emph);
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

