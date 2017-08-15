#include <cstdio>
#include <cmath>
#include <iostream>
#include <Eigen/Core>
#include "ppm2eigen.hh"
#include "enlarge.hh"
#include "fisheye.hh"
#include "edgedetect.hh"
#include "tilt.hh"
#include "scancontext.hh"

using namespace std;

void usage() {
  cout << "Usage: tools (enlarge|enlargeds|bump|detect|collect|tilt|lpoly|match) <input filename>.p[gp]m <output filename>.p[gp]m <args>?" << endl;
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
  else if(strcmp(argv[1], "detect") == 0)
    mode = 3;
  else if(strcmp(argv[1], "collect") == 0)
    mode = 4;
  else if(strcmp(argv[1], "tilt") == 0)
    mode = 6;
  else if(strcmp(argv[1], "lpoly") == 0)
    mode = 8;
  else if(strcmp(argv[1], "match") == 0)
    mode = 9;
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
  case 3:
    {
      edgedetect<float> detect;
      for(int i = 0; i < 3; i ++)
        data[i] = detect.detect(data[i], edgedetect<float>::DETECT_BOTH);
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
      const int M_TILT = 32;
      for(int i = 0; i < M_TILT; i ++) {
        Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> out[3];
        for(int j = 0; j < 3; j ++)
          out[j] = tilt.tilt(data[j], bump[0], i, M_TILT, .995);
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
      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> work(bump.rgb2l(data));
      data[0] = lf.getLowFreqImage(detect.detect(work, edgedetect<float>::COLLECT_BOTH));
      for(int i = 1; i < 3; i ++)
        data[i] = data[0];
      std::vector<Eigen::Matrix<float, 3, 1> > points(lf.getLowFreq(work));
      std::cout << "Handled points:" << std::endl;
      for(int i = 0; i < points.size(); i ++)
        std::cout << points[i].transpose() << std::endl;
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
      float thresh_para(.9995);
      float thresh_points(.025);
      float thresh_r(.125);
      float zr(.125);
      // bump to bump match.
      float r_max_theta(.01);
      int   div(120);
      int   nshow(6);
      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> mout[3];
      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> mbump;
      Eigen::Matrix<float, 3, 3> I3;
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
      for(int i = 0; i < shape0.size(); i ++) {
        shape0[i][2] *= zr;
        shape1[i][2] *= zr;
      }
      statmatch.init(shape0, thresh_para, thresh_points, thresh_r);
      std::vector<match_t<float> > matches(statmatch.match(shape1, div, r_max_theta));
      for(int n = 0; n < min(int(matches.size()), nshow); n ++) {
        Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> outs[3], outs2[3], outs3[3], outs4[3], outs5[3];
        tilter<float> tilt;
        reDig<float>  redig;
        tilt.initialize(zr * sqrt((data1[0].rows() * data1[0].cols()) / (data[0].rows() * data[0].cols())));
        cerr << "Writing " << n << " / " << matches.size() << "(" << float(1) / matches[n].rdepth << ", " << matches[n].rpoints << ", " << matches[n].ratio << ")" << endl;
        for(int idx = 0; idx < 3; idx ++)
          outs[idx] = tilt.tilt(showMatch<float>(mout[idx], shape1, matches[n].srcpoints), mbump, matches[n].rot, I3, matches[n].offset, matches[n].ratio, zero3);
        normalize<float>(outs, 1.);
        std::string outfile;
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-src.ppm");
        savep2or3<float>(outfile.c_str(), outs, false);
        
        for(int idx = 0; idx < 3; idx ++)
          outs2[idx] = showMatch<float>(data3[idx], shape0, matches[n].dstpoints);
        normalize<float>(outs2, 1.);
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-dst.ppm");
        savep2or3<float>(outfile.c_str(), outs2, false);
        
        for(int idx = 0; idx < 3; idx ++)
          outs3[idx] = outs[idx] + outs2[idx];
        normalize<float>(outs3, 1.);
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-match.ppm");
        savep2or3<float>(outfile.c_str(), outs3, false);
        
        for(int idx = 0; idx < 3; idx ++) {
          outs4[idx] = redig.emphasis(data3[idx], data[idx], shape0, shape1, matches[n], float(- 2));
          outs5[idx] = redig.emphasis(data3[idx], data[idx], shape0, shape1, matches[n], float(2));
        }
        normalize<float>(outs4, 1.);
        normalize<float>(outs5, 1.);
        outfile = std::string(argv[3]) + std::to_string(n + 1) + std::string("-emphasis-0.ppm");
        savep2or3<float>(outfile.c_str(), outs4, false);
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

