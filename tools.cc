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
      const int M_TILT = 4;
      tilt.initialize(16);
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
      lowFreq<float> lf;
      PseudoBump<float> bump;
      for(int i = 0; i < 3; i ++)
        data[i] = lf.getLowFreqImage(data[i]);
      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> work(bump.rgb2l(data));
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
      matchPartialPartial<float> statmatch;
      PseudoBump<float> bump;
      lowFreq<float> lf;
      std::vector<Eigen::Matrix<float, 3, 1> > shape0(lf.getLowFreq(bump.rgb2l(data), 20));
      std::vector<Eigen::Matrix<float, 3, 1> > shape1(lf.getLowFreq(bump.rgb2l(data1), 20));
      for(int i = 0; i < shape0.size(); i ++) {
        shape0[i][2] *= sqrt(float(data[0].rows())  * data[0].cols());
        shape1[i][2] *= sqrt(float(data1[0].rows()) * data1[0].cols());
      }
      statmatch.init(shape0, .85, .25);
      std::vector<match_t<float> > matches(statmatch.match(shape1, 20));
      for(int n = 0; n < min(int(matches.size()), 16); n ++) {
        Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> outs[3], zero2(data2[0].rows(), data2[0].cols());
        Eigen::Matrix<float, 3, 3> I3;
        for(int i = 0; i < 3; i ++)
          for(int j = 0; j < 3; j ++)
            I3(i, j) = (i == j ? float(1) : float(0));
        for(int i = 0; i < zero2.rows(); i ++)
          for(int j = 0; j < zero2.cols(); j ++)
            zero2(i, j) = float(0);
        tilter<float> tilt;
        tilt.initialize(16);
        cerr << "Writing " << n << " / " << matches.size() << endl;
        for(int idx = 0; idx < 3; idx ++)
          outs[idx] = tilt.tilt(showMatch<float>(data2[idx], shape1, matches[n].srcpoints), zero2, matches[n].rot, I3, matches[n].offset, matches[n].ratio);;
        //  outs[idx] = showMatch<float>(data2[idx], shape1, matches[n].srcpoints);
        normalize<float>(outs, 1.);
        std::string outfile;
        outfile = std::string(argv[3]) + std::string("-src-") + std::to_string(n + 1) + std::string(".ppm");
        savep2or3<float>(outfile.c_str(), outs, false);
        
        for(int idx = 0; idx < 3; idx ++)
          outs[idx] = showMatch<float>(data3[idx], shape0, matches[n].dstpoints);
        normalize<float>(outs, 1.);
        outfile = std::string(argv[3]) + std::string("-dst-") + std::to_string(n + 1) + std::string(".ppm");
        savep2or3<float>(outfile.c_str(), outs, false);
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

