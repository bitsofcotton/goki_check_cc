#if !defined(_OBJ2VECTOR_)

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using std::getline;
using std::istringstream;
using std::stringstream;
using std::vector;
using std::pair;
template <typename T> class match_t;
template <typename T> class tilter;

// N.B. delaunay trianglation algorithm best works but this isn't.
//      terrible inefficient and imcomplete temporary stub.
template <typename T> vector<Eigen::Matrix<int, 3, 1> > loadBumpSimpleMesh(const vector<Eigen::Matrix<T, 3, 1> >& dst, const vector<int>& dstpoints) {
  vector<Eigen::Matrix<int, 3, 1> > res;
  for(int i = 0; i < dstpoints.size(); i ++) {
    cerr << "mesh: " << i << "/" << dstpoints.size() << endl;
    for(int j = i + 1; j < dstpoints.size(); j ++)
      // when sorted order, brute force can be writed.
      for(int k = j + 1; k < dstpoints.size(); k ++) {
        bool flag(false);
        for(int l = 0; l < dstpoints.size(); l ++) {
          // make sure counter clockwise.
          int ii(i), jj(j), kk(k);
          Eigen::Matrix<T, 2, 1> p0, p1, p2;
          p0[0] = dst[dstpoints[ii]][0];
          p0[1] = dst[dstpoints[ii]][1];
          p1[0] = dst[dstpoints[jj]][0];
          p1[1] = dst[dstpoints[jj]][1];
          p2[0] = dst[dstpoints[kk]][0];
          p2[1] = dst[dstpoints[kk]][1];
          Eigen::Matrix<T, 2, 1> p10(p1 - p0), p20(p2 - p0);
          if(atan2(p10[1], p10[0]) - atan2(p20[1], p20[0]) > T(0)) {
            ii = j;
            jj = i;
          }
          Eigen::Matrix<T, 4, 4> dc;
          dc(0, 0) = dst[dstpoints[ii]][0];
          dc(0, 1) = dst[dstpoints[ii]][1];
          dc(1, 0) = dst[dstpoints[jj]][0];
          dc(1, 1) = dst[dstpoints[jj]][1];
          dc(2, 0) = dst[dstpoints[kk]][0];
          dc(2, 1) = dst[dstpoints[kk]][1];
          dc(3, 0) = dst[dstpoints[l]][0];
          dc(3, 1) = dst[dstpoints[l]][1];
          for(int ii = 0; ii < dc.rows(); ii ++) {
            dc(ii, 2) = dc(ii, 0) * dc(ii, 0) + dc(ii, 1) * dc(ii, 1);
            dc(ii, 3) = T(1);
          }
          if(dc.determinant() > T(0)) {
            flag = true;
            break;
          }
        }
        if(!flag) {
          Eigen::Matrix<int, 3, 1> buf;
          buf[0] = i;
          buf[1] = j;
          buf[2] = k;
          res.push_back(buf);
        }
      }
  }
  return res;
}

template <typename T> bool loadobj(vector<Eigen::Matrix<T, 3, 1> >& data, vector<Eigen::Matrix<int, 3, 1> >& polys, const char* filename) {
  ifstream input;
  input.open(filename);
  if(input.is_open()) {
    string work;
    while(getline(input, work)) {
      int i = 0;
      for(; i < work.size() && work[i] == ' '; i ++);
      if(i + 1 < work.size() && work[i] == 'v' && work[i + 1] == ' ') {
        stringstream sub(work.substr(i + 2, work.size() - (i + 2)));
        Eigen::Matrix<T, 3, 1> buf;
        sub >> buf[0];
        sub >> buf[1];
        sub >> buf[2];
        data.push_back(buf);
      } else if(i < work.size() && work[i] == 'f') {
        stringstream sub(work.substr(i + 1, work.size() - (i + 1)));
        Eigen::Matrix<int, 3, 1> wbuf;
        int  widx(0);
        bool flag(false);
        while(!sub.eof() && !sub.bad()) {
          sub >> wbuf[widx ++];
          if(sub.eof() || sub.bad())
            break;
          if(widx > 2)
            flag = true;
          widx %= 3;
          if(flag)
            polys.push_back(wbuf);
          sub.ignore(20, ' ');
        }
      }
    }
    for(int i = 0; i < polys.size(); i ++)
      for(int j = 0; j < polys[i].size(); j ++)
        // XXX fixme < or >, ++ or ;.
        polys[i][j] = ((polys[i][j] % (2 * data.size())) + 2 * data.size()) % data.size();
  } else {
    cerr << "Unable to open file: " << filename << endl;
    return false;
  }
  return true;
}

#define _OBJ2VECTOR_
#endif

