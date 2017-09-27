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
using std::ofstream;
using std::getline;
using std::istringstream;
using std::stringstream;
using std::vector;
using std::pair;
using std::floor;
template <typename T> class match_t;
template <typename T> class tilter;

template <typename T> bool clockwise(const Eigen::Matrix<T, 2, 1>& p0, const Eigen::Matrix<T, 2, 1>& p1, const Eigen::Matrix<T, 2, 1>& p2) {
  Eigen::Matrix<T, 2, 1> p10(p1 - p0), p20(p2 - p0);
  const T Pi(T(4) * atan2(T(1), T(1)));
  T theta((atan2(p10[1], p10[0]) - atan2(p20[1], p20[0]) + Pi) / (T(2) * Pi));
  theta -= floor(theta);
  theta += T(1);
  theta -= floor(theta);
  if(theta * T(2) * Pi - Pi < T(0))
    return true;
  return false;
}

// N.B. delaunay trianglation algorithm best works but this isn't.
//      terrible inefficient and imcomplete temporary stub.
template <typename T> vector<Eigen::Matrix<int, 3, 1> > loadBumpSimpleMesh(const vector<Eigen::Matrix<T, 3, 1> >& dst, const vector<int>& dstpoints) {
  vector<Eigen::Matrix<int, 3, 1> > res;
  for(int i = 0; i < dstpoints.size(); i ++) {
    cerr << "mesh: " << i << "/" << dstpoints.size() << endl;
    for(int j = i + 1; j < dstpoints.size(); j ++)
      for(int k = j + 1; k < dstpoints.size(); k ++) {
        int ii(i), jj(j), kk(k);
        Eigen::Matrix<T, 2, 1> p0, p1, p2;
        p0[0] = dst[dstpoints[ii]][0];
        p0[1] = dst[dstpoints[ii]][1];
        p1[0] = dst[dstpoints[jj]][0];
        p1[1] = dst[dstpoints[jj]][1];
        p2[0] = dst[dstpoints[kk]][0];
        p2[1] = dst[dstpoints[kk]][1];
        // make sure counter clockwise.
        if(clockwise(p0, p1, p2)) {
          jj = k;
          kk = j;
        }
        // XXX have a glitch, on the same line.
        if((p0[0] == p1[0] && p1[0] == p2[0]) ||
           (p0[1] == p1[1] && p1[1] == p2[1]))
          continue;
        bool flag(false);
        // brute force.
        for(int l = 0; l < dstpoints.size() && !flag; l ++) {
          Eigen::Matrix<T, 4, 4> dc;
          dc(0, 0) = dst[dstpoints[ii]][0];
          dc(0, 1) = dst[dstpoints[ii]][1];
          dc(1, 0) = dst[dstpoints[jj]][0];
          dc(1, 1) = dst[dstpoints[jj]][1];
          dc(2, 0) = dst[dstpoints[kk]][0];
          dc(2, 1) = dst[dstpoints[kk]][1];
          dc(3, 0) = dst[dstpoints[l]][0];
          dc(3, 1) = dst[dstpoints[l]][1];
          for(int m = 0; m < dc.rows(); m ++) {
            dc(m, 2) = dc(m, 0) * dc(m, 0) + dc(m, 1) * dc(m, 1);
            dc(m, 3) = T(1);
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

template <typename T> bool saveobj(const vector<Eigen::Matrix<T, 3, 1> >& data, const vector<Eigen::Matrix<int, 3, 1> >& polys, const char* filename) {
  ofstream output;
  output.open(filename, std::ios::out);
  if(output.is_open()) {
    for(int i = 0; i < data.size(); i ++)
      output << "v " << data[i][0] << " " << data[i][1] << " " << data[i][2] * T(8) << endl;
    for(int i = 0; i < polys.size(); i ++)
      output << "f " << polys[i][0] + 1 << " " << polys[i][1] + 1 << " " << polys[i][2] + 1 << endl;
    output.close();
  } else {
    cerr << "Unable to open file: " << filename << endl;
    return false;
  }
  return true;
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
          sub >> wbuf[widx];
          if(wbuf[widx] >= 0)
            wbuf[widx] --;
          widx ++;
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
        polys[i][j] = ((polys[i][j] % (2 * data.size())) + 2 * data.size()) % data.size();
    input.close();
  } else {
    cerr << "Unable to open file: " << filename << endl;
    return false;
  }
  return true;
}

#define _OBJ2VECTOR_
#endif

