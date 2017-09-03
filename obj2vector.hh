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
template <typename T> class match_t;
template <typename T> class tilter;

template <typename T> vector<Eigen::Matrix<int, 3, 1> > loadConvexHull(const vector<Eigen::Matrix<T, 3, 1> >& dst, const match_t<T>& match) {
  tilter<T> tilt;
  vector<Eigen::Matrix<int, 3, 1> > res;
  // convex hull of normal vector.
  for(int i = 0; i < match.dstpoints.size(); i ++) {
    cerr << "convexHull : " << i << "/" << match.dstpoints.size() << endl;
    for(int j = i + 1; j < match.dstpoints.size(); j ++)
      for(int k = max(i, j) + 1; k < match.dstpoints.size(); k ++) {
        Eigen::Matrix<T, 3, 1> n(tilt.solveN(dst[match.dstpoints[i]],
                                             dst[match.dstpoints[j]],
                                             dst[match.dstpoints[k]]));
        bool flag = true;
        for(int l = 1; l < match.dstpoints.size(); l ++)
          if(n.dot(dst[match.dstpoints[0]]) *
             n.dot(dst[match.dstpoints[l]]) < T(0)) {
            flag = false;
            break;
          } else if(n.dot(dst[match.dstpoints[l]]) < T(0))
            n = - n;
        if(!flag || n[2] <= T(0))
          continue;
        Eigen::Matrix<int, 3, 1> buf;
        buf[0] = i;
        buf[1] = j;
        buf[2] = k;
        res.push_back(buf);
      }
  }
  return res;
}

template <typename T> bool loadobj(vector<Eigen::Matrix<T, 3, 1> >& data, const char* filename) {
  ifstream input;
  input.open(filename);
  if(input.is_open()) {
    string work;
    while(getline(input, work)) {
      int i = 0;
      for(; i < work.size() && work[i] == ' '; i ++);
      if(i < work.size() && work[i] == 'v' &&
         i + 1 < work.size() && work[i + 1] != 'f' && work[i + 1] != 'n') {
        stringstream sub(work.substr(i + 2, work.size() - (i + 2)));
        Eigen::Matrix<T, 3, 1> buf;
        sub >> buf[0];
        sub >> buf[1];
        sub >> buf[2];
        data.push_back(buf);
      }
    }
  } else {
    cerr << "Unable to open file: " << filename << endl;
    return false;
  }
  return true;
 
}

#define _OBJ2VECTOR_
#endif

