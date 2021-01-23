/* BSD 3-Clause License:
 * Copyright (c) 2018, bitsofcotton.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 *    Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *    Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation or other materials provided with the distribution.
 *    Neither the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#if !defined(_FILEIO_GOKI_)

using std::max;
using std::min;
using std::cerr;
using std::endl;
using std::string;
using std::ifstream;
using std::getline;
using std::istringstream;
using std::stringstream;
using std::ofstream;
using std::vector;
using std::pair;
using std::make_pair;

template <typename T> class match_t;

template <typename T> class simpleFile {
public:
#if defined(_WITHOUT_EIGEN_)
  typedef SimpleMatrix<T> Mat;
  typedef SimpleMatrix<T> Mat3x3;
  typedef SimpleVector<T> Vec3;
  typedef SimpleVector<T> Vec;
  typedef SimpleVector<int> Veci3;
  typedef SimpleVector<int> Veci4;
#else
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
  typedef Eigen::Matrix<T,   3, 3> Mat3x3;
  typedef Eigen::Matrix<T,   3, 1> Vec3;
  typedef Eigen::Matrix<int, 3, 1> Veci3;
  typedef Eigen::Matrix<int, 4, 1> Veci4;
#endif
  
  inline bool whiteline(const std::string& s) {
    for(auto ss(s.begin()); ss < s.end(); ++ ss)
      if(! std::isspace(* ss) && *ss != '\n')
        return false;
    return true;
  }

  inline bool loadstub(ifstream& input, const int& nmax, const int& ncolor, Mat* datas) {
    int i = 0, j = 0, k = 0;
    char buf;
    int  work = 0;
    bool mode = false;
    while(input.get(buf) && j < datas[0].rows()) {
      if('0' <= buf && buf <= '9') {
        work *= 10;
        work += buf - '0';
        mode  = true;
        continue;
      } else if(mode) {
        mode = false;
        datas[k](j, i) = T(work) / T(nmax);
        work = 0;
        if(++ k >= ncolor) {
          if(++ i >= datas[0].cols()) {
            if(++ j >= datas[0].rows())
              return true;
            i = 0;
          }
          k = 0;
        }
      }
    }
    return true;
  }
  
  bool loadp2or3(Mat data[3], const char* filename) {
    string line;
    string line2;
    string line3;
    ifstream input;
    input.open(filename);
    if(input.is_open()) {
      try {
        getline(input, line);
        while((whiteline(line) || line[0] == '#') && getline(input, line) && !input.eof() && !input.bad()) ;
        getline(input, line2);
        while((whiteline(line2) || line2[0] == '#') && getline(input, line2) && !input.eof() && !input.bad()) ;
        getline(input, line3);
        while((whiteline(line3) || line3[0] == '#') && getline(input, line3) && !input.eof() && !input.bad()) ;
        istringstream iline2(line2);
        int w, h;
        iline2 >> w;
        iline2 >> h;
        if(line.size() < 2 || w <= 0 || h <= 0) {
          cerr << "unknown size." << endl;
          input.close();
          return false;
        } 
        istringstream iline3(line3);
        int nmax;
        iline3 >> nmax;
        if(line[0] == 'P') {
          if(line[1] == '2') {
            data[0] = Mat(h, w);
            for(int i = 0; i < h; i ++)
              for(int j = 0; j < w; j ++)
                data[0](i, j) = T(0);
            loadstub(input, nmax, 1, data);
            data[1] = data[2] = data[0];
          } else if(line[1] == '3') {
            for(int i = 0; i < 3; i ++)
              data[i] = Mat(h, w);
            for(int i = 0; i < h; i ++)
              for(int j = 0; j < w; j ++)
                data[0](i, j) = data[1](i, j) = data[2](i, j) = T(0);
            loadstub(input, nmax, 3, data);
          } else {
            cerr << "unknown file type." << endl;
            input.close();
            return false;
          }
        } else {
          cerr << "unknown file type." << endl;
          input.close();
          return false;
        }
      } catch (...) {
        cerr << "Exception while reading." << endl;
      }
      input.close();
    } else {
      cerr << "Unable to open file for read: " << filename << endl;
      return false;
    }
    return true;
  }

  bool savep2or3(const char* filename, Mat data[3], const bool& gray, const int& depth = 255) {
    ofstream output;
    output.open(filename);
    if(output.is_open()) {
      try {
        if(gray)
          output << "P2" << "\n";
        else
          output << "P3" << "\n";
        output << data[0].cols() << " " << data[0].rows() << "\n";
        output << depth << "\n";
        for(int i = 0; i < data[0].rows(); i ++)
          for(int j = 0; j < data[0].cols(); j ++)
            if(gray)
              output << int(data[0](i, j) * T(depth)) << "\n";
            else
              for(int k = 0; k < 3; k ++) {
                output << int(data[k](i, j) * T(depth)) << "\n";
              }
      } catch (...) {
        cerr << "An error has occured while writing file." << endl;
      }
      output.close();
    } else {
      cerr << "Unable to open file for write: " << filename << endl;
      return false;
    }
    return true;
  }

  bool saveobj(const vector<Vec3>& data, const T& Mw0, const T& Mh0, const vector<Veci3>& polys, const char* filename, const vector<vector<int> >& edges = vector<vector<int> >(), const T& addstand = T(0)) {
    ofstream output;
    output.open(filename, std::ios::out);
    if(output.is_open()) {
      int lfs(0);
      const T Mh(Mh0 / T(2));
      const T Mw(Mw0 / T(2));
      for(int fslash(0) ; filename[fslash]; fslash ++)
        if(filename[fslash] == '/')
          lfs = fslash;
      if(lfs) lfs ++;
      output << "mtllib " << &filename[lfs] << ".mtl" << endl;
      output << "usemtl material0" << endl;
      T lz(1);
      for(int i = 0; i < data.size(); i ++)
        lz = min(data[i][2], lz);
      if(addstand != T(0)) {
        for(int i = 0; i < data.size(); i ++)
          output << "v " << data[i][1] << " " << Mw - data[i][0] << " " << data[i][2] + addstand - lz << endl;
        for(int i = 0; i < data.size(); i ++)
          output << "v " << data[i][1] << " " << Mw - data[i][0] << " " << 0 << endl;
      } else {
        for(int i = 0; i < data.size(); i ++)
          output << "v " << data[i][1] << " " << - data[i][0] << " " << data[i][2] << endl;
      }
      for(int i = 0; i < data.size(); i ++)
        output << "vt " << data[i][1] / T(Mh) / T(2) << " " << T(1) - data[i][0] / T(Mw) / T(2) << endl;
      if(addstand != T(0))
        for(int i = 0; i < data.size(); i ++)
          output << "vt " << 0 <<" " << 0 << endl;
      // xchg with clockwise/counter clockwise.
      for(int i = 0; i < polys.size(); i ++) {
        const int i0(polys[i][0] + 1);
        const int i1(polys[i][1] + 1);
        const int i2(polys[i][2] + 1);
        if(i0 != i1 && i1 != i2 && i2 != i0) {
          output << "f " << i0 << "/" << i0 << "/" << i0;
          output << " "  << i1 << "/" << i1 << "/" << i1;
          output << " "  << i2 << "/" << i2 << "/" << i2 << endl;
        }
      }
      if(addstand != T(0)) {
        for(int i = 0; i < polys.size(); i ++) {
          const int i0(data.size() + polys[i][0] + 1);
          const int i1(data.size() + polys[i][2] + 1);
          const int i2(data.size() + polys[i][1] + 1);
          if(i0 != i1 && i1 != i2 && i2 != i0) {
            output << "f " << i0 << "/" << i0 << "/" << i0;
            output << " "  << i1 << "/" << i1 << "/" << i1;
            output << " "  << i2 << "/" << i2 << "/" << i2 << endl;
          }
        }
        assert(0 < edges.size());
        for(int ii = 0; ii < edges.size(); ii ++) if(edges[ii].size())
          for(int i0 = 0; i0 < edges[ii].size(); i0 ++) {
            const auto& outer(edges[ii]);
            const int   i1((i0 + 1) % edges[ii].size());
            const int   ii0(data.size() + outer[i0] + 1);
            const int   ii1(outer[i1] + 1);
            const int   ii2(data.size() + outer[i1] + 1);
            const int   ii3(outer[i0] + 1);
            assert(0 <= outer[i0] && outer[i0] < data.size());
            assert(0 <= outer[i1] && outer[i1] < data.size());
            if(ii0 != ii1 && ii1 != ii2 && ii2 != ii0) {
              output << "f " << ii0 << "/" << ii0 << "/" << ii0;
              output << " "  << ii1 << "/" << ii1 << "/" << ii1;
              output << " "  << ii2 << "/" << ii2 << "/" << ii2 << endl;
            }
            if(ii0 != ii3 && ii3 != ii1 && ii1 != ii0) {
              output << "f " << ii0 << "/" << ii0 << "/" << ii0;
              output << " "  << ii3 << "/" << ii3 << "/" << ii3;
              output << " "  << ii1 << "/" << ii1 << "/" << ii1 << endl;
            }
          }
      }
      output.close();
    } else {
      cerr << "Unable to open file for write: " << filename << endl;
      return false;
    }
    return true;
  }

  bool loadobj(vector<Vec3>& data, vector<Veci3>& polys, const char* filename) {
    ifstream input;
    input.open(filename);
    if(input.is_open()) {
      string work;
      while(getline(input, work) && !input.eof() && !input.bad()) {
        int i = 0;
        for( ; i < work.size() && work[i] == ' '; i ++);
        if(i + 1 < work.size() && work[i] == 'v' && work[i + 1] == ' ') {
          stringstream sub(work.substr(i + 2, work.size() - (i + 2)));
          Vec3 buf(3);
          sub >> buf[1];
          sub >> buf[0];
          sub >> buf[2];
          buf[0] = - buf[0];
          data.emplace_back(buf);
        } else if(i + 1 < work.size() && work[i] == 'f' && work[i + 1] == ' ') {
          stringstream sub(work.substr(i + 2, work.size() - (i + 2)));
          Veci3 wbuf(3);
          int  widx(0);
          bool flag(false);
          while(!sub.eof() && !sub.bad()) {
            sub >> wbuf[widx];
            if(wbuf[widx] >= 0)
              wbuf[widx] --;
            widx ++;
            if(widx > 2)
              flag = true;
            if(flag)
              polys.emplace_back(wbuf);
            widx %= 3;
            if(sub.eof() || sub.bad())
              break;
            sub.ignore(20, ' ');
          }
        }
      }
      for(int i = 0; i < polys.size(); i ++)
        for(int j = 0; j < polys[i].size(); j ++) {
          while(polys[i][j] < 0) polys[i][j] += data.size();
          polys[i][j] = polys[i][j] % data.size();
        }
      input.close();
    } else {
      cerr << "Unable to open file for read: " << filename << endl;
      return false;
    }
    return true;
  }
  
  bool saveMTL(const char* photo, const char* filename) {
    ofstream output;
    output.open(filename, std::ios::out);
    if(output.is_open()) {
      std::string pstr(photo);
      for(int i = 0; i < pstr.size(); i ++)
        if(pstr[pstr.size() - i - 1] == '.') {
          pstr = pstr.substr(0, max(0, int(pstr.size()) - i - 1)) + std::string(".ppm");
          break;
        }
      output << "newmtl material0" << std::endl;
      output << "Ka 1.000000 1.000000 1.000000" << std::endl;
      output << "Kd 1.000000 1.000000 1.000000" << std::endl;
      output << "Ks 0.000000 0.000000 0.000000" << std::endl;
      output << "illum 1" << std::endl;
      output << "map_Ka " << pstr << std::endl;
      output << "map_Kd " << pstr << std::endl << std::endl;
      output.close();
    } else {
      cerr << "Unable to open file for write: " << filename << endl;
      return false;
    }
    return true;
  }

  bool loaddat(const char* filename, std::string& header, vector<vector<T> >& data) {
    ifstream input;
    input.open(filename);
    if(input.is_open()) {
      string work;
      header = std::string("");
      data   = vector<vector<T> >();
      while(getline(input, work) && !input.eof() && !input.bad())
        if(whiteline(work))
          continue;
        else if(work[0] == ';')
          header += work + std::string("\n");
        else {
          std::stringstream ss(work);
          for(int i = 0, j = 0; ss.tellg() <= work.size(); j ++) {
            if(data.size() <= j)
              data.resize(j + 1, vector<T>());
            data[j].emplace_back(T(0));
            ss >> data[j][data[j].size() - 1];
          }
        }
      input.close();
    } else {
      cerr << "Unable to open file for read: " << filename << endl;
      return false;
    }
    return true;
  }
  
  bool savedat(const char* filename, std::string& header, vector<vector<T> >& data) {
    ofstream output;
    output.open(filename, std::ios::out);
    if(output.is_open()) {
      output << header;
      for(int i = 0; i < data[0].size(); i ++) {
        for(int j = 0; j < data.size(); j ++)
          output << (i < data[j].size() ? data[j][i] : T(0)) << " ";
        output << std::endl;
      }
      output.close();
    } else {
      cerr << "Unable to open file for write: " << filename << endl;
      return false;
    }
    return true;
  }
  
  bool loadcenterr(vector<Vec3>& center, vector<T>& r, const char* filename) {
    center = vector<Vec3>();
    r      = vector<T>();
    std::ifstream input;
    try {
      input.open(filename);
      std::string buf;
      while(getline(input, buf) && !input.eof() && !input.bad()) {
        std::stringstream sbuf(buf);
        typename simpleFile<num_t>::Vec3 work(3);
        sbuf >> work[0];
        sbuf >> work[1];
        sbuf >> work[2];
        center.emplace_back(work);
        num_t workr;
        sbuf >> workr;
        r.emplace_back(workr);
      }
      input.close();
    } catch(...) {
      std::cerr << "Something had occured when reading center - r txt." << std::endl;
      return false;
    }
    return center.size() == r.size();
  }
  
  bool savecenterr(const char* filename, const vector<Vec3>& center, const vector<T>& r) {
    ofstream output;
    output.open(filename, std::ios::out);
    if(output.is_open()) {
      assert(center.size() == r.size());
      for(int i = 0; i < center.size(); i ++) {
        assert(center[i].size() == 3);
        output << center[i][0] << " " << center[i][1] << " " << center[i][2] << " " << r[i] << std::endl;
      }
      output.close();
    } else {
      cerr << "Unable to open file for write: " << filename << endl;
      return false;
    }
    return true;
  }
};

#define _FILEIO_GOKI_
#endif

