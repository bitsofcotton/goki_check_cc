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

#if defined(_WITHOUT_EIGEN_)
#include "simplelin.hh"
#else
#include <Eigen/Core>
#endif

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#if defined(_WITH_GLTF2_)
#include <fx/gltf.h>
#endif

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
using std::sqrt;
using std::vector;
using std::pair;
using std::make_pair;
using std::atan2;

template <typename T> class match_t;

template <typename T> class simpleFile {
public:
#if defined(_WITHOUT_EIGEN_)
  typedef SimpleMatrix<T> Mat;
  typedef SimpleMatrix<T> Mat3x3;
  typedef SimpleVector<T> Vec3;
  typedef SimpleVector<int> Veci3;
  typedef SimpleVector<int> Veci4;
#else
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T,   3, 3> Mat3x3;
  typedef Eigen::Matrix<T,   3, 1> Vec3;
  typedef Eigen::Matrix<int, 3, 1> Veci3;
  typedef Eigen::Matrix<int, 4, 1> Veci4;
#endif

  bool loadstub(ifstream& input, const int& nmax, const int& ncolor, Mat* datas) {
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
        datas[k](j, i) = T(work) / nmax;
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
        while(line[0] == '#' && getline(input, line) && !input.eof() && !input.bad()) ;
        getline(input, line2);
        while(line2[0] == '#' && getline(input, line2) && !input.eof() && !input.bad()) ;
        getline(input, line3);
        while(line3[0] == '#' && getline(input, line3) && !input.eof() && !input.bad()) ;
        istringstream iline2(line2);
        int w, h;
        iline2 >> w;
        iline2 >> h;
        if(w <= 0 || h <= 0) {
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
            loadstub(input, nmax, 1, data);
          } else if(line[1] == '3') {
            for(int i = 0; i < 3; i ++)
              data[i] = Mat(h, w);
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
      } catch (const ifstream::failure& e) {
        cerr << "Exception while reading." << endl;
      } catch (const istringstream::failure& e) {
        cerr << "Exception while reading." << endl;
      }
      input.close();
    } else {
      cerr << "Unable to open file: " << filename << endl;
      return false;
    }
    return true;
  }

  bool savep2or3(const char* filename, Mat data[3], const bool& gray) {
    ofstream output;
    output.open(filename);
    if(output.is_open()) {
      try {
        if(gray)
          output << "P2" << "\n";
        else
          output << "P3" << "\n";
        output << data[0].cols() << " " << data[0].rows() << "\n";
        output << 255 << "\n";
        for(int i = 0; i < data[0].rows(); i ++)
          for(int j = 0; j < data[0].cols(); j ++)
            if(gray)
              output << int(data[0](i, j) * 255) << "\n";
            else
              for(int k = 0; k < 3; k ++)
                output << int(data[k](i, j) * 255) << "\n";
      } catch (const ofstream::failure& e) {
        cerr << "An error has occured while writing file." << endl;
      }
      output.close();
    } else {
      cerr << "Unable to open file: " << filename << endl;
      return false;
    }
    return true;
  }

  bool saveobj(const vector<Vec3>& data, const vector<Veci3>& polys, const char* filename, const bool& arout = false, const bool& addstand = false, const vector<vector<int> >& edges = vector<vector<int> >(), const T& zs = T(2), const T& aroffset = T(.2), const T& arrot = T(.025)) {
    ofstream output;
    output.open(filename, std::ios::out);
    if(output.is_open()) {
      double Mh(1), Mw(1), lz(0);
      int lfs(0);
      for(int fslash(0) ; filename[fslash]; fslash ++)
        if(filename[fslash] == '/')
          lfs = fslash;
      if(lfs) lfs ++;
      output << "mtllib " << &filename[lfs] << ".mtl" << endl;
      for(int i = 0; i < data.size(); i ++) {
        Mh = max(data[i][1], Mh);
        Mw = max(data[i][0], Mw);
        lz = min(- data[i][2], lz);
      }
      if(arout) {
        assert(!addstand);
        match_t<T> m;
        m.offset[0] += aroffset;
        m.offset[2] -= T(.5);
        const T Pi(T(4) * atan2(T(1), T(1)));
        const T theta(2. * Pi * T(aroffset < T(0) ? 0 : 1) / T(2));
        const T lpsi(Pi * arrot);
        Mat3x3 R0(3, 3);
        Mat3x3 R1(3, 3);
        R0(0, 0) =   cos(theta);
        R0(0, 1) = - sin(theta);
        R0(0, 2) = 0.;
        R0(1, 0) =   sin(theta);
        R0(1, 1) =   cos(theta);
        R0(1, 2) = 0.;
        R0(2, 0) = 0.;
        R0(2, 1) = 0.;
        R0(2, 2) = 1.;
        R1(0, 0) = 1.;
        R1(0, 1) = 0.;
        R1(0, 2) = 0.;
        R1(1, 0) = 0.;
        R1(1, 1) =   cos(lpsi);
        R1(1, 2) = - sin(lpsi);
        R1(2, 0) = 0.;
        R1(2, 1) =   sin(lpsi);
        R1(2, 2) =   cos(lpsi);
        m.rot    = R0.transpose() * R1 * R0;
        for(int i = 0; i < data.size(); i ++) {
          Vec3 workv(3);
          workv[0] =   (  Mw / T(2) - data[i][0]) / sqrt(Mh * Mh + Mw * Mw);
          workv[1] =   (- Mh / T(2) + data[i][1]) / sqrt(Mh * Mh + Mw * Mw);
          workv[2] = - (  lz / T(2) + data[i][2]) / sqrt(Mh * Mh + Mw * Mw);
          workv    = m.transform(workv);
          output << "v " << workv[1] << " " << workv[0] << " " << workv[2] << endl;
        }
      } else if(addstand) {
        for(int i = 0; i < data.size(); i ++)
          output << "v " << data[i][1] << " " << Mw - data[i][0] << " " << - data[i][2] + zs - lz << endl;
        for(int i = 0; i < data.size(); i ++)
          output << "v " << data[i][1] << " " << Mw - data[i][0] << " " << 0 << endl;
      } else
        for(int i = 0; i < data.size(); i ++)
          output << "v " << data[i][1] << " " << - data[i][0] << " " << - data[i][2] << endl;
      for(int i = 0; i < data.size(); i ++)
        output << "vt " << data[i][1] / Mh << " " << 1. - data[i][0] / Mw << endl;
      output << "usemtl material0" << endl;
      // xchg with clockwise/counter clockwise.
      for(int i = 0; i < polys.size(); i ++) {
        output << "f " << polys[i][0] + 1 << "/" << polys[i][0] + 1 << "/" << polys[i][0] + 1;
        output << " "  << polys[i][1] + 1 << "/" << polys[i][1] + 1 << "/" << polys[i][1] + 1;
        output << " "  << polys[i][2] + 1 << "/" << polys[i][2] + 1 << "/" << polys[i][2] + 1 << endl;
        if(addstand) {
          output << "f " << data.size() + polys[i][0] + 1;
          output << " "  << data.size() + polys[i][2] + 1;
          output << " "  << data.size() + polys[i][1] + 1 << endl;
        }
      }
      if(addstand && edges.size())
        cerr << edges.size() << "parts found." << endl;
        for(int ii = 0; ii < edges.size(); ii ++) if(edges[ii].size()) {
          const vector<int>& outer(edges[ii]);
          for(int i = 0; i < outer.size() - 1; i ++) {
            const int& i0(i);
            const int  i1(i + 1);
            output << "f " << data.size() + outer[i0] + 1;
            output << " "  << outer[i1] + 1;
            output << " "  << data.size() + outer[i1] + 1 << endl;
            output << "f " << data.size() + outer[i0] + 1;
            output << " "  << outer[i0] + 1;
            output << " "  << outer[i1] + 1 << endl;
          }
        }
      output.close();
    } else {
      cerr << "Unable to open file: " << filename << endl;
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
        for(; i < work.size() && work[i] == ' '; i ++);
        if(i + 1 < work.size() && work[i] == 'v' && work[i + 1] == ' ') {
          stringstream sub(work.substr(i + 2, work.size() - (i + 2)));
          Vec3 buf(3);
          sub >> buf[1];
          sub >> buf[0];
          sub >> buf[2];
          buf[0] = - buf[0];
          buf[2] = - buf[2];
          data.push_back(buf);
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
              polys.push_back(wbuf);
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
      cerr << "Unable to open file: " << filename << endl;
      return false;
    }
    return true;
  }
  
  // ( Please read the page before to compile )
  // thanks to: https://github.com/jessey-git/fx-gltf/
  bool loadglTF(vector<vector<Vec3> >& data, vector<vector<Veci3> >& polys, vector<Vec3>& center, vector<vector<Veci4> >& bone, const char* filename) {
#if defined(_WITH_GLTF2_)
    data   = vector<vector<Vec3> >();
    polys  = vector<vector<Veci3> >();
    center = vector<Vec3>();
    bone   = vector<vector<Veci4> >();
    const auto doc(fx::gltf::LoadFromText(filename));
    for(int nodeIndex = 0; nodeIndex < doc.nodes.size(); nodeIndex ++) {
      Vec3 lcenter(3);
      for(int i = 0; i < 3; i ++)
        lcenter[i] = doc.nodes[nodeIndex].translation[i];
      center.push_back(lcenter);
    }
    data.resize(center.size());
    polys.resize(center.size());
    bone.resize(center.size());
    for(int meshIndex = 0; meshIndex < doc.meshes.size(); meshIndex ++)
      for(int primitiveIndex = 0; primitiveIndex < doc.meshes[meshIndex].primitives.size(); primitiveIndex ++) {
        const auto& mesh(doc.meshes[meshIndex]);
        const auto& primitive(mesh.primitives[primitiveIndex]);
        uint8_t* m_vertexBuffer;
        int      n_vertexBuffer;
        uint8_t* m_jointBuffer;
        int      n_jointBuffer;
        bool found_pos(false);
        bool found_joint(false);
        fx::gltf::Accessor::ComponentType jointType;
        for (auto const & attrib : primitive.attributes)
          if (attrib.first == std::string("POSITION")) {
            const auto& bv(doc.bufferViews[doc.accessors[attrib.second].bufferView]);
            m_vertexBuffer = (uint8_t*)&doc.buffers[bv.buffer].data[static_cast<uint64_t>(bv.byteOffset) + doc.accessors[attrib.second].byteOffset];
            n_vertexBuffer = doc.accessors[attrib.second].count;
            assert(doc.accessors[attrib.second].componentType == fx::gltf::Accessor::ComponentType::Float);
            assert(doc.accessors[attrib.second].type == fx::gltf::Accessor::Type::Vec3);
            found_pos = true;
          } else if(attrib.first == std::string("JOINTS_0")) {
            const auto& bv(doc.bufferViews[doc.accessors[attrib.second].bufferView]);
            m_jointBuffer = (uint8_t*)&doc.buffers[bv.buffer].data[static_cast<uint64_t>(bv.byteOffset) + doc.accessors[attrib.second].byteOffset];
            n_jointBuffer = doc.accessors[attrib.second].count;
            jointType     = doc.accessors[attrib.second].componentType;
            assert(jointType == fx::gltf::Accessor::ComponentType::UnsignedByte ||
                   jointType == fx::gltf::Accessor::ComponentType::UnsignedShort);
            assert(doc.accessors[attrib.second].type == fx::gltf::Accessor::Type::Vec4);
            found_joint = true;
          }
        const auto bv(doc.bufferViews[doc.accessors[primitive.indices].bufferView]);
        const uint8_t* m_indexBuffer = (uint8_t*)&doc.buffers[bv.buffer].data[static_cast<uint64_t>(bv.byteOffset) + doc.accessors[primitive.indices].byteOffset];
        assert(doc.accessors[primitive.indices].componentType == fx::gltf::Accessor::ComponentType::UnsignedShort);
        assert(doc.accessors[primitive.indices].type == fx::gltf::Accessor::Type::Scalar);
        assert(found_pos && found_joint);
        for(int cIdx = 0; cIdx < center.size(); cIdx ++) {
          vector<Vec3>  vertices;
          vector<Veci3> indices;
          vector<Veci4> bones;
          Veci4         lbone(4);
          lbone[0] = lbone[1] = lbone[2] = lbone[3] = - 1;
          for(int i = 0; i < n_vertexBuffer; i ++) {
            bool idx_match(false);
            for(int j = 0; j < 4; j ++)
              if(jointType == fx::gltf::Accessor::ComponentType::UnsignedByte) {
                unsigned char c = *(unsigned char*)(&m_jointBuffer[4 * i + j]);
                lbone[j]   = static_cast<int>(c);
                assert(0 <= lbone[j] && lbone[j] < center.size());
                idx_match |= static_cast<int>(c) == cIdx;
              } else if(jointType == fx::gltf::Accessor::ComponentType::UnsignedShort) {
                unsigned short c = *(unsigned short*)(&m_jointBuffer[4 * 2 * i + 2 * j]);
                lbone[j]   = static_cast<int>(c);
                assert(0 <= lbone[j] && lbone[j] < center.size());
                idx_match |= static_cast<int>(c) == cIdx;
              }
            if(idx_match) {
              Vec3 work(3);
              for(int j = 0; j < 3; j ++)
                work[j] = T(*(float*)(&m_vertexBuffer[3 * 4 * i + 4 * j]));
              vertices.push_back(work);
              bones.push_back(lbone);
            }
          }
          for(int i = 2; i < vertices.size(); i ++) {
            // XXX fixme: delaunay?
            Veci3 work(3);
            work[0] = i - 2;
            work[1] = i - 1;
            work[2] = i;
            indices.push_back(work);
          }
          data[cIdx]  = vertices;
          polys[cIdx] = indices;
          bone[cIdx]  = bones;
        }
      }
    cerr << data.size() << "parts found." << endl;
    return true;
#else
    assert(0 && "Please compile with _WITH_GLTF2_");
    return false;
#endif
  }
};

#define _FILEIO_GOKI_
#endif

