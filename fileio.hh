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

template <typename T> class simpleFile {
public:
#if defined(_WITHOUT_EIGEN_)
  typedef SimpleMatrix<T> Mat;
  typedef SimpleMatrix<T> Mat3x3;
  typedef SimpleVector<T> Vec3;
  typedef SimpleVector<int> Veci3;
#else
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T,   3, 3> Mat3x3;
  typedef Eigen::Matrix<T,   3, 1> Vec3;
  typedef Eigen::Matrix<int, 3, 1> Veci3;
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

  bool saveobj(const vector<Vec3>& data, const vector<Veci3>& polys, const char* filename, const bool& arout = false, const bool& addstand = false, const vector<vector<int> >& edges = vector<vector<int> >(), const T& zs = T(2), const T& aroffset = T(.2)) {
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
        for(int i = 0; i < data.size(); i ++)
          output << "v " << aroffset + (data[i][1] - Mh / T(2)) / sqrt(Mh * Mh + Mw * Mw) << " " << (Mw / T(2) - data[i][0]) / sqrt(Mh * Mh + Mw * Mw) << " " << - (data[i][2] + lz / T(2)) / sqrt(Mh * Mh + Mw * Mw) - T(.5) << endl;
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
  
  // following 3 functions and 1 structure:
  //   thanks to: https://github.com/jessey-git/fx-gltf/
#if defined(_WITH_GLTF2_)
  struct BufferInfo {
    fx::gltf::Accessor const * Accessor;
    uint8_t const * Data;
    uint32_t DataStride;
    uint32_t TotalSize;
    bool HasData() const noexcept {
      return Data != nullptr;
    }
  };

  static uint32_t CalculateDataTypeSize(fx::gltf::Accessor const & accessor) noexcept {
        uint32_t elementSize = 0;
        switch (accessor.componentType)
        {
        case fx::gltf::Accessor::ComponentType::Byte:
        case fx::gltf::Accessor::ComponentType::UnsignedByte:
            elementSize = 1;
            break;
        case fx::gltf::Accessor::ComponentType::Short:
        case fx::gltf::Accessor::ComponentType::UnsignedShort:
            elementSize = 2;
            break;
        case fx::gltf::Accessor::ComponentType::Float:
        case fx::gltf::Accessor::ComponentType::UnsignedInt:
            elementSize = 4;
            break;
        }

        switch (accessor.type)
        {
        case fx::gltf::Accessor::Type::Mat2:
            return 4 * elementSize;
            break;
        case fx::gltf::Accessor::Type::Mat3:
            return 9 * elementSize;
            break;
        case fx::gltf::Accessor::Type::Mat4:
            return 16 * elementSize;
            break;
        case fx::gltf::Accessor::Type::Scalar:
            return elementSize;
            break;
        case fx::gltf::Accessor::Type::Vec2:
            return 2 * elementSize;
            break;
        case fx::gltf::Accessor::Type::Vec3:
            return 3 * elementSize;
            break;
        case fx::gltf::Accessor::Type::Vec4:
            return 4 * elementSize;
            break;
        }

        return 0;
  }
  
  struct BufferInfo GetData(fx::gltf::Document const & doc, fx::gltf::Accessor const & accessor) {
    fx::gltf::BufferView const & bufferView = doc.bufferViews[accessor.bufferView];
    fx::gltf::Buffer const & buffer = doc.buffers[bufferView.buffer];
    const uint32_t dataTypeSize = CalculateDataTypeSize(accessor);
    return BufferInfo{ &accessor, &buffer.data[static_cast<uint64_t>(bufferView.byteOffset) + accessor.byteOffset], dataTypeSize, accessor.count * dataTypeSize };
  }
#endif
  
  bool loadglTF(vector<vector<Vec3> >& data, vector<vector<Veci3> >& polys, vector<Vec3>& center, const char* filename) {
#if defined(_WITH_GLTF2_)
    data   = vector<vector<Vec3> >();
    polys  = vector<vector<Veci3> >();
    center = vector<Vec3>();
    fx::gltf::Document doc = fx::gltf::LoadFromText(filename);
    for(int nodeIndex = 0; nodeIndex < doc.nodes.size(); nodeIndex ++) {
      Vec3 lcenter(3);
      for(int i = 0; i < 3; i ++)
        lcenter[i] = doc.nodes[nodeIndex].translation[i];
      center.push_back(lcenter);
    }
    data.resize(center.size());
    polys.resize(center.size());
    for(int meshIndex = 0; meshIndex < doc.meshes.size(); meshIndex ++)
      for(int primitiveIndex = 0; primitiveIndex < doc.meshes[meshIndex].primitives.size(); primitiveIndex ++) {
        fx::gltf::Mesh const & mesh = doc.meshes[meshIndex];
        fx::gltf::Primitive const & primitive = mesh.primitives[primitiveIndex];
        BufferInfo m_vertexBuffer;
        BufferInfo m_jointBuffer;
        fx::gltf::Accessor::ComponentType jointType;
        for (auto const & attrib : primitive.attributes)
          if (attrib.first == "POSITION") {
            m_vertexBuffer = GetData(doc, doc.accessors[attrib.second]);
            assert(doc.accessors[attrib.second].componentType == fx::gltf::Accessor::ComponentType::Float);
            assert(doc.accessors[attrib.second].type == fx::gltf::Accessor::Type::Vec3);
          } else if(attrib.first == "JOINTS_0") {
            m_jointBuffer = GetData(doc, doc.accessors[attrib.second]);
            jointType     = doc.accessors[attrib.second].componentType;
            assert(doc.accessors[attrib.second].componentType == fx::gltf::Accessor::ComponentType::UnsignedByte || doc.accessors[attrib.second].componentType == fx::gltf::Accessor::ComponentType::UnsignedShort);
            assert(doc.accessors[attrib.second].type == fx::gltf::Accessor::Type::Vec4);
          }
        for(int cIdx = 0; cIdx < center.size(); cIdx ++) {
          vector<Vec3>  vertices;
          vector<Veci3> indices;
          for(int i = 0; i < m_vertexBuffer.DataStride; i ++) {
            bool idx_match(false);
            for(int j = 0; j < 4; j ++)
              if(jointType == fx::gltf::Accessor::ComponentType::UnsignedByte) {
                unsigned char c = *(unsigned char*)(&m_jointBuffer.Data[4 * i + j]);
                idx_match |= static_cast<int>(c) == cIdx;
              } else if(jointType == fx::gltf::Accessor::ComponentType::UnsignedShort) {
                unsigned short c = *(unsigned short*)(&m_jointBuffer.Data[4 * 2 * i + 2 * j]);
                idx_match |= static_cast<int>(c) == cIdx;
              }
            if(idx_match) {
              Vec3 work(3);
              for(int j = 0; j < 3; j ++)
                work[j] = *(float*)(&m_vertexBuffer.Data[4 * 3 * i + j]);
              vertices.push_back(work);
            }
          }
          for(int i = 2; i < vertices.size(); i ++) {
            Veci3 work(3);
            work[0] = i - 2;
            work[1] = i - 1;
            work[2] = i;
            indices.push_back(work);
          }
          data[cIdx]  = vertices;
          polys[cIdx] = indices;
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

