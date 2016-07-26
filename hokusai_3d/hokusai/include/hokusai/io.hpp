#ifndef HOKUSAI_IO_H
#define HOKUSAI_IO_H

#include <stdio.h>
#include <vector>
#include "strutil.hpp"
#include "common.hpp"
#include "system.hpp"

namespace hokusai
{

int countDigits(int number);
#define GET_LINE() do { if (!fgets(buf, 1024, f)) return false; } while (0)
#define COND_READ(cond, where, len) do { if ((cond) && !fread((void *)&(where), (len), 1, f)) return false; } while (0)
#define FPRINTF(...) do { if (fprintf(__VA_ARGS__) < 0) return false; } while (0)
#define FWRITE(ptr, size, nmemb, stream) do { if (fwrite((ptr), (size), (nmemb), (stream)) != (nmemb)) return false; } while (0)
#define LINE_IS(text) hokusai::begins_with(buf, text)

// Skip comments in an ASCII file (lines beginning with #)
void skip_comments(FILE *f);

// Tesselate an arbitrary n-gon.  Appends triangles to "tris".
void tess(const std::vector<Vec3r > &verts, const std::vector<int> &thisface, std::vector<Vec3i > &tris);

// Read an obj file
bool read_obj(FILE *f,  std::vector< Vec3r >& vertices,
              std::vector< Vec3r >& normals,
              std::vector< Vec3i >& triangles);


class HokusaiImporter
{
public:
    ~HokusaiImporter();
    HokusaiImporter(const std::string& fileName, System& system);
};

class HokusaiExporter
{
public:
    ~HokusaiExporter();
    HokusaiExporter(const std::string& fileName, const System& system);
};

class BlenderExporter
{
    public:
        ~BlenderExporter();
        BlenderExporter();
        BlenderExporter(const std::string& prefix, const std::string& path, const int& particleNumber, const int& frameNumber);
        void debugHeader(const std::string& fileName);
        void debugFrame(const std::string& fileName);
        void apply(const System& system);
        std::string& prefix();
        const std::string& prefix() const;
        std::string& path();
        const std::string& path() const;
        int& frameNumber();
        const int& frameNumber() const;
        int& particleNumber();
        const int& particleNumber() const;
        int& frameCounter();
        const int& frameCounter() const;
    private:
        std::string m_prefix;
        std::string m_path;
        int m_frameNumber;
        int m_frameCounter;
        int m_particleNumber;
};

} //namespace hokusai
#endif
