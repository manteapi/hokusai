#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <vector>
#include "Vec.hpp"
#include "strutil.hpp"

#define GET_LINE() do { if (!fgets(buf, 1024, f)) return false; } while (0)
#define COND_READ(cond, where, len) do { if ((cond) && !fread((void *)&(where), (len), 1, f)) return false; } while (0)
#define FPRINTF(...) do { if (fprintf(__VA_ARGS__) < 0) return false; } while (0)
#define FWRITE(ptr, size, nmemb, stream) do { if (fwrite((ptr), (size), (nmemb), (stream)) != (nmemb)) return false; } while (0)
#define LINE_IS(text) trimesh::begins_with(buf, text)

// Skip comments in an ASCII file (lines beginning with #)
void skip_comments(FILE *f);

// Tesselate an arbitrary n-gon.  Appends triangles to "tris".
void tess(const std::vector<Vec3<float> > &verts, const std::vector<int> &thisface, std::vector<Vec3<int> > &tris);

// Read an obj file
bool read_obj(FILE *f,  std::vector< Vec3<float> >& vertices, 
                        std::vector< Vec3<float> >& normals, 
                        std::vector< Vec3<int> >& triangles);

#endif
