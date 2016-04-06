#include "./../include/hokusai/io.hpp"

using namespace std;
namespace hokusai
{

// Skip comments in an ASCII file (lines beginning with #)
void skip_comments(FILE *f)
{
    int c;
    bool in_comment = false;
    while (1) {
        c = fgetc(f);
        if (c == EOF)
            return;
        if (in_comment) {
            if (c == '\n')
                in_comment = false;
        } else if (c == '#') {
            in_comment = true;
        } else if (!isspace(c)) {
            break;
        }
    }
    ungetc(c, f);
}

// Tesselate an arbitrary n-gon.  Appends triangles to "tris".
void tess(const vector<Vec3r > &verts, const vector<int> &thisface,
          vector<Vec3i > &tris)
{
    if (thisface.size() < 3)
        return;
    if (thisface.size() == 3) {
        tris.push_back(Vec3i(thisface[0],
                       thisface[1],
                thisface[2]));
        return;
    }
    if (thisface.size() == 4) {
        // Triangulate in the direction that
        // gives the shorter diagonal
        const Vec3r &p0 = verts[thisface[0]], &p1 = verts[thisface[1]];
        const Vec3r &p2 = verts[thisface[2]], &p3 = verts[thisface[3]];
        HReal d02 = (p0-p2).lengthSquared();
        HReal d13 = (p1-p3).lengthSquared();
        int i = (d02 < d13) ? 0 : 1;
        tris.push_back(Vec3i(thisface[i],
                                  thisface[(i+1)%4],
                       thisface[(i+2)%4]));
        tris.push_back(Vec3i(thisface[i],
                                  thisface[(i+2)%4],
                       thisface[(i+3)%4]));
        return;
    }

    // 5-gon or higher - just tesselate arbitrarily...
    for (size_t i = 2; i < thisface.size(); i++)
        tris.push_back(Vec3i(thisface[0],
                       thisface[i-1],
                thisface[i]));
}

// Read an obj file
bool read_obj(FILE *f, vector< Vec3r >& vertices, vector< Vec3r >& normals, vector< Vec3i >& triangles)
{


    vector<int> thisface;
    while (1) {
        skip_comments(f);
        if (feof(f))
            return true;
        char buf[1024];
        GET_LINE();
        if (LINE_IS("v ") || LINE_IS("v\t")) {
            HReal x, y, z;
            if (sscanf(buf+1, "%f %f %f", &x, &y, &z) != 3) {
                return false;
            }
            vertices.push_back(Vec3r(x,y,z));
        } else if (LINE_IS("vn ") || LINE_IS("vn\t")) {
            HReal x, y, z;
            if (sscanf(buf+2, "%f %f %f", &x, &y, &z) != 3) {
                return false;
            }
            normals.push_back(Vec3r(x,y,z));
        } else if (LINE_IS("f ") || LINE_IS("f\t") ||
                   LINE_IS("t ") || LINE_IS("t\t")) {
            thisface.clear();
            char *c = buf;
            while (1) {
                while (*c && *c != '\n' && !isspace(*c))
                    c++;
                while (*c && isspace(*c))
                    c++;
                int thisf;
                if (sscanf(c, " %d", &thisf) != 1)
                    break;
                if (thisf < 0)
                    thisf += vertices.size();
                else
                    thisf--;
                thisface.push_back(thisf);
            }
            tess(vertices, thisface, triangles);
        }
    }

    // XXX - FIXME
    // Right now, handling of normals is fragile: we assume that
    // if we have the same number of normals as vertices,
    // the file just uses per-vertex normals.  Otherwise, we can't
    // handle it.
    if (vertices.size() != normals.size())
        normals.clear();

    return true;
}
}
