#ifndef HOKUSAI_IO_H
#define HOKUSAI_IO_H

#include <vector>
#include <string>
#include "marchingSquare.hpp"

namespace hokusai
{

void exportOBJ(const std::string & filename, const std::vector<Edge> & edges);

void SVGExport(FILE* file, const Vec2d& center, const double& radius, const double& linewidth, std::string color = "#000000")
{
    fprintf(file, "<circle cx=\"%f\" cy=\"%f\" r=\"%f\" stroke=\"%s\" stroke-width=\"%f\" fill=\"%s\" />\n", center[0], center[1], radius, color.c_str(), linewidth, color.c_str()  );
}


void SVGExport(FILE* file, const Vec2d& p0, const Vec2d& p1, std::string color = "#000000")
{
    fprintf(file, "<g>\n");
    fprintf(file, "<path\n");
    fprintf(file, "style=\"fill:none;stroke:%s;stroke-width:%f;\"\n", color.c_str(), 2.0);
    fprintf(file, "d=\"M %f,%f %f,%f z\"\n", p0[0], p0[1], p1[0], p1[1]);
    fprintf(file, "/>\n");
    fprintf(file, "</g>\n");
}


void SVGExport(FILE* file, const std::vector<hokusai::Edge>& edges, int width, int height)
{
    fprintf(file, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n");

    fprintf(file, "<svg\n");
    fprintf(file, "width=\"%d\"\n",width);
    fprintf(file, "height=\"%d\"\n",height);
    fprintf(file, "version=\"1.1\">\n");

    for(const hokusai::Edge& e : edges)
    {
        SVGExport(file, e.p1, e.p2);
    }

    fprintf(file, "</svg>\n");
}

void SVGExport(FILE* file, const std::vector<Vec2d>& points, std::vector<std::string>& colors, std::vector<double>& radius, std::vector<double>& linewidth, int width, int height)
{
    fprintf(file, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n");

    fprintf(file, "<svg\n");
    fprintf(file, "width=\"%d\"\n",width);
    fprintf(file, "height=\"%d\"\n",height);
    fprintf(file, "version=\"1.1\">\n");

    for(size_t i=0; i<points.size(); ++i)
    {
        SVGExport(file, points[i], radius[i], linewidth[i], colors[i]);
    }

    fprintf(file, "</svg>\n");
}

void SVGExport(const std::string& filename, std::vector<Vec2d>& points, std::vector<std::string>& colors, std::vector<double>& radius, std::vector<double>& linewidth, int width, int height)
{
    Vec2d minBB, maxBB;
    for(size_t i=0; i<2; ++i)
    {
        minBB[i] = std::numeric_limits<double>::max();
        maxBB[i] = -std::numeric_limits<double>::max();
    }
    for(Vec2d& p : points)
    {
        for(size_t i=0; i<2; ++i)
        {
            minBB[i] = std::min(minBB[i], p[i]);
            maxBB[i] = std::max(maxBB[i], p[i]);
        }
    }

    double marginW = 0.01*width;
    double marginH = 0.01*height;
    std::array<double,3> scaleFactor = {1.0,1.0,1.0};
    scaleFactor[0] = (width-marginW)*(1.0/(maxBB[0]-minBB[0]));
    scaleFactor[1] = (height-marginH);
    for(Vec2d& p : points)
    {
        Vec2d tmp_point( (p[0]-minBB[0])*scaleFactor[0] + 0.5*marginW, -( (p[1]-minBB[1])*scaleFactor[1]  + 0.5*marginH) + height);
        p = tmp_point;
    }

    FILE *file;

    if((file=fopen(filename.c_str(),"w"))==NULL)
    {
        std::cout << "Unable to open : " << filename << std::endl;
        return;
    }

    std::cout << "Exporting '" << filename << "'..." << std::flush;
    SVGExport(file, points, colors, radius, linewidth, width, height);
    std::cout << " done." << std::endl;

    fclose(file);
}

void SVGExport(const std::string& filename, const Vec2d& minBB, const Vec2d& maxBB, const std::vector<hokusai::Edge>& edges, int width, int height)
{
    std::vector<hokusai::Edge> f_edge = edges;

    //Keep it
//    Vec2d minBB, maxBB;
//    for(size_t i=0; i<2; ++i)
//    {
//        minBB[i] = std::numeric_limits<double>::max();
//        maxBB[i] = -std::numeric_limits<double>::max();
//    }
//    for(hokusai::Edge& e : f_edge)
//    {
//        for(size_t i=0; i<2; ++i)
//        {
//            minBB[i] = std::min(minBB[i], e.p1[i]);
//            minBB[i] = std::min(minBB[i], e.p2[i]);

//            maxBB[i] = std::max(maxBB[i], e.p1[i]);
//            maxBB[i] = std::max(maxBB[i], e.p2[i]);
//        }
//    }

    double marginW = 0.01*width;
    double marginH = 0.01*height;
    std::array<double,3> scaleFactor = {1.0,1.0,1.0};
    scaleFactor[0] = (width-marginW)*(1.0/(maxBB[0]-minBB[0]));
    scaleFactor[1] = (height-marginH);
    for(hokusai::Edge& e : f_edge)
    {
        Vec2d p1( (e.p1[0]-minBB[0])*scaleFactor[0] + 0.5*marginW, -( (e.p1[1]-minBB[1])*scaleFactor[1]  + 0.5*marginH) + height);
        Vec2d p2( (e.p2[0]-minBB[0])*scaleFactor[0] + 0.5*marginW, -( (e.p2[1]-minBB[1])*scaleFactor[1]  + 0.5*marginH) + height);
        e = hokusai::Edge(p1,p2);
    }

    FILE *file;

    if((file=fopen(filename.c_str(),"w"))==NULL)
    {
        std::cout << "Unable to open : " << filename << std::endl;
        return;
    }

    std::cout << "Exporting '" << filename << "'..." << std::flush;
    SVGExport(file, f_edge, width, height);
    std::cout << " done." << std::endl;

    fclose(file);
}

}

#endif
