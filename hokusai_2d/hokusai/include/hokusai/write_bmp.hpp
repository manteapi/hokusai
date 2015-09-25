#ifndef HOKUSAI_WRITE_BMP_HPP
#define HOKUSAI_WRITE_BMP_HPP

namespace hokusai
{
void write_bmp( const char *filename, unsigned char *image, int width, int height, bool invertY );
}

#endif
