#ifndef HOKUSAI_HSL2RGB_HPP
#define HOKUSAI_HSL2RGB_HPP

#include "common.hpp"

namespace hokusai
{
inline unsigned char GetRValue(unsigned int color)
{
	return (unsigned char)((color>>16)&0xFF);
}

inline unsigned char GetGValue(unsigned int color)
{
	return (unsigned char)((color>>8)&0xFF);
}

inline unsigned char GetBValue(unsigned int color)
{
	return (unsigned char)(color&0xFF);
}

inline unsigned int RGB(unsigned char r,unsigned char g,unsigned char b)
{
	unsigned int color = ((unsigned int)r<<16) | ((unsigned int)g<<8) | b;
	return color;
}

void RGBtoHSL(unsigned int color,unsigned int& h, unsigned int& s, unsigned int& l);
unsigned int HSLtoRGB(const unsigned int& h, const unsigned int& s, const unsigned int& l);
unsigned int BrightenColor(unsigned int color,unsigned int amount);
unsigned int DarkenColor(unsigned int color,unsigned int amount);
}
#endif
