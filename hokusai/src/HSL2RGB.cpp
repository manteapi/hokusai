#include "../include/hokusai/HSL2RGB.hpp"
//

namespace hokusai
{

// This is a subfunction of HSLtoRGB
static void HSLtoRGB_Subfunction(unsigned int& c, const HReal& temp1, const HReal& temp2, const HReal& temp3)
{
    if((temp3 * 6) < 1)
        c = (unsigned int)((temp2 + (temp1 - temp2)*6*temp3)*100);
    else
        if((temp3 * 2) < 1)
            c = (unsigned int)(temp1*100);
        else
            if((temp3 * 3) < 2)
                c = (unsigned int)((temp2 + (temp1 - temp2)*(.66666 - temp3)*6)*100);
            else
                c = (unsigned int)(temp2*100);
    return;
}


// This function extracts the hue, saturation, and luminance from "color" 
// and places these values in h, s, and l respectively.
void RGBtoHSL(unsigned int color,unsigned int& h, unsigned int& s, unsigned int& l)
{
    unsigned int r = (unsigned int)GetRValue(color);
    unsigned int g = (unsigned int)GetGValue(color);
    unsigned int b = (unsigned int)GetBValue(color);

    HReal r_percent = ((HReal)r)/255;
    HReal g_percent = ((HReal)g)/255;
    HReal b_percent = ((HReal)b)/255;

    HReal max_color = 0;
    if((r_percent >= g_percent) && (r_percent >= b_percent))
    {
        max_color = r_percent;
    }
    if((g_percent >= r_percent) && (g_percent >= b_percent))
        max_color = g_percent;
    if((b_percent >= r_percent) && (b_percent >= g_percent))
        max_color = b_percent;

    HReal min_color = 0;
    if((r_percent <= g_percent) && (r_percent <= b_percent))
        min_color = r_percent;
    if((g_percent <= r_percent) && (g_percent <= b_percent))
        min_color = g_percent;
    if((b_percent <= r_percent) && (b_percent <= g_percent))
        min_color = b_percent;

    HReal L = 0;
    HReal S = 0;
    HReal H = 0;

    L = (max_color + min_color)/2;

    if(max_color == min_color)
    {
        S = 0;
        H = 0;
    }
    else
    {
        if(L < .50)
        {
            S = (max_color - min_color)/(max_color + min_color);
        }
        else
        {
            S = (max_color - min_color)/(2 - max_color - min_color);
        }
        if(max_color == r_percent)
        {
            H = (g_percent - b_percent)/(max_color - min_color);
        }
        if(max_color == g_percent)
        {
            H = 2 + (b_percent - r_percent)/(max_color - min_color);
        }
        if(max_color == b_percent)
        {
            H = 4 + (r_percent - g_percent)/(max_color - min_color);
        }
    }
    s = (unsigned int)(S*100);
    l = (unsigned int)(L*100);
    H = H*60;
    if(H < 0)
        H += 360;
    h = (unsigned int)H;
}

// This function converts the "color" object to the equivalent RGB values of
// the hue, saturation, and luminance passed as h, s, and l respectively
unsigned int HSLtoRGB(const unsigned int& h, const unsigned int& s, const unsigned int& l)
{
    unsigned int r = 0;
    unsigned int g = 0;
    unsigned int b = 0;

    HReal L = ((HReal)l)/100;
    HReal S = ((HReal)s)/100;
    HReal H = ((HReal)h)/360;

    if(s == 0)
    {
        r = l;
        g = l;
        b = l;
    }
    else
    {
        HReal temp1 = 0;
        if(L < .50)
        {
            temp1 = L*(1 + S);
        }
        else
        {
            temp1 = L + S - (L*S);
        }

        HReal temp2 = 2*L - temp1;

        HReal temp3 = 0;
        for(int i = 0 ; i < 3 ; i++)
        {
            switch(i)
            {
            case 0: // red
            {
                temp3 = H + .33333f;
                if(temp3 > 1)
                    temp3 -= 1;
                HSLtoRGB_Subfunction(r,temp1,temp2,temp3);
                break;
            }
            case 1: // green
            {
                temp3 = H;
                HSLtoRGB_Subfunction(g,temp1,temp2,temp3);
                break;
            }
            case 2: // blue
            {
                temp3 = H - .33333f;
                if(temp3 < 0)
                    temp3 += 1;
                HSLtoRGB_Subfunction(b,temp1,temp2,temp3);
                break;
            }
            default:
            {

            }
            }
        }
    }
    r = (unsigned int)((((HReal)r)/100)*255);
    g = (unsigned int)((((HReal)g)/100)*255);
    b = (unsigned int)((((HReal)b)/100)*255);
    return RGB(r,g,b);
}

unsigned int BrightenColor(unsigned int color,unsigned int amount)
{
    unsigned int h,s,l;

    RGBtoHSL(color,h,s,l);
    l+=amount;
    if ( l > 100 )
    {
        l = 100;
    }
    return HSLtoRGB(h,s,l);
}


unsigned int DarkenColor(unsigned int color,unsigned int amount)
{
    unsigned int h,s,l;

    RGBtoHSL(color,h,s,l);
    if ( amount >= l )
    {
        l = 0;
    }
    else
    {
        l-=amount;
    }
    return HSLtoRGB(h,s,l);
}
}
