#include "./../include/hokusai/kernel.hpp"

namespace hokusai
{

AkinciKernel::AkinciKernel()
{
    h = 0;
    m_v1 = 0;
    m_v2 = 0;
}

AkinciKernel::AkinciKernel( HReal _h)
{
    h = _h;
    m_v1 = 32.0/(M_PI*pow(h,9.0));
    m_v2 = pow(h,6.0)/64.0;
    adhesion = 0.007/pow(h,3.25);
}

AkinciKernel::AkinciKernel( const AkinciKernel& k)
{
    h = k.h;
    m_v1 = k.m_v1;
    m_v2 = k.m_v2;
    adhesion = k.adhesion;
}

AkinciKernel::~AkinciKernel(){}

HReal AkinciKernel::cohesionValue( const HReal r )
{
    HReal value=0;
    if( (2.0*r>h) && (r<=h) )
    {
        //value=m_v1*pow(h-r,3.0)*pow(r,3.0);
        value=m_v1*(h-r)*(h-r)*(h-r)*r*r*r;
    }
    else if( (r>0.0) && (2.0*r<=h) )
    {
        //value=m_v1*(2.0*pow(h-r,3.0)*pow(r,3.0)-m_v2);
        value=m_v1*(2.0*(h-r)*(h-r)*(h-r)*r*r*r-m_v2);
    }
    else
        value=0.0;
    return value;
}

HReal AkinciKernel::adhesionValue( const HReal r)
{
    HReal value=0;
    if( (2.0*r)>h && (r<=h) )
        value=adhesion*pow(-4.0*r*r/h + 6.0*r -2.0*h,1.0/4.0);
    else
        value=0.0;
    return value;
}

MonaghanKernel::MonaghanKernel()
{
    m_h = 0;
    m_invH = 0;
    m_v = 0;
    m_g = 0;
}

MonaghanKernel::MonaghanKernel( HReal _h )
{
    m_h = _h;
    m_invH = 1.0/m_h;
    m_v = 1.0/(4.0*M_PI*m_h*m_h*m_h);
    m_g = 1.0/(4.0*M_PI*m_h*m_h*m_h);
}

MonaghanKernel::MonaghanKernel( const MonaghanKernel& k)
{
    m_h = k.m_h;
    m_invH = k.m_invH;
    m_v = k.m_v;
    m_g = k.m_g;
}

MonaghanKernel::~MonaghanKernel(){}

HReal MonaghanKernel::monaghanValue( const Vec3r & r )
{
    HReal value = 0.0;
    HReal q = r.length()*m_invH;
    if( q >= 0 && q < 1 )
    {
        value = m_v*( (2-q)*(2-q)*(2-q) - 4.0f*(1-q)*(1-q)*(1-q));
    }
    else if ( q >=1 && q < 2 )
    {
        value = m_v*( (2-q)*(2-q)*(2-q) );
    }
    else
    {
        value = 0.0f;
    }
    return value;
}

void MonaghanKernel::monaghanGradient( const Vec3r& r, Vec3r& gradient )
{
    HReal dist = r.length();
    HReal q = dist*m_invH;
    gradient.fill(0.0);
    if( q >= 0 && q < 1 )
    {
        HReal scalar = -3.0f*(2-q)*(2-q);
        scalar += 12.0f*(1-q)*(1-q);
        gradient = (m_g*m_invH*scalar/dist)*r;
    }
    else if ( q >=1 && q < 2 )
    {
        HReal scalar = -3.0f*(2-q)*(2-q);
        gradient = (m_g*scalar*m_invH/dist)*r;
    }
}

BoundaryKernel::BoundaryKernel()
{
    h = 0;
    cs = 0;
}

BoundaryKernel::~BoundaryKernel(){}

BoundaryKernel::BoundaryKernel( HReal _h, HReal _cs )
{
    h = _h;
    cs = _cs;
}

BoundaryKernel::BoundaryKernel(const BoundaryKernel &k)
{
    h = k.h;
    cs = k.cs;
}

HReal BoundaryKernel::gamma( HReal distance )
{
    HReal q = distance / h;
    HReal coeff = 0.0;
    if( q > 0 && q <= 0.666666667f ) //2.0/3.0
    {
        coeff = 0.666666667f;
    }
    else if( q > 0.666666667f && q < 1 )
    {
        coeff = q*(2 - 1.5*q);
    }
    else if( q >= 1 && q < 2 )
    {
        coeff = 0.5*(2-q)*(2-q);
    }
    coeff *= (0.02*cs*cs/distance);
    return coeff;
}

} //namespace hokusai
