#ifndef HOKUSAI_KERNEL_HPP
#define HOKUSAI_KERNEL_HPP

#include "./../include/hokusai/common.hpp"

namespace hokusai
{

class AkinciKernel
{
public :
    AkinciKernel();
    AkinciKernel(HReal _h);
    AkinciKernel(const AkinciKernel& k);
    ~AkinciKernel();
    HReal cohesionValue( const HReal r);
    HReal adhesionValue(const HReal r);
    HReal h,m_v1,m_v2,adhesion;
};

//Monaghan 3D kernel
class MonaghanKernel
{
public :

    MonaghanKernel();
    MonaghanKernel( HReal _h );
    MonaghanKernel( const MonaghanKernel& k );
    ~MonaghanKernel();

public :

    HReal m_h;
    HReal m_invH;
    HReal m_v;
    HReal m_g;

public :

    HReal monaghanValue( const Vec3r& r );
    void monaghanGradient( const Vec3r& r, Vec3r& gradient );
};

class BoundaryKernel
{
public :

    BoundaryKernel();
    BoundaryKernel( HReal h, HReal cs );
    BoundaryKernel( const BoundaryKernel& k);
    ~BoundaryKernel();

public :

    HReal h;
    HReal cs;

public :

    HReal gamma( HReal distance );
};

}//namespace hokusai

#endif //HOKUSAI_KERNEL_HPP
