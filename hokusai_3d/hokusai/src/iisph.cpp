#include "./../include/hokusai/iisph.hpp"

namespace hokusai
{

IISPHParams::IISPHParams( const MonaghanKernel& mKernel, const AkinciKernel& aKernel, const BoundaryKernel& bKernel ) :
    m_pKernel(mKernel),
    m_aKernel(aKernel),
    m_bKernel(bKernel)
{
}

IISPHParams::IISPHParams()
{

}

IISPHParams::~IISPHParams()
{

}

AkinciKernel& IISPHParams::akinciKernel()
{
    return m_aKernel;
}

MonaghanKernel& IISPHParams::monaghanKernel()
{
    return m_pKernel;
}

BoundaryKernel& IISPHParams::boundaryKernel()
{
    return m_bKernel;
}

}//namespace hokusai
