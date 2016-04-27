#ifndef HOKUSAI_IISPH_HPP
#define HOKUSAI_IISPH_HPP

#include "kernel.hpp"

namespace hokusai
{

class IISPHParams
{
public:
    IISPHParams( const MonaghanKernel& mKernel, const AkinciKernel& aKernel, const BoundaryKernel& bKernel );
    IISPHParams();
    ~IISPHParams();

    AkinciKernel& akinciKernel();
    MonaghanKernel& monaghanKernel();
    BoundaryKernel& boundaryKernel();

private:

    AkinciKernel m_aKernel;
    MonaghanKernel m_pKernel;
    BoundaryKernel m_bKernel;
};

}//namespace hokusai

#endif //HOKUSAI_IISPH_HPP
