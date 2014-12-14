#ifndef PARTICLE_H
#define PARTICLE_H

#include "Vec.hpp"
#include <vector>

class Boundary
{
    typedef Vec3<double> Vec;

    public :
    double psi; //density number
    Vec x,v; ///position and velocity

    public :
    Boundary()
    {
        psi=0.0;
        x=Vec(0.0);
        v=Vec(0.0);
    }

    Boundary(Vec& _x, Vec _v, double _psi)
    {
        x=_x; 
        v=_v;
        psi=_psi;
    }

    Boundary(const Boundary& b)
    {
        x=b.x; 
        v=b.v;
        psi=b.psi;
    }

    ~Boundary(){}
};

class Particle
{
    typedef Vec3<double> Vec;

    public :
    double rho, rho_adv, rho_corr, p, p_l, previousP, aii;
    Vec x, v, v_adv, f_adv, f_p, dii_fluid, dii_boundary, sum_dij, n;
    Vec3<float> c;
    vector<int> fluidNeighbor;
    vector<int> boundaryNeighbor;
    
    public :
    Particle()
    {
        rho = rho_adv = rho_corr = p = p_l = previousP = aii = 0.0;
        x = v = v_adv = f_adv = f_p = dii_fluid = dii_boundary = sum_dij = n = Vec(0.0);
        c = Vec3<float>(0.0);
        fluidNeighbor.clear();
        boundaryNeighbor.clear();
    }

    Particle(Vec& _x, Vec& _v, Vec3<float> _c = Vec3<float>(0,0,1))
    {
        x = _x;
        v = _v;
        c = _c;

        rho = rho_adv = rho_corr = p = p_l = previousP = aii = 0.0;
        v_adv = f_adv = f_p = dii_fluid = dii_boundary = sum_dij = n = Vec(0.0);
        fluidNeighbor.clear();
        boundaryNeighbor.clear();
    }

    Particle(const Particle& _p)
    {
        rho=_p.rho; 
        rho_adv=_p.rho_adv; 
        rho_corr=_p.rho_corr;
        p=_p.p; 
        p_l=_p.p_l;
        previousP=_p.previousP;
        aii=_p.aii; 

        x = _p.x;
        v = _p.v;
        v_adv = _p.v_adv;
        f_adv = _p.f_adv;
        f_p = _p.f_p;
        dii_fluid=_p.dii_fluid;
        dii_boundary=_p.dii_boundary;
        sum_dij=_p.sum_dij;
        n = _p.n;

        c = _p.c;

        fluidNeighbor=_p.fluidNeighbor;
        boundaryNeighbor=_p.boundaryNeighbor;
    }

    ~Particle(){}
};

#endif
