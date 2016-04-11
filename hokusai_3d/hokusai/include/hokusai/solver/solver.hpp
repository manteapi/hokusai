/*
* Copyright 2015 Pierre-Luc Manteaux
*
*This file is part of Hokusai.
*
*Hokusai is free software: you can redistribute it and/or modify
*it under the terms of the GNU General Public License as published by
*the Free Software Foundation, either version 3 of the License, or
*(at your option) any later version.
*
*Hokusai is distributed in the hope that it will be useful,
*but WITHOUT ANY WARRANTY; without even the implied warranty of
*MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*GNU General Public License for more details.
*
*You should have received a copy of the GNU General Public License
*along with Hokusai.  If not, see <http://www.gnu.org/licenses/>.
*
* Author : Pierre-Luc Manteaux
* Contact : pierre-luc.manteaux@inria.fr
*/

#ifndef HOKUSAI_SOLVER_HPP
#define HOKUSAI_SOLVER_HPP

#include <iostream>
#include <fstream>
#include <ctime>
#include <memory>

#include <aljabr/AljabrCore>
#include "./../utils.hpp"
#include "./../common.hpp"
#include "./../particle.hpp"
#include "./../gridUtility.hpp"
#include "./../triMesh.hpp"
#include "./../sampler.hpp"
#include "./../particleSource.hpp"

namespace hokusai
{

class BaseSolver
{
public:
    ~BaseSolver(){}
    BaseSolver(){}
    virtual void init() = 0;
    virtual void solve() = 0;
};
typedef std::shared_ptr<BaseSolver> BaseSolverPtr;

class IISPHSolver : public BaseSolver
{
public:
    IISPHSolver(){}
    ~IISPHSolver(){}
    virtual void init(){}
    virtual void solve(){}

};
typedef std::shared_ptr<IISPHSolver> IISPHSolverPtr;

}


#endif // SOLVER_HPP
