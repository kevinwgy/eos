/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#pragma once
#include <stdarg.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <cassert>

using std::string;

/**************************
 * Utility functions
 **************************
*/
//--------------------------------------------------
//! MPI Rank 0 will print to stdout
void print(const char format[],...);
//--------------------------------------------------
//! MPI Rank 0 will print to stdout in red color
void print_error(const char format[],...);
//--------------------------------------------------
//! MPI Rank 0 will print to a file
void print(FILE* fd, const char format[],...);
//--------------------------------------------------
//! Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const string getCurrentDateTime();
//--------------------------------------------------
//! Call MPI_Finalize and exit (with error)
void exit_mpi();
//--------------------------------------------------
//! Check for NAN
template <class T>
inline int m2c_isnan(const T& t) {return (t != t);}
//--------------------------------------------------

//--------------------------------------------------
template<typename Functor>
void tabulate2Dfunction_uniform(Functor fun, double xmin, double xmax,
                                int Nx, double ymin, double ymax, int Ny,
                                std::vector<std::vector<double> > &result)
{
  //NOTE: fun must be of type "double fun(double x, double y)".

  //Check for the case of a 1D plot
  if(xmin==xmax || Nx==1) {
    result.resize(1);
    assert(Ny>=1);
    double dy = Ny==1 ? 0.0 : (ymax-ymin)/(Ny-1);
    result[0].reserve(Ny);
    for(int i=0; i<Ny; i++)
      result[0].push_back(fun(xmin, ymin+i*dy));
    return;
  }
  else if(ymin==ymax || Ny==1) {
    result.resize(1);
    assert(Nx>=1);
    double dx = Nx==1 ? 0.0 : (xmax-xmin)/(Nx-1);
    result[0].reserve(Nx);
    for(int i=0; i<Nx; i++)
      result[0].push_back(fun(xmin+i*dx, ymin));
    return;
  }

  // 2D
  assert(Nx>1 && Ny>1);
  double dx = (xmax-xmin)/(Nx-1);
  double dy = (ymax-ymin)/(Ny-1);
  result.resize(Ny);
  double y = ymin;
  for(auto&& res : result) {
    res.reserve(Nx);
    for(int i=0; i<Nx; i++)
      res.push_back(fun(xmin+i*dx, y));
    y += dy;
  }
}

//--------------------------------------------------
