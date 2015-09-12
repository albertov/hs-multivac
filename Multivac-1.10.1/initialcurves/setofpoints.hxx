// Copyright (C) 2002-2004 Vivien Mallet
//
// This file is part of Multivac library.
// Multivac library provides front-tracking algorithms.
//
// Multivac is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// Multivac is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License (file "license") for more details.
//
// For more information, please see the Multivac home page:
//     http://spacetown.free.fr/fronts/


#ifndef FILE_INITIALCURVES_SETOFPOINTS_HXX


#include "../errors.cxx"
#include <cstdio>
#include <cmath>


namespace Multivac
{



  /////////////
  // CCIRCLE //
  /////////////

  //! The initial curve is a circle.
  template <class T>
  class CSetOfPoints: public CInitialCurve<T>
  {


    /**************
     * ATTRIBUTES *
     **************/

  protected:
    Curve<T> Front;
    //! Is the bottom left corner outside the front?
    bool origin_out;


    /*****************************
     * CONSTRUCTORS & DESTRUCTOR *
     *****************************/

  public:

    CSetOfPoints();
    CSetOfPoints(Curve<T>& InitialFront);
    CSetOfPoints(string InitialFrontFile, int orientation, bool origin_out_);

    ~CSetOfPoints()  throw();


    /***********
     * METHODS *
     ***********/
    
  public:
  
    virtual void SetDistances(CMesh<T>& Mesh, CLevelSet<T>& Phi);
    virtual T GetDistance(T x, T y);

    virtual void GetProjectionOnMesh(T& x, T& y, CMesh<T>& Mesh);

    virtual void Save(string CurveFile) const;

  };  // CSetOfPoints.


}  // namespace Multivac.


#define FILE_INITIALCURVES_SETOFPOINTS_HXX
#endif
