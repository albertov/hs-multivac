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


#ifndef FILE_INITIALCURVES_THREECIRCLES_HXX


#include "../errors.cxx"
#include <cstdio>
#include <cmath>


namespace Multivac
{



  ///////////////////
  // CTHREECIRCLES //
  ///////////////////

  //! The initial curve is defined by three circles.
  template <class T>
  class CThreeCircles: public CInitialCurve<T>
  {


    /**************
     * ATTRIBUTES *
     **************/

  protected:

    T CenterX1;
    T CenterY1;
    T CenterX2;
    T CenterY2;
    T CenterX3;
    T CenterY3;

    T Radius1;
    T Radius2;
    T Radius3;


    /*****************************
     * CONSTRUCTORS & DESTRUCTOR *
     *****************************/

  public:

    CThreeCircles()  throw();
    CThreeCircles(T CenterX1_, T CenterY1_, T Radius1_,
		  T CenterX2_, T CenterY2_, T Radius2_,
		  T CenterX3_, T CenterY3_, T Radius3_,
		  bool reversed_ = false)  throw();

    ~CThreeCircles()  throw();


    /***********
     * METHODS *
     ***********/
    
  public:
  
    virtual void SetDistances(CMesh<T>& Mesh, CLevelSet<T>& Phi);
    virtual T GetDistance(T x, T y);

    virtual void GetProjectionOnMesh(T& x, T& y, CMesh<T>& Mesh);

    virtual void Save(string CurveFile) const;

  };  // CThreeCircles.


}  // namespace Multivac.


#define FILE_INITIALCURVES_THREECIRCLES_HXX
#endif
