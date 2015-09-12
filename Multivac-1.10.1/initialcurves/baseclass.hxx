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


#ifndef FILE_INITIALCURVES_BASECLASS_HXX


#include "../errors.cxx"
#include <cstdio>


namespace Multivac
{


  ///////////////////
  // CINITIALCURVE //
  ///////////////////

  //! Base class for initial curves.
  /*! Defines the initial curves interface.  All initial curves
    must be defined in the same way.
    \note
    This is an abstract class.
  */
  template <class T>
  class CInitialCurve
  {

    /**************
     * ATTRIBUTES *
     **************/

  protected:

    //! Current front representation.
    Curve<T> Curve;

    //! 'true' if outside and inside are swapped.
    bool reversed;


    /*****************************
     * CONSTRUCTORS & DESTRUCTOR *
     *****************************/

  public:

    CInitialCurve()  throw();

    virtual ~CInitialCurve()  throw();


    /***********
     * METHODS *
     ***********/

  public:
  
    virtual void SetDistances(CMesh<T>& Mesh, CLevelSet<T>& Phi) = 0;
    virtual T GetDistance(T x, T y) = 0;

    virtual void GetProjectionOnMesh(T& x, T& y, CMesh<T>& Mesh) = 0;

    virtual void Save(string CurveFile) const = 0;

  };  // CInitialCurve.


}  // namespace Multivac.


#define FILE_INITIALCURVES_BASECLASS_HXX
#endif
