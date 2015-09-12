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


#ifndef FILE_MESHES_BASECLASS_HXX


#include "../errors.cxx"
#include <cstdio>


namespace Multivac
{


  ///////////
  // CMESH //
  ///////////

  //! Base class for meshes.
  /*! Defines meshes interface.  All meshes must be defined in the same way.
    \note
    This is an abstract class.
  */
  template <class T>
  class CMesh
  {


    /**************
     * ATTRIBUTES *
     **************/

  protected:

    // Boundaries.
    //! Minimum abscissa.
    T Xmin;
    //! Maximum abscissa.
    T Xmax;
    //! Minimum ordinate.
    T Ymin;
    //! Maximum ordinate.
    T Ymax;

    //! Number of points along (x'x).
    int Nx;
    //! Space step along (x'x).
    T Delta_x;
    //! Number of points along (y'y).
    int Ny;
    //! Space step along (y'y).
    T Delta_y;


    /*****************************
     * CONSTRUCTORS & DESTRUCTOR *
     *****************************/

  public:

    CMesh()  throw();
    CMesh(T Xmin, T Xmax, T Ymin, T Ymax)  throw();

    virtual ~CMesh()  throw();


    /***********
     * METHODS *
     ***********/

  public:
  
    T GetXmin() const;
    T GetXmax() const;
    T GetYmin() const;
    T GetYmax() const;

    void SetXmin(T Xmin_);
    void SetXmax(T Xmax_);
    void SetYmin(T Ymin_);
    void SetYmax(T Ymax_);

    T GetDelta_x() const;
    void SetDelta_x(T Delta_x_);
    int GetNx() const;
    T GetDelta_y() const;
    void SetDelta_y(T Delta_y_);
    int GetNy() const;

    void GetClosestUpperPoint(T& x, T& y);

    virtual void Save(string XFile, string YFile) const = 0;
    virtual void SaveNonOrthogonalMesh(string PointsFile,
				       string EdgesFile,
				       string TrianglesFile) const = 0;

  };  // CMesh.


}  // namespace Multivac.


#define FILE_MESHES_BASECLASS_HXX
#endif
