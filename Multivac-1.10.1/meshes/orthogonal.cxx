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


#ifndef FILE_MESHES_ORTHOGONAL_CXX


#include "orthogonal.hxx"


namespace Multivac
{



  ///////////////////////////////
  // CONSTRUCTORS & DESTRUCTOR //
  ///////////////////////////////


  //! Default constructor.
  template <class T>
  COrthogonalMesh<T>::COrthogonalMesh()  throw()
  {

    this->Delta_x = 0.0;
    this->Delta_y = 0.0;

  }


  //! Main constructor.
  /*!
    
  */
  template <class T>
  COrthogonalMesh<T>::COrthogonalMesh(T Xmin_, T Xmax_, T Ymin_, T Ymax_,
				      int Nx_, int Ny_)  throw():
    CMesh<T>(Xmin_, Xmax_, Ymin_, Ymax_)
  {

    this->Nx = Nx_;
    this->Ny = Ny_;

    this->Delta_x = (this->Xmax - this->Xmin) / T(this->Nx - 1);
    this->Delta_y = (this->Ymax - this->Ymin) / T(this->Ny - 1);

  }


  //! Destructor.
  template <class T>
  COrthogonalMesh<T>::~COrthogonalMesh()  throw()
  {

  }



  /////////////
  // METHODS //
  /////////////


  //! Saves the mesh.
  /*! Saves the orthogonal mesh.  Abscissae and ordinates of grid points
    are saved.
    \param XFile files where abscissae will be saved.
    \param YFile files where ordinates will be saved.
  */
  template <class T>
  void COrthogonalMesh<T>::Save(string XFile, string YFile) const
  {

    int i;

    Vector<T> XGrid(this->Nx);
    T X = this->Xmin;
    Vector<T> YGrid(this->Ny);
    T Y = this->Ymin;

    for(i=0; i<this->Nx; i++)
      {
	XGrid(i) = X;
	X += this->Delta_x;
      }

    for(i=0; i<this->Ny; i++)
      {
	YGrid(i) = Y;
	Y += this->Delta_y;
      }

    XGrid.WriteText(XFile);
    YGrid.WriteText(YFile);

  }


  //! Saves non-orthogonal mesh.
  /*! Nothing is done because the current mesh is an orthogonal mesh.
    \param PointsFile files where mesh points will be saved.
    \param EdgesFile files where mesh edges will be saved.
    \param TrianglesFile files where mesh triangles will be saved.
  */
  template <class T>
  void COrthogonalMesh<T>::SaveNonOrthogonalMesh(string PointsFile,
						 string EdgesFile,
						 string TrianglesFile) const
  {

  }


}  // namespace Multivac.


#define FILE_MESHES_ORTHOGONAL_CXX
#endif
