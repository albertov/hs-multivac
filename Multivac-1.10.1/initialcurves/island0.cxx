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


#ifndef FILE_INITIALCURVE_ISLAND0_CXX


#include "island0.hxx"


namespace Multivac
{



  ///////////////////////////////
  // CONSTRUCTORS & DESTRUCTOR //
  ///////////////////////////////


  //! Default constructor.
  template <class T>
  CIsland0<T>::CIsland0()  throw()
  {

  }


  //! Main constructor.
  /*! Defines circles properties.
    \param CenterX1_ first circle center abscissa.
    \param CenterY1_ first circle center ordinate.
    \param Raduis1_ first circle radius.
    \param CenterX2_ second circle center abscissa.
    \param CenterY2_ second circle center ordinate.
    \param Raduis2_ second circle radius.
    \param CenterX3_ island center abscissa.
    \param CenterY3_ island center ordinate.
    \param Raduis3_ island radius.
  */
  template <class T>
  CIsland0<T>::CIsland0(T CenterX1_, T CenterY1_, T Radius1_,
			T CenterX2_, T CenterY2_, T Radius2_,
			T CenterX3_, T CenterY3_, T Radius3_,
			bool reversed_)  throw()
  {

    CenterX1 = CenterX1_;
    CenterY1 = CenterY1_;

    Radius1 = Radius1_;

    CenterX2 = CenterX2_;
    CenterY2 = CenterY2_;

    Radius2 = Radius2_;
 
    CenterX3 = CenterX3_;
    CenterY3 = CenterY3_;

    Radius3 = Radius3_;

    this->reversed = reversed_;

  }


  //! Destructor.
  template <class T>
  CIsland0<T>::~CIsland0()  throw()
  {

  }



  /////////////
  // METHODS //
  /////////////


  //! Updates the level set function Phi by setting distances to the front.
  /*! The distance outside the first (biggest) circle is positive,
    the distance inside the island is positive and the distance in between
    the two circles is negative.
    \param Mesh orthogonal mesh.
    \param Phi level set function defined on an orthogonal mesh.
    \exception Seldon::WrongCol attempt to reach a wrong column number.
    \exception Seldon::WrongRow attempt to reach a wrong row number.
  */
  template <class T>
  void CIsland0<T>::SetDistances(CMesh<T>& Mesh, CLevelSet<T>& Phi)
  {

    int i, j;
    T distance;
    T diffx, diffy, diffx0, diffy0;

    int Nx = Mesh.GetNx();
    int Ny = Mesh.GetNy();
    T Delta_x = Mesh.GetDelta_x();
    T Delta_y = Mesh.GetDelta_y();
    T X = Mesh.GetXmin();
    T Ymin = Mesh.GetYmin();
    T Y;

    for (i=0; i<Nx; i++)
      {
	Y = Ymin;
	for (j=0; j<Ny; j++)
	  {
	    diffx = CenterX1 - X;
	    diffy = CenterY1 - Y;
	    diffx0 = CenterX2 - X;
	    diffy0 = CenterY2 - Y;
	    if ( (distance =
		  min(sqrt( diffx * diffx + diffy * diffy ) - Radius1,
		      sqrt( diffx0 * diffx0 + diffy0 * diffy0 )
		      - Radius2) )
		 < 0)
	      {
		diffx = CenterX3 - X;
		diffy = CenterY3 - Y;
		Phi(i, j) = max( distance,
				 - sqrt( diffx * diffx + diffy * diffy )
				 + Radius3 );
	      }
	    else
	      Phi(i, j) = distance;
	    if (this->reversed)
	      Phi(i, j) = - Phi(i, j);
	    Y += Delta_y;
	  }
	X += Delta_x;
      }

  }


  //! Returns the distance from a given point to the front.
  /*!
    \param x point abscissa.
    \param y point ordinate.
    \return Distance from (x, y) to circles.
  */
  template <class T>
  inline T CIsland0<T>::GetDistance(T X, T Y)
  {
    
    T distance;

    T diffx = CenterX1 - X;
    T diffy = CenterY1 - Y;
    T diffx0 = CenterX2 - X;
    T diffy0 = CenterY2 - Y;
    if ( (distance =
	  min(sqrt( diffx * diffx + diffy * diffy ) - Radius1,
	      sqrt( diffx0 * diffx0 + diffy0 * diffy0 ) - Radius2) ) < 0)
      {
	diffx = CenterX3 - X;
	diffy = CenterY3 - Y;
	distance = max( distance,
			- sqrt( diffx * diffx + diffy * diffy ) + Radius3 );
      }

    if (this->reversed)
      distance = - distance;

    return distance;

  }

  
  //! Returns the closest mesh-point of the projection of (x, y) on the curve.
  /*!
    Let A = (x, y). Let B be the projection of A on the curve. Let C be
    the closest point to B that is on the mesh (included its vertices).
    On exit, C is returned through x and y: C = (x, y).
    \param x first coordinate of the point to be projected.
    \param y second coordinate of the point to be projected.
    \param Mesh mesh.
    \warning Undefined function.
  */
  template <class T>
  void CIsland0<T>::GetProjectionOnMesh(T& x, T& y, CMesh<T>& Mesh)
  {

    throw CError_Undefined
      ("void CIsland0<T>::GetProjectionOnMesh(T& x, T& y, CMesh<T>& Mesh)");

  }


  //! Saves front.
  /*!
    \param CurveFile the descriptor of the file where data are saved.
    \warning This function doesn't save any point because circles
    are analytically defined.
  */
  template <class T>
  void CIsland0<T>::Save(string CurveFile) const
  {

  }



}  // namespace Multivac.


#define FILE_INITIALCURVE_ISLAND0_CXX
#endif
