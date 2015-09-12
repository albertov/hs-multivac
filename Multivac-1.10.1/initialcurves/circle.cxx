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


#ifndef FILE_INITIALCURVE_CIRCLE_CXX


#include "circle.hxx"


namespace Multivac
{



  ///////////////////////////////
  // CONSTRUCTORS & DESTRUCTOR //
  ///////////////////////////////


  //! Default constructor.
  template <class T>
  CCircle<T>::CCircle()  throw()
  {

  }


  //! Main constructor.
  /*! Defines circle properties.
    \param CenterX_ center abscissa.
    \param CenterY_ center ordinate.
    \param Raduis_ circle radius.
  */
  template <class T>
  CCircle<T>::CCircle(T CenterX_, T CenterY_,
		      T Radius_, bool reversed_)  throw()
  {

    CenterX = CenterX_;
    CenterY = CenterY_;

    Radius = Radius_;

    this->reversed = reversed_;

  }


  //! Destructor.
  template <class T>
  CCircle<T>::~CCircle()  throw()
  {

  }



  /////////////
  // METHODS //
  /////////////


  //! Updates the level set function Phi by setting distances to the circle.
  /*!
    \param Mesh orthogonal mesh.
    \param Phi level set function defined on an orthogonal mesh.
    \exception Seldon::WrongCol attempt to reach a wrong column number.
    \exception Seldon::WrongRow attempt to reach a wrong row number.
  */
  template <class T>
  void CCircle<T>::SetDistances(CMesh<T>& Mesh, CLevelSet<T>& Phi)
  {

    int i, j;

    int Nx = Mesh.GetNx();
    int Ny = Mesh.GetNy();
    T Delta_x = Mesh.GetDelta_x();
    T Delta_y = Mesh.GetDelta_y();
    T X = Mesh.GetXmin();
    T Ymin = Mesh.GetYmin();
    T Y;
    T diffx, diffy;

    for (i=0; i<Nx; i++)
      {
	Y = Ymin;
	for (j=0; j<Ny; j++)
	  {
	    diffx = CenterX - X;
	    diffy = CenterY - Y;
	    if (this->reversed)
	      Phi(i, j) = Radius - sqrt( diffx * diffx + diffy * diffy );
	    else
	      Phi(i, j) = sqrt( diffx * diffx + diffy * diffy ) - Radius;
	    Y += Delta_y;
	  }
	X += Delta_x;
      }

  }


  //! Returns the distance from a given point to the circle.
  /*!
    \param x point abscissa.
    \param y point ordinate.
    \return Distance from (x, y) to the circle.
  */
  template <class T>
  inline T CCircle<T>::GetDistance(T x, T y)
  {

    T diffx = CenterX - x;
    T diffy = CenterY - y;
    if (this->reversed)
      return Radius - sqrt( diffx * diffx + diffy * diffy );
    else
      return sqrt( diffx * diffx + diffy * diffy ) - Radius;

  }

  
  //! Returns the closest mesh-point of the projection of (x, y)
  //! on the circle.
  /*!
    Let A = (x, y). Let B be the projection of A on the circle. Let C be
    the closest point to B that is on the mesh (included its vertices).
    On exit, C is returned through x and y: C = (x, y).
    \param x first coordinate of the point to be projected.
    \param y second coordinate of the point to be projected.
    \param Mesh mesh.
  */
  template <class T>
  void CCircle<T>::GetProjectionOnMesh(T& x, T& y, CMesh<T>& Mesh)
  {

    T x_C, y_C;

    T Delta_x = Mesh.GetDelta_x();
    T Delta_y = Mesh.GetDelta_y();

    T dir_x = x - CenterX;
    T dir_y = y - CenterY;
    T dist = sqrt(dir_x*dir_x + dir_y*dir_y);

    T x_B, y_B;
    
    // (x_B, y_B) is the projection of (x, y) on the circle.
    if (dist<exp(-20.))
      {
	x_B = CenterX + Radius;
	y_B = 0;
      }
    else
      {
	x_B = CenterX + Radius * dir_x / dist;
	y_B = CenterY + Radius * dir_y / dist;
      }

    T x_u, y_u;
    x_u = x_B; y_u = y_B;
    Mesh.GetClosestUpperPoint(x_u, y_u);

    T x_temp, y_temp, dist_temp;


    // M       N(x_u, y_u)
    // o-------o
    // |       |
    // |       |
    // |       |
    // o-------o O
    // P

    // For [MN].
    x_C = x_u; y_C = y_u;
    this->Intersection(x_C - Delta_y, y_C, x_C, y_C);
    dir_x = x - x_C;
    dir_y = y - y_C;
    dist = sqrt(dir_x*dir_x + dir_y*dir_y);
    

    // For [ON].
    x_temp = x_u; y_temp = y_u;
    this->Intersection(x_temp, y_temp - Delta_y, x_temp, y_temp);

    dir_x = x - x_temp;
    dir_y = y - y_temp;
    dist_temp = sqrt(dir_x*dir_x + dir_y*dir_y);
    if (dist_temp<dist)
      {
	dist = dist_temp;
	x_C = x_temp;
	y_C = y_temp;
      }
    

    // For [PM].
    x_temp = x_u - Delta_x; y_temp = y_u;
    this->Intersection(x_temp, y_temp - Delta_y, x_temp, y_temp);

    dir_x = x - x_temp;
    dir_y = y - y_temp;
    dist_temp = sqrt(dir_x*dir_x + dir_y*dir_y);
    if (dist_temp<dist)
      {
	dist = dist_temp;
	x_C = x_temp;
	y_C = y_temp;
      }

    // For [PO].
    x_temp = x_u; y_temp = y_u - Delta_y;
    this->Intersection(x_temp - Delta_x, y_temp, x_temp, y_temp);

    dir_x = x - x_temp;
    dir_y = y - y_temp;
    dist_temp = sqrt(dir_x*dir_x + dir_y*dir_y);
    if (dist_temp<dist)
      {
	dist = dist_temp;
	x_C = x_temp;
	y_C = y_temp;
      }

    x = x_C;
    y = y_C;

  }


  //! Saves the circle.
  /*!
    \param CurveFile the descriptor of the file where data are saved.
    \warning This function doesn't save any point because the circle
    is analytically defined.
  */
  template <class T>
  void CCircle<T>::Save(string CurveFile) const
  {

  }


  //! Returns the intersection between the circle and a segment.
  /*!
   */
  template <class T>
  bool CCircle<T>::Intersection(T x_A, T y_A, T& x_B, T& y_B) const
  {
    
    T dist_A, dist_B;
    T dir_x_A, dir_y_A;
    T dir_x_B, dir_y_B;

    dir_x_A = x_A - CenterX;
    dir_y_A = y_A - CenterY;

    dir_x_B = x_B - CenterX;
    dir_y_B = y_B - CenterY;

    dist_A = sqrt(dir_x_A*dir_x_A + dir_y_A*dir_y_A);
    dist_B = sqrt(dir_x_B*dir_x_B + dir_y_B*dir_y_B);

    if ((dist_A<Radius) && (dist_B<Radius))
      return false;

    if ((x_A==x_B) && (y_A==y_B))
      if (dist_A==Radius)
	return true;
      else
	return false;

    T slope = (y_B - y_A) / (x_B - x_A);
    T origin = y_A - slope * x_A;

    T a = 1.0 + slope*slope;
    T b = 2.0 * slope * origin;
    T c = origin*origin - Radius*Radius;

    T disc = sqrt(b*b - 4.0*a*c);
    T x_1, x_2;

    x_1 = (- b - disc) / (2.0 * a);
    x_2 = (- b + disc) / (2.0 * a);

    if ((x_2<x_A) && (x_2<x_B))
      return false;
    if ((x_1>x_A) && (x_1>x_B))
      return false;

    if ((x_2>x_B) && (x_2>x_A))
      {
	x_B = x_1;
	y_B = slope * x_1 + origin;
      }
    else if ((x_1<x_B) && (x_1<x_A))
      {
	x_B = x_2;
	y_B = slope * x_2 + origin;
      }
    else if (fabs(x_B-x_1)<fabs(x_B-x_2))
      {
	x_B = x_1;
	y_B = slope * x_1 + origin;
      }
    else
      {
	x_B = x_2;
	y_B = slope * x_2 + origin;
      }
    
    return true;

  }



}  // namespace Multivac.


#define FILE_INITIALCURVE_CIRCLE_CXX
#endif
