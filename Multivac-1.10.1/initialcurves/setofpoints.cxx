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


#ifndef FILE_INITIALCURVE_SETOFPOINTS_CXX


#include "setofpoints.hxx"


namespace Multivac
{



  ///////////////////////////////
  // CONSTRUCTORS & DESTRUCTOR //
  ///////////////////////////////


  //! Default constructor.
  template <class T>
  CSetOfPoints<T>::CSetOfPoints()
  {

  }


  //! Constructor.
  /*! Sets the front equal to 'InitialFront'.
    \param InitialFront the initial front.
   */
  template <class T>
  CSetOfPoints<T>::CSetOfPoints(Curve<T>& InitialFront)
  {

    Front.Copy(InitialFront);
    origin_out = true;

  }


  //! Main constructor.
  /*! Reads the file defining the front.
    \param InitialFrontFile file that stores points defining the front.
    \param orientation orientation of the front: -1 if the
    orientation is unknown, 0 if points are stored in
    the trigonometrical orientation and 1 if they are stored
    in the reverse orientation.
    \param origin_out_ position of the bottom left corner of the mesh.
    'true' if the corner is outside the front, 'false' if it is inside.
  */
  template <class T>
  CSetOfPoints<T>::CSetOfPoints(string InitialFrontFile,
				int orientation, bool origin_out_)
  {

    Front.ReadText(InitialFrontFile, orientation);
    origin_out = origin_out_;

  }


  //! Destructor.
  template <class T>
  CSetOfPoints<T>::~CSetOfPoints()  throw()
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
  void CSetOfPoints<T>::SetDistances(CMesh<T>& Mesh, CLevelSet<T>& Phi)
  {

    int i, j;

    int Nx = Mesh.GetNx();
    int Ny = Mesh.GetNy();
    T Delta_x = Mesh.GetDelta_x();
    T Delta_y = Mesh.GetDelta_y();
    T Xmin = Mesh.GetXmin();
    T Ymin = Mesh.GetYmin();
    T x, y;

    T dist_maj = Nx * Delta_x + Ny * Delta_y;

    List<Vector<T> >& points = Front.GetPoints();

    // Phi values.
    Matrix<T>& PhiValues = Phi.GetValues();

    PhiValues.Fill(Nx*Delta_x + Ny*Delta_y);

    // Temporary vectors.
    Vector<T> A(2), B(2);
    Vector<T> C(2), D(2);
    
    T norm, dist, u, vect;
    T tempx, tempy;
    
    /*** Slightly modifies the initial curve to properly handle distances ***/

    // Go to the first point of the front.
    points.GoToTheHead();

    x = (points.GetCurrentValue()(0) - Xmin) / Delta_x;
    y = (points.GetCurrentValue()(1) - Ymin) / Delta_y;

    // If the point is exactly on a mesh edge.
    if (floor(x) == x)
      // The point is slightly moved.
      if (points.GetCurrentValue()(0) <= Xmin)
	points.GetCurrent()->GetElement()(0) -= Delta_x * 1.e-5;
      else
	points.GetCurrent()->GetElement()(0) += Delta_x * 1.e-5;

    if (floor(y) == y)
      // The point is slightly moved.
      if (points.GetCurrentValue()(1) <= Ymin)
	points.GetCurrent()->GetElement()(1) -= Delta_y * 1.e-5;
      else
	points.GetCurrent()->GetElement()(1) += Delta_y * 1.e-5;

    // The same for all other points.
    while (points.GoToNext_StopAtTheTail())
      {
	x = (points.GetCurrentValue()(0) - Xmin) / Delta_x;
	y = (points.GetCurrentValue()(1) - Ymin) / Delta_y;

	// If the point is exactly on a mesh point.
	if (floor(x) == x)
	  // The point is slightly moved.
	  if (points.GetCurrentValue()(0) <= Xmin)
	    points.GetCurrent()->GetElement()(0) -= Delta_x * 1.e-5;
	  else
	    points.GetCurrent()->GetElement()(0) += Delta_x * 1.e-5;

	if (floor(y) == y)
	  // The point is slightly moved.
	  if (points.GetCurrentValue()(1) <= Ymin)
	    points.GetCurrent()->GetElement()(1) -= Delta_y * 1.e-5;
	  else
	    points.GetCurrent()->GetElement()(1) += Delta_y * 1.e-5;
      }

    /*** Sets the distances ***/

    // Go to the first point of the front.
    points.GoToTheHead();

    A.Copy(points.GetCurrentValue());

    bool out(true);

    x = Xmin;
    for (i=0; i<Nx; i++)
      {

	if (i==0)
	  out = origin_out;
	else
	  {
	    // Go to the first point of the front.
	    points.GoToTheHead();
	    A.Copy(points.GetCurrentValue());

	    out = (PhiValues(i-1, 0) > 0.0);
	    while (points.GoToNext_StopAtTheTail())
	      {
		// Gets the current point.
		B.Copy(points.GetCurrentValue());
		
		C(0) = x - Delta_x; C(1) = Ymin;
		D(0) = x; D(1) = Ymin;

		out = out ^ Intersect(A, B, C, D);

		A.Copy(B);
	      }
	  }

	if (out)
	  PhiValues(i, 0) = dist_maj;
	else
	  PhiValues(i, 0) = - dist_maj;

	y = Ymin + Delta_y;
	for (j=1; j<Ny; j++)
	  {
	    
	    // Go to the first point of the front.
	    points.GoToTheHead();
	    A.Copy(points.GetCurrentValue());

	    // For all points...
	    while (points.GoToNext_StopAtTheTail())
	      {
		// Gets the current point.
		B.Copy(points.GetCurrentValue());
		
		C(0) = x; C(1) = y - Delta_y;
		D(0) = x; D(1) = y;

		out = out ^ Intersect(A, B, C, D);

		A.Copy(B);
	      }
	    
	    if (out)
	      PhiValues(i, j) = dist_maj;
	    else
	      PhiValues(i, j) = - dist_maj;

	    y += Delta_y;
	  }

	x += Delta_x;

      }

    // Go to the first point of the front.
    points.GoToTheHead();
    
    A.Copy(points.GetCurrentValue());
    
    // For all points...
    while (points.GoToNext_StopAtTheTail())
      {
	// Gets the current point.
	B.Copy(points.GetCurrentValue());
	
	norm = DistanceBetween(A, B);
	norm *= norm;
	
	// If 'A' and 'B' are very close, ignore 'B'.
	if (norm > exp(-30.) )
	  {
	    
	    x = Xmin;
	    for (i=0; i<Nx; i++)
	      {
		y = Ymin;
		for (j=0; j<Ny; j++)
		  {
		    
		    tempx = x - A(0);
		    tempy = y - A(1);
		    
		    vect = (B(0)-A(0)) * tempy - (B(1)-A(1)) * tempx;

		    if ( fabs(vect) != 0 )
		      {
			
			u = ( tempx * (B(0)-A(0)) + tempy * (B(1) - A(1)) )
			  / norm;
			
			if (u <= 0)
			  {
			    if ( fabs(PhiValues(i, j)) >
				 ( dist = sqrt(tempx*tempx + tempy*tempy) ) )
			      PhiValues(i, j) = sign(PhiValues(i, j)) * dist;
			  }
			else if (u >= 1)
			  {
			    tempx = x - B(0);
			    tempy = y - B(1);
			    if ( fabs(PhiValues(i, j)) >
				 ( dist = sqrt(tempx*tempx + tempy*tempy) ) )
			      PhiValues(i, j) = sign(PhiValues(i, j)) * dist;
			  }
			else
			  {
			    tempx = x - A(0) - u * (B(0) - A(0));
			    tempy = y - A(1) - u * (B(1) - A(1));
			    if ( fabs(PhiValues(i, j)) >
				 ( dist = sqrt(tempx*tempx + tempy*tempy) ) )
			      PhiValues(i, j) = sign(PhiValues(i, j)) * dist;
			  }
			
		      }
		    y += Delta_y;
		  }
		x += Delta_x;
	      }

	  }

	// 'A' <- 'B' and 'B' is going to be the next point on the front.
	A.Copy(B);

      }
  
  }


  //! Returns the distance from a given point to the circle.
  /*!
    \param x point abscissa.
    \param y point ordinate.
    \return Distance from (x, y) to the circle.
  */
  template <class T>
  inline T CSetOfPoints<T>::GetDistance(T x, T y)
  {

    return Front.GetSignedDistanceToTheCurve(x, y);

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
  void CSetOfPoints<T>::GetProjectionOnMesh(T& x, T& y, CMesh<T>& Mesh)
  {

    throw CError_Undefined(string("void CSetOfPoints<T>::GetProjectionOnMesh")
			   + "(T& x, T& y, CMesh<T>& Mesh)");

  }


  //! Saves the circle.
  /*!
    \param CurveFile the descriptor of the file where data are saved.
    \warning This function doesn't save any point because the circle
    is analytically defined.
  */
  template <class T>
  void CSetOfPoints<T>::Save(string CurveFile) const
  {

  }



}  // namespace Multivac.


#define FILE_INITIALCURVE_SETOFPOINTS_CXX
#endif
