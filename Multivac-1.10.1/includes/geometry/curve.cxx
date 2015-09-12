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


#ifndef FILE_CURVE_CXX


#include "curve.hxx"
#include <cstdio>


namespace Multivac
{


  /***************
   * CONSTRUCTOR *
   ***************/
  

  //! Default constructor.
  template <class T>
  Curve<T>::Curve()  throw()
  {

    n_=0;
    orientation_ = -1;

  }


  /**************
   * DESTRUCTOR *
   **************/
  

  //! Destructor.
  template <class T>
  Curve<T>::~Curve()  throw()
  {

    n_=0;

  }
  
  
  /*****************
   * BASIC METHODS *
   *****************/
  
  //! Builds a front (a sequence of sorted points).
  /*! A set of points is taken and sorted. The result is a sorted list
    of points. The first point is an arbitrary point.  The (i+1)-th point
    is the closest point to the i-th point among points whose index is
    strictly greater than i.
    \param points initial list of points.
    \note 'points' is the input and the output.
  */
  template <class T>
  void Curve<T>::
  InitWithUnsortedPoints_ClosestPointMethod(List<Vector<T> >& points)
  {

    points_.ClearAll();
    n_=0;

    Cell<Vector<T> >* closest;
    T dist, dist0, diff, diff0;

    if (!points.IsEmpty())
      {
	points.GoToTheHead();
	points_.AddAtTheEnd(points.RemoveCurrent());
	n_++;
      }

    while (!points.IsEmpty())
      {
	points.GoToTheHead();
	closest = points.GetCurrent();
	diff = (points_.GetTailValue())(0) - (points.GetCurrentValue())(0);
	diff0 = (points_.GetTailValue())(1) - (points.GetCurrentValue())(1);
	dist = sqrt( diff * diff + diff0 * diff0 );

	while (points.GoToNext_StopAtTheTail())
	  {
	    diff = (points_.GetTailValue())(0)
	      - (points.GetCurrentValue())(0);
	    diff0 = (points_.GetTailValue())(1)
	      - (points.GetCurrentValue())(1);
	    dist0 = sqrt( diff * diff + diff0 * diff0 );
	    if (dist0<dist)
	      {
		closest = points.GetCurrent();
		dist = dist0;
	      }
	  }

	points_.AddAtTheEnd(points.Remove(closest));
	n_++;

      }

  }


  //! Adds a point in the curve.
  /*! Adds the point at the end of the list. The curve is not sorted.
    \param point point to be added.
  */
  template <class T>
  inline void Curve<T>::AddPoint(Vector<T>& point)
  {

    if (n_>=0)
      n_++;

    points_.AddAtTheEnd(point);

  }


  //! Copies a curve.
  /*! \param C curve to be copied.
   */
  template <class T>
  inline void Curve<T>::Copy(Curve<T>& C)
  {

    int i = 0;
    List<Vector<T> >& Cpoints = C.GetPoints();

    points_.ClearAll();
    Cpoints.GoToTheHead();

    if (!Cpoints.IsEmpty())
      {
	points_.AddAtTheEnd(Cpoints.GetCurrentValue());
	i++;
      }

    while (Cpoints.GoToNext_StopAtTheTail())
      {
	points_.AddAtTheEnd(Cpoints.GetCurrentValue());
	i++;
      }

    n_ = i;

  }


  //! Returns the number of points.
  /*! \return Number of points.
    \note The function returns -1 is the number of points is unknown.
  */
  template <class T>
  inline int Curve<T>::GetNbPoints() const
  {

    return n_;

  }


  //! Returns list that stores points.
  /*! \return The list that stores points.
   */
  template <class T>
  List<Vector<T> >& Curve<T>::GetPoints()
  {

    return points_;

  }


  //! Returns orientation.
  /*! \return -1 if the orientation is unknown, 0 if points are
    stored in the trigonometrical orientation and 1 if they are
    stored in the reverse orientation.
  */
  template <class T>
  int Curve<T>::GetOrientation() const
  {

    return orientation_;

  }


  //! Sets the orientation.
  /*! \param orientation the new orientation that is:
    -1 if the orientation is unknown, 0 if points are
    stored in the trigonometrical orientation and 1 if they are
    stored in the reverse orientation.
    \return The previous orientation.
  */
  template <class T>
  int Curve<T>::SetOrientation(int orientation)
  {

    int temp = orientation_;

    orientation_ = orientation;

    return temp;

  }


  //! Removes all points.
  template <class T>
  void Curve<T>::ClearAll()
  {

    n_=0;
    points_.ClearAll();

  }


  //! Returns the distance from a given point to the curve.
  /*! \param x point abscissa.
    \param y point ordinate.
    \return Distance from (x, y) to the curve.
    \note The curve is considered to be a set of segments and
    to be closed.
  */
  template <class T>
  T Curve<T>::GetDistanceToTheCurve(T x, T y)
  {
    
    // Temporary vectors.
    Vector<T> A(2);
    Vector<T> B(2);
    
    T dist, u, res;
    T tempx, tempy;
    
    // Go to the first point of the front.
    points_.GoToTheHead();

    // Stores the first point in 'A' and adds this point at the tail.
    // Adding the first point at the tail is convenient in order to
    // compute distances to the segment between the first point and the
    // last point of the front.
    if (!points_.IsEmpty())
      {
	A.Copy(points_.GetCurrentValue());
	points_.AddAtTheEnd(A);
	tempx = x - A(0);
	tempy = y - A(1);
	res = sqrt( tempx*tempx + tempy*tempy );
      }
    else
      return T(-1);

    // For all points...
    while (points_.GoToNext_StopAtTheTail())
      {

	// Gets the current point.
	B.Copy(points_.GetCurrentValue());

	dist = DistanceBetween(A, B);

	// If 'A' and 'B' are very close, just computes
	// the distance to those two points.
	if (dist < exp(-30.) )
	  {
	    
	    tempx = x - A(0);
	    tempy = y - A(1);
	    res = min(res, sqrt( tempx*tempx + tempy*tempy ) );
	    tempx = x - B(0);
	    tempy = y - B(1);
	    res = min(res, sqrt( tempx*tempx + tempy*tempy ) );
	    
	  }
	else  // 'A' and 'B' are not very close.
	  {
	    
	    tempx = x - A(0);
	    tempy = y - A(1);

	    u = ( tempx * (B(0)-A(0)) + tempy * (B(1) - A(1)) )
	      / (dist * dist);
	    
	    if (u < 0)
	      res = min( res, sqrt( tempx*tempx + tempy*tempy ) );
	    else if (u > 1)
	      {
		tempx = x - B(0);
		tempy = y - B(1);
		res = min( res, sqrt( tempx*tempx + tempy*tempy ) );
	      }
	    else
	      {
		tempx = x - A(0) - u * (B(0) - A(0));
		tempy = y - A(1) - u * (B(1) - A(1));
		res = min( res, sqrt( tempx*tempx + tempy*tempy ) );
	      }
	    
	  }
	
	// 'A' <- 'B' and 'B' is going to be the next point on the front.
	A.Copy(B);
	
      }
  
    // Removes the last point (which is currently the first point --
    // see before the previous 'while' loop).
    points_.GoToTheTail();
    points_.DeleteCurrent();

    return res;
    
  }


  //! Projects a given point on the curve.
  /*! Projects a given point on the curve and returns the distance
    to the curve.  The parameter X is set to the projection.
    \param X point to be projected.
    \return Distance from (x, y) to the curve.
    \note The curve is considered to be a set of segments and
    to be closed.
  */
  template <class T>
  T Curve<T>::GetProjection(Vector<T>& X)
  {
    
    // Temporary vectors.
    Vector<T> A(2);
    Vector<T> B(2);
    
    T dist, u, res, res0;
    T tempx, tempy;
    
    T x = X(0);
    T y = X(1);

    // Go to the first point of the front.
    points_.GoToTheHead();

    // Stores the first point in 'A' and adds this point at the tail.
    // Adding the first point at the tail is convenient in order to
    // compute distances to the segment between the first point and the
    // last point of the front.
    if (!points_.IsEmpty())
      {
	A.Copy(points_.GetCurrentValue());
	points_.AddAtTheEnd(A);
	tempx = x - A(0);
	tempy = y - A(1);
	res = sqrt( tempx*tempx + tempy*tempy );
	X.Copy(A);
      }
    else
      return T(-1);

    // For all points...
    while (points_.GoToNext_StopAtTheTail())
      {

	// Gets the current point.
	B.Copy(points_.GetCurrentValue());

	dist = DistanceBetween(A, B);

	// If 'A' and 'B' are very close, just computes
	// the distance to those two points.
	if (dist < exp(-30.) )
	  {
	    
	    tempx = x - A(0);
	    tempy = y - A(1);
	    res0 = sqrt( tempx*tempx + tempy*tempy );
	    if (res0 < res)
	      {
		res = res0;
		X.Copy(A);
	      }
	    tempx = x - B(0);
	    tempy = y - B(1);
	    res0 = sqrt( tempx*tempx + tempy*tempy );
	    if (res0 < res)
	      {
		res = res0;
		X.Copy(B);
	      }
	    
	  }
	else  // 'A' and 'B' are not very close.
	  {
	    
	    tempx = x - A(0);
	    tempy = y - A(1);

	    u = ( tempx * (B(0)-A(0)) + tempy * (B(1) - A(1)) )
	      / (dist * dist);
	    
	    if (u < 0)
	      {
		res0 = sqrt( tempx*tempx + tempy*tempy );
		if (res0 < res)
		  {
		    res = res0;
		    X.Copy(A);
		  }
	      }
	    else if (u > 1)
	      {
		tempx = x - B(0);
		tempy = y - B(1);
		res0 = sqrt( tempx*tempx + tempy*tempy );
		if (res0 < res)
		  {
		    res = res0;
		    X.Copy(B);
		  }
	      }
	    else
	      {
		tempx = x - A(0) - u * (B(0) - A(0));
		tempy = y - A(1) - u * (B(1) - A(1));
		res0 = sqrt( tempx*tempx + tempy*tempy );
		if (res0 < res)
		  {
		    res = res0;
		    X(0) = A(0) + u * (B(0) - A(0));
		    X(1) = A(1) + u * (B(1) - A(1));
		  }
	      }
	    
	  }
	
	// 'A' <- 'B' and 'B' is going to be the next point on the front.
	A.Copy(B);
	
      }
  
    // Removes the last point (which is currently the first point --
    // see before the previous 'while' loop).
    points_.GoToTheTail();
    points_.DeleteCurrent();

    return res;
    
  }


  //! Returns the signed distance from a given point to the curve.
  /*! \param x point abscissa.
    \param y point ordinate.
    \return Signed distance from (x, y) to the curve.
    \note The curve is considered to be a set of segments and
    to be closed.
  */
  template <class T>
  T Curve<T>::GetSignedDistanceToTheCurve(T x, T y)
  {
    
    T limit = exp(-25.);

    // Temporary vectors.
    Vector<T> A(2), A0(2);
    Vector<T> B(2), B0(2);
    
    T dist, u, u_old, res, res0;
    T tempx, tempy;
    
    // Go to the first point of the front.
    points_.GoToTheHead();

    if (!points_.IsEmpty())
      {
	A.Copy(points_.GetCurrentValue());
	tempx = x - A(0);
	tempy = y - A(1);
	res = sqrt( tempx*tempx + tempy*tempy );
	u_old = 0;
	A0.Copy(A);
	(points_.GetCurrent())->GetPrevious()->GetElement(B0);
      }
    else
      return T(-1);

    // For all points...
    while (points_.GoToNext_StopAtTheTail())
      {

	// Gets the current point.
	B.Copy(points_.GetCurrentValue());

	dist = DistanceBetween(A, B);

	// If 'A' and 'B' are very close, just computes
	// the distance to those two points.
	if (dist < exp(-30.) )
	  {
	    
	    tempx = x - A(0);
	    tempy = y - A(1);
	    res0 = sqrt( tempx*tempx + tempy*tempy );
	    if (res0 < res)
	      {
		res = res0;
		A0.Copy(A);
		B0.Copy(B);
	      }
	    tempx = x - B(0);
	    tempy = y - B(1);
	    res0 = sqrt( tempx*tempx + tempy*tempy );
	    if (res0 < res)
	      {
		res = res0;
		A0.Copy(A);
		B0.Copy(B);
	      }
	    u_old = 0;
	    
	  }
	else  // 'A' and 'B' are not very close.
	  {
	    
	    tempx = x - A(0);
	    tempy = y - A(1);

	    u = ( tempx * (B(0)-A(0)) + tempy * (B(1) - A(1)) )
	      / (dist * dist);
	    
	    if (u <= 0)
	      {
		res0 = sqrt( tempx*tempx + tempy*tempy );
		if ( (res0 <= res)
		     && ( (res - res0 > limit)
			  || (fabs(u-0.5) < fabs(u_old-0.5))) )
		  {
		    u_old = u;
		    res = res0;
		    A0.Copy(A);
		    B0.Copy(B);
		  }
	      }
	    else if (u >= 1)
	      {
		tempx = x - B(0);
		tempy = y - B(1);
		res0 = sqrt( tempx*tempx + tempy*tempy );
		if ( (res0 <= res)
		     && ( (res - res0 > limit)
			  || (fabs(u-0.5) < fabs(u_old-0.5))) )
		  {
		    u_old = u;
		    res = res0;
		    A0.Copy(A);
		    B0.Copy(B);
		  }
	      }
	    else
	      {
		tempx = x - A(0) - u * (B(0) - A(0));
		tempy = y - A(1) - u * (B(1) - A(1));
		res0 = sqrt( tempx*tempx + tempy*tempy );
		if ( (res0 <= res)
		     && ( (res - res0 > limit)
			  || (fabs(u-0.5) < fabs(u_old-0.5))) )
		  {
		    u_old = u;
		    res = res0;
		    A0.Copy(A);
		    B0.Copy(B);
		  }
	      }
	    
	  }

	// 'A' <- 'B' and 'B' is going to be the next point on the front.
	A.Copy(B);
	
      }
  
    if ( ((B0(0)-A0(0)) * (y-A0(1))) - ((B0(1)-A0(1)) * (x-A0(0))) > 0.0 )
      res = - res;

    if (orientation_ == 1)
      res = -res;

    return res;
    
  }


  /*****************
   * I/O FUNCTIONS *
   *****************/


  //! Saves the curve.
  /*!
    \param f file where points will be saved.
  */
  template <class T>
  void Curve<T>::WriteText(ofstream& f)
  {

    if (!points_.IsEmpty())
      {
	points_.GoToTheHead();
	
	(points_.GetCurrentValue()).WriteText(f);
	f << endl;

	while (points_.GoToNext_StopAtTheTail())
	  {
	    (points_.GetCurrentValue()).WriteText(f);
	    f << endl;
	  }
      }

  }


  //! Reads curve coordinates stored in a text file.
  /*! Reads curve coordinates stored in a text file. For each
    point of the curve, coordinates are stored on a single line.
    Each line of the input file represents a single point.
    The order of points in the curve is the same as the order
    in the file.
    \param f file name.
    \param orientation orientation of the curve: -1 if the
    orientation is unknown, 0 if points are stored in
    the trigonometrical orientation and 1 if they are stored
    in the reverse orientation.
  */
  template <class T>
  void Curve<T>::ReadText(string file_name, int orientation)
  {
    
    int i;
    string line_buffer;
    
    ifstream file;
    file.open(file_name.c_str(), ifstream::in);

    if (!file.good())
      throw IOError("Curve::ReadText(string, int)",
		    "Unable to open file \"" + file_name + "\".");
    
    Vector<T> point(2);
    int length = 0;
    
    points_.ClearAll();
    n_ = 0;
    orientation_ = orientation;
    
    i = 0;
    while ( (file.good()) && (length == i) )
      {
	
	getline(file, line_buffer);
	
	istringstream line(line_buffer);
	line.flags(istringstream::skipws);
    
	i = 0;
	while ( line.good() )
	  {
	    line >> point(i);
	    i++;
	  }
	
	if ( (length == 0) && (i != 0) )
	  length = i;
	
	if ( i == length )
	  {
	    points_.AddAtTheEnd(point);
	    n_++;
	  }
	
      }
    
  }


}  // namespace Multivac.


#define FILE_CURVE_CXX
#endif
