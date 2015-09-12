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


#ifndef FILE_CURVE_HXX

#include <cstdio>
#include <cmath>

#include <string>
#include <sstream>
#include <fstream>

using std::string;
using std::istringstream;
using std::ifstream;
using std::getline;


namespace Multivac
{


  //! Distance between two points (in 2D).
  /*!
    \param A coordinates of the first point (vector).
    \param B coordinates of the second point (vector).
    \return The distance between A and B.
  */
  template <class T>
  T DistanceBetween(Vector<T>& A, Vector<T>& B)
  {

    T tempx = A(0) - B(0);
    T tempy = A(1) - B(1);

    return sqrt( tempx*tempx + tempy*tempy );

  }


  //! Checks whether two segments intersect and have a different slope.
  /*! Checks whether ]AB] (or [AB[ if (AB) slope is less than (CD) slope)
    intersects ]CD].
    \param A origin of the first segment.
    \param B end of the first segment.
    \param C origin of the second segment.
    \param D end of the second segment.
    \return 'true' if ]AB] (or [AB[ if (AB) slope is less than (CD) slope)
    and ]CD] intersect and have a different slope, 'false' otherwise.
  */
  template <class T>
  bool Intersect(Vector<T>& A, Vector<T>& B,
		 Vector<T>& C, Vector<T>& D)
  {

    T slopeAB, slopeCD;

    slopeAB = B(0) - A(0);

    if (B(0) == D(0) && B(1) == D(1))
      return true;
    else if ( slopeAB == 0.0 )
      {
	slopeCD = D(0) - C(0);
	if ( slopeCD == 0.0 )
	  return false;
	slopeCD = (D(1) - C(1)) / slopeCD;
	T y = slopeCD * (A(0) - C(0)) + C(1);
	return ( ( (y > A(1) && y <= B(1)) || (y < A(1) && y >= B(1))
		   || (y == A(1) && y == B(1)) )
		 && ( (C(0) > A(0) && D(0) <= A(0))
		      || (C(0) < A(0) && D(0) >= A(0)) ) );
      }
    else
      {
	slopeAB = (B(1) - A(1)) / slopeAB;
	slopeCD = D(0) - C(0);
	if ( slopeCD == 0.0 )
	  {
	    T y = slopeAB * (C(0) - A(0)) + A(1);
	    return ( ( (y > C(1) && y <= D(1)) || (y < C(1) && y >= D(1))
		       || (y == C(1) && y == D(1)) )
		     && ( (D(0) > A(0) && D(0) <= B(0))
			  || (D(0) < A(0) && D(0) >= B(0)) ) );
	  }
	else
	  {
	    slopeCD = (D(1) - C(1)) / slopeCD;
	    if ( slopeAB - slopeCD == 0.0 )
	      return false;
	    T x = ( A(1) - slopeAB * A(0) - C(1) + slopeCD * C(0) )
	      / (slopeCD - slopeAB);
	    return ( ( (x > C(0) && x <= D(0)) || (x < C(0) && x >= D(0))
		       || (x == C(0) && x == D(0)) )
		     && ( (x > A(0) && x <= B(0))
			  || (x < A(0) && x >= B(0)) ) );
	  }
      }

  }



  /////////////////
  //  C U R V E  //
  /////////////////


  template <class T>
  class Curve
  {

    /************************
     * TYPEDEF DECLARATIONS *
     ************************/

  public:
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;


    /**************
     * ATTRIBUTES *
     **************/

  protected:
    //! Number of points.
    int n_;
    //! Points.
    List<Vector<T> > points_;
    //! Orientation (-1: unkonw orientation, 0: trigonometrical
    //! orientation, 1: reverse orientation).
    int orientation_;


    /***********
     * METHODS *
     ***********/

  public:
    // Constructors.
    Curve()  throw();
  
    // Destructor.
    ~Curve()  throw();

    // Basic methods.
    void InitWithUnsortedPoints_ClosestPointMethod(List<Vector<T> >& points);
    void AddPoint(Vector<T>& point);

    void Copy(Curve<T>& C);

    int GetNbPoints() const;
    List<Vector<T> >& GetPoints();

    int GetOrientation() const;
    int SetOrientation(int orientation);

    void ClearAll();

    T GetDistanceToTheCurve(T x, T y);
    T GetProjection(Vector<T>& X);
    T GetSignedDistanceToTheCurve(T x, T y);

    // Input/output methods.
    void WriteText(ofstream& f);
    void ReadText(string file_name, int orientation = -1);

  };  // Curve.


}  // namespace Multivac.


#define FILE_CURVE_HXX
#endif
