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


#ifndef FILE_INITIALIZER_BASECLASS_CXX


#include "baseclass.hxx"
#include <cstdio>


namespace Multivac
{



  ///////////////////////////////
  // CONSTRUCTORS & DESTRUCTOR //
  ///////////////////////////////


  //! Default constructor.
  template <class T>
  CInitializer<T>::CInitializer()  throw():
    LastCurveUpdate(-1),
    LastCurveUpdateForDisplay(-1)
  {

  }


  //! Destructor.
  template <class T>
  CInitializer<T>::~CInitializer()  throw()
  {

  }



  /////////////
  // METHODS //
  /////////////


  //! Returns the current stored front.
  template <class T>
  Curve<T>& CInitializer<T>::GetFront()
  {

    return Front;

  }


  //! Saves current stored front.
  /*!
    \param CurvesFile file in which the front will be stored.
    \param CurveLengthsFile file in which the number of points (on the front)
    will be stored.
    \note
    The curve contains points of the front.
    Those points are saved in an ASCII file, in
    two columns:
    \par 0.31500 1.24779
    \par 1.88779 1.22890
    \par ...
    \warning The saved curve is the curve currently stored by the initializer.
    No update is computed in this routine.
  */
  template <class T>
  void CInitializer<T>::Save(string CurvesFile, string CurveLengthsFile)
  {

    ofstream CurvesStream(CurvesFile.c_str(), ofstream::app);
    // Saves the current curves which contains points of the front.
    Front.WriteText(CurvesStream);
    CurvesStream.close();

    // Saves the number of points.
    // The use of the Seldon::Vector<int> N is convenient.
    Vector<int> N(1);
    N(0) = Front.GetNbPoints();
    ofstream CurveLengthsStream(CurveLengthsFile.c_str(), ofstream::app);
    N.WriteText(CurveLengthsStream);
    CurveLengthsStream << "\n";
    CurveLengthsStream.close();

  }

  

}  // namespace Multivac.


#define FILE_INITIALIZER_BASECLASS_CXX
#endif
