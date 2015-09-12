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


#ifndef FILE_SPEEDFUNCTIONS_BASECLASS_CXX


#include "baseclass.hxx"
#include <cstdio>


namespace Multivac
{



  ///////////////////////////////
  // CONSTRUCTORS & DESTRUCTOR //
  ///////////////////////////////


  //! Default constructor.
  template <class T>
  CSpeedFunction<T>::CSpeedFunction()  throw()
  {

    dependence_position = false;
    dependence_time = false;
    dependence_normal = false;
    dependence_curvature = false;

  }


  //! Destructor.
  template <class T>
  CSpeedFunction<T>::~CSpeedFunction()  throw()
  {

  }



  /////////////
  // METHODS //
  /////////////



  //! Does the speed function depend upon the position?
  /*!
    \return 'true' if the velocity depends upon the position,
    'false' otherwise.
  */
  template <class T>
  inline bool CSpeedFunction<T>::IsPositionDependent()  const
  {

    return dependence_position;

  }


  //! Does the speed function depend upon the time?
  /*!
    \return 'true' if the velocity depends upon the time,
    'false' otherwise.
  */
  template <class T>
  inline bool CSpeedFunction<T>::IsTimeDependent()  const
  {

    return dependence_time;

  }


  //! Does the speed function depend upon the normal?
  /*!
    \return 'true' if the velocity depends upon the normal,
    'false' otherwise.
  */
  template <class T>
  inline bool CSpeedFunction<T>::IsNormalDependent()  const
  {

    return dependence_normal;

  }


  //! Does the speed function depend upon the curvature?
  /*!
    \return 'true' if the velocity depends upon the curvature,
    'false' otherwise.
  */
  template <class T>
  inline bool CSpeedFunction<T>::IsCurvatureDependent()  const
  {

    return dependence_curvature;

  }


  //! Returns the matrix that stores speed rates on grid points.
  /*! The matrix is returned by reference.
    \return a reference to the matrix that stores speed rates on grid points.
  */
  template <class T>
  inline Matrix<T>& CSpeedFunction<T>::GetValues()
  {

    return Values;

  }


  //! Returns the speed rate at a given grid point.
  /*!
    \param i row index.
    \param j column index.
    \return The speed rate at (i, j).
  */
  template <class T>
  inline T CSpeedFunction<T>::operator() (int i, int j) const
  {

    return Values(i, j);

  }


  //! Saves current speed rates.
  /*!
    \param FFile files where rates will be saved.
  */
  template <class T>
  void CSpeedFunction<T>::Save(string FFile) const
  {

    // Saves the matrix 'Values' in FFile.
    Values.WriteText(FFile);

  }


}  // namespace Multivac.


#define FILE_SPEEDFUNCTIONS_BASECLASS_CXX
#endif
