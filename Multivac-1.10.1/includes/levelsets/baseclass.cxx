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


#ifndef FILE_LEVELSETS_BASECLASS_CXX


#include "baseclass.hxx"
#include <cstdio>


namespace Multivac
{



  ///////////////////////////////
  // CONSTRUCTORS & DESTRUCTOR //
  ///////////////////////////////


  //! Default constructor.
  template <class T>
  CLevelSet<T>::CLevelSet()  throw()
  {

  }


  //! Destructor.
  template <class T>
  CLevelSet<T>::~CLevelSet()  throw()
  {

  }



  /////////////
  // METHODS //
  /////////////


  //! Reallocates memory to store level set function values
  //! on an orthogonal mesh.
  /*!
    \param i number of grid point along (x'x).
    \param j number of grid point along (y'y).
    \exception Seldon::NoMemory there is not enough available memory.
  */
  template <class T>
  inline void CLevelSet<T>::Reallocate(int i, int j)
  {

    Values.Reallocate(i, j);

  }


  //! Returns the matrix that stores level set values on grid points.
  /*! The matrix is returned by reference.
    \return a reference to the matrix that stores level set values
    on grid points.
  */
  template <class T>
  inline Matrix<T>& CLevelSet<T>::GetValues()
  {

    return Values;

  }


  //! Access operator to get level set function values
  //! on grid points (orthogonal mesh).
  /*!
    \param i grid point index (along (x'x)).
    \param j grid point index (along (y'y)).
    \exception Seldon::WrongCol attempt to reach a wrong column number.
    \exception Seldon::WrongRow attempt to reach a wrong row number.
    \return Level set function value on (i, j).
  */
  template <class T>
  inline T& CLevelSet<T>::operator() (int i, int j)
  {

    return Values(i, j);

  }


}  // namespace Multivac.


#define FILE_LEVELSETS_BASECLASS_CXX
#endif
