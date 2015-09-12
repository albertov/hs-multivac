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


#ifndef FILE_LEVELSETS_ORTHOGONAL_CXX


#include "orthogonal.hxx"


namespace Multivac
{



  ///////////////////////////////
  // CONSTRUCTORS & DESTRUCTOR //
  ///////////////////////////////


  //! Default constructor.
  template <class T>
  COrthogonalLevelSet<T>::COrthogonalLevelSet()  throw()
  {

  }


  //! Destructor.
  template <class T>
  COrthogonalLevelSet<T>::~COrthogonalLevelSet()  throw()
  {

  }



  /////////////
  // METHODS //
  /////////////


  //! Inits the level set.
  /*! The matrix 'Values' is allocated.
    \param Mesh orthogonal mesh.
  */
  template <class T>
  void COrthogonalLevelSet<T>::Init(CMesh<T>& Mesh)
  {

    this->Values.Reallocate(Mesh.GetNx(), Mesh.GetNy());

  }


  //! Saves current level set function.
  /*!
    \param PhiFile files where values will be saved.
  */
  template <class T>
  void COrthogonalLevelSet<T>::Save(string PhiFile) const
  {

    // Saves the matrix 'Values' in PhiFile.
    this->Values.WriteText(PhiFile);

  }


}  // namespace Multivac.


#define FILE_LEVELSETS_ORTHOGONAL_CXX
#endif

