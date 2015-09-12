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


#ifndef FILE_LEVELSETS_BASECLASS_HXX


#include "../errors.cxx"
#include <cstdio>


namespace Multivac
{


  ///////////////
  // CLEVELSET //
  ///////////////

  //! Base class for level set functions.
  /*! Defines the level set interface.  All level set functions
    must be defined in the same way.
    \note
    This is an abstract class.
  */
  template <class T>
  class CLevelSet
  {


    /**************
     * ATTRIBUTES *
     **************/

  protected:

    //! Level set values on grid points (for an orthogonal mesh).
    Matrix<T> Values;


    /*****************************
     * CONSTRUCTORS & DESTRUCTOR *
     *****************************/

  public:

    CLevelSet()  throw();

    virtual ~CLevelSet()  throw();


    /***********
     * METHODS *
     ***********/

  public:
  
    virtual void Init(CMesh<T>& Mesh) = 0;
    virtual void Reallocate(int i, int j);
    
    Matrix<T>& GetValues();
    virtual T& operator() (int i, int j);

    virtual void Save(string PhiFile) const = 0;

  };  // CLevelSet.


}  //  namespace Multivac.


#define FILE_LEVELSETS_BASECLASS_HXX
#endif
