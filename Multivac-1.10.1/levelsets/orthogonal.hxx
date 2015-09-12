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


#ifndef FILE_LEVELSETS_ORTHOGONAL_HXX


#include "../errors.cxx"
#include <cstdio>


namespace Multivac
{


  /////////////////////////
  // CORTHOGONALLEVELSET //
  /////////////////////////

  //! Level set function defined on an orthogonal mesh.
  template <class T>
  class COrthogonalLevelSet: public CLevelSet<T>
  {


    /**************
     * ATTRIBUTES *
     **************/

  protected:


    /*****************************
     * CONSTRUCTORS & DESTRUCTOR *
     *****************************/

  public:

    COrthogonalLevelSet()  throw();

    ~COrthogonalLevelSet()  throw();


    /***********
     * METHODS *
     ***********/

  public:

    virtual void Init(CMesh<T>& Mesh);

    virtual void Save(string PhiFile) const;

  };  // COrthogonal.


}  //  namespace Multivac.


#define FILE_LEVELSETS_ORTHOGONAL_HXX
#endif
