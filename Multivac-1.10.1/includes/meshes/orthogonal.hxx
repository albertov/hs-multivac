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


#ifndef FILE_MESHES_ORTHOGONAL_HXX


#include "../errors.cxx"
#include <cstdio>


namespace Multivac
{


  /////////////////////
  // CORTHOGONALMESH //
  /////////////////////

  //! Orthogonal mesh.
  template <class T>
  class COrthogonalMesh: public CMesh<T>
  {


    /**************
     * ATTRIBUTES *
     **************/

  protected:


    /*****************************
     * CONSTRUCTORS & DESTRUCTOR *
     *****************************/
    
  public:

    COrthogonalMesh()  throw();
    COrthogonalMesh(T Xmin, T Xmax, T Ymin, T Ymax,
		    int Nx, int Ny)  throw();

    ~COrthogonalMesh()  throw();


    /***********
     * METHODS *
     ***********/

  public:
  
    virtual void Save(string XFile, string YFile) const;
    virtual void SaveNonOrthogonalMesh(string PointsFile,
				       string EdgesFile,
				       string TrianglesFile) const;

  };  // COrthogonalMesh.


}  // namespace Multivac.


#define FILE_MESHES_ORTHOGONAL_HXX
#endif
