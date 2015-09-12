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


#ifndef FILE_SAVER_SAVEATTHEEND_HXX


#include "../errors.cxx"
#include <cstdio>


namespace Multivac
{


  //////////////////
  // CCURVESSAVER //
  //////////////////

  //! This saver stores curves the last mesh, the last speed function
  //! and the last level set function.
  /*!
    Nothing is saved until the end of the last iteration.  Then, the mesh,
    the speed function and the level set function are saved.
  */
  template <class T>
  class CSaveAtTheEnd: public CSaver<T>
  {


    /**************
     * ATTRIBUTES *
     **************/

  protected:


    /*****************************
     * CONSTRUCTORS & DESTRUCTOR *
     *****************************/

  public:

    CSaveAtTheEnd()  throw();
    CSaveAtTheEnd(string TimeFileName, string CurvesFileName,
		  string CurveLengthsFileName, string PhiFileName,
		  string FFileName, string XFileName, string YFileName,
		  string PointsFileName, string EdgesFileName,
		  string TrianglesFileName, int Period_);

    ~CSaveAtTheEnd()  throw();


    /***********
     * METHODS *
     ***********/
    
  public:
  
    virtual void SaveAtTheBeginning(CMesh<T>& Mesh,
				    CSpeedFunction<T>& F,
				    CLevelSet<T>& Phi,
				    CInitializer<T>& Initializer);
    virtual void SaveAtCurrentIteration(CMesh<T>& Mesh,
					CSpeedFunction<T>& F,
					CLevelSet<T>& Phi,
					T time, int iter,
					CInitializer<T>& Initializer);
    virtual void SaveAtTheEnd(CMesh<T>& Mesh,
			      CSpeedFunction<T>& F,
			      CLevelSet<T>& Phi,
			      Vector<T>& time,
			      int iter,
			      CInitializer<T>& Initializer);

  };  // CCurvesSaver.


}  // namespace Multivac.


#define FILE_SAVER_SAVEATTHEEND_HXX
#endif
