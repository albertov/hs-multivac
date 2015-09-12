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


#ifndef FILE_SAVER_CURVESSAVER_HXX


#include "../errors.cxx"
#include <cstdio>


namespace Multivac
{


  //////////////////
  // CCURVESSAVER //
  //////////////////

  //! This saver stores curves with a given period.
  /*!
    If the current iteration is a multiple of the given period,
    then this saver stores the curves (stored in the initializer and probably
    built from the current level set function). If the curve is not built yet,
    the initializer curve builder is called (member function
    'BuildCurveForDisplay'). In addition to that, the mesh is saved
    at the beginning. At the end, the speed function, the level set are saved
    (as they are at this time). Finally, iteration times (of the whole
    simulation) are saved.
  */
  template <class T>
  class CCurvesSaver: public CSaver<T>
  {


    /**************
     * ATTRIBUTES *
     **************/

  protected:


    /*****************************
     * CONSTRUCTORS & DESTRUCTOR *
     *****************************/

  public:

    CCurvesSaver()  throw();
    CCurvesSaver(string TimeFileName, string CurvesFileName,
		 string CurveLengthsFileName, string PhiFileName,
		 string FFileName, string XFileName, string YFileName,
		 string PointsFileName, string EdgesFileName,
		 string TrianglesFileName, int Period_);

    virtual ~CCurvesSaver()  throw();


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
			      Vector<T>& time, int iter,
			      CInitializer<T>& Initializer);

  };  // CCurvesSaver.


}  // namespace Multivac.


#define FILE_SAVER_CURVESSAVER_HXX
#endif
