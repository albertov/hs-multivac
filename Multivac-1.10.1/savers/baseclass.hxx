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


#ifndef FILE_SAVER_BASECLASS_HXX


#include "../errors.cxx"
#include <cstdio>
#include <string>
#include <fstream>


namespace Multivac
{


  ////////////
  // CSAVER //
  ////////////

  //! Base class for savers (which save data).
  /*! Defines the savers interface.  All savers must be defined
    in the same way.
    \note
    This is an abstract class.
  */
  template <class T>
  class CSaver
  {

    /**************
     * ATTRIBUTES *
     **************/

  protected:
    //! Saves iteration times.
    string TimeFile;

    //! Saves curves (points abscissa and ordinate).
    //! All curves are saved in the same file.
    string CurvesFile;
    //! Saves curves lengths.
    string CurveLengthsFile;
    //! Saves level set function(s).
    string PhiFile;
    //! Saves speed function(s).
    string FFile;

    //! Saves grid abscissae.
    string XFile;
    //! Saves grid ordinates.
    string YFile;
    //! Saves mesh points.
    string PointsFile;
    //! Saves mesh edges.
    string EdgesFile;
    //! Saves mesh triangles.
    string TrianglesFile;

    //! The savers will save data if the current iteration is
    //! a multiple of 'Period'.
    int Period;
    //! Stores the last iteration when data were saved.
    int LastSaved;


    /*****************************
     * CONSTRUCTORS & DESTRUCTOR *
     *****************************/

  public:

    CSaver()  throw();
    CSaver(string TimeFileName,
	   string CurvesFileName, string CurveLengthsName,
	   string PhiFileName, string FFileName, string XFileName,
	   string YFileName, string PointsFileName, string EdgesFileName,
	   string TrianglesFileName, int Period_);

    virtual ~CSaver()  throw();


    /***********
     * METHODS *
     ***********/

  public:
  
    virtual void SaveAtTheBeginning(CMesh<T>& Mesh,
				    CSpeedFunction<T>& F,
				    CLevelSet<T>& Phi,
				    CInitializer<T>& Initializer) = 0;
    virtual void SaveAtCurrentIteration(CMesh<T>& Mesh,
					CSpeedFunction<T>& F,
					CLevelSet<T>& Phi,
					T time, int iter,
					CInitializer<T>& Initializer) = 0;
    virtual void SaveAtTheEnd(CMesh<T>& Mesh, CSpeedFunction<T>& F,
			      CLevelSet<T>& Phi, Vector<T>& time,
			      int iter, CInitializer<T>& Initializer) = 0;

  private:
    void CheckAndInit(string FileName);

  };  // CSaver.

}  // namespace Multivac.


#define FILE_SAVER_BASECLASS_HXX
#endif
