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


#ifndef FILE_SAVER_BASECLASS_CXX


#include "baseclass.hxx"
#include <cstdio>


namespace Multivac
{



  ///////////////////////////////
  // CONSTRUCTORS & DESTRUCTOR //
  ///////////////////////////////


  //! Default constructor.
  template <class T>
  CSaver<T>::CSaver()  throw():
    TimeFile(NULL), CurvesFile(NULL),
    CurveLengthsFile(NULL), PhiFile(NULL),
    FFile(NULL), XFile(NULL), YFile(NULL),
    PointsFile(NULL), EdgesFile(NULL),
    TrianglesFile(NULL), Period(-1), LastSaved(-1)
  {

  }


  //! Main contructor.
  /*!
    \param TimeFileName file in which iteration times will be stored.
    \param CurvesFileName file in which fronts will be stored.
    \param CurveLengthsFileName file in which numbers of points (on fronts)
    will be stored.
    \param PhiFileName file in which the level set function will be stored.
    \param FFileName file in which the speed function will be stored.
    \param XFileName file in which grid point abscissae will be stored.
    \param YFileName file in which grid point ordinates will be stored.
    \param PointsFileName file in which mesh points will be stored.
    \param EdgesFileName file in which mesh edges will be stored.
    \param TrianglesFileName file in which mesh triangles will be stored.
    \param Period_ Data will be saved if the current iteration is
    a multiple of 'Period_'.
  */
  template <class T>
  CSaver<T>::
  CSaver(string TimeFileName,
	 string CurvesFileName, string CurveLengthsFileName,
	 string PhiFileName, string FFileName, string XFileName,
	 string YFileName, string PointsFileName, string EdgesFileName,
	 string TrianglesFileName, int Period_):
    LastSaved(-1)
  {

    if (Period_<=0)
      Period = 1;
    else
      Period = Period_;
    
    TimeFile = TimeFileName;
    CheckAndInit(TimeFile);

    CurvesFile = CurvesFileName;
    CheckAndInit(CurvesFile);

    CurveLengthsFile = CurveLengthsFileName;
    CheckAndInit(CurveLengthsFile);

    PhiFile = PhiFileName;
    CheckAndInit(PhiFile);

    FFile = FFileName;
    CheckAndInit(FFile);

    XFile = XFileName;
    CheckAndInit(XFile);

    YFile = YFileName;
    CheckAndInit(YFile);

    PointsFile = PointsFileName;
    CheckAndInit(PointsFile);

    EdgesFile = EdgesFileName;
    CheckAndInit(EdgesFile);

    TrianglesFile = TrianglesFileName;
    CheckAndInit(TrianglesFile);

  }


  //! Destructor.
  template <class T>
  CSaver<T>::~CSaver()  throw()
  {

  }


  //! Checks and inits a file.
  /*! Checks whether a file is ready for output operations
    and clears it (i.e. it is empty on exit).
    \param FileName file name.
  */
  template <class T>
  void CSaver<T>::CheckAndInit(string FileName)
  {
    ofstream FileStream(FileName.c_str(), ofstream::out);
#ifdef MULTIVAC_CHECK_IO
    //! Checks whether the file is ready.
    if (!FileStream.good())
      {
	string comment = "Unable to write in file \"" + FileName + "\"";
	throw Multivac::CError_FileIO("CSaver::CheckAndInit(string)",
				      comment.c_str());
      }
#endif
    FileStream.close();
  }


}  // namespace Multivac.


#define FILE_SAVER_BASECLASS_CXX
#endif
