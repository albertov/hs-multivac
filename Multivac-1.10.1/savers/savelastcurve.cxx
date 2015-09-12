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


#ifndef FILE_SAVER_SAVELASTCURVE_CXX


#include "savelastcurve.hxx"


namespace Multivac
{



  ///////////////////////////////
  // CONSTRUCTORS & DESTRUCTOR //
  ///////////////////////////////


  //! Default constructor.
  template <class T>
  CSaveLastCurve<T>::CSaveLastCurve()  throw()
  {

  }


  //! Main contructor.
  /*! This constructor opens all files.
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
  CSaveLastCurve<T>::
  CSaveLastCurve(string TimeFileName, string CurvesFileName,
		 string CurveLengthsFileName, string PhiFileName,
		 string FFileName, string XFileName, string YFileName,
		 string PointsFileName, string EdgesFileName,
		 string TrianglesFileName, int Period_):
    CSaver<T>(TimeFileName, CurvesFileName,
	      CurveLengthsFileName, PhiFileName, FFileName, XFileName,
	      YFileName, PointsFileName, EdgesFileName,
	      TrianglesFileName, Period_)
  {

  }


  //! Destructor.
  template <class T>
  CSaveLastCurve<T>::~CSaveLastCurve()  throw()
  {

  }



  /////////////
  // METHODS //
  /////////////


  //! Empty function.
  /*! \param Mesh initial mesh.
    \param F initial speed function.
    \param Phi initial level set.
    \param Initializer initializer.
  */
  template <class T>
  void CSaveLastCurve<T>::SaveAtTheBeginning(CMesh<T>& Mesh,
					     CSpeedFunction<T>& F,
					     CLevelSet<T>& Phi,
					     CInitializer<T>& Initializer)
  {

  }


  //! Empty function.
  /*! \param Mesh initial mesh.
    \param F initial speed function.
    \param Phi initial level set.
    \param time current time.
    \param iter current iteration.
    \param Initializer initializer.
  */
  template <class T>
  void CSaveLastCurve<T>::SaveAtCurrentIteration(CMesh<T>& Mesh,
						 CSpeedFunction<T>& F,
						 CLevelSet<T>& Phi,
						 T time, int iter,
						 CInitializer<T>& Initializer)
  {

  }


  //! Saves the last mesh, the last speed function, the last
  //! level set function and the last curve.
  /*! \param Mesh initial mesh.
    \param F initial speed function.
    \param Phi initial level set.
    \param Time vector that constains iteration times
    (of the whole simulation).
    \param iter current iteration.
    \param Initializer initializer.
    \exception Seldon::IOError an error occured during I/O operations.
  */
  template <class T>
  void CSaveLastCurve<T>::SaveAtTheEnd(CMesh<T>& Mesh,
				       CSpeedFunction<T>& F,
				       CLevelSet<T>& Phi,
				       Vector<T>& Time,
				       int iter,
				       CInitializer<T>& Initializer)
  {

    Mesh.Save(this->XFile, this->YFile);
    F.Save(this->FFile);
    Phi.Save(this->PhiFile);
    Time.WriteText(this->TimeFile);

    // The curve may need to be built before being saved.
    Initializer.BuildCurveForDisplay(iter, Mesh, Phi);
    Initializer.Save(this->CurvesFile, this->CurveLengthsFile);

  }



}  // namespace Multivac.


#define FILE_SAVER_SAVELASTCURVE_CXX
#endif
