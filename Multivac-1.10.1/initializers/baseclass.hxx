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


#ifndef FILE_INITIALIZER_BASECLASS_HXX


#include "../errors.cxx"
#include <cstdio>


namespace Multivac
{


  //////////////////
  // CINITIALIZER //
  //////////////////

  //! Base class for initializers.
  /*! Defines the initializers interface.  All initializers must be defined
    in the same way.
    \note
    This is an abstract class.
  */
  template <class T>
  class CInitializer
  {


    /**************
     * ATTRIBUTES *
     **************/

  protected:

    //! The front may be saved in this curve.
    Curve<T> Front;

    //! Stores the last iteration when the front was built
    //! on updating purpose.
    int LastCurveUpdate;
    //! Stores the last iteration when the front was built
    //! on displaying purpose.
    int LastCurveUpdateForDisplay;


    /*****************************
     * CONSTRUCTORS & DESTRUCTOR *
     *****************************/

  public:

    CInitializer()  throw();

    virtual ~CInitializer()  throw();


    /***********
     * METHODS *
     ***********/

  public:

    virtual bool IsNarrowBand() const = 0;
    virtual bool IsFastMarching() const = 0;
  
    virtual void FirstInitMesh(CMesh<T>& Mesh) const = 0;
    virtual void FirstInitInitialCurve(CMesh<T>& Mesh,
				       CInitialCurve<T>& InitialCurve) const
      = 0;
    virtual void FirstInitPhiAndF(CMesh<T>& Mesh,
				  CInitialCurve<T>& InitialCurve,
				  CLevelSet<T>& Phi, CSpeedFunction<T>& F,
				  CUpdater<T>& Updater) = 0;

    virtual void InitMesh(int iter, CMesh<T>& Mesh, CLevelSet<T>& Phi,
			  CSpeedFunction<T>& F, CUpdater<T>& Updater,
			  T CurrentTime) const = 0;
    virtual void InitPhiAndF(int iter, CMesh<T>& Mesh, CLevelSet<T>& Phi,
			     CSpeedFunction<T>& F, CUpdater<T>& Updater,
			     T CurrentTime) = 0;

    virtual void BuildCurveForDisplay(int iter, CMesh<T>& Mesh,
				      CLevelSet<T>& Phi) = 0;

    Curve<T>& GetFront();

    virtual void Save(string CurvesFile, string CurveLengthsFile);

  };  // CInitializer.


}  //  namespace Multivac.


#define FILE_INITIALIZER_BASECLASS_HXX
#endif
