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


#ifndef FILE_INITIALIZER_FASTMARCHINGNEVERINIT_HXX


#include "../errors.cxx"
#include <cstdio>


namespace Multivac
{


  ////////////////////////////
  // CFASTMARCHINGNEVERINIT //
  ////////////////////////////

  //! This initializer only provides first initializations.
  //! It doesn't provide any reinitialization or update.
  /*!
    The mesh is initialized.  The level set is initialized according
    to the initial curve position. Then, the speed function is
    initialized according to the initial level set function.
    The reinitialization only calls the speed function updater.
    \note
    This initializer is designed for the fast marching method.
  */
  template <class T>
  class CFastMarchingNeverInit: public CInitializer<T>
  {


    /**************
     * ATTRIBUTES *
     **************/

  protected:


    /*****************************
     * CONSTRUCTORS & DESTRUCTOR *
     *****************************/

  public:

    CFastMarchingNeverInit()  throw();

    ~CFastMarchingNeverInit()  throw();


    /***********
     * METHODS *
     ***********/
    
  public:
  
    virtual bool IsNarrowBand() const;
    virtual bool IsFastMarching() const;
  
    virtual void FirstInitMesh(CMesh<T>& Mesh) const;
    virtual void FirstInitInitialCurve(CMesh<T>& Mesh,
				       CInitialCurve<T>& InitialCurve) const;
    virtual void FirstInitPhiAndF(CMesh<T>& Mesh,
				  CInitialCurve<T>& InitialCurve,
				  CLevelSet<T>& Phi,
				  CSpeedFunction<T>& F,
				  CUpdater<T>& Updater);

    virtual void InitMesh(int iter, CMesh<T>& Mesh,
			  CLevelSet<T>& Phi, CSpeedFunction<T>& F,
			  CUpdater<T>& Updater, T CurrentTime) const;
    virtual void InitPhiAndF(int iter, CMesh<T>& Mesh,
			     CLevelSet<T>& Phi, CSpeedFunction<T>& F,
			     CUpdater<T>& Updater, T CurrentTime);

    virtual void BuildCurveForDisplay(int iter, CMesh<T>& Mesh,
				      CLevelSet<T>& Phi);

  };  // CCFastMarchingNeverInit.


}  // namespace Multivac.


#define FILE_INITIALIZER_FASTMARCHINGNEVERINIT_HXX
#endif
