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


#ifndef FILE_INITIALIZER_NARROWBANDNEVERINIT_HXX


#include "../errors.cxx"
#include <cstdio>


namespace Multivac
{


  ////////////////
  // CNEVERINIT //
  ////////////////

  //! This initializer can be used to perform initializations and
  //! reinitializations of the level set function and the tube
  //! involved in the narrow band level set method.
  /*!
    First, the mesh is initialized.  The level set is initialized according
    to the initial curve position.  Then, the speed function is
    initialized according to the initial level set function.
    Reinitializations build the tube and init the corresponding level set
    function.  As for the speed function, reinitializations only call the
    speed function updater.
    \note
    This initializer is designed for the narrow band
    level set method and for an orthogonal mesh.
  */
  template <class T>
  class CNarrowBandNeverInit: public CInitializer<T>
  {


    /**************
     * ATTRIBUTES *
     **************/

  protected:


    /*****************************
     * CONSTRUCTORS & DESTRUCTOR *
     *****************************/

  public:

    CNarrowBandNeverInit()  throw();

    ~CNarrowBandNeverInit()  throw();


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
    virtual void InitPhi(int iter, CMesh<T>& Mesh,
			 CLevelSet<T>& Phi, CUpdater<T>& Updater);
    virtual void InitF(CMesh<T>& Mesh, Vector<List<Vector<int> >, Vect_Full,
		       NewAlloc<List<Vector<int> > > >& Tube,
		       CLevelSet<T>& Phi, T CurrentTime,
		       CSpeedFunction<T>& F);

    virtual void BuildCurveForDisplay(int iter, CMesh<T>& Mesh,
				      CLevelSet<T>& Phi);

  };  // CCNarrowBandNeverInit.


}  // namespace Multivac.


#define FILE_INITIALIZER_NARROWBANDNEVERINIT_HXX
#endif
