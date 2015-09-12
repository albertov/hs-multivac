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


#ifndef FILE_UPDATER_NARROWBANDFIRSTORDERENGQUISTOSHER_HXX


#include "../errors.cxx"
#include <cstdio>


namespace Multivac
{


  ////////////////////////////////////////
  // CNARROWBANDFIRSTORDERENGQUISTOSHER //
  ////////////////////////////////////////

  //! This updater uses the first order Engquist-Osher scheme.
  //! The time integration is performed by an explicit Euler scheme.
  /*! \note This updater is designed for the narrow band level set method.
   */
  template <class T>
  class CNarrowBandFirstOrderEngquistOsher: public CUpdater<T>
  {


    /**************
     * ATTRIBUTES *
     **************/

  protected:


    /*****************************
     * CONSTRUCTORS & DESTRUCTOR *
     *****************************/

  public:

    CNarrowBandFirstOrderEngquistOsher()  throw();
    CNarrowBandFirstOrderEngquistOsher(int TubeSemiWidth_,
				       int BarrierWidth_,
				       int OutSpaceWidth_)  throw();

    ~CNarrowBandFirstOrderEngquistOsher()  throw();


    /***********
     * METHODS *
     ***********/
    
  public:
  
    virtual bool IsNarrowBand() const;
    virtual bool IsFastMarching() const;

    virtual void Init(CMesh<T>& Mesh, CLevelSet<T>& Phi);
    virtual void UpdateLevelSet(T Delta_t,
				CMesh<T>& Mesh,
				CSpeedFunction<T>& F,
				CLevelSet<T>& Phi,
				T CurrentTime);

  };  // CNarrowBandFirstOrderEngquistOsher.


}  // namespace Multivac.


#define FILE_UPDATER_NARROWBANDFIRSTORDERENGQUISTOSHER_HXX
#endif
