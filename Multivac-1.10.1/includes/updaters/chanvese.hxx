// Copyright (C) 2006 Vivien Mallet
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


#ifndef FILE_UPDATER_CHANVESE_HXX


#include "../errors.cxx"
#include <cstdio>


namespace Multivac
{


  ///////////////
  // CCHANVESE //
  ///////////////

  //! This updater implements the Chan-Vese algorithm for image segmentation.
  /*! \warning The speed function to be used in combination with this updater
    is supposed to provide the image intensity instead of the front velocity.
  */
  template <class T>
  class CChanVese: public CUpdater<T>
  {


    /**************
     * ATTRIBUTES *
     **************/

  protected:
    
    //! inside_ weight of the inside integral.
    T inside;
    //! outside_ weight of the outside integral.
    T outside;
    //! mu_ weight of the penalty term on the length.
    T mu;
    //! nu_ additional velocity.
    T nu;

    /*! \brief Dirac_threshold_ scale factor applied to the smoothed Dirac
      function. It should be above 1. */
    T Dirac_threshold;
    T pi_threshold;
    T two_threshold;


    /*****************************
     * CONSTRUCTORS & DESTRUCTOR *
     *****************************/

  public:

    CChanVese()  throw();
    CChanVese(int TubeSemiWidth_,
	      int BarrierWidth_,
	      int OutSpaceWidth_,
	      T inside_, T outside_,
	      T mu_, T nu_, T Dirac_threshold_)  throw();

    ~CChanVese()  throw();


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

    virtual T Dirac(T phi);

  };  // CChanVese.


}  // namespace Multivac.


#define FILE_UPDATER_CHANVESE_HXX
#endif
