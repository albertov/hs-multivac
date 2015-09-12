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


#ifndef FILE_UPDATER_FASTMARCHINGFIRSTORDERENGQUISTOSHER_HXX


#include "../errors.cxx"
#include <cstdio>


namespace Multivac
{


  //////////////////////////////////////////
  // CFASTMARCHINGFIRSTORDERENGQUISTOSHER //
  //////////////////////////////////////////

  //! This updater uses the first order Engquist Osher scheme.
  /*! \note This updater is designed for the fast marching method.
   */
  template <class T>
  class CFastMarchingFirstOrderEngquistOsher: public CUpdater<T>
  {


    /***********************
     * TYPEDEF DECLARATION *
     ***********************/

  public:

    typedef CFastMarchingNode<T> heap_node_value_type;
    typedef CFastMarchingNode<T>& heap_node_reference;
    typedef const CFastMarchingNode<T>& heap_node_const_reference;
    typedef CFastMarchingNode<T>* heap_node_pointer;
    typedef const CFastMarchingNode<T>* heap_node_const_pointer;


    /**************
     * ATTRIBUTES *
     **************/

  protected:


    /*****************************
     * CONSTRUCTORS & DESTRUCTOR *
     *****************************/

  public:

    CFastMarchingFirstOrderEngquistOsher()  throw();
    CFastMarchingFirstOrderEngquistOsher(T TMax_)  throw();

    ~CFastMarchingFirstOrderEngquistOsher()  throw();


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

    virtual bool KeepOnWorking() const;

  private:
    
    T  RefreshTime(CMesh<T>& Mesh, CSpeedFunction<T>& F,
		   CLevelSet<T>& Phi, int i, int j);

  };  // CFastMarchingFirstOrderEngquistOsher.


}  // namespace Multivac.


#define FILE_UPDATER_FASTMARCHINGFIRSTORDERENGQUISTOSHER_HXX
#endif
