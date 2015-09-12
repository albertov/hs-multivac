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


#ifndef FILE_UPDATER_BASECLASS_HXX


#include "../errors.cxx"
#include <cstdio>


namespace Multivac
{


  ////////////////////////
  // CLASS DECLARATIONS //
  ////////////////////////

  template <class T>
  class CFastMarchingNode;


  //////////////
  // CUPDATER //
  //////////////

  //! Base class for updaters which updates the level set function
  //! at each iteration.
  /*! Defines the updater interface.  All initial curves must be defined
    in the same way.  An updater provides a numerical scheme (in space)
    and an integrator (in time).
    \note
    This is an abstract class.
  */
  template <class T>
  class CUpdater
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

    /**** For all methods ****/

    //! Temporary values may be stored in 'Temp' in order to perform
    //! computations.
    Matrix<T> Temp;

    //! Number of ghost cells.
    int offset;

    //! 'NeedSpeedUpdateFlag' indicates whether speed rates should be
    //! updated (by calling the speed class updater).
    bool NeedSpeedUpdateFlag;

    //! 'NeedInitializationFlag' indicates whether reinitialization should
    //! be performed.
    bool NeedInitializationFlag;

    /**** For the full matrix level set method ****/


    /**** For the narrow band level set method ****/

    //! Tube points are stored in 'Tube'.  It defines the tube that is used
    //! in the narrow band level set method.
    Vector<List<Vector<int> >, Vect_Full,
	   NewAlloc<List<Vector<int> > > > Tube;
    //! Tube semi width: number of cells on each side of the front.
    int TubeSemiWidth;

    //! Stores points of the barrier.  In the narrow band level set method, if
    //! a 'Barrier' point is reached by the front, then the tube has to be
    //! rebuilt.
    Matrix<int> Barrier;
    //! Barrier width: number of cells along the barrier width.
    int BarrierWidth;

    //! Stores points close to the tube edges.  In the narrow band level set
    //! method, if a 'OutSpace' point is reached by the front,
    //! then an unstability
    //! occurred.  'OutSpace' points should not be reached by the front.
    Matrix<int> OutSpace;
    //! "Outspace" width: number of cells along the "outspace" width.
    int OutSpaceWidth;

    //! In the fast marching method, 'TrialPoints' stores points whose arrival
    //! time are not known yet.  But, temporary arrival times have been
    //! computed for those points yet and are stored in 'TrialPoints', with
    //! point coordinates.
    ArrayHeap<CFastMarchingNode<T> > TrialPoints;
    //! 'PointersToNodes' element (i, j) points to the corresponding node in
    //! 'TrialPoints' (if it exists).
    Matrix<int> PointersToNodes;
    //! For the fast marching method, 'Tmax' is a majorant of arrival times.
    T TMax;


    /*****************************
     * CONSTRUCTORS & DESTRUCTOR *
     *****************************/

  public:

    CUpdater()  throw();

    virtual ~CUpdater()  throw();


    /***********
     * METHODS *
     ***********/

  public:
  
    /**** For all methods ****/

    virtual bool IsNarrowBand() const = 0;
    virtual bool IsFastMarching() const = 0;

    virtual void Init(CMesh<T>& Mesh, CLevelSet<T>& Phi) = 0;
    virtual void UpdateLevelSet(T Delta_t, CMesh<T>& Mesh,
				CSpeedFunction<T>& F,
				CLevelSet<T>& Phi,
				T CurrentTime) = 0;

    virtual bool NeedSpeedUpdate() const;
    virtual bool NeedInitialization() const;

    virtual Matrix<T>& GetTemp();

    int GetOffset();

    /**** For full matrix level set methods ****/


    /**** For narrow band level set methods ****/

    Vector<List<Vector<int> >, Vect_Full,
	   NewAlloc<List<Vector<int> > > >& GetTube();
    int GetTubeSemiWidth() const;

    Matrix<int>& GetBarrier();
    int GetBarrierWidth() const;
    Matrix<int>& GetOutSpace();
    int GetOutSpaceWidth() const;
    
    /**** For fast marching methods ****/

    virtual bool KeepOnWorking() const;

    ArrayHeap<CFastMarchingNode<T> >& GetTrialPoints();
    Matrix<int>& GetPointersToNodes();
    
    T GetTMax();

  };  // CUpdater.


}  // namespace Multivac.


#define FILE_UPDATER_BASECLASS_HXX
#endif
