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


#ifndef FILE_UPDATER_BASECLASS_CXX


#include "baseclass.hxx"
#include <cstdio>


namespace Multivac
{



  ///////////////////////////////
  // CONSTRUCTORS & DESTRUCTOR //
  ///////////////////////////////


  //! Default constructor.
  template <class T>
  CUpdater<T>::CUpdater()  throw()
  {

  }


  //! Destructor.
  template <class T>
  CUpdater<T>::~CUpdater()  throw()
  {

  }



  /////////////
  // METHODS //
  /////////////

  
  //! Should speed rates be updated?
  /*! If 'NeedSpeedUpdate' returns 'true', then the speed class updater should
    be called so as to update the speed function.
    \return A boolean that indicates whether or not speed rates should be
    updated.
  */
  template <class T>
  inline bool CUpdater<T>::NeedSpeedUpdate() const
  {

    return NeedSpeedUpdateFlag;

  }


  //! Should any initialization be performed?
  /*! If 'NeedInitialization' returns 'true', then the initializer
    should reinit the method (level set function, speed function, mesh, ...).
    \return A boolean that indicates whether or not calculations need
    a reinitialization.
  */
  template <class T>
  inline bool CUpdater<T>::NeedInitialization() const
  {

    return NeedInitializationFlag;

  }


  //! Returns matrix 'Temp' that stores temporary values.
  /*! \return A reference to matrix 'Temp' which stored temporary values.
   */
  template <class T>
  inline Matrix<T>& CUpdater<T>::GetTemp()
  {

    return Temp;

  }


  //! Returns the number of ghost cells.
  /*! Returns the number of ghost layers to be added around the.domain.
    \return The number of ghost cells.
  */
  template <class T>
  inline int CUpdater<T>::GetOffset()
  {

    return offset;

  }


  //! Returns the tube.
  /*! \return A reference to the vector of lists of vectors that stores
    tube points.
  */
  template <class T>
  inline Vector<List<Vector<int> >, Vect_Full,
		NewAlloc<List<Vector<int> > > >& CUpdater<T>::GetTube()
  {

    return Tube;

  }


  //! Returns the tube semi width.
  /*! \return Number of cells on each side of the front.
   */
  template <class T>
  inline int CUpdater<T>::GetTubeSemiWidth() const
  {

    return TubeSemiWidth;

  }


  //! Returns barrier points.
  /*! If one of 'Barrier' points is reached by the front, then the tube
    has to be rebuilt.
    \return A reference to the matrix 'Barrier' which contains barrier points.
  */
  template <class T>
  inline Matrix<int>& CUpdater<T>::GetBarrier()
  {

    return Barrier;

  }


  //! Returns the barrier width.
  /*!  If one of 'Barrier' points is reached by the front, then the tube has
    to be rebuilt.
    \return Number of cells along the barrier width.
  */
  template <class T>
  inline int CUpdater<T>::GetBarrierWidth() const
  {

    return BarrierWidth;

  }


  //! Returns "outspace" points.
  /*! The "outspace" is the set of points that should never be reached by the
    front.  If the front reachs one of these points, an instability occurred.
    \return A reference to the matrix 'OutSpace' which contains "outspace"
    points.
  */
  template <class T>
  inline Matrix<int>& CUpdater<T>::GetOutSpace()
  {

    return OutSpace;

  }


  //! Returns the "outspace" width.
  /*! The "outspace" is the set of points that should never be reached by the
    front.  If the front reachs one of these points, an instability occurred.
    \return Number of cells along the "outspace" width.
  */
  template <class T>
  inline int CUpdater<T>::GetOutSpaceWidth() const
  {

    return OutSpaceWidth;

  }


  //! Should computations be continued?
  /*! If 'KeepOnWorking' returns 'true', then computations should be
    continued.
    \return A boolean that indicates whether or not the program should keep
    on computing.
  */
  template <class T>
  inline bool CUpdater<T>::KeepOnWorking() const
  {

    return true;

  }


  //! Returns the heap of trial points.
  /*! \return A reference to the heap (binary tree) 'TrialPoints' which
    stores trial points and their arrival times.
  */
  template <class T>
  inline ArrayHeap<CFastMarchingNode<T> >& CUpdater<T>::GetTrialPoints()
  {

    return TrialPoints;

  }
  

  //! Returns pointers to nodes of 'TrialPoints'.
  /*! \return A reference to the matrix 'PointersToNodes' which stores
    pointers to nodes of 'TrialPoints'.
  */
  template <class T>
  inline Matrix<int>& CUpdater<T>::GetPointersToNodes()
  {

    return PointersToNodes;

  }
 

  //! Returns the given majorant of arrival times.
  /*! This majorant was set when the updater was constructed.  This majorant
    was given by the user.
    \return the given majorant of arrival times.
  */
  template <class T>
  inline T CUpdater<T>::GetTMax()
  {

    return TMax;

  }
 

}  // namespace Multivac.


#define FILE_UPDATER_BASECLASS_CXX
#endif
