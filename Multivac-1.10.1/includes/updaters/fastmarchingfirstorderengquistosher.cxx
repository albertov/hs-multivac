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


#ifndef FILE_UPDATER_FASTMARCHINGFIRSTORDERENGQUISTOSHER_CXX


#include "fastmarchingfirstorderengquistosher.hxx"


namespace Multivac
{



  ///////////////////////////////
  // CONSTRUCTORS & DESTRUCTOR //
  ///////////////////////////////


  //! Default constructor.
  template <class T>
  CFastMarchingFirstOrderEngquistOsher<T>::
  CFastMarchingFirstOrderEngquistOsher()  throw()
  {

  }


  //! Main constructor.
  /*! \param Tmax_ a majorant of arrival times.
   */
  template <class T>
  CFastMarchingFirstOrderEngquistOsher<T>::
  CFastMarchingFirstOrderEngquistOsher(T TMax_)  throw()
  {

    this->NeedSpeedUpdateFlag = true;
    this->NeedInitializationFlag = false;
    
    this->TMax = TMax_;

  }


  //! Destructor.
  template <class T>
  CFastMarchingFirstOrderEngquistOsher<T>::
  ~CFastMarchingFirstOrderEngquistOsher()  throw()
  {

  }



  /////////////
  // METHODS //
  /////////////


  //! Should this updater be used for the narrow band level set method?
  /*! \return 'false' because the updater cannot be used for the narrow band
    level set method.
  */
  template <class T>
  inline bool CFastMarchingFirstOrderEngquistOsher<T>::IsNarrowBand() const
  {

    return false;

  }


  //! Should this updater be used for the fast marching method?
  /*! \return 'true' because the updater can be used for the fast marching
    method.
  */
  template <class T>
  inline bool CFastMarchingFirstOrderEngquistOsher<T>::IsFastMarching() const
  {

    return true;

  }


  //! Inits the updater.
  /*! The matrix 'PointersToNodes' is allocated and filled with -1.
    If this->PointersToNodes(i, j) is -1, the grid point (i, j) is not
    in trial points.
    \param Mesh orthogonal mesh.
    \param Phi level set function defined on an orthogonal mesh.
  */
  template <class T>
  inline void CFastMarchingFirstOrderEngquistOsher<T>::Init(CMesh<T>& Mesh,
							    CLevelSet<T>& Phi)
  {

    this->PointersToNodes.Reallocate(Mesh.GetNx(), Mesh.GetNy());
    this->PointersToNodes.Fill(-1);

  }


  //! Updates the level set function Phi.
  /*! This function updates the level set function Phi.  The trial point
    with the smallest time arrival time is removed from trial points set
    (and becomes a "known" point).  Then, neighbors arrival times are updated
    (if needed).  Updates are performed with the first order Engquist-Osher
    scheme.
    \param Delta_t time step.
    \param Mesh orthogonal mesh.
    \param F speed function defined on Mesh.
    \param Phi level set function defined on Mesh.
    \param CurrentTime current time.
  */
  template <class T>
  inline void CFastMarchingFirstOrderEngquistOsher<T>::
  UpdateLevelSet(T Delta_t, CMesh<T>& Mesh, CSpeedFunction<T>& F,
		 CLevelSet<T>& Phi, T CurrentTime)
  {

    int i, j, index;

    // Phi values.
    Matrix<T>& PhiValues = Phi.GetValues();

    T ArrivalTime;
    // A pointer to a node that will be added to the heap (trial points).
    heap_node_pointer NewTrialPoint;

    int Nx = Mesh.GetNx();
    int Ny = Mesh.GetNy();

    // The root (or head) of the sorted binary tree (or heap) is the trial
    // point with the smallest arrival time.
    i = this->TrialPoints.GetRoot()->GetX();
    j = this->TrialPoints.GetRoot()->GetY();

    // The point is removed from the trial points set.
    this->TrialPoints.DeleteRoot(this->PointersToNodes);


    // Four neighbors may be updated.
    // A neighbor is updated if (1) it is a grid point and if
    // (2) it's arrival time is not known (and definitively set).
    //
    // If the neighbor point is in trial points set, then, the
    // following operations are performed:
    //   1- update the new arrival time (RefreshTime);
    //   2- update the corresponding node (in the heap) with this
    //      new arrival time;
    //   3- sort the heap and update this->PointersToNodes.
    //
    // If the neighbor point is unknown (not in trial points set and
    // hasn't a set arrival time yet), then, the following operations
    // are performed:
    //   1- compute the arrival time at that point;
    //   2- add this point to trial points set (in the heap);
    //   3- sort the heap and update this->PointersToNodes.

    // East.
    if (i < Nx - 1)
      if (this->PointersToNodes(i + 1, j)>=0)
	{
	  ArrivalTime = RefreshTime(Mesh, F, Phi, i + 1, j);
	  index = this->PointersToNodes(i + 1, j);
	  this->TrialPoints(index)->SetValue(ArrivalTime);
	  this->TrialPoints.MoveUp(index, this->PointersToNodes);
	}
      else if  (PhiValues(i + 1, j)==this->TMax)
	{
	  ArrivalTime = RefreshTime(Mesh, F, Phi, i + 1, j);
	  NewTrialPoint = new heap_node_value_type(ArrivalTime, i + 1, j);
	  this->TrialPoints.Add(NewTrialPoint, this->PointersToNodes);
	}

    // West.
    if (i > 0)
      if (this->PointersToNodes(i - 1, j)>=0)
	{
	  ArrivalTime = RefreshTime(Mesh, F, Phi, i - 1, j);
	  index = this->PointersToNodes(i - 1, j);
	  this->TrialPoints(index)->SetValue(ArrivalTime);
	  this->TrialPoints.MoveUp(index, this->PointersToNodes);
	}
      else if  (PhiValues(i - 1, j)==this->TMax)
	{
	  ArrivalTime = RefreshTime(Mesh, F, Phi, i - 1, j);
	  NewTrialPoint = new heap_node_value_type(ArrivalTime, i - 1, j);
	  this->TrialPoints.Add(NewTrialPoint, this->PointersToNodes);
	}

    // North.
    if (j < Ny - 1)
      if (this->PointersToNodes(i, j + 1)>=0)
	{
	  ArrivalTime = RefreshTime(Mesh, F, Phi, i, j + 1);
	  index = this->PointersToNodes(i, j + 1);
	  this->TrialPoints(index)->SetValue(ArrivalTime);
	  this->TrialPoints.MoveUp(index, this->PointersToNodes);
	}
      else if  (PhiValues(i, j + 1)==this->TMax)
	{
	  ArrivalTime = RefreshTime(Mesh, F, Phi, i, j + 1);
	  NewTrialPoint = new heap_node_value_type(ArrivalTime, i, j + 1);
	  this->TrialPoints.Add(NewTrialPoint, this->PointersToNodes);
	}

    // South.
    if (j > 0)
      if (this->PointersToNodes(i, j - 1)>=0)
	{
	  ArrivalTime = RefreshTime(Mesh, F, Phi, i, j - 1);
	  index = this->PointersToNodes(i, j - 1);
	  this->TrialPoints(index)->SetValue(ArrivalTime);
	  this->TrialPoints.MoveUp(index, this->PointersToNodes);
	}
      else if  (PhiValues(i, j - 1)==this->TMax)
	{
	  ArrivalTime = RefreshTime(Mesh, F, Phi, i, j - 1);
	  NewTrialPoint = new heap_node_value_type(ArrivalTime, i, j - 1);
	  this->TrialPoints.Add(NewTrialPoint, this->PointersToNodes);
	}

  }


  //! Should computations be continued?
  /*! If 'KeepOnWorking' returns 'true', then computations should be
    continued. Computations will be continued if there are still trial points
    and if the current smallest arrival time is strictly less than TMax
    (majorant of arrival times).
    \return A boolean that indicates whether or not the program should keep
    on computing.
  */
  template <class T>
  inline bool CFastMarchingFirstOrderEngquistOsher<T>::KeepOnWorking() const
  {

    return ( !(this->TrialPoints.IsEmpty())
	     && (this->TrialPoints.GetRoot()->GetValue() < this->TMax) );

  }



  /******************
   * PRIVATE METHOD *
   ******************/


  //! Compute the arrival time at point (i, j).
  /*! The arrival time at point (i, j) is computed according to the current
    level set function.  That implies that neighbors of the point (i, j) have
    the right values.  The level set is updated with the computed arrival
    time which is returned.
    \param Mesh an orthogonal mesh.
    \param F speed function defined on Mesh.
    \param Phi level set function defined on Mesh.
    \param i first index of the point.
    \param j second index of the point.
    \return The arrival time at point (i, j).
  */
  template <class T>
  inline T CFastMarchingFirstOrderEngquistOsher<T>::
  RefreshTime(CMesh<T>& Mesh, CSpeedFunction<T>& F, CLevelSet<T>& Phi,
	      int i, int j)
  {

    // Phi values.
    Matrix<T>& PhiValues = Phi.GetValues();

    // F values.
    Matrix<T>& FValues = F.GetValues();

    int Nx = Mesh.GetNx();
    int Ny = Mesh.GetNy();

    // North, South, West, East and two temporary values.
    T N, S, W, E;
    T diff, time;
  
    // First, we set neighbor points arrival times.  If no arrival time
    // is available (because the point is outside the domain or because
    // no arrival time has been computed at this point yet), then the arrival
    // time is set to TMax (which will lead to the same behavior as with
    // an infinite arrival time).

    if ( (j<Ny-1) && (this->PointersToNodes(i, j + 1)<0.0)
	 && (fabs(PhiValues(i, j + 1))!=this->TMax) )
      N = PhiValues(i, j + 1);
    else
      N = this->TMax;

    if ( (j>0) && (this->PointersToNodes(i, j - 1)<0.0)
	 && (fabs(PhiValues(i, j - 1))!=this->TMax) )
      S = PhiValues(i, j - 1);
    else
      S = this->TMax;

    if ( (i<Nx-1) && (this->PointersToNodes(i + 1, j)<0.0)
	 && (fabs(PhiValues(i + 1, j))!=this->TMax) )
      E = PhiValues(i + 1, j);
    else
      E = this->TMax;

    if ( (i>0) && (this->PointersToNodes(i - 1, j)<0.0)
	 && (fabs(PhiValues(i - 1, j))!=this->TMax) )
      W = PhiValues(i - 1, j);
    else
      W = this->TMax;

    // If only one direction (north-south or west-east) has to be considered.
    if ( (W==this->TMax) && (E==this->TMax) )
      {
	if (S<N)
	  PhiValues(i, j) = min(PhiValues(i, j),
				PhiValues(i, j - 1)
				+ Mesh.GetDelta_y() / FValues(i, j));
	else
	  PhiValues(i, j) = min(PhiValues(i, j),
				PhiValues(i ,j + 1)
				+ Mesh.GetDelta_y() / FValues(i, j));
      }
    else if ( (N==this->TMax) && (S==this->TMax) )
      {
	if (E<W)
	  PhiValues(i,j) = min(PhiValues(i, j),
			       PhiValues(i + 1, j)
			       + Mesh.GetDelta_x() / FValues(i, j));
	else
	  PhiValues(i,j) = min(PhiValues(i, j),
			       PhiValues(i - 1, j)
			       + Mesh.GetDelta_x() / FValues(i, j));
      }
    else  // The two directions (north-south and west-east) are relevent.
      {
	// Chooses the smallest arrival time in each direction.
	N = N<S?N:S;
	W = W<E?W:E;
	// Updates the arrival time at point (i, j).
	if ( (N<W) && (N + Mesh.GetDelta_y() / FValues(i, j)<W) )
	  PhiValues(i, j) = min(PhiValues(i, j),
				N + Mesh.GetDelta_y() / FValues(i, j));
	else if ( (W<N) && (W + Mesh.GetDelta_x() / FValues(i, j)<N) )
	  PhiValues(i, j) = min(PhiValues(i, j),
				W + Mesh.GetDelta_x() / FValues(i, j));
	else
	  {
	    diff = N - W;
	    time = Mesh.GetDelta_x() / FValues(i, j);
	    PhiValues(i, j) = min(PhiValues(i, j),
				  ( N + W
				    + sqrt( - diff * diff
					    + 2.0 * time * time ) ) / 2.0 );
	  }
      }

    return PhiValues(i, j);

  }


}  // namespace Multivac.


#define FILE_UPDATER_FASTMARCHINGFIRSTORDERENGQUISTOSHER_CXX
#endif
