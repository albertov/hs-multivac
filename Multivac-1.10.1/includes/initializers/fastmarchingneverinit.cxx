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


#ifndef FILE_INITIALIZER_FASTMARCHINGNEVERINIT_CXX


#include "fastmarchingneverinit.hxx"


namespace Multivac
{



  ///////////////////////////////
  // CONSTRUCTORS & DESTRUCTOR //
  ///////////////////////////////


  //! Default constructor.
  template <class T>
  CFastMarchingNeverInit<T>::CFastMarchingNeverInit()  throw()
  {

  }


  //! Destructor.
  template <class T>
  CFastMarchingNeverInit<T>::~CFastMarchingNeverInit()  throw()
  {

  }



  /////////////
  // METHODS //
  /////////////


  //! Should this initializer be used for the narrow band level set method?
  /*!
    This initializer is dedicated to the fast marching method and is
    therefore incompatible with the narrow band level method.
    \return 'false'.
  */
  template <class T>
  inline bool CFastMarchingNeverInit<T>::IsNarrowBand() const
  {

    return false;

  }


  //! Should this initializer be used for the fast marching method?
  /*!
    This initializer is dedicated to the fast marching method.
    \return 'true'.
  */
  template <class T>
  inline bool CFastMarchingNeverInit<T>::IsFastMarching() const
  {

    return true;

  }


  //! First initialization of the mesh.
  /*! Nothing is done.
    \param Mesh mesh.
  */
  template <class T>
  void CFastMarchingNeverInit<T>::FirstInitMesh(CMesh<T>& Mesh) const
  {

  }


  //! First initialization of the initial curve.
  /*! Nothing is done.
    \param Mesh mesh.
    \param InitialCurve initial curve.
  */
  template <class T>
  void CFastMarchingNeverInit<T>::
  FirstInitInitialCurve(CMesh<T>& Mesh,
			CInitialCurve<T>& InitialCurve) const
  {

  }


  //! First initialization of the level set function and the speed function.
  /*! The level set function is allocated and filled with TMax.  A few
    arrival times on computed on point close to the front.  They are trial
    points. The speed function is initialized.  Finally, the updater
    is initialized.
    \param Mesh mesh.
    \param InitialCurve initial curve.
    \param Phi level set function.
    \param F speed function.
    \param Updater updater.
  */
  template <class T>
  void CFastMarchingNeverInit<T>::
  FirstInitPhiAndF(CMesh<T>& Mesh, CInitialCurve<T>& InitialCurve,
		   CLevelSet<T>& Phi, CSpeedFunction<T>& F,
		   CUpdater<T>& Updater)
  {

    int i, j;

    // Phi values.
    Matrix<T>& PhiValues = Phi.GetValues();

    T Delta_x = Mesh.GetDelta_x();
    T Delta_y = Mesh.GetDelta_y();
    T Xmin = Mesh.GetXmin();
    T Ymin = Mesh.GetYmin();
    T X = Mesh.GetXmin();
    T Y = Mesh.GetYmin();
    int Nx = Mesh.GetNx();
    int Ny = Mesh.GetNy();

    /**** Initializes the speed function ****/
    
#ifdef MULTIVAC_CHECK_COMPATIBILITY
    // Checks whether dependencies are compatible with CFastMarchingNeverInit.
    // Actually, the fast marching method does not allow dependencies
    // upon the normal or the curvature.
    if (F.IsNormalDependent() || F.IsCurvatureDependent())
      throw CError_Incompatibility(string("CFastMarchingNeverInit<T>::")
				   + "FirstInitPhiAndF",
				   string("The speedfunction should not ")
				   + string("depend on the normal to the ")
				   + string("front or on the curvature.\n")
				   + string("   This limitation is inherent ")
				   + string("in the fast marching method. ")
				   + string("Use the narrow band level set ")
				   + "method instead");
#endif
    
    // Memory allocation.
    F.Init(Mesh);

    // Builds the speed function.
    Matrix<T>& FValues = F.GetValues();
    
    X = Xmin; Y = Ymin;
    for (i=0; i<Nx; i++)
      {
	for (j=0; j<Ny; j++)
	  {
	    FValues(i, j) = F(X, Y, 0.0);
	    Y += Delta_y;
	  }
	X += Delta_x;
	Y = Ymin;
      }

    /**** First initialization of the level set functionn ****/

    // Basic initialization of the level set function
    // (memory allocation).
    Phi.Init(Mesh);

    // Builds the level set.
    InitialCurve.SetDistances(Mesh, Phi);

    /**** Initializes first arrival times (trial points set) ****/

    // Inits the updater.
    Updater.Init(Mesh, Phi);

    // 'TrialPoints' will store all trial points in a sorted heap.
    ArrayHeap<typename CUpdater<T>::heap_node_value_type>& TrialPoints
      = Updater.GetTrialPoints();
    // New trial point to be added in 'TrialPoints'.
    typename CUpdater<T>::heap_node_value_type* NewTrialPoint;
    // Pointers to heap nodes.  'PointersToNodes(i, j)' gives the index
    // of the node (correponding to the point (i, j)) in the heap.  If the
    // point (i, j) is not in the heap, then 'PointersToNodes(i, j)' is
    // set to -1.
    Matrix<int>& PointersToNodes = Updater.GetPointersToNodes();
    PointersToNodes.Fill(-1);
    // Majorant of all arrival times.
    T TMax = Updater.GetTMax();

    // Allocates 'TrialPoints'.
    TrialPoints.Reallocate(0);

    // For all points on the grid...
    for (i=1; i<Nx-1; i++)
      for (j=1; j<Ny-1; j++)
	// If the current point is outside the front but close to the front.
	if ( (PhiValues(i, j) > 0.0)
	     && ( (PhiValues(i, j) * PhiValues(i-1, j) < 0.0)
		  || (PhiValues(i, j) * PhiValues(i+1, j) < 0.0)
		  || (PhiValues(i, j) * PhiValues(i, j-1) < 0.0)
		  || (PhiValues(i, j) * PhiValues(i, j+1) < 0.0)
		  )
	     )
	  {
	    // Sets arrival time at the current point.
	    PhiValues(i, j) = PhiValues(i, j)
	      / FValues(i, j);
	    // Adds the current point to the trial points set.
	    NewTrialPoint =
	      new typename CUpdater<T>::heap_node_value_type(PhiValues(i, j),
							     i, j);
	    TrialPoints.Add(NewTrialPoint, PointersToNodes);
	  }
	else if (PhiValues(i, j) > 0.0)
	  PhiValues(i, j) = TMax;
	else
	  PhiValues(i, j) = - TMax;

    // WARNING: it is assumed that the distance between domain edges
    // and the front is at least one cell.
    i = 0;
    for (j=0; j<Ny; j++)
      PhiValues(i, j) = TMax;

    i = Nx - 1;
    for (j=0; j<Ny; j++)
      PhiValues(i, j) = TMax;

    j = 0;
    for (i=0; i<Nx; i++)
      PhiValues(i, j) = TMax;

    j = Ny - 1;
    for (i=0; i<Nx; i++)
      PhiValues(i, j) = TMax;

  }


  //! Updates (reinitialization) the mesh.
  /*! Nothing is done.
    \param iter current iteration.
    \param Mesh mesh.
    \param Phi level set function.
    \param F speed function.
    \param Updater updater.
    \param CurrentTime current time.
  */
  template <class T>
  void CFastMarchingNeverInit<T>::InitMesh(int iter,
					   CMesh<T>& Mesh,
					   CLevelSet<T>& Phi,
					   CSpeedFunction<T>& F,
					   CUpdater<T>& Updater,
					   T CurrentTime) const
  {

  }


  //! Updates (reinitialization) the level set and the speed function.
  /*! Only the speed function is updated, using the updater supplied
    by the speed function.
    \param iter current iteration.
    \param Mesh mesh.
    \param Phi level set function.
    \param F speed function.
    \param Updater updater.
    \param CurrentTime current time.
  */
  template <class T>
  void CFastMarchingNeverInit<T>::InitPhiAndF(int iter,
					      CMesh<T>& Mesh,
					      CLevelSet<T>& Phi,
					      CSpeedFunction<T>& F,
					      CUpdater<T>& Updater,
					      T CurrentTime)
  {
    
    if (Updater.NeedSpeedUpdate() && F.IsTimeDependent())
      // then updates the speed function.
      {

	int i, j;
	
	T Delta_x = Mesh.GetDelta_x();
	T Delta_y = Mesh.GetDelta_y();
	T Xmin = Mesh.GetXmin();
	T Ymin = Mesh.GetYmin();
	T X(Xmin);
	T Y(Ymin);
	T Nx = Mesh.GetNx();
	T Ny = Mesh.GetNy();

	// Memory allocation.
	F.Init(Mesh);
	
	// Builds the speed function.
	Matrix<T>& FValues = F.GetValues();
	
	X = Xmin; Y = Ymin;
	for (i=0; i<Nx; i++)
	  {
	    for (j=0; j<Ny; j++)
	      {
		FValues(i, j) = F(X, Y, CurrentTime);
		Y += Delta_y;
	      }
	    X += Delta_x;
	    Y = Ymin;
	  }
	
      }

  }


  //! Builds the front on display purpose.
  /*! This function is empty because of the fact that the fast marching
    is a stationary method.
    \param iter current iteration.
    \param Mesh mesh.
    \param Phi level set function.
  */
  template <class T>
  void CFastMarchingNeverInit<T>::BuildCurveForDisplay(int iter,
						       CMesh<T>& Mesh,
						       CLevelSet<T>& Phi)
  {

  }


}  // namespace Multivac.


#define FILE_INITIALIZER_FASTMARCHINGNEVERINIT_CXX
#endif
