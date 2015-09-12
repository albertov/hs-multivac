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


#ifndef FILE_INITIALIZER_NARROWBANDNEVERINIT_CXX


#include "narrowbandneverinit.hxx"


namespace Multivac
{



  ///////////////////////////////
  // CONSTRUCTORS & DESTRUCTOR //
  ///////////////////////////////


  //! Default constructor.
  template <class T>
  CNarrowBandNeverInit<T>::CNarrowBandNeverInit()  throw()
  {

  }


  //! Destructor.
  template <class T>
  CNarrowBandNeverInit<T>::~CNarrowBandNeverInit()  throw()
  {

  }



  /////////////
  // METHODS //
  /////////////


  //! Should this initializer be used for the narrow band level set method?
  /*!
    This initializer is dedicated to the narrow band level set method.
    \return 'true'.
  */
  template <class T>
  inline bool CNarrowBandNeverInit<T>::IsNarrowBand() const
  {

    return true;

  }


  //! Should this initializer be used for the fast marching method?
  /*!
    This initializer is dedicated to the narrow band level set method and is
    therefore incompatible with the fast marching method.
    \return 'false'.
  */
  template <class T>
  inline bool CNarrowBandNeverInit<T>::IsFastMarching() const
  {

    return false;

  }


  //! First initialization of the mesh.
  /*! Nothing is done.
    \param Mesh mesh.
  */
  template <class T>
  void CNarrowBandNeverInit<T>::FirstInitMesh(CMesh<T>& Mesh) const
  {

  }


  //! First initialization of the initial curve.
  /*! Nothing is done.
    \param Mesh mesh.
    \param InitialCurve initial curve.
  */
  template <class T>
  void CNarrowBandNeverInit<T>::
  FirstInitInitialCurve(CMesh<T>& Mesh, CInitialCurve<T>& InitialCurve) const
  {

  }


  //! First initialization of the level set function and the speed function.
  /*! The level set function is allocated and initialized according to
    distance to the initial curve.  The speed function is initialized as well.
    The updater is initialized.  Then, the tube is built.
    \param Mesh orthogonal mesh.
    \param InitialCurve initial curve.
    \param Phi level set function.
    \param F speed function.
    \param Updater updater.
  */
  template <class T>
  void CNarrowBandNeverInit<T>::
  FirstInitPhiAndF(CMesh<T>& Mesh, CInitialCurve<T>& InitialCurve,
		   CLevelSet<T>& Phi, CSpeedFunction<T>& F,
		   CUpdater<T>& Updater)
  {
    
    int i, j;

    // Mesh properties.
    int Nx = Mesh.GetNx();
    int Ny = Mesh.GetNy();
    T Xmin = Mesh.GetXmin();
    T Ymin = Mesh.GetYmin();
    T Delta_x = Mesh.GetDelta_x();
    T Delta_y = Mesh.GetDelta_y();
    T X(Xmin), Y(Ymin);
    
    // Other values.
    T TubeSemiWidth = Delta_x * T(Updater.GetTubeSemiWidth());
    Vector<List<Vector<int> >, Vect_Full, NewAlloc<List<Vector<int> > > >&
      Tube = Updater.GetTube();

    // Basic initialization of the level set function
    // (memory allocation).
    Phi.Init(Mesh);

    // Basic initialization of the updater (memory
    // allocation).
    Updater.Init(Mesh, Phi);
    
#ifdef MULTIVAC_CHECK_COMPATIBILITY
    // Checks whether dependencies are compatible with CNarrowBandNeverInit.
    // No dependence upon the curvature is allowed.
    if (F.IsCurvatureDependent())
      throw CError_Incompatibility(string("CNarrowBandNeverInit<T>::")
				   + "FirstInitPhiAndF",
				   string("The speedfunction should not ")
				   + string("depend on the curvature.\n")
				   + string("   Use another Initializer ")
				   + "instead, e.g. CNarrowBandExtension");
#endif
    
    // Basic initialization of the speed function
    // (memory allocation).
    F.Init(Mesh);

    // Builds the level set.
    InitialCurve.SetDistances(Mesh, Phi);
  
    // Phi values.
    Matrix<T>& PhiValues = Phi.GetValues();
    // Set distances to the curve (in the whole grid).
    for (i=0; i<Nx; i++)
      for (j=0; j<Ny; j++)
	// The level set function should not have values (in absolute value)
	// above the tube semi width.
	if (fabs(PhiValues(i, j)) > TubeSemiWidth)
	  PhiValues(i, j) = sign(PhiValues(i,j)) * TubeSemiWidth;
    
    Tube.Reallocate(Nx);
    InitPhi(0, Mesh, Phi, Updater);
    
    // Builds the speed function.
    if (!F.IsTimeDependent()
	&& !F.IsNormalDependent()
	&& !F.IsCurvatureDependent())
      {

	Matrix<T>& FValues = F.GetValues();

	X = Xmin;
	Y = Ymin;

	for (i=0; i<Nx; i++)
	  {
	    for (j=0; j<Ny; j++)
	      {
		FValues(i, j) = F(X, Y, 0.0, 0.0, 0.0, 0.0);
		Y += Delta_y;
	      }
	    X += Delta_x;
	    Y = Ymin;
	  }

      }

    else
      
      InitF(Mesh, Updater.GetTube(), Phi, 0.0, F);
    
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
  void CNarrowBandNeverInit<T>::InitMesh(int iter,
					 CMesh<T>& Mesh,
					 CLevelSet<T>& Phi,
					 CSpeedFunction<T>& F,
					 CUpdater<T>& Updater,
					 T CurrentTime) const
  {

  }


  //! Updates (reinitialization) the level set and the speed function.
  /*! The front is reconstruted.  Then, a new tube is generated and
    the level set is defined in it.
    \param iter current iteration.
    \param Mesh mesh.
    \param Phi level set function.
    \param F speed function.
    \param Updater updater.
    \param CurrentTime current time.
  */
  template <class T>
  void CNarrowBandNeverInit<T>::InitPhiAndF(int iter,
					    CMesh<T>& Mesh,
					    CLevelSet<T>& Phi,
					    CSpeedFunction<T>& F,
					    CUpdater<T>& Updater,
					    T CurrentTime)
  {
  
    if (Updater.NeedInitialization())
      InitPhi(iter, Mesh, Phi, Updater);

    if ((Updater.NeedInitialization() && F.IsPositionDependent())
	|| F.IsTimeDependent()
	|| F.IsNormalDependent()
	|| F.IsCurvatureDependent())
      // Updates the speed function.
      InitF(Mesh, Updater.GetTube(), Phi, CurrentTime, F);
    
  }


  //! Updates (reinitialization) the level set only.
  /*! The front is reconstruted.  Then, a new tube is generated and
    the level set is defined in it.
    \param iter current iteration.
    \param Mesh mesh.
    \param Phi level set function.
    \param Updater updater.
  */
  template <class T>
  void CNarrowBandNeverInit<T>::InitPhi(int iter,
					CMesh<T>& Mesh,
					CLevelSet<T>& Phi,
					CUpdater<T>& Updater)
  {
    
    int i, j;

    // Phi values.
    Matrix<T>& PhiValues = Phi.GetValues();

    // Mesh properties.
    T Delta_x = Mesh.GetDelta_x();
    T Delta_y = Mesh.GetDelta_y();

    // Mesh properties.
    int Nx = Mesh.GetNx();
    int Ny = Mesh.GetNy();
    int offset = Updater.GetOffset();
    int NX = Nx + 2 * offset;
    int NY = Ny + 2 * offset;

    // Mesh properties.
    T Xmin = Mesh.GetXmin();
    T Ymin = Mesh.GetYmin();
    T X = Mesh.GetXmin();
    T Y = Mesh.GetYmin();

    /**** Builds the front ****/

    // Temporary list of points (from which the front is reconstructed).
    List<Vector<T> > Points;
    // Temporary vector.
    Vector<T> Point(2);

    // All points on the front are stored in the list 'Points'.
    // A point on the front lies where the level set function sign changes.
    for (i=0; i<Nx-1; i++)
      for (j=0; j<Ny-1; j++)
	{
	  if (PhiValues(i, j) * PhiValues(i+1, j) < 0)
	    {
	      Point(0) = Xmin
		+ (T(i) + PhiValues(i,j) / (PhiValues(i,j)
					    - PhiValues(i+1,j))) * Delta_x;
	      Point(1) = Ymin + T(j) * Delta_y;
	      Points.AddAtTheEnd(Point);
	    }
	  if (PhiValues(i, j) * PhiValues(i, j+1) < 0)
	    {
	      Point(0) = Xmin + T(i) * Delta_x;
	      Point(1) = Ymin
		+ (T(j) + PhiValues(i, j) / (PhiValues(i, j)
					     - PhiValues(i, j+1))) * Delta_y;
	      Points.AddAtTheEnd(Point);
	    }
	}

    // Edges.
    i = Nx - 1;
    for (j=0; j<Ny-1; j++)
      if (PhiValues(i, j) * PhiValues(i, j+1) < 0)
	{
	  Point(0) = Xmin + i * Delta_x;
	  Point(1) = Ymin
	    + (j + PhiValues(i, j) / (PhiValues(i, j)
				      - PhiValues(i, j+1))) * Delta_y;
	  Points.AddAtTheEnd(Point);
	}

    // Edges.
    j = Ny - 1;
    for (i=0; i<Nx-1; i++)
      if (PhiValues(i, j) * PhiValues(i+1, j) < 0)
	{
	  Point(0) = Xmin
	    + (i + PhiValues(i, j) / (PhiValues(i, j)
				      - PhiValues(i+1, j))) * Delta_x;
	  Point(1) = Ymin + j * Delta_y;
	  Points.AddAtTheEnd(Point);
	}

    // The curve is built by sorting points lying on the front.
    this->Front.InitWithUnsortedPoints_ClosestPointMethod(Points);

    // Now, the last time the front was reconstructed is
    // the current iteration.
    this->LastCurveUpdate = iter;
    // This reconstruction is fine for display purposes.
    this->LastCurveUpdateForDisplay = iter;

    /**** Computes distances to the front ****/

    // Temporary values used by the updater.  This temporary matrix
    // has been allocated because the updater needs at least one such matrix.
    // Then, 'Temp' can be used.  And it should be used to save memory.
    Matrix<T>& Temp = Updater.GetTemp();
    Temp.Zero();

    for (i=0; i<Nx; i++)
      for (j=0; j<Ny; j++)
	PhiValues(i, j) =
	  sign(PhiValues(i, j)) * Delta_x * Updater.GetTubeSemiWidth();

    // Number of points that will be in 'Barrier' and 'OutSpace'.
    int CountBarrierElements = 0;
    int CountOutSpaceElements = 0;

    // Tube, barrier and "outspace" dimensions.  It is needed to put points
    // in the right set.
    T TubeSemiWidth = Delta_x * T(Updater.GetTubeSemiWidth());
    T BarrierDistance = Delta_x * T(Updater.GetTubeSemiWidth()
				    - Updater.GetBarrierWidth());
    T OutSpaceDistance = Delta_x * T(Updater.GetTubeSemiWidth()
				     - Updater.GetOutSpaceWidth());

    // Gets front points.
    List<Vector<T> >& FrontPoints = this->Front.GetPoints();

    // Temporary vectors.
    Vector<T> A(2);
    Vector<T> B(2);

    // The maximum distance between two connected points on the front.
    T LimitDistance = sqrt( Delta_x * Delta_x + Delta_y * Delta_y );

    int imin, imax, jmin, jmax;
    T dist, u;
    T res, tempx, tempy, tempx0, tempy0;

    // Go to the first point of the front.
    FrontPoints.GoToTheHead();

    // Stores the first point in 'A' and adds this point at the tail.
    // Adding the first point at the tail is convenient in order to
    // compute distances to the segment between the first point and the
    // last point of the front.
    if (!FrontPoints.IsEmpty())
      {
	A.Copy(FrontPoints.GetCurrentValue());
	FrontPoints.AddAtTheEnd(A);
      }

    // For all points...
    while (FrontPoints.GoToNext_StopAtTheTail())
      {

	// Gets the current point.
	B.Copy(FrontPoints.GetCurrentValue());

	// If 'A' and 'B' are connected.
	if ( (dist = DistanceBetween(A, B)) < LimitDistance)
	  {

	    // Rectangle boundaries where distances to the front
	    // will be computed.
	    imin = max( 0,
			int(ceil( (min(A(0),B(0)) - Xmin - TubeSemiWidth)
				  / Delta_x ) )
			);
	    imax = min( Nx - 1,
			int(floor( (max(A(0),B(0)) - Xmin + TubeSemiWidth)
				   / Delta_x ) )
			);
	    jmin = max( 0,
			int(ceil( (min(A(1),B(1)) - Ymin - TubeSemiWidth)
				  / Delta_y ) )
			);
	    jmax = min( Ny - 1,
			int(floor( (max(A(1),B(1)) - Ymin + TubeSemiWidth)
				   / Delta_y ) )
			);

	    // Corner.
	    X = Xmin + imin * Delta_x;
	    Y = Ymin + jmin * Delta_y;

	    // In the rectangle, distances to the front are updated.
	    for (i=imin; i<=imax; i++)
	      {

		for (j=jmin; j<=jmax; j++)
		  {

		    // If 'A' and 'B' are very close, just computes
		    // the distance to those two points.
		    if (dist < exp(-30.) )
		      {

			if (PhiValues(i, j) > 0.0)
			  {
			    tempx = X - A(0);
			    tempy = Y - A(1);
			    tempx0 = X - B(0);
			    tempy0 = Y - B(1);
			    res = sqrt( min( tempx*tempx + tempy*tempy,
					     tempx0*tempx0
					     + tempy0*tempy0 ) );
			    PhiValues(i, j) = min(PhiValues(i, j), res);
			  }
			else
			  {
			    tempx = X - A(0);
			    tempy = Y - A(1);
			    tempx0 = X - B(0);
			    tempy0 = Y - B(1);
			    res = sqrt( min( tempx*tempx + tempy*tempy,
					     tempx0*tempx0
					     + tempy0*tempy0 ) );
			    PhiValues(i, j) = max(PhiValues(i, j), - res);
			  }

		      }
		    else  // 'A' and 'B' are not very close.
		      {

			u = ( (X - A(0)) * (B(0)-A(0)) +
			      (Y - A(1)) * (B(1) - A(1)) ) / (dist*dist);

			if (PhiValues(i, j) > 0.0)
			  {
			    if (u < 0)
			      {
				tempx = X - A(0);
				tempy = Y - A(1);
				PhiValues(i, j) = min(PhiValues(i, j),
						      sqrt(tempx*tempx
							   + tempy*tempy));
			      }
			    else if (u > 1)
			      {
				tempx = X - B(0);
				tempy = Y - B(1);
				PhiValues(i, j) = min(PhiValues(i, j),
						      sqrt(tempx*tempx
							   + tempy*tempy));
			      }
			    else
			      {
				tempx = X - A(0) - u * (B(0) - A(0));
				tempy = Y - A(1) - u * (B(1) - A(1));
				PhiValues(i, j) = min(PhiValues(i, j),
						      sqrt(tempx*tempx
							   + tempy*tempy));
			      }

			  }
			else  // (PhiValues(i, j) < 0.0).
			  {

			    if (u < 0)
			      {
				tempx = X - A(0);
				tempy = Y - A(1);
				PhiValues(i, j) = max(PhiValues(i, j),
						      - sqrt(tempx*tempx
							     + tempy*tempy));
			      }
			    else if (u > 1)
			      {
				tempx = X - B(0);
				tempy = Y - B(1);
				PhiValues(i, j) = max(PhiValues(i, j),
						      - sqrt(tempx*tempx
							     + tempy*tempy));
			      }
			    else
			      {
				tempx = X - A(0) - u * (B(0) - A(0));
				tempy = Y - A(1) - u * (B(1) - A(1));
				PhiValues(i, j) = max(PhiValues(i, j),
						      - sqrt(tempx*tempx
							     + tempy*tempy));
			      }

			  }

		      }

		    // To count 'Barrier' and 'OutSpace' element.
		    if (Temp(i, j) == 2.0)
		      CountBarrierElements--;
		    else if (Temp(i, j) == 3.0)
		      CountOutSpaceElements--;
		    
		    // Determines where the point lies (outside the tube,
		    // in the barrier, ...).
		    if (fabs(PhiValues(i, j)) < TubeSemiWidth)
		      if (fabs(PhiValues(i, j)) >= OutSpaceDistance)
			Temp(i, j) = 3.0;
		      else if (fabs(PhiValues(i, j)) >= BarrierDistance)
			Temp(i, j) = 2.0;
		      else
			Temp(i, j) = 1.0;

		    // To count 'Barrier' and 'OutSpace' element.
		    if (Temp(i, j) == 2.0)
		      CountBarrierElements++;
		    else if (Temp(i, j) == 3.0)
		      CountOutSpaceElements++;
		
		    Y += Delta_y;
		
		  }

		X += Delta_x;
		Y = Ymin + jmin * Delta_y;
	    
	      }

	  }
      
	// 'A' <- 'B' and 'B' is going to be the next point on the front.
	A.Copy(B);

      }

    // Checks, for all points in the narrow band, if there are
    // enough neighbors within the narrow band along x and along y.

    // Number of neighbors.
    int count;
    int prev(0), ind0, ind1;

    // Along y.
    for (i=0; i<Nx; i++)
      {
	count = 0;
	for (j=0; j<Ny; j++)
	  if (Temp(i, j)>0)  // In the narrow band.
	    {
	      if (count == 0)
		prev = j;
	      count++;
	    }
	  else if ( (count != 0) && (count <= offset) )
	    // More neighbors are required.
	    {
	      ind0 = prev;
	      ind1 = j-1;
	      // The narrow band is expanded to include more points.
	      while ( (ind1 - ind0 < offset)
		      && ( (ind1 != Ny-1) || (ind0 != 0) ) )
		{
		  if (PhiValues(i, ind0) < PhiValues(i, ind1))
		    ind0--;
		  else
		    ind1++;
		  if (ind0 == -1)
		    {
		      ind0 = 0;
		      ind1++;
		    }
		  if (ind1 == Ny)
		    {
		      ind1 = Ny-1;
		      ind0--;
		    }
		  // The new point is added in the 'OutSpace'.
		  if (Temp(i, ind0) != 3.0)
		    CountOutSpaceElements++;
		  Temp(i, ind0) = 3;
		  if (Temp(i, ind1) != 3.0)
		    CountOutSpaceElements++;
		  Temp(i, ind1) = 3;
		  while ( (ind0 != 0) && (Temp(i, ind0) != 0) )
		    ind0--;
		  while ( (ind1 != Ny-1) && (Temp(i, ind1) != 0) )
		    ind1++;
		}
	      if (ind1 == j)
		count++;
	      else
		count = 0;
	    }
	  else
	    count = 0;
      }

    // Along x.
    for (j=0; j<Ny; j++)
      {
	count = 0;
	for (i=0; i<Nx; i++)
	  if (Temp(i, j)>0)  // In the narrow band.
	    {
	      if (count == 0)
		prev = i;
	      count++;
	    }
	  else if ( (count != 0) && (count <= offset) )
	    // More neighbors are required.
	    {
	      ind0 = prev;
	      ind1 = i-1;
	      // The narrow band is expanded to include more points.
	      while ( (ind1 - ind0 < offset)
		      && ( (ind1 != Nx-1) || (ind0 != 0) ) )
		{
		  if (PhiValues(ind0, j) < PhiValues(ind1, j))
		    ind0--;
		  else
		    ind1++;
		  if (ind0 == -1)
		    {
		      ind0 = 0;
		      ind1++;
		    }
		  if (ind1 == Nx)
		    {
		      ind1 = Nx-1;
		      ind0--;
		    }
		  // The new point is added in the 'OutSpace'.
		  if (Temp(ind0, j) != 3.0)
		    CountOutSpaceElements++;
		  Temp(ind0, j) = 3;
		  if (Temp(ind1, j) != 3.0)
		    CountOutSpaceElements++;
		  Temp(ind1, j) = 3;
		  while ( (ind0 != 0) && (Temp(ind0, j) != 0) )
		    ind0--;
		  while ( (ind1 != Nx-1) && (Temp(ind1, j) != 0) )
		    ind1++;
		}
	      if (ind1 == i)
		count++;
	      else
		count = 0;
	    }
	  else
	    count = 0;
      }

    for (i=0; i<Nx; i++)
      for (j=0; j<Ny; j++)
	if (Temp(i, j)==0)
	  PhiValues(i, j) = sign(PhiValues(i,j)) * TubeSemiWidth;

    // Removes the last point (which is currently the first point --
    // see before the previous 'while' loop).
    FrontPoints.GoToTheTail();
    FrontPoints.DeleteCurrent();

    /**** Builds the tube ****/

    // Gets the tube address.
    Vector<List<Vector<int> >, Vect_Full, NewAlloc<List<Vector<int> > > >&
      Tube = Updater.GetTube();

    // Gets 'Barrier' and 'OutSpace' addresses.
    Matrix<int>& Barrier = Updater.GetBarrier();
    Matrix<int>& OutSpace = Updater.GetOutSpace();

    // Allocations.
    Barrier.Reallocate(CountBarrierElements, 2);
    OutSpace.Reallocate(CountOutSpaceElements, 2);
  
    // Temporary vector.
    Vector<int> Interval(2);
  
    // Index counters.
    CountBarrierElements = 0;
    CountOutSpaceElements = 0;

    // All grid points are "scanned" to build 'Barrier' and 'OutSpace'.
    for (i=0; i<Nx; i++)
      for (j=0; j<Ny; j++)
	if (Temp(i, j) != 0.0)  // If the point is not outside.
	  {
	    if (Temp(i, j) == 2.0)
	      {
		Barrier(CountBarrierElements, 0) = i;
		Barrier(CountBarrierElements, 1) = j;
		CountBarrierElements++;
	      }
	    else if (Temp(i, j) == 3.0)
	      {
		OutSpace(CountOutSpaceElements, 0) = i;
		OutSpace(CountOutSpaceElements, 1) = j;
		CountOutSpaceElements++;
	      }
	  }

    // Temporary index.
    int lastj;

    // Builds the tube.
    for (i=0; i<Nx; i++)
      {

	Tube(i).ClearAll();

	lastj = -1;

	for (j=0; j<Ny; j++)
	  // If the current point is outside the tube and if
	  // previous points were inside.
	  if ((Temp(i, j) == 0.0) && (lastj != -1))
	    {
	      Interval(0) = lastj;
	      Interval(1) = j - 1;
	      Tube(i).AddAtTheEnd(Interval);
	      lastj = -1;
	    }
	// If the current point is in the tube and if
	// previous points were outside.
	  else if ((Temp(i, j) != 0.0) && (lastj == -1))
	    lastj = j;
	
	// If we reached the domain edge and if
	// previous points were inside the tube.
	if (lastj != -1)
	  {
	    Interval(0) = lastj;
	    Interval(1) = Ny - 1;
	    Tube(i).AddAtTheEnd(Interval);
	  }

      }

    // Fills 'Temp' with 'Phi' values (the updater needs it).
    for (i=0; i<Nx; i++)
      for (j=0; j<Ny; j++)
	Temp(i + offset, j + offset) = PhiValues(i, j);

    // Ghost cells.
    for (i=offset; i<NX-offset; i++)
      for (j=0; j<offset; j++)
	{
	  Temp(i, j) = Temp(i, offset);
	  Temp(i, Ny+offset+j) = Temp(i, Ny-1+offset);
	}
    for (i=0; i<offset; i++)
      for (j=offset; j<NY-offset; j++)
	{
	  Temp(i, j) = Temp(offset, j);
	  Temp(Nx+offset+i, j) = Temp(Nx-1+offset, j);
	}

  }


  //! Updates (reinitialization) the speed function.
  /*! Updates (reinitialization) on a tube the speed function
    according to a given levet set function.
    \param Mesh mesh.
    \param Tube tube in which values are updated.
    \param Phi level set function.
    \param CurrentTime current time.
  */
  template <class T>
  void CNarrowBandNeverInit<T>::
  InitF(CMesh<T>& Mesh,
	Vector<List<Vector<int> >, Vect_Full,
	NewAlloc<List<Vector<int> > > >& Tube,
	CLevelSet<T>& Phi, T CurrentTime, CSpeedFunction<T>& F)
  {
  
    int i, j, jmin, jmax, jm(0), jp(0);
    bool FirstLoop;

    // Speed values.
    Matrix<T>& Values = F.GetValues();

    // Phi values.
    Matrix<T>& PhiValues = Phi.GetValues();

    // Mesh properties.
    int Nx = Mesh.GetNx();
    int Ny = Mesh.GetNy();
    int NX = Nx + 1;
    int NY = Ny + 1;

    // Mesh properties.
    T Delta_x = Mesh.GetDelta_x();
    T Delta_y = Mesh.GetDelta_y();

    // Mesh properties.
    T Xmin = Mesh.GetXmin();
    T Ymin = Mesh.GetYmin();

    // Finite differences.
    T Dmx, Dpx, Dmy, Dpy;
    // Squared.
    T Dmx2, Dpx2, Dmy2, Dpy2;

    // Normal vector coordinates.
    T NormalX, NormalY;
    // Norm of the normal vector.
    T norm;
    // theta: normal angle.
    // costheta = cos(theta) and sintheta = sin(theta).
    T costheta = 0.0;
    T sintheta = 0.0;
    // Squared.
    T cos2, sin2;

    // Curvature.
    T curvature = 0.0;


    // Computes bounds on partial derivatives of the Hamiltonian.
    for (i=1; i<NX; i++)
      // If the tube intersects the (i-1)th "row" of the mesh.
      if (!Tube(i-1).IsEmpty())
	{
	    
	  // Go to the first Vector<int>.
	  Tube(i-1).GoToTheHead();
	    
	  FirstLoop = true;
	  
	  // 'Tube(i-1)' stores vectors of integers.
	  // Those vectors are intervals:
	  // (15, 20) is stored if points (i-1, j-1) are
	  // inside the tube where 15 <= j <= 20.
	  // The following loop performs calculations in all intervals.
	  while (FirstLoop || Tube(i-1).GoToNext_StopAtTheTail())
	    {
		
	      FirstLoop = false;

	      if (F.IsNormalDependent())
		{
		  
		  if ( (i==1) || (Tube(i-2).IsEmpty()) )
		    jm = NY;
		  else
		    {
		      Tube(i-2).GoToTheHead();
		      jm = (Tube(i-2).GetCurrentValue())(0) + 1;
		    }
		  
		  if ( (i==Nx) || (Tube(i).IsEmpty()) )
		    jp = NY;
		  else
		    {
		      Tube(i).GoToTheHead();
		      jp = (Tube(i).GetCurrentValue())(0) + 1;
		    }
		  
		}

	      jmin = (Tube(i-1).GetCurrentValue())(0) + 1;
	      jmax = (Tube(i-1).GetCurrentValue())(1) + 1;

	      for (j=jmin; j<=jmax; j++)
		{
		  
		  if (F.IsNormalDependent())
		    {
		      if (jm<j)
			while (jm<j)
			  {
			    if (jm!=(Tube(i-2).GetCurrentValue())(1) + 1)
			      jm++;
			    else if (Tube(i-2).GoToNext_StopAtTheTail())
			      jm = (Tube(i-2).GetCurrentValue())(0) + 1;
			    else
			      jm = NY;
			  }
		      
		      if (jp<j)
			while (jp<j)
			  {
			    if (jp!=( Tube(i).GetCurrentValue())(1) + 1)
			      jp++;
			    else if (Tube(i).GoToNext_StopAtTheTail())
			      jp = (Tube(i).GetCurrentValue())(0) + 1;
			    else
			      jp = NY;
			  }
		      
		      if (jm == j)
			{
			  Dmx = ( PhiValues(i-1, j-1) - PhiValues(i-2, j-1) )
			    / Delta_x;
			  if (jp == j)
			    Dpx = ( PhiValues(i, j-1) - PhiValues(i-1, j-1) )
			      / Delta_x;
			  else
			    Dpx = Dmx;
			}
		      else
			{
			  if (jp == j)
			    Dpx = ( PhiValues(i, j-1) - PhiValues(i-1, j-1) )
			      / Delta_x;
			  else
			    Dpx = 0.0;
			  Dmx = Dpx;
			}

		      if (j > jmin)
			{
			  Dmy = ( PhiValues(i-1, j-1) - PhiValues(i-1, j-2) )
			    / Delta_y;
			  if (j < jmax)
			    Dpy = ( PhiValues(i-1, j) - PhiValues(i-1, j-1) )
			      / Delta_y;
			  else
			    Dpy = Dmy;
			}
		      else
			{
			  if (j < jmax)
			    Dpy = ( PhiValues(i-1, j) - PhiValues(i-1, j-1) )
			      / Delta_y;
			  else
			    Dpy = 0.0;
			  Dmy = Dpy;
			}
		      
		      // To speed up calculations.
		      Dmx2 = Dmx * Dmx;
		      Dpx2 = Dpx * Dpx;
		      Dmy2 = Dmy * Dmy;
		      Dpy2 = Dpy * Dpy;
		      
		      NormalX = 0.0;
		      NormalY = 0.0;
		      
		      norm = sqrt( Dpx2 + Dpy2 );
		      // If 'norm' is significant.
		      if ( norm > _limit)
			{
			  NormalX += Dpx / norm;
			  NormalY += Dpy / norm;
			}
		      
		      norm = sqrt( Dmx2 + Dpy2 );
		      // If 'norm' is significant.
		      if ( norm > _limit)
			{
			  NormalX += Dmx / norm;
			  NormalY += Dpy / norm;
			}
		      
		      norm = sqrt( Dpx2 + Dmy2 );
		      // If 'norm' is significant.
		      if ( norm > _limit)
			{
			  NormalX += Dpx / norm;
			  NormalY += Dmy / norm;
			}
		      
		      norm = sqrt( Dmx2 + Dmy2 );
		      // If 'norm' is significant.
		      if ( norm > _limit)
			{
			  NormalX += Dmx / norm;
			  NormalY += Dmy / norm;
			}
		      
		      // Norm of the normal vector.
		      norm = sqrt( NormalX * NormalX + NormalY * NormalY );
		      
		      if (norm > _limit)
			{
			  costheta = NormalX / norm;
			  sintheta = NormalY / norm;
			}
		      else
			{
			  costheta = 1.0;
			  sintheta = 0.0;
			}

		      sin2 = sintheta * sintheta;
		      cos2 = costheta * costheta;
		      
		    }

		  Values(i-1, j-1) = F(Xmin + Delta_x * (i-1),
				       Ymin + Delta_y * (j-1),
				       CurrentTime,
				       costheta, sintheta,
				       curvature);
		  
		}
	      
	    }
	  
	}
  
  }


  //! Builds the front on display purpose.
  /*! A curve is built to approximate the front.  This curve contains
    points that should be used to display the front.
    \param iter current iteration.
    \param Mesh mesh.
    \param Phi level set function.
  */
  template <class T>
  void CNarrowBandNeverInit<T>::BuildCurveForDisplay(int iter,
						     CMesh<T>& Mesh,
						     CLevelSet<T>& Phi)
  {

    // If the curve has already been built, no need to building
    // the front again.
    if (iter == this->LastCurveUpdateForDisplay)
      return;

    // Points will store all points of the curve and will be used
    // to build the curve (by sorting points).
    List<Vector<T> > Points;
    // Used to store points that are going to be added to Points.
    Vector<T> X(2);
    
    int i, j;

    // Phi values.
    Matrix<T>& PhiValues = Phi.GetValues();

    T Delta_x = Mesh.GetDelta_x();
    T Delta_y = Mesh.GetDelta_y();
    int Nx = Mesh.GetNx();
    int Ny = Mesh.GetNy();
    T Xmin = Mesh.GetXmin();
    T Ymin = Mesh.GetYmin();
    
    // Initialization of the list of points.
    Points.ClearAll();

    // Searching for points on the front.
    // If the sign of the level set function Phi changes from one grid
    // point to the other, then a point on the front lies between these two
    // grid points.  An approximation of the coordinates of this point is
    // computed and stored in Points.
    for (i=0; i<Nx-1; i++)
      for (j=0; j<Ny-1; j++)
	{
	  if (PhiValues(i, j) * PhiValues(i+1, j) < 0)
	    {
	      X(0) = Xmin
		+ ( T(i) + PhiValues(i,j)
		    / (PhiValues(i,j) - PhiValues(i+1,j)) ) * Delta_x;
	      X(1) = Ymin + T(j) * Delta_y;
	      Points.AddAtTheEnd(X);
	    }
	  if (PhiValues(i, j) * PhiValues(i, j+1) < 0)
	    {
	      X(0) = Xmin + T(i) * Delta_x;
	      X(1) = Ymin
		+ (T(j) + PhiValues(i, j)
		   / (PhiValues(i, j) - PhiValues(i, j+1)) ) * Delta_y;
	      Points.AddAtTheEnd(X);
	    }
	}

    // The last column.
    j = Ny - 1;
    for (i=0; i<Nx-1; i++)
      if (PhiValues(i, j) * PhiValues(i+1, j) < 0)
	{
	  X(0) = Xmin + (i + PhiValues(i, j)
			 / (PhiValues(i, j) - PhiValues(i+1, j))) * Delta_x;
	  X(1) = Ymin + j * Delta_y;
	  Points.AddAtTheEnd(X);
	}

    // The last row.
    i = Nx - 1;
    for (j=0; j<Ny-1; j++)
      if (PhiValues(i, j) * PhiValues(i, j+1) < 0)
	{
	  X(0) = Xmin + i * Delta_x;
	  X(1) = Ymin + (j + PhiValues(i, j)
			 / (PhiValues(i, j) - PhiValues(i, j+1))) * Delta_y;
	  Points.AddAtTheEnd(X);
	}

    // The list of points ('Points') is sorted to get the curve.
    this->Front.InitWithUnsortedPoints_ClosestPointMethod(Points);

    // The curve has been built at iteration 'iter'.
    this->LastCurveUpdateForDisplay = iter;
    this->LastCurveUpdate = iter;

  }


}  // namespace Multivac.


#define FILE_INITIALIZER_NARROWBANDNEVERINIT_CXX
#endif
