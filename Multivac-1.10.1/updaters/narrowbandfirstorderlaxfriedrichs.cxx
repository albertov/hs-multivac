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


#ifndef FILE_UPDATER_NARROWBANDFIRSTORDERLAXFRIEDRICHS_CXX


#include "narrowbandfirstorderlaxfriedrichs.hxx"


namespace Multivac
{



  ///////////////////////////////
  // CONSTRUCTORS & DESTRUCTOR //
  ///////////////////////////////


  //! Default constructor.
  template <class T>
  CNarrowBandFirstOrderLaxFriedrichs<T>::
  CNarrowBandFirstOrderLaxFriedrichs()  throw()
  {
    this->offset = 1;
  }


  //! Main constructor.
  /*! Tube semi width, barrier width and "outspace" width are set.
    \param TubeSemiWidth_ tube semi width, namely the number of cells
    on each side of the front.
    \param BarrierWidth_ barrier width, namely the number of cells along
    the barrier width.
    \param OutSpaceWidth_ "outspace" width, namely the number of cells along
    the "outspace" width.
  */
  template <class T>
  CNarrowBandFirstOrderLaxFriedrichs<T>::
  CNarrowBandFirstOrderLaxFriedrichs(int TubeSemiWidth_,
				     int BarrierWidth_,
				     int OutSpaceWidth_)  throw()
  {

    this->offset = 1;

    this->TubeSemiWidth = TubeSemiWidth_;
    this->BarrierWidth = BarrierWidth_;
    this->OutSpaceWidth = OutSpaceWidth_;

    this->NeedSpeedUpdateFlag = true;
    this->NeedInitializationFlag = false;

    this->TMax = 0.;

  }


  //! Destructor.
  template <class T>
  CNarrowBandFirstOrderLaxFriedrichs<T>::
  ~CNarrowBandFirstOrderLaxFriedrichs()  throw()
  {

  }



  /////////////
  // METHODS //
  /////////////


  //! Should this updater be used for the narrow band level set method?
  /*! \return 'true' because the updater can be used for the narrow band
    level set method.
  */
  template <class T>
  inline bool CNarrowBandFirstOrderLaxFriedrichs<T>::IsNarrowBand() const
  {

    return true;

  }


  //! Should this updater be used for the fast marching method?
  /*! \return 'false' because the updater cannot be used for the fast
    marching method.
  */
  template <class T>
  inline bool CNarrowBandFirstOrderLaxFriedrichs<T>::IsFastMarching() const
  {

    return false;

  }


  //! Inits the updater.
  /*! The matrix 'Temp' is allocated.
    \param Mesh orthogonal mesh.
    \param Phi level set function defined on an orthogonal mesh.
  */
  template <class T>
  void CNarrowBandFirstOrderLaxFriedrichs<T>::Init(CMesh<T>& Mesh,
						   CLevelSet<T>& Phi)
  {

    // 'Temp' contains more points than the mesh because it
    // contains ghost cells.
    this->Temp.Reallocate(Mesh.GetNx() + 2, Mesh.GetNy() + 2);

  }


  //! Updates the level set function Phi.
  /*! This function updates the level set function Phi on the tube,
    according to the speed function F.  The Lax-Friedrichs scheme is the
    space scheme and the time integration is performed by the Euler explicit
    method.
    \param Delta_t time step.
    \param Mesh orthogonal mesh.
    \param F speed function defined on Mesh.
    \param Phi level set function defined on Mesh.
    \param CurrentTime current time.
  */
  template <class T>
  void CNarrowBandFirstOrderLaxFriedrichs<T>::
  UpdateLevelSet(T Delta_t, CMesh<T>& Mesh, CSpeedFunction<T>& F,
		 CLevelSet<T>& Phi, T CurrentTime)
  {

    int i, j, jm, jp;
    int jmin, jmax;
    T sum_x, sum_y;

    // Phi values.
    Matrix<T>& PhiValues = Phi.GetValues();

    // F values.
    Matrix<T>& FValues = F.GetValues();

    T Delta_x = Mesh.GetDelta_x();
    T Delta_y = Mesh.GetDelta_y();

    int Nx = Mesh.GetNx();
    int Ny = Mesh.GetNy();
    int NX = Mesh.GetNx() + 1;
    int NY = Mesh.GetNy() + 1;

    T FMax(0);

    T Dx, Dmx, Dpx, Dy, Dmy, Dpy;
    T DxMax = - 10.0 * Delta_t / Delta_x;
    T DxMin = 10.0 * Delta_t / Delta_x;
    T DyMax = - 10.0 * Delta_t / Delta_y;
    T DyMin = 10.0 * Delta_t / Delta_y;
    T norm;

    T alphax = 0.0;
    T alphay = 0.0;

    T Speed, CentralValue;

    bool FirstLoop;

    // Computes bounds on partial derivatives of the Hamiltonian.
    for (i=1; i<NX; i++)
      // If the tube intersects the (i-1)th "row" of the mesh.
      if (!this->Tube(i-1).IsEmpty())
	{
	    
	  // Go to the first Vector<int>.
	  this->Tube(i-1).GoToTheHead();
	    
	  FirstLoop = true;
	  
	  // 'Tube(i-1)' stores vectors of integers.
	  // Those vectors are intervals:
	  // (15, 20) is stored if points (i-1, j-1) are inside
	  // the tube where 15 <= j <= 20.
	  // The following loop performs calculations in all intervals.
	  while (FirstLoop || this->Tube(i-1).GoToNext_StopAtTheTail())
	    {
		
	      FirstLoop = false;

	      if ( (i==1) || (this->Tube(i-2).IsEmpty()) )
		jm = NY;
	      else
		{
		  this->Tube(i-2).GoToTheHead();
		  jm = (this->Tube(i-2).GetCurrentValue())(0) + 1;
		}

	      if ( (i==Nx) || (this->Tube(i).IsEmpty()) )
		jp = NY;
	      else
		{
		  this->Tube(i).GoToTheHead();
		  jp = (this->Tube(i).GetCurrentValue())(0) + 1;
		}

	      jmin = (this->Tube(i-1).GetCurrentValue())(0) + 1;
	      jmax = (this->Tube(i-1).GetCurrentValue())(1) + 1;

	      for (j=jmin; j<=jmax; j++)
		{
		  
		  FMax = max(FMax, fabs(FValues(i-1, j-1)));

		  if (jm<j)
		    while (jm<j)
		      {
			if (jm!=(this->Tube(i-2).GetCurrentValue())(1) + 1)
			  jm++;
			else if (this->Tube(i-2).GoToNext_StopAtTheTail())
			  jm = (this->Tube(i-2).GetCurrentValue())(0) + 1;
			else
			  jm = NY;
		      }

		  if (jp<j)
		    while (jp<j)
		      {
			if (jp!=( this->Tube(i).GetCurrentValue())(1) + 1)
			  jp++;
			else if (this->Tube(i).GoToNext_StopAtTheTail())
			  jp = (this->Tube(i).GetCurrentValue())(0) + 1;
			else
			  jp = NY;
		      }

		  if ((jm == j) && (jp == j))
		    {
		      Dx = ( this->Temp(i+1, j) - this->Temp(i-1, j) )
			/ (2.0 * Delta_x);
		      if (Dx < DxMin)
			DxMin = Dx;
		      if (Dx > DxMax)
			DxMax = Dx;
		    }
		  else if (jm == j)
		    {
		      Dx = ( this->Temp(i, j) - this->Temp(i-1, j) )
			/ Delta_x;
		      if (Dx < DxMin)
			DxMin = Dx;
		      if (Dx > DxMax)
			DxMax = Dx;
		    }
		  else if (jp == j)
		    {
		      Dx = ( this->Temp(i+1, j) - this->Temp(i, j) )
			/ Delta_x;
		      if (Dx < DxMin)
			DxMin = Dx;
		      if (Dx > DxMax)
			DxMax = Dx;
		    }

		  if ((j > jmin) && (j < jmax))
		    {
		      Dy = ( this->Temp(i, j+1) - this->Temp(i, j-1) )
			/ (2.0 * Delta_y);
		      if (Dy < DyMin)
			DyMin = Dy;
		      if (Dy > DyMax)
			DyMax = Dy;
		    }
		  else if (j > jmin)
		    {
		      Dy = ( this->Temp(i, j) - this->Temp(i, j-1) )
			/ Delta_y;
		      if (Dy < DyMin)
			DyMin = Dy;
		      if (Dy > DyMax)
			DyMax = Dy;
		    }
		  else if (j < jmax)
		    {
		      Dy = ( this->Temp(i, j+1) - this->Temp(i, j) )
			/ Delta_y;
		      if (Dy < DyMin)
			DyMin = Dy;
		      if (Dy > DyMax)
			DyMax = Dy;
		    }

		}
	      
	    }
	  
	}

    norm = sqrt( max(DxMax*DxMax, DxMin*DxMin)
		 + max(DyMax*DyMax, DyMin*DyMin) );

    // alphax and alphay...
    alphax = FMax + F.GetMaxF1(DxMin, DxMax, DyMin, DyMax, norm);
    alphay = FMax + F.GetMaxF2(DxMin, DxMax, DyMin, DyMax, norm);

#ifdef MULTIVAC_CHECK_INSTABILITY
    if (Delta_t * alphax > Delta_x)
      throw CError_Instability(string("CNarrowBandFirstOrderLaxFriedrichs::")
			       + "UpdateLevelSet",
			       string("The CFL is not satisfied along x. ")
			       + "Delta_t = " + Talos::to_str(Delta_t)
			       + "; Delta_x = " + Talos::to_str(Delta_x)
			       + "; Delta_x / Delta_t = "
			       + Talos::to_str(Delta_x / Delta_t)
			       + "; alpha_x = " + Talos::to_str(alphax));
    
    if (Delta_t * alphay > Delta_y)
      throw CError_Instability(string("CNarrowBandFirstOrderLaxFriedrichs::")
			       + "UpdateLevelSet",
			       string("The CFL is not satisfied along y. ")
			       + "Delta_t = " + Talos::to_str(Delta_t)
			       + "; Delta_y = " + Talos::to_str(Delta_y)
			       + "; Delta_y / Delta_t = "
			       + Talos::to_str(Delta_y / Delta_t)
			       + "; alpha_y = " + Talos::to_str(alphay));
#endif

    // Updates are performed on tube points.
    for (i=0; i<Nx; i++)
      // If the tube intersects the (i-1)th "row" of the mesh.
      if (!this->Tube(i).IsEmpty())
	{
	    
	  // Go to the first Vector<int>.
	  this->Tube(i).GoToTheHead();
	    
	  FirstLoop = true;
	  
	  // 'Tube(i-1)' stores vectors of integers.
	  // Those vectors are intervals:
	  // (15, 20) is stored if points (i-1, j-1) are inside
	  // the tube where 15 <= j <= 20.
	  // The following loop performs calculations in all intervals.
	  while (FirstLoop || this->Tube(i).GoToNext_StopAtTheTail())
	    {
		
	      FirstLoop = false;
		
	      for (j=(this->Tube(i).GetCurrentValue())(0);
		   j<=(this->Tube(i).GetCurrentValue())(1); j++)
		{
		  
		  Speed = FValues(i, j);

		  // Central Value involved in the 5 point numerical scheme.
		  CentralValue = this->Temp(i+1, j+1);

		  Dmx = ( CentralValue - this->Temp(i, j+1) ) / Delta_x;
		  Dpx = ( this->Temp(i+2, j+1) - CentralValue ) / Delta_x;
		  Dmy = ( CentralValue - this->Temp(i+1, j) ) / Delta_y;
		  Dpy = ( this->Temp(i+1, j+2) - CentralValue ) / Delta_y;

		  sum_x = Dmx + Dpx;
		  sum_y = Dmy + Dpy;

		  PhiValues(i, j) -=
		    Delta_t * (
			       Speed * sqrt( sum_x * sum_x
					     + sum_y * sum_y )
			       - ( alphax * ( Dpx - Dmx )
				   + alphay * ( Dpy - Dmy ) )
			       ) / 2.0;
		  
		}
	      
	    }
	  
	}


    this->NeedInitializationFlag = false;

    // If a 'Barrier' point sign has changed, the front has reached
    //  the barrier; then, a reinitialization should be performed.
    for (i=0; (i<this->Barrier.GetM()) && (!this->NeedInitializationFlag);
	 i++)
      this->NeedInitializationFlag = (this->Temp(this->Barrier(i,0) + 1,
						 this->Barrier(i,1) + 1)
				      * PhiValues(this->Barrier(i,0),
						  this->Barrier(i,1)) < 0);
    
#ifdef MULTIVAC_CHECK_INSTABILITY
    // If a 'OutSpace' point sign has changed, the front has reached
    // the outspace; then, an instability has occurred.
    for (i=0; i<this->OutSpace.GetM(); i++)
      if (this->Temp(this->OutSpace(i,0) + 1, this->OutSpace(i,1) + 1) *
	  PhiValues(this->OutSpace(i,0), this->OutSpace(i,1)) < 0)
  	throw CError_Instability(string("CNarrowBandFirstOrderLaxFriedrichs")
				 + "::UpdateLevelSet",
				 "A point of the outspace has been reached");
#endif
    
  
    // 'Temp' is updated for next calculations (inside the tube and
    // on ghost cells).
    for (i=0; i<Nx; i++)
      if (!this->Tube(i).IsEmpty())
	{
	  
	  this->Tube(i).GoToTheHead();
	  
	  FirstLoop = true;
	  
	  while (FirstLoop || this->Tube(i).GoToNext_StopAtTheTail())
	    {
	      
	      FirstLoop = false;
	      
	      for (j=(this->Tube(i).GetCurrentValue())(0);
		   j<=(this->Tube(i).GetCurrentValue())(1); j++)
		this->Temp(i+1, j+1) = PhiValues(i, j);
	    }
	  
	}
    for (i=1; i<NX; i++)
      {
	this->Temp(i, 0) = PhiValues(i - 1, 0);
	this->Temp(i, NY) = PhiValues(i - 1, Ny - 1);
      }
    for (j=1; j<NY; j++)
      {
	this->Temp(0, j) = PhiValues(0, j - 1);
	this->Temp(NX, j) = PhiValues(Nx - 1, j - 1);
      }

  }


}  // namespace Multivac.


#define FILE_UPDATER_NARROWBANDFIRSTORDERLAXFRIEDRICHS_CXX
#endif
