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


#ifndef FILE_UPDATER_NARROWBANDENO2ENGQUISTOSHER_CXX


#include "narrowbandeno2engquistosher.hxx"


namespace Multivac
{



  ///////////////////////////////
  // CONSTRUCTORS & DESTRUCTOR //
  ///////////////////////////////


  //! Default constructor.
  template <class T>
  CNarrowBandEno2EngquistOsher<T>::CNarrowBandEno2EngquistOsher()  throw()
  {
    this->offset = 2;
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
  CNarrowBandEno2EngquistOsher<T>::
  CNarrowBandEno2EngquistOsher(int TubeSemiWidth_,
			       int BarrierWidth_,
			       int OutSpaceWidth_)  throw()
  {

    this->offset = 2;

    this->TubeSemiWidth = TubeSemiWidth_;
    this->BarrierWidth = BarrierWidth_;
    this->OutSpaceWidth = OutSpaceWidth_;

    this->NeedSpeedUpdateFlag = true;
    this->NeedInitializationFlag = false;

    this->TMax = 0.;

  }


  //! Destructor.
  template <class T>
  CNarrowBandEno2EngquistOsher<T>::~CNarrowBandEno2EngquistOsher()  throw()
  {

  }



  /////////////
  // METHODS //
  /////////////


  //! Should this updater be used for the narrow band level set method?
  /*! \return 'true' because the updater can be used for the narrow band level
    set method.
  */
  template <class T>
  inline bool CNarrowBandEno2EngquistOsher<T>::IsNarrowBand() const
  {

    return true;

  }


  //! Should this updater be used for the fast marching method?
  /*! \return 'false' because the updater cannot be used for the fast marching
    method.
  */
  template <class T>
  inline bool CNarrowBandEno2EngquistOsher<T>::IsFastMarching() const
  {

    return false;

  }


  //! Inits the updater.
  /*! The matrix 'Temp' is allocated.
    \param Mesh orthogonal mesh.
    \param Phi level set function defined on an orthogonal mesh.
  */
  template <class T>
  void CNarrowBandEno2EngquistOsher<T>::Init(CMesh<T>& Mesh,
					     CLevelSet<T>& Phi)
  {

    // 'Temp' contains more points than the mesh because it
    // contains ghost cells.
    this->Temp.Reallocate(Mesh.GetNx() + 4, Mesh.GetNy() + 4);

  }


  //! Updates the level set function Phi.
  /*! This function updates the level set function Phi on the tube,
    according to the speed function F.  The Engquist-Osher scheme is the space
    scheme and the time integration is performed by the Euler explicit method.
    \param Delta_t time step.
    \param Mesh orthogonal mesh.
    \param F speed function defined on Mesh.
    \param Phi level set function defined on Mesh.
    \param CurrentTime current time.
  */
  template <class T>
  void CNarrowBandEno2EngquistOsher<T>::UpdateLevelSet(T Delta_t,
						       CMesh<T>& Mesh,
						       CSpeedFunction<T>& F,
						       CLevelSet<T>& Phi,
						       T CurrentTime)
  {

    int i, j;

    // Phi values.
    Matrix<T>& PhiValues = Phi.GetValues();

    // F values.
    Matrix<T>& FValues = F.GetValues();

    T Delta_x = Mesh.GetDelta_x();
    T DiscrX = 1.0 / (Delta_x * Delta_x);
    T Delta_y = Mesh.GetDelta_y();
    T DiscrY = 1.0 / (Delta_y * Delta_y);

    int Nx = Mesh.GetNx();
    int Ny = Mesh.GetNy();
    int NX = Mesh.GetNx() + 4;
    int NY = Mesh.GetNy() + 4;

    T Speed, CentralValue;
    T Dmx, Dpx, Dmy, Dpy;
    T amx, apx, bx;
    T amy, apy, by;

    bool FirstLoop;
    
    // Updates are performed on tube points.
    for (i=0; i<Nx; i++)
      // If the tube intersects the (i-1)th "row" of the mesh.
      if (!this->Tube(i).IsEmpty())
	{
	    
	  // Go to the first Vector<int>.
	  this->Tube(i).GoToTheHead();
	    
	  FirstLoop = true;
	  
	  // 'Tube(i-1)' stores vectors of integers.  Those vectors are
	  // intervals:
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
		  
		  CentralValue = this->Temp(i+2, j+2);
		  
		  apx = this->Temp(i+4, j+2) - 2.0 * this->Temp(i+3, j+2)
		    + CentralValue;
		  bx = this->Temp(i+3, j+2) - 2.0 * CentralValue
		    + this->Temp(i+1, j+2);
		  amx = CentralValue - 2.0 * this->Temp(i+1, j+2)
		    + this->Temp(i, j+2);
		  
		  apy = this->Temp(i+2, j+4) - 2.0 * this->Temp(i+2, j+3)
		    + CentralValue;
		  by = this->Temp(i+2, j+3) - 2.0 * CentralValue
		    + this->Temp(i+2, j+1);
		  amy = CentralValue - 2.0 * this->Temp(i+2, j+1)
		    + this->Temp(i+2, j);

		  if (Speed > 0)
		    {
		      if (abs(amx)>abs(bx))
			Dmx = max(CentralValue - this->Temp(i+1,j+2)
				  + bx / 2.0, 0.0);
		      else
			Dmx = max(CentralValue - this->Temp(i+1,j+2)
				  + amx / 2.0, 0.0);
		      if (abs(bx)>abs(apx))
			Dpx = min(this->Temp(i+3,j+2) - CentralValue
				  - apx / 2.0, 0.0);
		      else
			Dpx = min(this->Temp(i+3,j+2) - CentralValue
				  - bx / 2.0, 0.0);

		      if (abs(amy)>abs(by))
			Dmy = max(CentralValue - this->Temp(i+2,j+1)
				  + by / 2.0, 0.0);
		      else
			Dmy = max(CentralValue - this->Temp(i+2,j+1)
				  + amy / 2.0, 0.0);
		      if (abs(by)>abs(apy))
			Dpy = min(this->Temp(i+2,j+3) - CentralValue
				  - apy / 2.0, 0.0);
		      else
			Dpy = min(this->Temp(i+2,j+3) - CentralValue
				  - by / 2.0, 0.0);
		      
		      PhiValues(i, j) -=
			Delta_t * Speed
			* sqrt( DiscrX * ( Dmx * Dmx + Dpx * Dpx  )
				+ DiscrY * ( Dmy * Dmy + Dpy * Dpy) );
		    }
		  else
		    {
		      if (abs(amx)>abs(bx))
			Dmx = min(CentralValue - this->Temp(i+1,j+2)
				  + bx / 2.0, 0.0);
		      else
			Dmx = min(CentralValue - this->Temp(i+1,j+2)
				  + amx / 2.0, 0.0);
		      if (abs(bx)>abs(apx))
			Dpx = max(this->Temp(i+3,j+2) - CentralValue
				  - apx / 2.0, 0.0);
		      else
			Dpx = max(this->Temp(i+3,j+2) - CentralValue
				  - bx / 2.0, 0.0);

		      if (abs(amy)>abs(by))
			Dmy = min(CentralValue - this->Temp(i+2,j+1)
				  + by / 2.0, 0.0);
		      else
			Dmy = min(CentralValue - this->Temp(i+2,j+1)
				  + amy / 2.0, 0.0);
		      if (abs(by)>abs(apy))
			Dpy = max(this->Temp(i+2,j+3) - CentralValue
				  - apy / 2.0, 0.0);
		      else
			Dpy = max(this->Temp(i+2,j+3) - CentralValue
				  - by / 2.0, 0.0);
		      
		      PhiValues(i, j) -=
			Delta_t * Speed
			* sqrt( DiscrX * ( Dmx * Dmx + Dpx * Dpx  )
				+ DiscrY * ( Dmy * Dmy + Dpy * Dpy) );
		    }
		  
		}
	  
	    }

	}

    this->NeedInitializationFlag = false;

    // If a 'Barrier' point sign has changed, the front has reached
    // the barrier; then, a reinitialization should be performed.
    for (i=0; (i<this->Barrier.GetM()) && (!this->NeedInitializationFlag);
	 i++)
      this->NeedInitializationFlag = (this->Temp(this->Barrier(i,0) + 2,
						 this->Barrier(i,1) + 2)
				      * PhiValues(this->Barrier(i,0),
						  this->Barrier(i,1)) < 0);
    
#ifdef MULTIVAC_CHECK_INSTABILITY
    // If a 'OutSpace' point sign has changed, the front has reached
    // the outspace; then, an instability has occurred.
    for (i=0; i<this->OutSpace.GetM(); i++)
      if (this->Temp(this->OutSpace(i,0) + 2, this->OutSpace(i,1) + 2) *
	  PhiValues(this->OutSpace(i,0), this->OutSpace(i,1)) < 0)
	throw CError_Instability(string("CNarrowBandEno2EngquistOsher::")
				 + "UpdateLevelSet",
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
		this->Temp(i+2, j+2) = PhiValues(i, j);
	    }
	  
	}
    for (i=2; i<NX-4; i++)
      {
	this->Temp(i, 1) = PhiValues(i - 2, 0);
	this->Temp(i, NY-2) = PhiValues(i - 2, Ny - 2);
      }
    for (j=2; j<NY-4; j++)
      {
	this->Temp(1, j) = PhiValues(0, j - 2);
	this->Temp(NX-2, j) = PhiValues(Nx - 2, j - 2);
      }

  }


}  // namespace Multivac.


#define FILE_UPDATER_NARROWBANDENO2ENGQUISTOSHER_CXX
#endif
