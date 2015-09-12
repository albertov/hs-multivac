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


#ifndef FILE_UPDATER_CHANVESE_CXX


#include "chanvese.hxx"


namespace Multivac
{



  ///////////////////////////////
  // CONSTRUCTORS & DESTRUCTOR //
  ///////////////////////////////


  //! Default constructor.
  template <class T>
  CChanVese<T>::CChanVese()  throw()
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
    \param inside_ weight of the inside integral.
    \param outside_ weight of the outside integral.
    \param mu_ weight of the penalty term on the length.
    \param nu_ additional velocity.
    \param Dirac_threshold_ scale factor applied to the smoothed Dirac
    function. It should be above 1.
  */
  template <class T>
  CChanVese<T>::CChanVese(int TubeSemiWidth_,
			  int BarrierWidth_,
			  int OutSpaceWidth_,
			  T inside_, T outside_,
			  T mu_, T nu_, T Dirac_threshold_)  throw()
  {

    this->offset = 2;

    this->TubeSemiWidth = TubeSemiWidth_;
    this->BarrierWidth = BarrierWidth_;
    this->OutSpaceWidth = OutSpaceWidth_;

    this->NeedSpeedUpdateFlag = true;
    this->NeedInitializationFlag = false;

    this->TMax = 0.;

    mu = mu_;
    nu = nu_;
    inside = inside_;
    outside = outside_;

    Dirac_threshold = Dirac_threshold_;
    pi_threshold = 3.1415926535898 / Dirac_threshold;
    two_threshold = 2. * Dirac_threshold;

  }


  //! Destructor.
  template <class T>
  CChanVese<T>::~CChanVese()  throw()
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
  inline bool CChanVese<T>::IsNarrowBand() const
  {

    return true;

  }


  //! Should this updater be used for the fast marching method?
  /*! \return 'false' because the updater cannot be used for the fast
    marching method.
  */
  template <class T>
  inline bool CChanVese<T>::IsFastMarching() const
  {

    return false;

  }


  //! Inits the updater.
  /*! The matrix 'Temp' is allocated.
    \param Mesh orthogonal mesh.
    \param Phi level set function defined on an orthogonal mesh.
  */
  template <class T>
  void CChanVese<T>::Init(CMesh<T>& Mesh,
			  CLevelSet<T>& Phi)
  {

    // 'Temp' contains more points than the mesh because it
    // contains ghost cells.
    this->Temp.Reallocate(Mesh.GetNx() + 4, Mesh.GetNy() + 4);

    // Rescales the Dirac threshold onto the grid.
    Dirac_threshold *= min(Mesh.GetDelta_x(), Mesh.GetDelta_y());

  }


  //! Updates the level set function Phi.
  /*! This function updates the level set function Phi on the tube,
    according to the Chan-Vese algorithm.
    \param Delta_t time step.
    \param Mesh orthogonal mesh.
    \param F speed function defined on Mesh.
    \param Phi level set function defined on Mesh.
    \param CurrentTime current time.
  */
  template <class T>
  void CChanVese<T>::
  UpdateLevelSet(T Delta_t, CMesh<T>& Mesh, CSpeedFunction<T>& F,
		 CLevelSet<T>& Phi, T CurrentTime)
  {

    int i, j;

    // Phi values.
    Matrix<T>& PhiValues = Phi.GetValues();

    // F values.
    Matrix<T>& FValues = F.GetValues();

    T Delta_x = Mesh.GetDelta_x();
    T Delta_y = Mesh.GetDelta_y();

    int Nx = Mesh.GetNx();
    int Ny = Mesh.GetNy();
    int NX = Mesh.GetNx() + 2;
    int NY = Mesh.GetNy() + 2;

    T Dpx, Dcx, Dpy, Dcy;

    T Speed, CentralValue;

    /*** Computes inside and outside integrals ***/

    // Inside and outside integrals.
    T c1(0.), c2(0.);
    T count1(0.), count2(0.);

    for (i=0; i<Nx; i++)
      for (j=0; j<Ny; j++)
	if (PhiValues(i, j) > 0.)
	  {
	    c1 += FValues(i, j);
	    count1 += 1;
	  }
	else
	  {
	    c2 += FValues(i, j);
	    count2 += 1;
	  }
    c1 /= count1;
    c2 /= count2;

    /*** Updates are performed on tube points ***/

    T grad_cx, grad_cy;
    T grad_mx, grad_my;
    T grad_x, grad_y;
    T diff1, diff2;
    bool FirstLoop;

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

		  if (mu != 0.)
		    {ERR(Ici);
		    /*** Central derivatives ***/
		  
		    CentralValue = this->Temp(i+2, j+2);

		    Dpx = ( this->Temp(i+3, j+2) - CentralValue ) / Delta_x;
		    Dcx = ( this->Temp(i+3, j+2) -  this->Temp(i+1, j+2) )
		      / (2. * Delta_x);
		    Dpy = ( this->Temp(i+2, j+3) - CentralValue ) / Delta_y;
		    Dcy = ( this->Temp(i+2, j+3) -  this->Temp(i+2, j+1) )
		      / (2. * Delta_y);

		    grad_cx = sqrt(Dpx * Dpx + Dcy * Dcy);
		    if (grad_cx != 0.)
		      grad_cx = Dpx / grad_cx;
		    grad_cy = sqrt(Dcx * Dcx + Dpy * Dpy);
		    if (grad_cy != 0.)
		      grad_cy = Dpy / grad_cy;

		    /*** Along x ***/

		    CentralValue = this->Temp(i+1, j+2);

		    Dpx = ( this->Temp(i+2, j+2) - CentralValue ) / Delta_x;
		    Dcy = ( this->Temp(i+1, j+3) -  this->Temp(i+1, j+1) )
		      / (2. * Delta_y);

		    grad_mx = sqrt(Dpx * Dpx + Dcy * Dcy);
		    if (grad_mx != 0.)
		      grad_mx = Dpx / grad_mx;

		    /*** Along y ***/
		  
		    CentralValue = this->Temp(i+2, j+1);

		    Dcx = ( this->Temp(i+3, j+1) -  this->Temp(i+1, j+1) )
		      / (2. * Delta_x);
		    Dpy = ( this->Temp(i+2, j+1) - CentralValue ) / Delta_y;

		    grad_my = sqrt(Dcx * Dcx + Dpy * Dpy);
		    if (grad_my != 0.)
		      grad_my = Dpy / grad_my;

		    /*** Total ***/

		    grad_x = 0.;
		    if (grad_cx - grad_mx != 0)
		      grad_x = mu / (Delta_x * Delta_x)
			* (grad_cx - grad_mx);

		    grad_y = 0.;
		    if (grad_cy - grad_my != 0)
		      grad_y = mu / (Delta_y * Delta_y)
			* (grad_cy - grad_my);
		    }
		  else
		    {
		      grad_x = 0.;
		      grad_y = 0.;
		    }
		  
		  diff1 = Speed - c1;
		  diff2 = Speed - c2;
		  PhiValues(i, j) += Delta_t * Dirac(PhiValues(i, j))
		    * (grad_x + grad_y - nu - inside * diff1 * diff1
		       + outside * diff2 * diff2);
		}
	      
	    }
	  
	}

    this->NeedInitializationFlag = false;

    // If a 'Barrier' point sign has changed, the front has reached
    //  the barrier; then, a reinitialization should be performed.
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
  	throw CError_Instability(string("CChanVese") + "::UpdateLevelSet",
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
    for (i=2; i<NX; i++)
      {
	this->Temp(i, 0) = PhiValues(i - 2, 0);
	this->Temp(i, 1) = PhiValues(i - 2, 0);
	this->Temp(i, NY) = PhiValues(i - 2, Ny - 1);
	this->Temp(i, NY + 1) = PhiValues(i - 2, Ny - 1);
      }
    for (j=2; j<NY; j++)
      {
	this->Temp(0, j) = PhiValues(0, j - 2);
	this->Temp(1, j) = PhiValues(0, j - 2);
	this->Temp(NX, j) = PhiValues(Nx - 1, j - 2);
	this->Temp(NX + 1, j) = PhiValues(Nx - 1, j - 2);
      }

  }


  //! Dirac function.
  /*! Discretized Dirac function.
    \param phi input value.
  */
  template <class T>
  inline T CChanVese<T>::Dirac(T phi)
  {
    if (phi < -Dirac_threshold || phi > Dirac_threshold)
      return 0.;
    else
      return (1. + cos(phi * pi_threshold))
	/ two_threshold;
  }


}  // namespace Multivac.


#define FILE_UPDATER_CHANVESE_CXX
#endif
