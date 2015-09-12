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


#ifndef FILE_SPEEDFUNCTIONS_PIECEWISECONSTANTSPEED_CXX


#include "piecewiseconstantspeed.hxx"


namespace Multivac
{



  ///////////////////////////////
  // CONSTRUCTORS & DESTRUCTOR //
  ///////////////////////////////


  //! Default constructor.
  /*! \warning Speed rates are set to 1 and the ordinate limit to 0.
   */
  template <class T>
  CPiecewiseConstantSpeed<T>::CPiecewiseConstantSpeed()  throw()
  {

    this->dependence_position = true;
    this->dependence_time = false;
    this->dependence_normal = false;
    this->dependence_curvature = false;

    SpeedRate = T(1);
    SpeedRate0 = T(1);
    Limit = 0;
    LimitIndex = -1;

  }


  //! Main constructor.
  /*! Initializes the object with a given constant speed rate.
    \param SpeedRate_ constant speed rate in the upper band.
    \param SpeedRate0_ constant speed rate in the lower band.
    \param Limit_ ordinate of the limit between the two bands.
  */
  template <class T>
  CPiecewiseConstantSpeed<T>::
  CPiecewiseConstantSpeed(T SpeedRate_, T SpeedRate0_, T Limit_)  throw()
  {

    this->dependence_position = true;
    this->dependence_time = false;
    this->dependence_normal = false;
    this->dependence_curvature = false;

    SpeedRate = SpeedRate_;
    SpeedRate0 = SpeedRate0_;

    Limit = Limit_;
    LimitIndex = -1;

  }


  //! Destructor.
  template <class T>
  CPiecewiseConstantSpeed<T>::~CPiecewiseConstantSpeed()  throw()
  {

  }



  /////////////
  // METHODS //
  /////////////


  //! Initializes speed function: memory allocation.
  /*! Allocates memory to store speed rates.
    \param Mesh orthogonal mesh.
  */
  template <class T>
  inline void CPiecewiseConstantSpeed<T>::Init(CMesh<T>& Mesh)
  {
    
    // Memory allocation.
    this->Values.Reallocate(Mesh.GetNx(), Mesh.GetNy());
    this->Values.Fill(0.0);

  }


  //! Returns speed rate at some point.
  /*!
    \param x abscissa.
    \param y ordinate.
    \return The speed rate on (x, y).
  */
  template <class T>
  inline T CPiecewiseConstantSpeed<T>::operator() (T x, T y, T time) const
  {
    
    return (y<Limit?SpeedRate0:SpeedRate);

  }


  //! Returns speed rate at some point.
  /*!
    \param x abscissa.
    \param y ordinate.
    \param time time.
    \param nx normal fisrt-coordinate.
    \param ny normal second-coordinate.
    \param curvature curvature.
    \return The constant speed rate.
  */
  template <class T>
  inline T CPiecewiseConstantSpeed<T>::
  operator() (T x, T y, T time, T nx, T ny, T curvature) const
  {

    return (y<Limit?SpeedRate0:SpeedRate);

  }


  //! Returns an upper bound of the maximum of the first derivative of F
  //! with respect to Phi_x, multiplied by |\nabla Phi|_2.
  /*! Returns an upper bound of the maximum (absolute value) of the first
    derivative of F with respect to Phi_x, where Phi_x = d(Phi)/dx,
    multiplied by |\nabla Phi|_2, i.e. max | F_{Phi_x} |\nabla Phi|_2 |.
    Returns this bound on [DxMin, DxMax] x [DyMin, DyMax].
    \param DxMin minimum of Phi_x.
    \param DxMax maximum of Phi_x.
    \param DyMin minimum of Phi_y.
    \param DyMax maximum of Phi_y.
    \param norm2 an upper bound of |\nabla Phi|_2.
    \return An upper bound of the maximum (absolute value) of the first
    derivative of F with respect to Phi_x, multiplied by |\nabla Phi|_2.
  */
  template <class T>
  inline T CPiecewiseConstantSpeed<T>::GetMaxF1(T DxMin, T DxMax,
						T DyMin, T DyMax,
						T norm2) const
  {

    return 0.0;

  }


  //! Returns an upper bound of the maximum of the first derivative of F
  //! with respect to Phi_y, multiplied by |\nabla Phi|_2.
  /*! Returns an upper bound of the maximum (absolute value) of the first
    derivative of F with respect to Phi_y, where Phi_y = d(Phi)/dy,
    multiplied by |\nabla Phi|_2, i.e. max | F_{Phi_y} |\nabla Phi|_2 |.
    Returns this bound on [DxMin, DxMax] x [DyMin, DyMax].
    \param DxMin minimum of Phi_x.
    \param DxMax maximum of Phi_x.
    \param DyMin minimum of Phi_y.
    \param DyMax maximum of Phi_y.
    \param norm2 an upper bound of |\nabla Phi|_2.
    \return An upper bound of the maximum (absolute value) of the first
    derivative of F with respect to Phi_y, multiplied by |\nabla Phi|_2.
  */
  template <class T>
  inline T CPiecewiseConstantSpeed<T>::GetMaxF2(T DxMin, T DxMax,
						T DyMin, T DyMax,
						T norm2) const
  {

    return 0.0;

  }


  //! Returns speed rate on grid points.
  /*!
    \param i row index.
    \param j column index.
    \return The speed rate.
    \exception Multivac::CError_OutOfDomain attempt to compute outside
    the domain.
  */
  template <class T>
  inline T CPiecewiseConstantSpeed<T>::operator() (int i, int j) const
  {
    
#ifdef MULTIVAC_CHECK_CALCULATIONS_BOUNDARIES
    if (i < 0)
      throw CError_OutOfDomain(string("CPiecewiseConstantSpeed<T>::")
			       + "operator() (int, int)",
			       "Row index is strictly negative");
    if (i >= this->Values.GetM())
      throw CError_OutOfDomain(string("CPiecewiseConstantSpeed<T>::")
			       + "operator() (int, int)",
			       "Row index is too high");
    if (j < 0)
      throw CError_OutOfDomain(string("CPiecewiseConstantSpeed<T>::")
			       + "operator() (int, int)",
			       "Column index is strictly negative");
    if (j >= this->Values.GetN())
      throw CError_OutOfDomain(string("CPiecewiseConstantSpeed<T>::")
			       + "operator() (int, int)",
			       "Column index is too high");
#endif

    return (j<LimitIndex?SpeedRate0:SpeedRate);

  }


  /*! Returns speed rate and derivatives at some point.
    \param x abscissa.
    \param y ordinate.
    \param nx normal fisrt-coordinate.
    \param ny normal second-coordinate.
    \param t date.
    \param dFdp partial derivative of F with respect to p.
    \param dFdx partial derivative of F with respect to x.
    \param dFdy partial derivative of F with respect to y.
    \param dFdnx partial derivative of F with respect to nx.
    \param dFdny partial derivative of F with respect to ny.
    \return The speed rate on (x, y).
  */
  template <class T>
  inline T CPiecewiseConstantSpeed<T>::
  GetDerivatives(T x, T y, T nx, T ny, T t,
		 T& dFdp, T& dFdx, T& dFdy,
		 T& dFdnx, T& dFdny) const
  {
    
    dFdp = 1.0;
    dFdx = 0.0; dFdy = 0.0;
    dFdnx = 0.0; dFdny = 0.0;

    return (y<Limit?SpeedRate0:SpeedRate);

  }


  /*! Returns speed rate and second derivatives at some point.
    \param x abscissa.
    \param y ordinate.
    \param nx normal fisrt-coordinate.
    \param ny normal second-coordinate.
    \param t date.
    \return The speed rate on (x, y).
  */
  template <class T>
  inline T CPiecewiseConstantSpeed<T>::
  Get2ndDerivatives(T x, T y, T nx, T ny, T t,
		    T& dFdpdp, T& dFdpdx, T& dFdpdy,
		    T& dFdpdnx, T& dFdpdny,
		    T& dFdxdx, T& dFdxdy,
		    T& dFdxdnx, T& dFdxdny,
		    T& dFdydy, T& dFdydnx,
		    T& dFdydny, T& dFdnxdnx,
		    T& dFdnxdny, T& dFdnydny) const
  {

    dFdpdp = 0.0;
    dFdpdx = 0.0;
    dFdpdy = 0.0;
    dFdpdnx = 0.0;
    dFdpdny = 0.0;
    dFdxdx = 0.0;
    dFdxdy = 0.0;
    dFdxdnx = 0.0;
    dFdxdny = 0.0;
    dFdydy = 0.0;
    dFdydnx = 0.0;
    dFdydny = 0.0;
    dFdnxdnx = 0.0;
    dFdnxdny = 0.0;
    dFdnydny = 0.0;

    return (y<Limit?SpeedRate0:SpeedRate);

  }


}  // namespace Multivac.


#define FILE_SPEEDFUNCTIONS_PIECEWISECONSTANTSPEED_CXX
#endif
