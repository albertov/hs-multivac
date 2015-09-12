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


#ifndef FILE_SPEEDFUNCTIONS_FIREMODEL_CXX

#include "firemodel.hxx"

#ifndef _limit
#define _limit 0.000001
#endif

namespace Multivac
{



  ///////////////////////////////
  // CONSTRUCTORS & DESTRUCTOR //
  ///////////////////////////////


  //! Default constructor.
  /*! \warning Nothing is defined.
   */
  template <class T>
  CFireModel<T>::CFireModel()  throw()
  {

    this->dependence_position = false;
    this->dependence_time = false;
    this->dependence_normal = true;
    this->dependence_curvature = false;

  }


  //! Main constructor.
  /*! Initializes the object with given parameters.
    \param U_ magnitude of wind velocity (for v_f).
    \param m_ parameter (for v_f).
    \param c_1_ parameter (for v_f).
    \param epsilon_0_ parameter (for v_f, epsilon and beta).
    \param a_ parameter (for beta).
    \param b_ parameter (for beta).
    \param epsilon_1_ parameter (for epsilon).
    \param CFL_ CFL number.
    \param coeff_ correction to the CFL number.
  */
  template <class T>
  CFireModel<T>::CFireModel(T U_, T m_, T c_1_, T epsilon_0_,
			    T a_, T b_, T epsilon_1_)  throw()
  {

    this->dependence_position = false;
    this->dependence_time = false;
    this->dependence_normal = true;
    this->dependence_curvature = false;

    U = U_;

    m = m_;
    c_1 = c_1_;

    epsilon_0 = epsilon_0_;

    a = a_;
    b = b_;
    epsilon_1 = epsilon_1_;

  }


  //! Destructor.
  template <class T>
  CFireModel<T>::~CFireModel()  throw()
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
  inline void CFireModel<T>::Init(CMesh<T>& Mesh)
  {
    
    // Memory allocation.
    this->Values.Reallocate(Mesh.GetNx(), Mesh.GetNy());
    this->Values.Fill(0.0);

  }


  //! Returns speed rate at some point.
  /*!
    \param x abscissa.
    \param y ordinate.
    \return The constant speed rate.
    \warning Undefined function (because the speed function
    cannot be computed from only x and y).
  */
  template <class T>
  inline T CFireModel<T>::operator() (T x, T y, T time) const
  {
    
    throw CError_Undefined
      ("T CFireModel<T>::operator() (T x, T y, T time) const",
       "The velocity depends on the normal");

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
  inline T CFireModel<T>::operator() (T x, T y, T time,
				      T nx, T ny, T curvature) const
  {
    
    return Model(nx, ny, time);

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
  inline T CFireModel<T>::GetMaxF1(T DxMin, T DxMax,
				   T DyMin, T DyMax,
				   T norm2) const
  {
    
    return c_1 * (m + 1.) * sqrt(U);
    
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
  inline T CFireModel<T>::GetMaxF2(T DxMin, T DxMax,
				   T DyMin, T DyMax,
				   T norm2) const
  {
    
    return c_1 * (m + 1.) * sqrt(U);
    
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
    \return The speed rate at (x, y).
  */
  template <class T>
  inline T CFireModel<T>::GetDerivatives(T x, T y, T nx, T ny, T t,
					 T& dFdp, T& dFdx, T& dFdy,
					 T& dFdnx, T& dFdny) const
  {
    
    // For U.
    if (nx>0.0)
      dFdp = c_1 / (2.0 * sqrt(U)) * pow(nx, 1.5);
    else
      dFdp = 0.0;

    dFdx = 0.0;
    dFdy = 0.0;

    dFdnx = 0.0; dFdny = 0.0;

    return ( Model(nx, ny, 0.0) );

  }


  /*! Returns speed rate and second derivatives at some point.
    \param x abscissa.
    \param y ordinate.
    \param nx normal fisrt-coordinate.
    \param ny normal second-coordinate.
    \param t date.
    \return The speed rate at (x, y).
  */
  template <class T>
  inline T CFireModel<T>::Get2ndDerivatives(T x, T y, T nx, T ny, T t,
					    T& dFdpdp, T& dFdpdx, T& dFdpdy,
					    T& dFdpdnx, T& dFdpdny,
					    T& dFdxdx, T& dFdxdy,
					    T& dFdxdnx, T& dFdxdny,
					    T& dFdydy, T& dFdydnx,
					    T& dFdydny, T& dFdnxdnx,
					    T& dFdnxdny, T& dFdnydny) const
  {

    // For U.
    if (nx>0.0)
      dFdpdp = - c_1 / (4.0 * U * sqrt(U)) * pow(nx, 1.5);
    else
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

    return ( Model(nx, ny, 0.0) );

  }


  /*******************
   * PRIVATE METHODS *
   *******************/


  //! Returns the speed rate as prescribed by the model.
  /*!
    \param costheta cos(theta) where theta is the angle between the normal
    to the front and the wind direction.
    \param sintheta sin(theta) where theta is the angle between the normal
    to the front and the wind direction.
    \param CurrentTime time.
    \return The speed rate for in direction (costheta, sintheta)
    at time CurrentTime.
  */
  template <class T>
  inline T CFireModel<T>::Model(T costheta, T sintheta, T CurrentTime) const
  {

    // Speed function:
    // F(theta) = epsilon_0 * sin2 + a * U * sin2 * Exp(-b * U * sin2)
    //            + G(theta)
    // where G(theta) is defined as follows:
    // If |theta|<pi/2 then G(theta) = epsilon_0 * cos2
    //                                 + c_1 * srqt(U) * costheta^(m)
    // else G(theta) = epsilon_0 * cos2 * Exp(-epsilon_1 * U * cos2)
    
    T cos2 = costheta * costheta;
    T sin2 = sintheta * sintheta;

    if (costheta < 0)
      return ( epsilon_0 * sin2 + a * U * sin2 * exp(-b * U * sin2)
	       + epsilon_0 * cos2 * exp(-epsilon_1 * U * cos2) );
    else if (m!=1.5)
      return ( epsilon_0 + c_1 * sqrt(U) * pow(costheta, m)
	       + a * U * sin2 * exp(-b * U * sin2) );
    else
      return ( epsilon_0 + c_1 * sqrt(U * costheta) * costheta
	       + a * U * sin2 * exp(-b * U * sin2) );

  }


}  // namespace Multivac.


#define FILE_SPEEDFUNCTIONS_FIREMODEL_CXX
#endif
