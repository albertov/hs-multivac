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


#ifndef FILE_SPEEDFUNCTIONS_FIREMODEL_HXX


#include "../errors.cxx"
#include <cstdio>


namespace Multivac
{


  ////////////////
  // CFIREMODEL //
  ////////////////

  //! The speed rate is given by the model proposed Fendell and Wolff.
  /*!
    The speed rate is defined as follows:
    \par
    \par theta is the angle between the normal to the front and
    the wind direction.
    \par If |theta| < pi / 2
    \par F(U, theta) = v_f(U cos(theta)) + beta(U sin^mu(theta))
    \par else
    \par F(U, theta) = beta(U sin^mu(theta)) + epsilon(U cos^2(theta)).
    \par
    \par where v_f defines the speed rate at the head of the front,
    \par beta defines the speed rate at the flank, and
    \par epsilon is the speed rate at the rear.
  */
  template <class T>
  class CFireModel: public CSpeedFunction<T>
  {


    /**************
     * ATTRIBUTES *
     **************/

  protected:

    //! Magnitude of the wind velocity.
    T U;

    //! Parameter (for v_f).
    T m;
    //! Parameter (for v_f).
    T c_1;

    //! Parameter (for v_f, epsilon and beta).
    T epsilon_0;
    //! Parameter (for beta).
    T a;
    //! Parameter (for beta).
    T b;
    //! Parameter (for epsilon).
    T epsilon_1;


    /*****************************
     * CONSTRUCTORS & DESTRUCTOR *
     *****************************/

  public:

    CFireModel()  throw();
    CFireModel(T U_, T m_, T c_1_, T epsilon_0_,
	       T a_, T b_, T epsilon_1_)  throw();

    ~CFireModel()  throw();


    /***********
     * METHODS *
     ***********/

  public:
  
    virtual void Init(CMesh<T>& Mesh);

    virtual inline T operator() (T x, T y, T time) const;
    virtual inline T operator() (T x, T y, T time,
				 T nx, T ny, T curvature) const;

    virtual T GetMaxF1(T Xmin, T Xmax, T Ymin, T Ymax, T norm2) const;
    virtual T GetMaxF2(T Xmin, T Xmax, T Ymin, T Ymax, T norm2) const;

    virtual T GetDerivatives(T x, T y, T nx, T ny, T t,
			     T& dFdp, T& dFdx, T& dFdy,
			     T& dFdnx, T& dFdny) const;
    virtual T Get2ndDerivatives(T x, T y, T nx, T ny, T t,
				T& dFdpdp, T& dFdpdx, T& dFdpdy,
				T& dFdpdnx, T& dFdpdny,
				T& dFdxdx, T& dFdxdy,
				T& dFdxdnx, T& dFdxdny,
				T& dFdydy, T& dFdydnx,
				T& dFdydny, T& dFdnxdnx,
				T& dFdnxdny, T& dFdnydny) const;

    // For optimization.
    void SetRate(T new_parameter)
    {
      U = new_parameter;
    }

  private:
    
    T Model(T, T, T) const;

  };  // CFireModel.


}  // namespace Multivac.


#define FILE_MESHES_FIREMODEL_HXX
#endif
