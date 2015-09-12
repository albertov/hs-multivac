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


#ifndef FILE_SPEEDFUNCTIONS_CONSTANTSPEED_HXX


#include "../errors.cxx"
#include <cstdio>


namespace Multivac
{


  ////////////////////
  // CCONSTANTSPEED //
  ////////////////////

  //! The speed rate is constant in time and in space.
  template <class T>
  class CConstantSpeed: public CSpeedFunction<T>
  {


    /**************
     * ATTRIBUTES *
     **************/

  protected:

    //! The constant speed rate.
    T SpeedRate;


    /*****************************
     * CONSTRUCTORS & DESTRUCTOR *
     *****************************/

  public:

    CConstantSpeed()  throw();
    CConstantSpeed(T SpeedRate_)  throw();

    ~CConstantSpeed()  throw();


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

    virtual inline T operator() (int i, int j) const;

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

    void SetRate(T rate);
    void SetParameters(const Vector<T>& parameters);

  };  // CConstantSpeed.


}  // namespace Multivac.


#define FILE_MESHES_CONSTANTSPEED_HXX
#endif
