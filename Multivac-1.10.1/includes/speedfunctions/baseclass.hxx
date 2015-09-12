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


#ifndef FILE_SPEEDFUNCTIONS_BASECLASS_HXX


#include "../errors.cxx"
#include <cstdio>


namespace Multivac
{


  ////////////////////
  // CSPEEDFUNCTION //
  ////////////////////

  //! Base class for speed functions.
  /*! Defines the speed functions interface.  All speed functions must
    be defined in the same way.
    \note
    This is an abstract class.
  */
  template <class T>
  class CSpeedFunction
  {


    /**************
     * ATTRIBUTES *
     **************/

  protected:

    //! Speed rates on grid points.
    Matrix<T> Values;

    //! Does the speed function depend upon the position?
    bool dependence_position;
    //! Does the speed function depend upon the time?
    bool dependence_time;
    //! Does the speed function depend upon the normal?
    bool dependence_normal;
    //! Does the speed function depend upon the curvature?
    bool dependence_curvature;


    /*****************************
     * CONSTRUCTORS & DESTRUCTOR *
     *****************************/

  public:

    CSpeedFunction()  throw();

    virtual ~CSpeedFunction()  throw();


    /***********
     * METHODS *
     ***********/

  public:
  
    bool IsPositionDependent() const;
    bool IsTimeDependent() const;
    bool IsNormalDependent() const;
    bool IsCurvatureDependent() const;

    virtual void Init(CMesh<T>& Mesh) = 0;

    Matrix<T>& GetValues();
    virtual T operator() (T x, T y, T time) const = 0;
    virtual T operator() (T x, T y, T time,
			  T nx, T ny, T curvature) const = 0;

    virtual T GetMaxF1(T Xmin, T Xmax, T Ymin, T Ymax, T norm2) const = 0;
    virtual T GetMaxF2(T Xmin, T Xmax, T Ymin, T Ymax, T norm2) const = 0;

    virtual T operator() (int i, int j) const;

    virtual T GetDerivatives(T x, T y, T nx, T ny, T t,
			     T& dFdp, T& dFdx, T& dFdy,
			     T& dFdnx, T& dFdny) const = 0;
    virtual T Get2ndDerivatives(T x, T y, T nx, T ny, T t,
				T& dFdpdp, T& dFdpdx, T& dFdpdy,
				T& dFdpdnx, T& dFdpdny,
				T& dFdxdx, T& dFdxdy,
				T& dFdxdnx, T& dFdxdny,
				T& dFdydy, T& dFdydnx,
				T& dFdydny, T& dFdnxdnx,
				T& dFdnxdny, T& dFdnydny) const = 0;

    virtual void Save(string FFile) const;

  };  // CSpeedFunction.


}  // namespace Multivac.


#define FILE_MESHES_BASECLASS_HXX
#endif
