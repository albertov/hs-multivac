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


#ifndef FILE_SPEEDFUNCTIONS_GRADIENTSEGMENTATION_HXX


#include <cstdio>


namespace Multivac
{


  ///////////////////////////
  // CGRADIENTSEGMENTATION //
  ///////////////////////////

  //! This class provides a velocity in order to performe image segmentation.
  template <class T>
  class CGradientSegmentation: public CSpeedFunction<T>
  {


    /**************
     * ATTRIBUTES *
     **************/

		
  protected:
    //! Lower-left corner abscissa.
    T x_min;
    //! Step along x.
    T delta_x;
    //! Number of points along x.
    int Nx;
    //! Lower-left corner ordinate.
    T y_min;
    //! Step along y.
    T delta_y;
    //! Number of points along y.
    int Ny;

    //! Path to image gradient.
    string image_gradient_file;
    //! Curvature scale factor.
    T epsilon_c;
    //! Velocity scale.
    T scale;
    //! Velocity slow-down scale.
    T slow_down_scale;
    
  public:
    Matrix<T> image_gradient;


    /*****************************
     * CONSTRUCTORS & DESTRUCTOR *
     *****************************/

  public:

    CGradientSegmentation()  throw();
    CGradientSegmentation(T x_min_, T delta_x_, int Nx_,
			  T y_min_, T delta_y_, int Ny_,
			  string image_gradient_file_,
			  T epsilon_c_, T scale_,
			  T slow_down_scale_)  throw();
    ~CGradientSegmentation()  throw();


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

  };  // CGradientSegmentation.


}  // namespace Multivac.


#define FILE_MESHES_GRADIENTSEGMENTATION_HXX
#endif
