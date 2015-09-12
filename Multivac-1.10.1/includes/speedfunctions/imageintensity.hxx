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


#ifndef FILE_SPEEDFUNCTIONS_IMAGEINTENSITY_HXX


#include <cstdio>


namespace Multivac
{


  ///////////////////////////
  // CSEGMENTATIONVELOCITY //
  ///////////////////////////

  //! This class provides image intensities.
  template <class T>
  class CImageIntensity: public CSpeedFunction<T>
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

    //! Path to image intensity.
    string image_intensity_file;
    
    //! Intensity threshold.
    T threshold;
    
  public:
    Matrix<T> image_intensity;


    /*****************************
     * CONSTRUCTORS & DESTRUCTOR *
     *****************************/

  public:

    CImageIntensity()  throw();
    CImageIntensity(T x_min_, T delta_x_, int Nx_,
		    T y_min_, T delta_y_, int Ny_,
		    string image_intensity_file_,
		    T threshold_)  throw();
    ~CImageIntensity()  throw();


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

  };  // CImageIntensity.


}  // namespace Multivac.


#define FILE_MESHES_IMAGEINTENSITY_HXX
#endif
