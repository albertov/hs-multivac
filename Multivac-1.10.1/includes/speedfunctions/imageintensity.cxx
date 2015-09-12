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


#ifndef FILE_SPEEDFUNCTIONS_IMAGEINTENSITY_CXX

#include "imageintensity.hxx"


namespace Multivac
{



  ///////////////////////////////
  // CONSTRUCTORS & DESTRUCTOR //
  ///////////////////////////////


  //! Default constructor.
  /*! \warning Nothing is defined.
   */
  template <class T>
  CImageIntensity<T>::CImageIntensity()  throw()
  {
    this->dependence_position = true;
    this->dependence_time = false;
    this->dependence_normal = false;
    this->dependence_curvature = false;
  }
  
  //! Main constructor.
  template <class T>
  CImageIntensity<T>
  ::CImageIntensity(T x_min_, T delta_x_, int Nx_,
		    T y_min_, T delta_y_, int Ny_,
		    string image_intensity_file_,
		    T threshold_)  throw()
  {
    this->dependence_position = true;
    this->dependence_time = false;
    this->dependence_normal = false;
    this->dependence_curvature = false;
    
    x_min = x_min_;
    delta_x = delta_x_;
    Nx = Nx_;
    y_min = y_min_;
    delta_y = delta_y_;
    Ny = Ny_;

    image_intensity_file = image_intensity_file_;

    threshold = threshold_;
  }



  //! Destructor.
  template <class T>
  CImageIntensity<T>::~CImageIntensity()  throw()
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
  inline void CImageIntensity<T>::Init(CMesh<T>& Mesh)
  {
    int i, j;

    // Reads image intensity.
    int size = Nx * Ny;
    double* data = new double[size];
    ifstream data_stream(image_intensity_file.c_str());
    i = 0;
    while (data_stream.good() && i < size)
      {
	data_stream.read(reinterpret_cast<char*>(&data[i]), sizeof(double));
	i++;
      }
    if (i != size)
      throw string("Not enough bytes in \"")
	+ image_intensity_file + "\".";

    // Puts it in matrix.
    image_intensity.Reallocate(Nx, Ny);
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++)
	if (T(data[j * Nx + i]) < threshold)
	  image_intensity(i, j) = 0.;
	else
	  image_intensity(i, j) = 1.;
    
    delete[] data;

    // Memory allocation.
    this->Values.Reallocate(Mesh.GetNx(), Mesh.GetNy());
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
  inline T CImageIntensity<T>::operator() (T x, T y, T time) const
  {
    int i = int((x - x_min) / delta_x);
    int j = Ny - 1 - int((y - y_min) / delta_y);

    i = max(0, i);
    i = min(Nx - 1, i);

    j = max(0, j);
    j = min(Ny - 1, j);

    return image_intensity(i, j);
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
  inline T CImageIntensity<T>::operator() (T x, T y, T time,
					   T nx, T ny, T curvature) const
  {
    throw CError_Undefined
      ("T CImageIntensity<T>::operator() const");
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
  inline T CImageIntensity<T>::GetMaxF1(T DxMin, T DxMax,
					T DyMin, T DyMax,
					T norm2) const
  {
    throw CError_Undefined
      ("T CImageIntensity<T>::GetMaxF1 const");
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
  inline T CImageIntensity<T>::GetMaxF2(T DxMin, T DxMax,
					T DyMin, T DyMax,
					T norm2) const
  {
    throw CError_Undefined
      ("T CImageIntensity<T>::GetMaxF2 const");
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
  inline T CImageIntensity<T>::GetDerivatives(T x, T y, T nx, T ny, T t,
					      T& dFdp, T& dFdx, T& dFdy,
					      T& dFdnx, T& dFdny) const
  {
    throw CError_Undefined
      ("T CImageIntensity<T>::GetDerivatives const");
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
  inline T CImageIntensity<T>
  ::Get2ndDerivatives(T x, T y, T nx, T ny, T t,
		      T& dFdpdp, T& dFdpdx, T& dFdpdy,
		      T& dFdpdnx, T& dFdpdny,
		      T& dFdxdx, T& dFdxdy,
		      T& dFdxdnx, T& dFdxdny,
		      T& dFdydy, T& dFdydnx,
		      T& dFdydny, T& dFdnxdnx,
		      T& dFdnxdny, T& dFdnydny) const
  {
    throw CError_Undefined
      ("T CImageIntensity<T>::Get2ndDerivatives const");
  }


}  // namespace Multivac.


#define FILE_SPEEDFUNCTIONS_IMAGEINTENSITY_CXX
#endif
