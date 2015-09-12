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


#ifndef FILE_MESHES_BASECLASS_CXX


#include "baseclass.hxx"
#include <cstdio>


namespace Multivac
{



  ///////////////////////////////
  // CONSTRUCTORS & DESTRUCTOR //
  ///////////////////////////////


  //! Default constructor.
  template <class T>
  CMesh<T>::CMesh()  throw()
  {

    Xmin = 0.0;
    Xmax = 0.0;
    Ymin = 0.0;
    Ymax = 0.0;

    Delta_x = 0;
    Nx = 0;
    Delta_y = 0;
    Ny = 0;

  }


  //! Main constructor.
  /*!
    \param Xmin_ minimum abscissa.
    \param Xmax_ maximum abscissa.
    \param Ymin_ minimum ordinate.
    \param Ymax_ maximum ordinate.
  */
  template <class T>
  CMesh<T>::CMesh(T Xmin_, T Xmax_, T Ymin_, T Ymax_)  throw()
  {

    Xmin = Xmin_;
    Xmax = Xmax_;
    Ymin = Ymin_;
    Ymax = Ymax_;

    Delta_x = 0;
    Nx = 0;
    Delta_y = 0;
    Ny = 0;

  }


  //! Destructor.
  template <class T>
  CMesh<T>::~CMesh()  throw()
  {

  }



  /////////////
  // METHODS //
  /////////////


  //! Returns the minimum abscissa.
  /*!
    \return Minimum abscissa.
  */
  template <class T>
  inline T CMesh<T>::GetXmin() const
  {

    return Xmin;

  }


  //! Sets the minimum abscissa.
  /*!
    \param Xmin_ the new value of Xmin.
  */
  template <class T>
  inline void CMesh<T>::SetXmin(T Xmin_)
  {

    Xmin = Xmin_;

  }


  //! Returns the maximum abscissa.
  /*!
    \return Maximum abscissa.
  */
  template <class T>
  inline T CMesh<T>::GetXmax() const
  {

    return Xmax;

  }


  //! Sets the maximum abscissa.
  /*!
    \param XMax_ the new value of Xmax.
  */
  template <class T>
  inline void CMesh<T>::SetXmax(T Xmax_)
  {

    Xmax = Xmax_;

  }


  //! Returns the minimum ordinate.
  /*!
    \return Minimum ordinate.
  */
  template <class T>
  inline T CMesh<T>::GetYmin() const
  {

    return Ymin;

  }


  //! Sets the minimum ordinate.
  /*!
    \param Ymin_ the new value of Ymin.
  */
  template <class T>
  inline void CMesh<T>::SetYmin(T Ymin_)
  {

    Ymin = Ymin_;

  }


  //! Returns the maximum ordinate.
  /*!
    \return Maximum ordinate.
  */
  template <class T>
  inline T CMesh<T>::GetYmax() const
  {

    return Ymax;

  }


  //! Sets the maximum ordinate.
  /*!
    \param YMax_ the new value of Ymax.
  */
  template <class T>
  inline void CMesh<T>::SetYmax(T Ymax_)
  {

    Ymax = Ymax_;

  }


  //! Returns the spacestep along (x'x).
  /*!
    \return the spacestep along (x'x).
  */
  template <class T>
  inline T CMesh<T>::GetDelta_x() const
  {

    return Delta_x;

  }


  //! Sets the spacestep along (x'x).
  /*!
    \param the new spacestep along (x'x).
  */
  template <class T>
  inline void CMesh<T>::SetDelta_x(T Delta_x_)
  {

    Delta_x = Delta_x_;

  }


  //! Returns the number of points along (x'x).
  /*!
    \return the number of points along (x'x).
  */
  template <class T>
  inline int CMesh<T>::GetNx() const
  {

    return Nx;

  }


  //! Returns the space step along (y'y).
  /*!
    \return the space step along (y'y).
  */
  template <class T>
  inline T CMesh<T>::GetDelta_y() const
  {

    return Delta_y;

  }


  //! Sets the spacestep along (y'y).
  /*!
    \param the new spacestep along (y'y).
  */
  template <class T>
  inline void CMesh<T>::SetDelta_y(T Delta_y_)
  {

    Delta_y = Delta_y_;

  }


  //! Returns the number of points along (y'y).
  /*!
    \return the number of points along (y'y).
  */
  template <class T>
  inline int CMesh<T>::GetNy() const
  {

    return Ny;

  }


  //! Returns the closest "upper" mesh-point to (x, y).
  /*!
    Let B = (x_B, y_B) be the closest mesh-point to A = (x, y) such that
    (x_B > x) and (y_B > y). On exit, (x, y) = B.
  */
  template <class T>
  inline void CMesh<T>::GetClosestUpperPoint(T& x, T& y)
  {

    T pos_real;

    pos_real = (x - Xmin) / Delta_x;
    if (pos_real != int(pos_real))
      x = (int(pos_real)+1) * Delta_x;

    pos_real = (y - Xmin) / Delta_y;
    if (pos_real != int(pos_real))
      y = (int(pos_real)+1) * Delta_y;

  }


}  // namespace Multivac.


#define FILE_MESHES_BASECLASS_CXX
#endif
