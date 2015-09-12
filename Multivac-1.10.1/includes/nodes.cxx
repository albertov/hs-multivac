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


#ifndef FILE_NODES_CXX


#include "nodes.hxx"


namespace Multivac
{


  ///////////////////////////////
  // CONSTRUCTORS & DESTRUCTOR //
  ///////////////////////////////


  template <class T>
  CFastMarchingNode<T>::CFastMarchingNode()  throw()
  {
    
  }
  
   
  template <class T>
  CFastMarchingNode<T>::CFastMarchingNode(T time_, int X_, int Y_)  throw():
    time(time_), X(X_), Y(Y_)
  {

  }


  template <class T>
  CFastMarchingNode<T>::~CFastMarchingNode()  throw()
  {

  }



  /////////////
  // METHODS //
  /////////////

  
  template <class T>
  void CFastMarchingNode<T>::Init(T time_, int X_, int Y_)
  {
    time=time_;
    X=X_;
    Y=Y_;
  }


  template <class T>
  typename CFastMarchingNode<T>::value_type
  CFastMarchingNode<T>::GetValue() const
  {
    return time;
  }


  template <class T>
  int CFastMarchingNode<T>::GetX() const
  {
    return X;
  }


  template <class T>
  int CFastMarchingNode<T>::GetY() const
  {
    return Y;
  }


  template <class T>
  void CFastMarchingNode<T>::SetValue(value_type time_)
  {
    time = time_;
  }


}  // namespace Multivac.


#define FILE_NODES_CXX
#endif
