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


#ifndef FILE_NODES_HXX


#include "errors.cxx"
#include <cstdio>


namespace Multivac
{

  //////////////////////
  // FASTMARCHINGNODE //
  //////////////////////


  template <class T>
  class CFastMarchingNode
  {


    /************************
     * TYPEDEF DECLARATIONS *
     ************************/

  public:

    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;


    /**************
     * ATTRIBUTES *
     **************/

  protected:

    T time;
    int X;
    int Y;


    /*****************************
     * CONSTRUCTORS & DESTRUCTOR *
     *****************************/

  public:

    CFastMarchingNode()  throw();
    CFastMarchingNode(T time_, int X_, int Y_)  throw();

    ~CFastMarchingNode()  throw();


    /***********
     * METHODS *
     ***********/

  public:
    
    void Init(T time_, int X_, int Y_);
    value_type GetValue() const;
    int GetX() const;
    int GetY() const;
    void SetValue(value_type time_);
    

  };  // CFastMarchingNode.


}  //  namespace Multivac.


#define FILE_NODES_HXX
#endif
