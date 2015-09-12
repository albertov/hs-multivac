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


#ifndef FILE_ARRAYHEAP_HXX


#include "errors.cxx"

#include <iostream>
using std::cout;
using std::endl;
#include <cstdlib>
#include <ctime>
#include <cstdio>


namespace Multivac
{

  
  ///////////////
  // ARRAYHEAP //
  ///////////////

  template <class T>
  class ArrayHeap
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
    //! Vector storing pointers to nodes.
    Vector<T*> Nodes;
    //! Number of nodes.
    int NbNodes;


    /***********
     * METHODS *
     ***********/

  public:
    //Constructors.
    ArrayHeap()  throw();
    ArrayHeap(int depth)  throw();

    // Destructor.
    ~ArrayHeap()  throw();

    // Initializations.
    void Reallocate(int length);
    void Resize(int length, int lastelement = 0);

    int Add(T* X);
    int Add(T* X, Matrix<int>& Pointers);
    int MoveUp(int XIndex);
    int MoveUp(int XIndex, Matrix<int>& Pointers);

    void DeleteRoot();
    void DeleteRoot(Matrix<int>& Pointers);

    // Convenient functions.
    T* operator() (int i);

    T* GetRoot();
    T* GetRoot() const;

    int GetNbNodes() const;
    bool IsEmpty() const;

  };  // ArrayHead.


}  // namespace Multivac.


#define FILE_ARRAYHEAP_HXX
#endif
