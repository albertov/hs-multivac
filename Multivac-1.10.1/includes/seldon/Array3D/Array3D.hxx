// Copyright (C) 2001-2004 Vivien Mallet
//
// This file is part of Seldon library.
// Seldon library provides matrices and vectors structures for
// linear algebra.
// 
// Seldon is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// Seldon is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License (file "license") for more details.
//
// For more information, please see the Seldon home page:
//     http://spacetown.free.fr/lib/seldon/

// To be included by Seldon.hxx


#ifndef SELDON_FILE_ARRAY3D_HXX

#include "../Common/Common.hxx"
#include "../Common/Errors.cxx"
#include "../Common/Allocator.hxx"

namespace Seldon
{

  
  //! 3D array.
  /*!
    This class implements 3D arrays.
  */
  template <class T, class Allocator = SELDON_DEFAULT_ALLOCATOR<T> >
  class Array3D
  {
    // typdef declarations.
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;

    // Static attributes.
  protected:
    static Allocator array3D_allocator_;

    // Attributes.
  protected:
    // Length along dimension #1.
    int length1_;
    // Length along dimension #2.
    int length2_;
    // Length along dimension #3.
    int length3_;

    // Size of a slice (i.e. length1_ by length2_).
    int length23_;

    // Pointer to stored elements.
    pointer data_;

    // Methods.
  public:
    // Constructors.
    Array3D();
    Array3D(int i, int j, int k);
  
    // Destructor.
    ~Array3D();

    // Basic methods.
    int GetLength1() const;
    int GetLength2() const;
    int GetLength3() const;
    int GetSize() const;
    int GetDataSize() const;
    pointer GetData() const;

    // Memory management.
    void Reallocate(int i, int j, int k);

    // Element access and affectation.
    reference operator() (int i, int j, int k);
    const_reference operator() (int i, int j, int k) const;
    void Copy(const Array3D<T, Allocator>& A);

    // Convenient functions.
    void Zero();
    void Fill();
    template <class T0>
    void Fill(const T0& x);
    void FillRand();
    void Print() const;

  };


  // 3D array allocator.
  template <class T, class Allocator>
  Allocator Array3D<T, Allocator>::array3D_allocator_;


  template <class T, class Allocator>
  ostream& operator << (ostream& out,
			const Array3D<T, Allocator>& A);


} // namespace Seldon.


#define SELDON_FILE_ARRAY3D_HXX
#endif
