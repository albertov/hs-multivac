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

#ifndef SELDON_FILE_MATRIX_HERMPACKED_HXX

#include "../Common/Common.hxx"
#include "../Common/Properties.hxx"
#include "../Common/Storage.hxx"
#include "../Common/Errors.cxx"
#include "../Common/Allocator.hxx"

namespace Seldon
{
  

  //! Hermitian packed matrix class.
  template <class T, class Prop, class Storage,
	    class Allocator = SELDON_DEFAULT_ALLOCATOR<T> >
  class Matrix_HermPacked: public Spacetown, public Matrix_Base<T, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;

    // Attributes.
  protected:

    // Methods.
  public:
    // Constructor.
    Matrix_HermPacked();
    Matrix_HermPacked(int i, int j = 0);

    // Destructor.
    ~Matrix_HermPacked();
    void Clear();

    // Basic methods.
    int GetDataSize() const;

    // Memory management.
    void Reallocate(int i, int j);
    void SetData(int i, int j, pointer data);

    // Element access and affectation.
    value_type operator() (int i, int j);
    value_type operator() (int i, int j) const;
    reference Val(int i, int j);
    const_reference Val(int i, int j) const;
    reference operator[] (int i);
    const_reference operator[] (int i) const;
    Matrix_HermPacked<T, Prop, Storage, Allocator>&
    operator= (const Matrix_HermPacked<T, Prop, Storage, Allocator>& A);
    void Copy(const Matrix_HermPacked<T, Prop, Storage, Allocator>& A);

    // Convenient functions.
    void Zero();
    void SetIdentity();
    void Fill();
    template <class T0>
    void Fill(const T0& x);
    template <class T0>
    Matrix_HermPacked<T, Prop, Storage, Allocator>& operator= (const T0& x);
    void FillRand();
    void Print() const;
    void Print(int a, int b, int m, int n) const;
    void Print(int l) const;

    // Norms.
    value_type GetNormInf() const;

    // Input/output functions.
    void Write(string FileName) const;
    void Write(ofstream& FileStream) const;
    void WriteText(string FileName) const;
    void WriteText(ofstream& FileStream) const;
    void Read(string FileName);
    void Read(ifstream& FileStream);

  };


  //! Column-major hermitian packed matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColHermPacked, Allocator>:
    public Matrix_HermPacked<T, Prop, ColHermPacked, Allocator>
  {
  public:
    Matrix();
    Matrix(int i, int j = 0);

    template <class T0>
    Matrix<T, Prop, ColHermPacked, Allocator>& operator= (const T0& x);
  };


  //! Row-major hermitian packed matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowHermPacked, Allocator>:
    public Matrix_HermPacked<T, Prop, RowHermPacked, Allocator>
  {
  public:
    Matrix();
    Matrix(int i, int j = 0);

    template <class T0>
    Matrix<T, Prop, RowHermPacked, Allocator>& operator= (const T0& x);
  };


} // namespace Seldon.

#define SELDON_FILE_MATRIX_HERMPACKED_HXX
#endif
