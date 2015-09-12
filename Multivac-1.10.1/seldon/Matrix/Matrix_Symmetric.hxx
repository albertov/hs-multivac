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

#ifndef SELDON_FILE_MATRIX_SYMMETRIC_HXX

#include "../Common/Common.hxx"
#include "../Common/Properties.hxx"
#include "../Common/Storage.hxx"
#include "../Common/Errors.cxx"
#include "../Common/Allocator.hxx"

namespace Seldon
{
  

  //! Symmetric matrix stored in a full matrix.
  template <class T, class Prop, class Storage,
	    class Allocator = SELDON_DEFAULT_ALLOCATOR<T> >
  class Matrix_Symmetric: public Spacetown, public Matrix_Base<T, Allocator>
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
    pointer* me_;

    // Methods.
  public:
    // Constructor.
    Matrix_Symmetric();
    Matrix_Symmetric(int i, int j);

    // Destructor.
    ~Matrix_Symmetric();
    void Clear();

    // Basic methods.
    int GetDataSize() const;

    // Memory management.
    void Reallocate(int i, int j);
    void SetData(int i, int j, pointer data);
    void Nullify();

    // Element access and affectation.
    value_type operator() (int i, int j);
    value_type operator() (int i, int j) const;
    const_reference Val(int i, int j) const;
    reference Val(int i, int j);
    reference operator[] (int i);
    const_reference operator[] (int i) const;
    Matrix_Symmetric<T, Prop, Storage, Allocator>&
    operator= (const Matrix_Symmetric<T, Prop, Storage, Allocator>& A);
    void Copy(const Matrix_Symmetric<T, Prop, Storage, Allocator>& A);

    // Convenient functions.
    void Zero();
    void SetIdentity();
    void Fill();
    template <class T0>
    void Fill(const T0& x);
    template <class T0>
    Matrix_Symmetric<T, Prop, Storage, Allocator>& operator= (const T0& x);
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


  //! Column-major symmetric full-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColSym, Allocator>:
    public Matrix_Symmetric<T, Prop, ColSym, Allocator>
  {
  public:
    Matrix()  throw();
    Matrix(int i, int j);

    template <class T0>
    Matrix<T, Prop, ColSym, Allocator>& operator= (const T0& x);
  };


  //! Row-major symmetric full-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowSym, Allocator>:
    public Matrix_Symmetric<T, Prop, RowSym, Allocator>
  {
  public:
    Matrix()  throw();
    Matrix(int i, int j);

    template <class T0>
    Matrix<T, Prop, RowSym, Allocator>& operator= (const T0& x);
  };


} // namespace Seldon.

#define SELDON_FILE_MATRIX_SYMMETRIC_HXX
#endif
