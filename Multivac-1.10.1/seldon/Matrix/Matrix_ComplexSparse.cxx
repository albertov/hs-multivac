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

#ifndef SELDON_FILE_MATRIX_COMPLEXSPARSE_CXX

#include "Matrix_ComplexSparse.hxx"

namespace Seldon
{


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_ComplexSparse<T, Prop, Storage, Allocator>
  ::Matrix_ComplexSparse(): Matrix_Base<T, Allocator>()
  {
    real_nz_ = 0;
    imag_nz_ = 0;
    real_ptr_ = NULL;
    imag_ptr_ = NULL;
    real_ind_ = NULL;
    imag_ind_ = NULL;
  }


  //! Constructor.
  /*!
    Builds an empty i by j sparse matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_ComplexSparse<T, Prop, Storage, Allocator>
  ::Matrix_ComplexSparse(int i, int j): Matrix_Base<T, Allocator>(i, j)
  {
    real_nz_ = 0;
    imag_nz_ = 0;
    real_ptr_ = NULL;
    imag_ptr_ = NULL;
    real_ind_ = NULL;
    imag_ind_ = NULL;
  }


  //! Constructor.
  /*! Builds a sparse matrix of size i by j , with real_nz
    non-zero elements in the real part of the matrix and imag_nz
    non-zero elements in the imaginary part of the matrix.
    \param i number of rows.
    \param j number of columns.
    \param real_nz number of non-zero elements in the real part.
    \param imag_nz number of non-zero elements in the imaginary part.
    \note Matrix values are not initialized. Indices of non-zero entries
    are not initialized either.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_ComplexSparse<T, Prop, Storage, Allocator>::
  Matrix_ComplexSparse(int i, int j, int real_nz, int imag_nz):
    Matrix_Base<T, Allocator>(i, j)
  {
    this->real_nz_ = real_nz;
    this->imag_nz_ = imag_nz;

#ifdef SELDON_CHECK_DIMENSIONS
    if ( (static_cast<long int>(real_nz_-1) / static_cast<long int>(j)
	  >= static_cast<long int>(i)) ||
	 (static_cast<long int>(imag_nz_-1) / static_cast<long int>(j)
	  >= static_cast<long int>(i)) )
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_ComplexSparse::")
		       + "Matrix_ComplexSparse(int, int, int, int)",
		       string("There are more values (") + to_str(real_nz)
		       + " values for the real part and " + to_str(imag_nz)
		       + string(" values for the imaginary part) than")
		       + " elements in the matrix (" + to_str(i) + " by "
		       + to_str(j) + ").");
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	real_ptr_ = reinterpret_cast<int*>( calloc(Storage::GetFirst(i, j)+1,
						   sizeof(int)) );

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_ptr_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	imag_ptr_ = 0;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_ptr_ == NULL && i != 0 && j != 0)
      throw NoMemory(string("Matrix_ComplexSparse::")
		     + "Matrix_ComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * (Storage::GetFirst(i, j)+1))
		     + " bytes to store " + to_str(Storage::GetFirst(i, j)+1)
		     + string(" row or column start indices (for the real")
		     + " part), for a "
		     + to_str(i) + " by " + to_str(j) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	imag_ptr_ = reinterpret_cast<int*>( calloc(Storage::GetFirst(i, j)+1,
						   sizeof(int)) );

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	free(real_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (imag_ptr_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	free(real_ptr_);
	real_ptr_ = 0;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (imag_ptr_ == NULL && i != 0 && j != 0)
      throw NoMemory(string("Matrix_ComplexSparse::")
		     + "Matrix_ComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * (Storage::GetFirst(i, j)+1))
		     + " bytes to store " + to_str(Storage::GetFirst(i, j)+1)
		     + string(" row or column start indices (for the")
		     + string(" imaginary part), for a ")
		     + to_str(i) + " by " + to_str(j) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	real_ind_ = reinterpret_cast<int*>( calloc(real_nz_, sizeof(int)) );
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	free(real_ptr_);
	free(imag_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_ind_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	free(real_ptr_);
	free(imag_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_ind_ == NULL && i != 0 && j != 0)
      throw NoMemory(string("Matrix_ComplexSparse::")
		     + "Matrix_ComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * real_nz)
		     + " bytes to store " + to_str(real_nz)
		     + " row or column indices (real part), for a "
		     + to_str(i) + " by " + to_str(j) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	imag_ind_ = reinterpret_cast<int*>( calloc(imag_nz_, sizeof(int)) );
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	free(real_ptr_);
	free(imag_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	free(imag_ind_);
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_ind_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	free(real_ptr_);
	free(imag_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	free(imag_ind_);
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (imag_ind_ == NULL && i != 0 && j != 0)
      throw NoMemory(string("Matrix_ComplexSparse::")
		     + "Matrix_ComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * imag_nz)
		     + " bytes to store " + to_str(imag_nz)
		     + " row or column indices (imaginary part), for a "
		     + to_str(i) + " by " + to_str(j) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->real_data_ = this->allocator_.allocate(real_nz_, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	free(real_ptr_);
	free(imag_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	free(real_ind_);
	free(imag_ind_);
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_data_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	free(real_ptr_);
	free(imag_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	free(real_ind_);
	free(imag_ind_);
	real_ind_ = NULL;
	imag_ind_ = NULL;
	imag_data_ = NULL;
      }
    if (real_data_ == NULL && i != 0 && j != 0)
      throw NoMemory(string("Matrix_ComplexSparse::")
		     + "Matrix_ComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * real_nz)
		     + " bytes to store " + to_str(real_nz)
		     + " values (real part), for a "
		     + to_str(i) + " by " + to_str(j) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->imag_data_ = this->allocator_.allocate(imag_nz_, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	free(real_ptr_);
	free(imag_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	free(real_ind_);
	free(imag_ind_);
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->allocator_.deallocate(this->real_data_, real_nz_);
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_data_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	free(real_ptr_);
	free(imag_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	free(real_ind_);
	free(imag_ind_);
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->allocator_.deallocate(this->real_data_, real_nz_);
	real_data_ = NULL;
      }
    if (imag_data_ == NULL && i != 0 && j != 0)
      throw NoMemory(string("Matrix_ComplexSparse::")
		     + "Matrix_ComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * imag_nz)
		     + " bytes to store " + to_str(imag_nz)
		     + " values (imaginary part), for a "
		     + to_str(i) + " by " + to_str(j) + " matrix.");
#endif

  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix with non-zero values and indices
    provided by 'real_values' (values of the real part), 'real_ptr'
    (pointers for the real part), 'real_ind' (indices for the real part),
    'imag_values' (values of the imaginary part), 'imag_ptr'
    (pointers for the imaginary part) and 'imag_ind' (indices for the
    imaginary part). Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param real_values values of non-zero entries for the real part.
    \param real_ptr row or column start indices for the real part.
    \param real_ind row or column indices for the real part.
    \param imag_values values of non-zero entries for the imaginary part.
    \param imag_ptr row or column start indices for the imaginary part.
    \param imag_ind row or column indices for the imaginary part.
    \warning Input vectors 'real_values', 'real_ptr' and 'real_ind',
    'imag_values', 'imag_ptr' and 'imag_ind' are empty on exit.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  inline Matrix_ComplexSparse<T, Prop, Storage, Allocator>::
  Matrix_ComplexSparse(int i, int j,
		       Vector<T, Storage0, Allocator0>& real_values,
		       Vector<int, Storage1, Allocator1>& real_ptr,
		       Vector<int, Storage2, Allocator2>& real_ind,
		       Vector<T, Storage0, Allocator0>& imag_values,
		       Vector<int, Storage1, Allocator1>& imag_ptr,
		       Vector<int, Storage2, Allocator2>& imag_ind): 
    Matrix_Base<T, Allocator>(i, j)
  {
    
    real_nz_ = real_values.GetLength();
    imag_nz_ = imag_values.GetLength();
    
#ifdef SELDON_CHECK_DIMENSIONS
    // Checks whether vector sizes are acceptable.
    
    if (real_ind.GetLength() != real_nz_)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_ComplexSparse::")
		       + string("Matrix_ComplexSparse(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("There are ") + to_str(real_nz_)
		       + " values (real part) but "
		       + to_str(real_ind.GetLength())
		       + " row or column indices.");
      }

    if (imag_ind.GetLength() != imag_nz_)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_ComplexSparse::")
		       + string("Matrix_ComplexSparse(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("There are ") + to_str(imag_nz_)
		       + " values (imaginary part) but "
		       + to_str(imag_ind.GetLength())
		       + " row or column indices.");
      }

    if (real_ptr.GetLength()-1 != Storage::GetFirst(i, j))
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_ComplexSparse::")
		       + string("Matrix_ComplexSparse(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("The vector of start indices (real part)")
		       + " contains " + to_str(real_ptr.GetLength()-1)
		       + string(" row or column start indices (plus the")
		       + " number of non-zero entries) but there are "
		       + to_str(Storage::GetFirst(i, j))
		       + " rows or columns ("
		       + to_str(i) + " by " + to_str(j) + " matrix).");
      }

    if (imag_ptr.GetLength()-1 != Storage::GetFirst(i, j))
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_ComplexSparse::")
		       + string("Matrix_ComplexSparse(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("The vector of start indices (imaginary part)")
		       + " contains " + to_str(imag_ptr.GetLength()-1)
		       + string(" row or column start indices (plus the")
		       + " number of non-zero entries) but there are "
		       + to_str(Storage::GetFirst(i, j))
		       + " rows or columns ("
		       + to_str(i) + " by " + to_str(j) + " matrix).");
      }

    if ( (static_cast<long int>(real_nz_-1) / static_cast<long int>(j)
	  >= static_cast<long int>(i)) ||
	 (static_cast<long int>(imag_nz_-1) / static_cast<long int>(j)
	  >= static_cast<long int>(i)) )
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_ComplexSparse::")
		       + string("Matrix_ComplexSparse(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("There are more values (")
		       + to_str(real_values.GetLength())
		       + " values for the real part and "
		       + to_str(real_values.GetLength()) + string(" values")
		       + string(" for the imaginary part) than elements")
		       + " in the matrix ("
		       + to_str(i) + " by " + to_str(j) + ").");
      }
#endif

    this->real_ptr_ = real_ptr.GetData();
    this->imag_ptr_ = imag_ptr.GetData();
    this->real_ind_ = real_ind.GetData();
    this->imag_ind_ = imag_ind.GetData();
    this->real_data_ = real_values.GetData();
    this->imag_data_ = imag_values.GetData();

    real_ptr.Nullify();
    imag_ptr.Nullify();
    real_ind.Nullify();
    imag_ind.Nullify();
    real_values.Nullify();
    imag_values.Nullify();
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_ComplexSparse<T, Prop, Storage, Allocator>
  ::~Matrix_ComplexSparse()
  {
    this->m_ = 0;
    this->n_ = 0;

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	if (real_ptr_ != NULL)
	  {
	    free(real_ptr_);
	    real_ptr_ = NULL;
	  }
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	real_ptr_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	if (imag_ptr_ != NULL)
	  {
	    free(imag_ptr_);
	    imag_ptr_ = NULL;
	  }
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	imag_ptr_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	if (real_ind_ != NULL)
	  {
	    free(real_ind_);
	    real_ind_ = NULL;
	  }
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	real_ind_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	if (imag_ind_ != NULL)
	  {
	    free(imag_ind_);
	    imag_ind_ = NULL;
	  }
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	imag_ind_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	if (this->real_data_ != NULL)
	  {
	    this->allocator_.deallocate(this->real_data_, real_nz_);
	    this->real_data_ = NULL;
	  }
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->real_nz_ = 0;
	this->real_data_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	if (this->imag_data_ != NULL)
	  {
	    this->allocator_.deallocate(this->imag_data_, imag_nz_);
	    this->imag_data_ = NULL;
	  }
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->imag_nz_ = 0;
	this->imag_data_ = NULL;
      }
#endif

    this->real_nz_ = 0;
    this->imag_nz_ = 0;
  }
  

  //! Clears the matrix.
  /*! This methods is equivalent to the destructor. On exit, the matrix
    is empty (0x0).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ComplexSparse<T, Prop, Storage, Allocator>::Clear()
  {
    this->~Matrix_ComplexSparse();
  }


  /*********************
   * MEMORY MANAGEMENT *
   *********************/


  //! Redefines the matrix.
  /*! It clears the matrix and sets it to a new matrix defined by
    'real_values' (values of the real part), 'real_ptr'
    (pointers for the real part), 'real_ind' (indices for the real part),
    'imag_values' (values of the imaginary part), 'imag_ptr'
    (pointers for the imaginary part) and 'imag_ind' (indices for the
    imaginary part). Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param real_values values of non-zero entries for the real part.
    \param real_ptr row or column start indices for the real part.
    \param real_ind row or column indices for the real part.
    \param imag_values values of non-zero entries for the imaginary part.
    \param imag_ptr row or column start indices for the imaginary part.
    \param imag_ind row or column indices for the imaginary part.
    \warning Input vectors 'real_values', 'real_ptr' and 'real_ind',
    'imag_values', 'imag_ptr' and 'imag_ind' are empty on exit.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  void Matrix_ComplexSparse<T, Prop, Storage, Allocator>::
  SetData(int i, int j,
	  Vector<T, Storage0, Allocator0>& real_values,
	  Vector<int, Storage1, Allocator1>& real_ptr,
	  Vector<int, Storage2, Allocator2>& real_ind,
	  Vector<T, Storage0, Allocator0>& imag_values,
	  Vector<int, Storage1, Allocator1>& imag_ptr,
	  Vector<int, Storage2, Allocator2>& imag_ind)
  {
    
    this->m_ = i;
    this->n_ = j;
    real_nz_ = real_values.GetLength();
    imag_nz_ = imag_values.GetLength();
    
#ifdef SELDON_CHECK_DIMENSIONS
    // Checks whether vector sizes are acceptable.
    
    if (real_ind.GetLength() != real_nz_)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_ComplexSparse::SetData(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("There are ") + to_str(real_nz_)
		       + " values (real part) but "
		       + to_str(real_ind.GetLength())
		       + " row or column indices.");
      }

    if (imag_ind.GetLength() != imag_nz_)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_ComplexSparse::SetData(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("There are ") + to_str(imag_nz_)
		       + " values (imaginary part) but "
		       + to_str(imag_ind.GetLength())
		       + " row or column indices.");
      }

    if (real_ptr.GetLength()-1 != Storage::GetFirst(i, j))
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_ComplexSparse::SetData(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("The vector of start indices (real part)")
		       + " contains " + to_str(real_ptr.GetLength()-1)
		       + string(" row or column start indices (plus the")
		       + " number of non-zero entries) but there are "
		       + to_str(Storage::GetFirst(i, j))
		       + " rows or columns ("
		       + to_str(i) + " by " + to_str(j) + " matrix).");
      }

    if (imag_ptr.GetLength()-1 != Storage::GetFirst(i, j))
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_ComplexSparse::SetData(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("The vector of start indices (imaginary part)")
		       + " contains " + to_str(imag_ptr.GetLength()-1)
		       + string(" row or column start indices (plus the")
		       + " number of non-zero entries) but there are "
		       + to_str(Storage::GetFirst(i, j))
		       + " rows or columns ("
		       + to_str(i) + " by " + to_str(j) + " matrix).");
      }

    if ( (static_cast<long int>(real_nz_-1) / static_cast<long int>(j)
	  >= static_cast<long int>(i)) ||
	 (static_cast<long int>(imag_nz_-1) / static_cast<long int>(j)
	  >= static_cast<long int>(i)) )
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_ComplexSparse::SetData(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("There are more values (")
		       + to_str(real_values.GetLength())
		       + " values for the real part and "
		       + to_str(real_values.GetLength()) + string(" values")
		       + string(" for the imaginary part) than elements")
		       + " in the matrix ("
		       + to_str(i) + " by " + to_str(j) + ").");
      }
#endif

    this->real_ptr_ = real_ptr.GetData();
    this->imag_ptr_ = imag_ptr.GetData();
    this->real_ind_ = real_ind.GetData();
    this->imag_ind_ = imag_ind.GetData();
    this->real_data_ = real_values.GetData();
    this->imag_data_ = imag_values.GetData();

    real_ptr.Nullify();
    imag_ptr.Nullify();
    real_ind.Nullify();
    imag_ind.Nullify();
    real_values.Nullify();
    imag_values.Nullify();
  }

  
  //! Redefines the matrix.
  /*! It clears the matrix and sets it to a new matrix defined by arrays
    'real_values' (values of the real part), 'real_ptr'
    (pointers for the real part), 'real_ind' (indices for the real part),
    'imag_values' (values of the imaginary part), 'imag_ptr'
    (pointers for the imaginary part) and 'imag_ind' (indices for the
    imaginary part).
    \param i number of rows.
    \param j number of columns.
    \param real_nz number of non-zero entries (real part).
    \param real_values values of non-zero entries for the real part.
    \param real_ptr row or column start indices for the real part.
    \param real_ind row or column indices for the real part.
    \param imag_nz number of non-zero entries (imaginary part).
    \param imag_values values of non-zero entries for the imaginary part.
    \param imag_ptr row or column start indices for the imaginary part.
    \param imag_ind row or column indices for the imaginary part.
    \warning On exit, arrays 'real_values', 'real_ptr', 'real_ind',
    'imag_values', 'imag_ptr' and 'imag_ind' are managed by the matrix.
    For example, it means that the destructor will release those arrays;
    therefore, the user mustn't release those arrays.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ComplexSparse<T, Prop, Storage, Allocator>::
  SetData(int i, int j, int real_nz,
	  typename Matrix_ComplexSparse<T, Prop, Storage, Allocator>
	  ::pointer real_values,
	  int* real_ptr, int* real_ind, int imag_nz,
	  typename Matrix_ComplexSparse<T, Prop, Storage, Allocator>
	  ::pointer imag_values,
	  int* imag_ptr, int* imag_ind)
  {
    this->Clear();

    this->m_ = i;
    this->n_ = j;

    this->real_nz_ = real_nz;
    this->imag_nz_ = imag_nz;

    real_data_ = real_values;
    imag_data_ = imag_values;
    real_ind_ = real_ind;
    imag_ind_ = imag_ind;
    real_ptr_ = real_ptr;
    imag_ptr_ = imag_ptr;
  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/
  

  //! Returns the number of elements stored in memory.
  /*!
    Returns the number of elements stored in memory, i.e.
    the cumulated number of non-zero entries of both the real and
    the imaginary part.
    \return The number of elements stored in memory.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int Matrix_ComplexSparse<T, Prop, Storage, Allocator>::GetDataSize() const
  {
    return real_nz_ + imag_nz_;
  }


  //! Returns (row or column) start indices for the real part.
  /*!
    Returns the array ('ptr_') of start indices for the real part.
    \return The array of start indices for the real part.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int* Matrix_ComplexSparse<T, Prop, Storage, Allocator>::GetRealPtr() const
  {
    return real_ptr_;
  }


  //! Returns (row or column) start indices for the imaginary part.
  /*!
    Returns the array ('ptr_') of start indices for the imaginary part.
    \return The array of start indices for the imaginary part.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int* Matrix_ComplexSparse<T, Prop, Storage, Allocator>::GetImagPtr() const
  {
    return imag_ptr_;
  }


  //! Returns (row or column) indices of non-zero entries for the real part.
  /*!
    Returns the array ('ind_') of (row or column) indices
    of non-zero entries for the real part. This array defines non-zero
    entries indices if coupled with (column or row) start indices.
    \return The array of (row or column) indices of
    non-zero entries for the real part.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int* Matrix_ComplexSparse<T, Prop, Storage, Allocator>::GetRealInd() const
  {
    return real_ind_;
  }


  //! Returns (row or column) indices of non-zero entries
  //! for the imaginary part.
  /*!
    Returns the array ('ind_') of (row or column) indices
    of non-zero entries for the imaginary part. This array defines non-zero
    entries indices if coupled with (column or row) start indices.
    \return The array of (row or column) indices of
    non-zero entries for the imaginary part.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int* Matrix_ComplexSparse<T, Prop, Storage, Allocator>::GetImagInd() const
  {
    return imag_ind_;
  }


  //! Returns the length of the array of start indices for the real part.
  /*!
    \return The length of the array of start indices for the real part.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int Matrix_ComplexSparse<T, Prop, Storage, Allocator>
  ::GetRealPtrSize() const
  {
    return (Storage::GetFirst(this->m_, this->n_) + 1);
  }

  
  //! Returns the length of the array of start indices for the imaginary part.
  /*!
    \return The length of the array of start indices for the imaginary part.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int Matrix_ComplexSparse<T, Prop, Storage, Allocator>
  ::GetImagPtrSize() const
  {
    return (Storage::GetFirst(this->m_, this->n_) + 1);
  }

  
  //! Returns the length of the array of (column or row) indices
  //! for the real part.
  /*!
    Returns the length of the array ('ind_') of (row or column) indices
    of non-zero entries for the real part. This array defines non-zero entries
    indices if coupled with (column or row) start indices.
    \return The length of the array of (column or row) indices
    for the real part.
    \note The length of the array of (column or row) indices is the
    number of non-zero entries.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int Matrix_ComplexSparse<T, Prop, Storage, Allocator>
  ::GetRealIndSize() const
  {
    return real_nz_;
  }


  //! Returns the length of the array of (column or row) indices
  //! for the imaginary part.
  /*!
    Returns the length of the array ('ind_') of (row or column) indices
    of non-zero entries for the imaginary part. This array defines non-zero
    entries indices if coupled with (column or row) start indices.
    \return The length of the array of (column or row) indices
    for the imaginary part.
    \note The length of the array of (column or row) indices is the
    number of non-zero entries.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int Matrix_ComplexSparse<T, Prop, Storage, Allocator>
  ::GetImagIndSize() const
  {
    return imag_nz_;
  }


  //! Returns the array of values of the real part.
  /*!
    \return The array 'real_data_' of values of the real part..
  */
  template <class T, class Prop, class Storage, class Allocator>
  T* Matrix_ComplexSparse<T, Prop, Storage, Allocator>::GetRealData() const
  {
    return real_data_;
  }


  //! Returns the array of values of the imaginary part.
  /*!
    \return The array 'imag_data_' of values of the imaginary part..
  */
  template <class T, class Prop, class Storage, class Allocator>
  T* Matrix_ComplexSparse<T, Prop, Storage, Allocator>::GetImagData() const
  {
    return imag_data_;
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline complex<typename Matrix_ComplexSparse<T, Prop, Storage, Allocator>
  ::value_type>
  Matrix_ComplexSparse<T, Prop, Storage, Allocator>
  ::operator() (int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDARIES
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_ComplexSparse::operator()",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_ComplexSparse::operator()",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
#endif

    int real_k, imag_k, l;
    int real_a, real_b;
    int imag_a, imag_b;

    real_a = real_ptr_[Storage::GetFirst(i, j)];
    real_b = real_ptr_[Storage::GetFirst(i, j) + 1];

    imag_a = imag_ptr_[Storage::GetFirst(i, j)];
    imag_b = imag_ptr_[Storage::GetFirst(i, j) + 1];
    
    if (real_a != real_b)
      {
	l = Storage::GetSecond(i, j);
	for (real_k = real_a;
	     (real_k <real_b-1) && (real_ind_[real_k] < l);
	     real_k++);
	if (imag_a != imag_b)
	  {
	    for (imag_k = imag_a;
		 (imag_k < imag_b-1) && (imag_ind_[imag_k] < l);
		 imag_k++);
	    if (real_ind_[real_k] == l)
	      {
		if (imag_ind_[imag_k] == l)
		  return complex<T>(real_data_[real_k], imag_data_[imag_k]);
		else
		  return complex<T>(real_data_[real_k], T(0));
	      }
	    else
	      if (imag_ind_[imag_k] == l)
		return complex<T>(T(0), imag_data_[imag_k]);
	      else
		return complex<T>(T(0), T(0));
	  }
	else
	  {
	    if (real_ind_[real_k] == l)
	      return complex<T>(real_data_[real_k], T(0));
	    else
	      return complex<T>(T(0), T(0));
	  }
      }
    else
      {
	if (imag_a != imag_b)
	  {
	    l = Storage::GetSecond(i, j);
	    for (imag_k = imag_a;
		 (imag_k < imag_b-1) && (imag_ind_[imag_k] < l);
		 imag_k++);
	    if (imag_ind_[imag_k] == l)
	      return complex<T>(T(0), imag_data_[imag_k]);
	    else
	      return complex<T>(T(0), T(0));
	  }
	else
	  return complex<T>(T(0), T(0));
      }

  }


  /************************
   * CONVENIENT FUNCTIONS *
   ************************/


  //! Displays the matrix on the standard output.
  /*!
    Displays elements on the standard output, in text format.
    Each row is displayed on a single line and elements of
    a row are delimited by tabulations.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ComplexSparse<T, Prop, Storage, Allocator>::Print() const
  {
    for (int i = 0; i < this->m_; i++)
      {
	for (int j = 0; j < this->n_; j++)
	  cout << (*this)(i, j) << "\t";
	cout << endl;
      }
  }


 
  //////////////////////////////
  // MATRIX<COLCOMPLEXSPARSE> //
  //////////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/

  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, ColComplexSparse, Allocator>::Matrix()  throw():
    Matrix_ComplexSparse<T, Prop, ColComplexSparse, Allocator>()
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix with real_nz and imag_nz non-zero elements for
    the real part and the imaginary part respectively.
    \param i number of rows.
    \param j number of columns.
    \param real_nz number of non-zero elements for the real part.
    \param imag_nz number of non-zero elements for the imaginary part.
    \note Matrix values are not initialized.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, ColComplexSparse, Allocator>::Matrix(int i, int j,
						       int real_nz,
						       int imag_nz):
    Matrix_ComplexSparse<T, Prop, ColComplexSparse, Allocator>(i, j,
							       real_nz,
							       imag_nz)
  {
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix with non-zero values and indices
    provided by 'real_values' (values of the real part), 'real_ptr'
    (pointers for the real part), 'real_ind' (indices for the real part),
    'imag_values' (values of the imaginary part), 'imag_ptr'
    (pointers for the imaginary part) and 'imag_ind' (indices for the
    imaginary part). Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param real_values values of non-zero entries for the real part.
    \param real_ptr row or column start indices for the real part.
    \param real_ind row or column indices for the real part.
    \param imag_values values of non-zero entries for the imaginary part.
    \param imag_ptr row or column start indices for the imaginary part.
    \param imag_ind row or column indices for the imaginary part.
    \warning Input vectors 'real_values', 'real_ptr' and 'real_ind',
    'imag_values', 'imag_ptr' and 'imag_ind' are empty on exit.
  */
  template <class T, class Prop, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  Matrix<T, Prop, ColComplexSparse, Allocator>::
  Matrix(int i, int j,
	 Vector<T, Storage0, Allocator0>& real_values,
	 Vector<int, Storage1, Allocator1>& real_ptr,
	 Vector<int, Storage2, Allocator2>& real_ind,
	 Vector<T, Storage0, Allocator0>& imag_values,
	 Vector<int, Storage1, Allocator1>& imag_ptr,
	 Vector<int, Storage2, Allocator2>& imag_ind):
    Matrix_ComplexSparse<T, Prop, ColComplexSparse, Allocator>(i, j,
							       real_values,
							       real_ptr,
							       real_ind,
							       imag_values,
							       imag_ptr,
							       imag_ind)
  {
  }



  //////////////////////////////
  // MATRIX<ROWCOMPLEXSPARSE> //
  //////////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/

  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, RowComplexSparse, Allocator>::Matrix()  throw():
    Matrix_ComplexSparse<T, Prop, RowComplexSparse, Allocator>()
  {
  }


  /*! Builds a i by j matrix with real_nz and imag_nz non-zero elements for
    the real part and the imaginary part respectively.
    \param i number of rows.
    \param j number of columns.
    \param real_nz number of non-zero elements for the real part.
    \param imag_nz number of non-zero elements for the imaginary part.
    \note Matrix values are not initialized.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, RowComplexSparse, Allocator>::Matrix(int i, int j,
						       int real_nz,
						       int imag_nz):
    Matrix_ComplexSparse<T, Prop, RowComplexSparse, Allocator>(i, j,
							       real_nz,
							       imag_nz)
  {
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix with non-zero values and indices
    provided by 'real_values' (values of the real part), 'real_ptr'
    (pointers for the real part), 'real_ind' (indices for the real part),
    'imag_values' (values of the imaginary part), 'imag_ptr'
    (pointers for the imaginary part) and 'imag_ind' (indices for the
    imaginary part). Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param real_values values of non-zero entries for the real part.
    \param real_ptr row or column start indices for the real part.
    \param real_ind row or column indices for the real part.
    \param imag_values values of non-zero entries for the imaginary part.
    \param imag_ptr row or column start indices for the imaginary part.
    \param imag_ind row or column indices for the imaginary part.
    \warning Input vectors 'real_values', 'real_ptr' and 'real_ind',
    'imag_values', 'imag_ptr' and 'imag_ind' are empty on exit.
  */
  template <class T, class Prop, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  Matrix<T, Prop, RowComplexSparse, Allocator>::
  Matrix(int i, int j,
	 Vector<T, Storage0, Allocator0>& real_values,
	 Vector<int, Storage1, Allocator1>& real_ptr,
	 Vector<int, Storage2, Allocator2>& real_ind,
	 Vector<T, Storage0, Allocator0>& imag_values,
	 Vector<int, Storage1, Allocator1>& imag_ptr,
	 Vector<int, Storage2, Allocator2>& imag_ind):
    Matrix_ComplexSparse<T, Prop, RowComplexSparse, Allocator>(i, j,
							       real_values,
							       real_ptr,
							       real_ind,
							       imag_values,
							       imag_ptr,
							       imag_ind)
  {
  }


} // namespace Seldon.

#define SELDON_FILE_MATRIX_COMPLEXSPARSE_CXX
#endif
