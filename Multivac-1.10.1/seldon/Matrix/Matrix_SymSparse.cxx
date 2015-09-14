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

#ifndef SELDON_FILE_MATRIX_SYMSPARSE_CXX

#include "Matrix_SymSparse.hxx"

namespace Seldon
{


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Builds a empty matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_SymSparse<T, Prop, Storage, Allocator>::Matrix_SymSparse():
    Matrix_Base<T, Allocator>()
  {
    nz_ = 0;
    ptr_ = NULL;
    ind_ = NULL;
  }


  //! Constructor.
  /*!
    Builds an empty i by j sparse matrix.
    \param i number of rows.
    \param j number of columns.
    \note 'j' is assumed to be equal to 'i' and is therefore discarded.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::Matrix_SymSparse(int i, int j): Matrix_Base<T, Allocator>(i, i)
  {
    nz_ = 0;
    ptr_ = NULL;
    ind_ = NULL;
  }


  //! Constructor.
  /*! Builds a sparse matrix of size i x j , with nz non-zero (stored)
    elements.
    \param i number of rows.
    \param j number of columns.
    \param nz number of non-zero elements that are stored.
    \note Matrix values are not initialized. Indices of non-zero entries
    are not initialized either.
    \note 'j' is assumed to be equal to 'i' and is therefore discarded.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_SymSparse<T, Prop, Storage, Allocator>::
  Matrix_SymSparse(int i, int j, int nz):
    Matrix_Base<T, Allocator>(i, i)
  {
    this->nz_ = nz;

#ifdef SELDON_CHECK_DIMENSIONS
    if (static_cast<long int>(2 * nz_ - 2) / static_cast<long int>(i + 1)
	>= static_cast<long int>(i))
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim("Matrix_SymSparse::Matrix_SymSparse(int, int, int)",
		       string("There are more values (") + to_str(nz)
		       + string(" values) than elements in the upper")
		       + " part of the matrix ("
		       + to_str(i) + " by " + to_str(i) + ").");
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	ptr_ = reinterpret_cast<int*>( calloc(i + 1, sizeof(int)) );

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
      }
    if (ptr_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ind_ = NULL;
	this->data_ = NULL;
      }
    if (ptr_ == NULL && i != 0)
      throw NoMemory("Matrix_SymSparse::Matrix_SymSparse(int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * (i+1) ) + " bytes to store "
		     + to_str(i+1) + " row or column start indices, for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	ind_ = reinterpret_cast<int*>( calloc(nz_, sizeof(int)) );
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	free(ptr_);
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
      }
    if (ind_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	free(ptr_);
	ptr_ = NULL;
	this->data_ = NULL;
      }
    if (ind_ == NULL && i != 0)
      throw NoMemory("Matrix_SymSparse::Matrix_SymSparse(int, int, int)",
		     string("Unable to allocate ") + to_str(sizeof(int) * nz)
		     + " bytes to store " + to_str(nz)
		     + " row or column indices, for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->data_ = this->allocator_.allocate(nz_, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	free(ptr_);
	ptr_ = NULL;
	free(ind_);
	ind_ = NULL;
	this->data_ = NULL;
      }
    if (this->data_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	free(ptr_);
	ptr_ = NULL;
	free(ind_);
	ind_ = NULL;
      }
    if (this->data_ == NULL && i != 0)
      throw NoMemory("Matrix_SymSparse::Matrix_SymSparse(int, int, int)",
		     string("Unable to allocate ") + to_str(sizeof(int) * nz)
		     + " bytes to store " + to_str(nz) + " values, for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix with non-zero values and indices
    provided by 'values' (values), 'ptr' (pointers) and 'ind' (indices).
    Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param values values of non-zero entries.
    \param ptr row or column start indices.
    \param ind row or column indices.
    \warning Input vectors 'values', 'ptr' and 'ind' are empty on exit.
    Moreover 'j' is assumed to be equal to i so that 'j' is discarded.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  inline Matrix_SymSparse<T, Prop, Storage, Allocator>::
  Matrix_SymSparse(int i, int j,
		   Vector<T, Storage0, Allocator0>& values,
		   Vector<int, Storage1, Allocator1>& ptr,
		   Vector<int, Storage2, Allocator2>& ind):
    Matrix_Base<T, Allocator>(i, j)
  {
    nz_ = values.GetLength();
    
#ifdef SELDON_CHECK_DIMENSIONS
    // Checks whether vector sizes are acceptable.
    
    if (ind.GetLength() != nz_)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim(string("Matrix_SymSparse::Matrix_SymSparse(int, int, ")
		       + "const Vector&, const Vector&, const Vector&)",
		       string("There are ") + to_str(nz_) + " values but "
		       + to_str(ind.GetLength()) + " row or column indices.");
      }

    if (ptr.GetLength() - 1 != i)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim(string("Matrix_SymSparse::Matrix_SymSparse(int, int, ")
		       + "const Vector&, const Vector&, const Vector&)",
		       string("The vector of start indices contains ")
		       + to_str(ptr.GetLength()-1) + string(" row or column ")
		       + string("start  indices (plus the number of non-zero")
		       + " entries) but there are " + to_str(i)
		       + " rows or columns ("
		       + to_str(i) + " by " + to_str(i) + " matrix).");
      }

    if (static_cast<long int>(2 * nz_ - 2) / static_cast<long int>(i + 1)
	>= static_cast<long int>(i))
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim(string("Matrix_SymSparse::Matrix_SymSparse(int, int, ")
		       + "const Vector&, const Vector&, const Vector&)",
		       string("There are more values (")
		       + to_str(values.GetLength())
		       + " values) than elements in the matrix ("
		       + to_str(i) + " by " + to_str(i) + ").");
      }
#endif

    this->ptr_ = ptr.GetData();
    this->ind_ = ind.GetData();
    this->data_ = values.GetData();

    ptr.Nullify();
    ind.Nullify();
    values.Nullify();
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_SymSparse<T, Prop, Storage, Allocator>::~Matrix_SymSparse()
  {
    this->m_ = 0;
    this->n_ = 0;

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	if (ptr_ != NULL)
	  {
	    free(ptr_);
	    ptr_ = NULL;
	  }
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	ptr_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	if (ind_ != NULL)
	  {
	    free(ind_);
	    ind_ = NULL;
	  }
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	ind_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	if (this->data_ != NULL)
	  {
	    this->allocator_.deallocate(this->data_, nz_);
	    this->data_ = NULL;
	  }
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->nz_ = 0;
	this->data_ = NULL;
      }
#endif

    this->nz_ = 0;
  }
  

  //! Clears the matrix.
  /*! This methods is equivalent to the destructor. On exit,
    the matrix is empty (0x0).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymSparse<T, Prop, Storage, Allocator>::Clear()
  {
    this->~Matrix_SymSparse();
  }


  /*********************
   * MEMORY MANAGEMENT *
   *********************/


  //! Redefines the matrix.
  /*! It clears the matrix and sets it to a new matrix defined by
    'values' (values), 'ptr' (pointers) and 'ind' (indices).
    Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param values values of non-zero entries.
    \param ptr row or column start indices.
    \param ind row or column indices.
    \warning Input vectors 'values', 'ptr' and 'ind' are empty on exit.
    Moreover 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  void Matrix_SymSparse<T, Prop, Storage, Allocator>::
  SetData(int i, int j,
	  Vector<T, Storage0, Allocator0>& values,
	  Vector<int, Storage1, Allocator1>& ptr,
	  Vector<int, Storage2, Allocator2>& ind)
  {
    this->m_ = i;
    this->n_ = i;
    this->nz_ = values.GetLength();
    
#ifdef SELDON_CHECK_DIMENSIONS
    // Checks whether vector sizes are acceptable.
    
    if (ind.GetLength() != nz_)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim(string("Matrix_SymSparse::Reallocate(int, int, ")
		       + "const Vector&, const Vector&, const Vector&)",
		       string("There are ") + to_str(nz_) + " values but "
		       + to_str(ind.GetLength()) + " row or column indices.");
      }

    if (ptr.GetLength()-1 != i)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim(string("Matrix_SymSparse::Reallocate(int, int, ")
		       + "const Vector&, const Vector&, const Vector&)",
		       string("The vector of start indices contains ")
		       + to_str(ptr.GetLength()-1) + string(" row or column")
		       + string(" start indices (plus the number of")
		       + string(" non-zero entries) but there are ")
		       + to_str(i) + " rows or columns ("
		       + to_str(i) + " by " + to_str(i) + " matrix).");
      }

    if (static_cast<long int>(2 * nz_ - 2) / static_cast<long int>(i + 1)
	>= static_cast<long int>(i))
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim(string("Matrix_SymSparse::Reallocate(int, int, ")
		       + "const Vector&, const Vector&, const Vector&)",
		       string("There are more values (")
		       + to_str(values.GetLength())
		       + " values) than elements in the matrix ("
		       + to_str(i) + " by " + to_str(i) + ").");
      }
#endif

    this->ptr_ = ptr.GetData();
    this->ind_ = ind.GetData();
    this->data_ = values.GetData();

    ptr.Nullify();
    ind.Nullify();
    values.Nullify();
  }

  
  //! Redefines the matrix.
  /*! It clears the matrix and sets it to a new matrix defined by arrays
    'values' (values), 'ptr' (pointers) and 'ind' (indices).
    \param i number of rows.
    \param j number of columns.
    \param nz number of non-zero entries that are stored.
    \param values values of non-zero entries.
    \param ptr row or column start indices.
    \param ind row or column indices.
    \warning On exit, arrays 'values', 'ptr' and 'ind' are managed by the matrix.
    For example, it means that the destructor will released those arrays;
    therefore, the user mustn't release those arrays.
    Moreover 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::SetData(int i, int j, int nz,
	    typename Matrix_SymSparse<T, Prop, Storage, Allocator>
	    ::pointer values,
	    int* ptr,
	    int* ind)
  {
    this->Clear();

    this->m_ = i;
    this->n_ = i;

    this->nz_ = nz;

    this->data_ = values;
    ind_ = ind;
    ptr_ = ptr;
  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/
  
  
  //! Returns the number of elements stored in memory.
  /*!
    \return The number of elements stored in memory (the number
    of non-zero elements).
  */
  template <class T, class Prop, class Storage, class Allocator>
  int Matrix_SymSparse<T, Prop, Storage, Allocator>::GetDataSize() const
  {
    return nz_;
  }


  //! Returns (row or column) start indices.
  /*!
    Returns the array ('ptr_') of start indices.
    \return The array of start indices.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int* Matrix_SymSparse<T, Prop, Storage, Allocator>::GetPtr() const
  {
    return ptr_;
  }


  //! Returns (row or column) indices of non-zero entries.
  /*!
    Returns the array ('ind_') of (row or column) indices
    of non-zero entries. This array defines non-zero entries
    indices if coupled with (column or row) start indices.
    \return The array of (row or column) indices of
    non-zero entries.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int* Matrix_SymSparse<T, Prop, Storage, Allocator>::GetInd() const
  {
    return ind_;
  }


  //! Returns the length of the array of start indices.
  /*!
    \return The length of the array of start indices.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int Matrix_SymSparse<T, Prop, Storage, Allocator>::GetPtrSize() const
  {
    return (this->m_ + 1);
  }

  
  //! Returns the length of the array of (column or row) indices.
  /*!
    Returns the length of the array ('ind_') of (row or column) indices
    of non-zero entries (that are stored). This array defines non-zero
    entries indices (that are stored) if coupled with (column or row)
    start indices.
    \return The length of the array of (column or row) indices.
    \note The length of the array of (column or row) indices is the
    number of non-zero entries that are stored.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int Matrix_SymSparse<T, Prop, Storage, Allocator>::GetIndSize() const
  {
    return nz_;
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
  inline typename Matrix_SymSparse<T, Prop, Storage, Allocator>::value_type
  Matrix_SymSparse<T, Prop, Storage, Allocator>::operator() (int i,
							     int j) const
  {

#ifdef SELDON_CHECK_BOUNDARIES
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_SymSparse::operator()",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_SymSparse::operator()",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
#endif

    int k, l;
    int a, b;

    // Only the upper part is stored.
    if (i>j)
      {
	l = i;
	i = j;
	j = l;
      }

    a = ptr_[Storage::GetFirst(i, j)];
    b = ptr_[Storage::GetFirst(i, j) + 1];

    if (a == b)
      return T(0);

    l = Storage::GetSecond(i, j);

    for (k = a; (k<b-1) && (ind_[k]<l); k++);

    if (ind_[k] == l)
      return this->data_[k];
    else
      return T(0);
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
  void Matrix_SymSparse<T, Prop, Storage, Allocator>::Print() const
  {
    for (int i = 0; i < this->m_; i++)
      {
	for (int j = 0; j < this->n_; j++)
	  cerr << (*this)(i, j) << "\t";
	cerr << endl;
      }
  }


 
  //////////////////////////
  // MATRIX<COLSYMSPARSE> //
  //////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/

  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, ColSymSparse, Allocator>::Matrix()  throw():
    Matrix_SymSparse<T, Prop, ColSymSparse, Allocator>()
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix with nz non-zero (stored) elements.
    \param i number of rows.
    \param j number of columns.
    \param nz number of non-zero elements that are stored.
    \note Matrix values are not initialized.
    \warning 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, ColSymSparse, Allocator>::Matrix(int i, int j, int nz):
    Matrix_SymSparse<T, Prop, ColSymSparse, Allocator>(i, j, nz)
  {
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix with non-zero values and indices
    provided by 'values' (values), 'ptr' (pointers) and 'ind' (indices).
    Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param values values of non-zero entries.
    \param ptr row start indices.
    \param ind column indices.
    \warning Input vectors 'values', 'ptr' and 'ind' are empty on exit.
    Moreover 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  Matrix<T, Prop, ColSymSparse, Allocator>::
  Matrix(int i, int j,
	 Vector<T, Storage0, Allocator0>& values,
	 Vector<int, Storage1, Allocator1>& ptr,
	 Vector<int, Storage2, Allocator2>& ind):
    Matrix_SymSparse<T, Prop, ColSymSparse, Allocator>(i, j, values, ptr, ind)
  {
  }



  //////////////////////////
  // MATRIX<ROWSYMSPARSE> //
  //////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, RowSymSparse, Allocator>::Matrix()  throw():
    Matrix_SymSparse<T, Prop, RowSymSparse, Allocator>()
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix with nz non-zero (stored) elements.
    \param i number of rows.
    \param j number of columns.
    \param nz number of non-zero elements that are stored.
    \note Matrix values are not initialized.
    \warning 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, RowSymSparse, Allocator>::Matrix(int i, int j, int nz):
    Matrix_SymSparse<T, Prop, RowSymSparse, Allocator>(i, j, nz)
  {
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix with non-zero values and indices
    provided by 'values' (values), 'ptr' (pointers) and 'ind' (indices).
    Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param values values of non-zero entries.
    \param ptr column start indices.
    \param ind row indices.
    \warning Input vectors 'values', 'ptr' and 'ind' are empty on exit.
    Moreover 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  Matrix<T, Prop, RowSymSparse, Allocator>::
  Matrix(int i, int j,
	 Vector<T, Storage0, Allocator0>& values,
	 Vector<int, Storage1, Allocator1>& ptr,
	 Vector<int, Storage2, Allocator2>& ind):
    Matrix_SymSparse<T, Prop, RowSymSparse, Allocator>(i, j, values, ptr, ind)
  {
  }


} // namespace Seldon.

#define SELDON_FILE_MATRIX_SYMSPARSE_CXX
#endif
