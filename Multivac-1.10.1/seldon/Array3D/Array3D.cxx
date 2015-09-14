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

#ifndef SELDON_FILE_ARRAY3D_CXX

#include "Array3D.hxx"

namespace Seldon
{


  /****************
   * CONSTRUCTORS *
   ****************/
  

  //! Default constructor.
  /*!
    On exit, the array is an empty 0x0x0 3D array.
  */
  template <class T, class Allocator>
  inline Array3D<T, Allocator>::Array3D()
  {
    length1_ = 0;
    length2_ = 0;
    length3_ = 0;

    length23_ = 0;

    data_ = NULL;
  }


  //! Main constructor.
  /*! Builds a i x j x k 3D array, but data is not initialized.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
  */
  template <class T, class Allocator>
  inline Array3D<T, Allocator>::Array3D(int i, int j, int k)
  {
    length1_ = i;
    length2_ = j;
    length3_ = k;

    length23_ = length2_ * length3_;

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->data_ = array3D_allocator_.allocate(i*j*k, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->length1_ = 0;
	this->length2_ = 0;
	this->length3_ = 0;
	this->length23_ = 0;
	this->data_ = NULL;
      }
    if (data_ == NULL)
      {
	this->length1_ = 0;
	this->length2_ = 0;
	this->length3_ = 0;
	this->length23_ = 0;
      }
    if (data_ == NULL && i != 0 && j != 0 && k != 0)
      throw NoMemory("Array3D::Array3D(int, int, int)",
		     string("Unable to allocate memory for an array of size ")
		     + to_str(static_cast<long int>(i)
			      * static_cast<long int>(j)
			      * static_cast<long int>(k)
			      * static_cast<long int>(sizeof(T)))
		     + " bytes (" + to_str(i) + " x " + to_str(j)
		     + " x " + to_str(k) + " elements).");
#endif

  }


  /**************
   * DESTRUCTOR *
   **************/
  

  //! Destructor.
  template <class T, class Allocator>
  inline Array3D<T, Allocator>::~Array3D()
  {
    this->length1_ = 0;
    this->length2_ = 0;
    this->length3_ = 0;
    this->length23_ = 0;
    
#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	if (this->data_ != NULL)
	  {
	    array3D_allocator_.deallocate(this->data_,
					  this->length1_ * this->length23_);
	    this->data_ = NULL;
	  }
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->data_ = NULL;
      }
#endif
    
  }
  
  
  /*********************
   * MEMORY MANAGEMENT *
   *********************/


  //! Reallocates memory to resize the 3D array.
  /*!
    On exit, the array is a i x j x k 3D array.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \warning Depending on your allocator, data may be lost.
  */
  template <class T, class Allocator>
  inline void Array3D<T, Allocator>::Reallocate(int i, int j, int k)
  {
    if (i != this->length1_ || j != this->length2_ || k != length3_)
      {
	this->length1_ = i;
	this->length2_ = j;
	this->length3_ = k;

	this->length23_ = j * k;

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif
	    
	    this->data_ = 
	      reinterpret_cast<pointer>(array3D_allocator_.reallocate(this->data_,
								      i * j * k,
								      this));

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    this->length1_ = 0;
	    this->length2_ = 0;
	    this->length3_ = 0;
	    this->length23_ = 0;
	    this->data_ = NULL;
	  }
	if (this->data_ == NULL)
	  {
	    this->length1_ = 0;
	    this->length2_ = 0;
	    this->length3_ = 0;
	    this->length23_ = 0;
	  }
	if (data_ == NULL && i != 0 && j != 0 && k != 0)
	  throw NoMemory("Array3D::Reallocate(int, int, int)",
			 string("Unable to reallocate memory")
			 + " for an array of size "
			 + to_str(static_cast<long int>(i)
				  * static_cast<long int>(j)
				  * static_cast<long int>(k)
				  * static_cast<long int>(sizeof(T)))
			 + " bytes (" + to_str(i) + " x " + to_str(j)
			 + " x " + to_str(k) + " elements).");
#endif

      }
  }


  /*****************
   * BASIC METHODS *
   *****************/
  

  //! Returns the length in dimension #1.
  /*!
    \return The length in dimension #1.
  */
  template <class T, class Allocator>
  int Array3D<T, Allocator>::GetLength1() const
  {
    return this->length1_;
  }


  //! Returns the length in dimension #2.
  /*!
    \return The length in dimension #2.
  */
  template <class T, class Allocator>
  int Array3D<T, Allocator>::GetLength2() const
  {
    return this->length2_;
  }


  //! Returns the length in dimension #3.
  /*!
    \return The length in dimension #3.
  */
  template <class T, class Allocator>
  int Array3D<T, Allocator>::GetLength3() const
  {
    return this->length3_;
  }


  //! Returns the number of elements in the 3D array.
  /*!
    Returns the number of elements stored by the 3D array, i.e.
    the product of the lengths in the three dimensions.
    \return The number of elements in the 3D array.
  */
  template <class T, class Allocator>
  int Array3D<T, Allocator>::GetSize() const
  {
    return this->length1_ * this->length23_;
  }


  //! Returns the number of elements stored in memory.
  /*!
    Returns the number of elements stored in memory by
    the array, i.e. the product of lengths in the three
    dimensions.
    \return The number of elements stored in the array.
  */
  template <class T, class Allocator>
  int Array3D<T, Allocator>::GetDataSize() const
  {
    return this->length1_ * this->length23_;
  }


  //! Returns a pointer to the data array.
  /*!
    Returns a pointer to data, i.e. the data array 'data_' which stores the
    values.
    \return A pointer to the data array.
  */
  template <class T, class Allocator>
  typename Array3D<T, Allocator>::pointer Array3D<T, Allocator>
  ::GetData() const
  {
    return data_;
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


  //! Access operator.
  /*!
    Returns the value of element (i, j, k).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \return Element (i, j, k) of the 3D array.
  */
  template <class T, class Allocator>
  inline typename Array3D<T, Allocator>::reference
  Array3D<T, Allocator>::operator() (int i, int j, int k)
  {

#ifdef SELDON_CHECK_BOUNDARIES
    if (i < 0 || i >= this->length1_)
      throw WrongIndex("Array3D::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(this->length1_-1) + "], but is equal to "
		       + to_str(i) + ".");
    if (j < 0 || j >= this->length2_)
      throw WrongIndex("Array3D::operator()",
		       string("Index along dimension #2 should be in [0, ")
		       + to_str(this->length2_-1) + "], but is equal to "
		       + to_str(j) + ".");
    if (k < 0 || k >= this->length3_)
      throw WrongIndex("Array3D::operator()",
		       string("Index along dimension #3 should be in [0, ")
		       + to_str(this->length3_-1) + "], but is equal to "
		       + to_str(k) + ".");
#endif

    return data_[i * length23_ + j * length3_ + k];
  }

 
  //! Access operator.
  /*!
    Returns the value of element (i, j, k).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \return Element (i, j, k) of the 3D array.
  */
  template <class T, class Allocator>
  inline typename Array3D<T, Allocator>::const_reference
  Array3D<T, Allocator>::operator() (int i, int j, int k) const
  {
    
#ifdef SELDON_CHECK_BOUNDARIES
    if (i < 0 || i >= this->length1_)
      throw WrongIndex("Array3D::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(this->length1_-1) + "], but is equal to "
		       + to_str(i) + ".");
    if (j < 0 || j >= this->length2_)
      throw WrongIndex("Array3D::operator()",
		       string("Index along dimension #2 should be in [0, ")
		       + to_str(this->length2_-1) + "], but is equal to "
		       + to_str(j) + ".");
    if (k < 0 || k >= this->length3_)
      throw WrongIndex("Array3D::operator()",
		       string("Index along dimension #3 should be in [0, ")
		       + to_str(this->length3_-1) + "], but is equal to "
		       + to_str(k) + ".");
#endif

    return data_[i*length23_ + j*length3_ + k];
  }


  //! Duplicates a 3D array.
  /*!
    \param A 3D array to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline void Array3D<T, Allocator>::Copy(const Array3D<T, Allocator>& A)
  {
    this->Reallocate(A.GetLength1(), A.GetLength2(), A.GetLength3());

    array3D_allocator_.memorycpy(this->data_, A.GetData(),
				 this->GetDataSize());
  }


  /************************
   * CONVENIENT FUNCTIONS *
   ************************/


  //! Sets all elements to zero.
  /*!
    \warning It fills the memory with zeros. If the 3D array stores complex
    structures, use 'Fill' instead.
  */
  template <class T, class Allocator>
  void Array3D<T, Allocator>::Zero()
  {
    array3D_allocator_.memoryset(this->data_, char(0),
				 this->GetDataSize()*sizeof(value_type));
  }


  //! Fills the array.
  /*!
    On exit, the 3D array is filled with 1, 2, 3, 4, ... The order of
    those numbers depends on the storage.
  */
  template <class T, class Allocator>
  void Array3D<T, Allocator>::Fill()
  {
    for (int i = 0; i < this->GetDataSize(); i++)
      this->data_[i] = i;
  }


  //! Fills the 3D array with a given value.
  /*!
    On exit, the 3D array is filled with 'x'.
    \param x the value to fill the 3D array with.
  */
  template <class T, class Allocator>
  template <class T0>
  void Array3D<T, Allocator>::Fill(const T0& x)
  {
    for (int i = 0; i < this->GetDataSize(); i++)
      this->data_[i] = x;
  }


  //! Fills the 3D array randomly.
  /*!
    On exit, the 3D array is filled with random values.
  */
  template <class T, class Allocator>
  void Array3D<T, Allocator>::FillRand()
  {
    srand(time(NULL));
    for (int i = 0; i < this->GetDataSize(); i++)
      this->data_[i] = rand();
  }


  //! Displays the array on the standard output.
  /*!
    Displays elements on the standard output, in text format.
  */
  template <class T, class Allocator>
  void Array3D<T, Allocator>::Print() const
  {
    int i, j, k;

    for (i = 0; i < this->GetLength1(); i++)
      {
	for (j = 0; j < this->GetLength2(); j++)
	  {
	    for (k = 0; k < this->GetLength3(); k++)
	      cerr << (*this)(i, j, k) << '\t';
	    cerr << endl;
	  }
	cerr << endl;
      }
  }


  //! operator<< overloaded for a 3D array.
  /*!
    \param out output stream.
    \param A the 3D array.
    \return The updated stream.
  */
  template <class T, class Allocator>
  ostream& operator << (ostream& out,
			const Array3D<T, Allocator>& A)
  {
    int i, j, k;

    for (i = 0; i < A.GetLength1(); i++)
      {
	for (j = 0; j < A.GetLength2(); j++)
	  {
	    for (k = 0; k < A.GetLength3(); k++)
	      out << A(i, j, k) << '\t';
	    out << endl;
	  }
	out << endl;
      }

    return out;
  }


} // namespace Seldon.

#define SELDON_FILE_ARRAY3D_CXX
#endif
