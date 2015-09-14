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

#ifndef SELDON_FILE_MATRIX_HERMPACKED_CXX

#include "Matrix_HermPacked.hxx"

namespace Seldon
{


  /****************
   * CONSTRUCTORS *
   ****************/
  

  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_HermPacked<T, Prop, Storage, Allocator>::Matrix_HermPacked():
    Matrix_Base<T, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j hermitian matrix in packed form.
    \param i number of rows.
    \param j number of columns.
    \note 'j' is assumed to be equal to 'i' and is therefore discarded.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_HermPacked<T, Prop, Storage, Allocator>
  ::Matrix_HermPacked(int i, int j):
    Matrix_Base<T, Allocator>(i, i)
  {

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->data_ = this->allocator_.allocate((i * (i + 1)) / 2, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	this->data_ = NULL;
	return;
      }
    if (this->data_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	return;
      }
#endif

  }

  
  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_HermPacked<T, Prop, Storage, Allocator>::~Matrix_HermPacked()
  {
    
#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	if (this->data_ != NULL)
	  {
	    this->allocator_.deallocate(this->data_,
					(this->m_ * (this->m_ + 1)) / 2);
	    this->data_ = NULL;
	  }
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	this->data_ = NULL;
      }
#endif

  }


  //! Clears the matrix.
  /*!
    Destructs the matrix.
    \warning On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_HermPacked<T, Prop, Storage, Allocator>::Clear()
  {
    this->~Matrix_HermPacked();
  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of elements stored in memory.
  /*!
    \return The number of elements stored in memory.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int Matrix_HermPacked<T, Prop, Storage, Allocator>::GetDataSize() const
  {
    return (this->m_ * (this->m_ + 1)) / 2;
  }


  /*********************
   * MEMORY MANAGEMENT *
   *********************/


  //! Reallocates memory to resize the matrix.
  /*!
    On exit, the matrix is a i x j matrix.
    \param i new number of rows.
    \param j new number of columns.
    \warning Depending on your allocator, data may be lost.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_HermPacked<T, Prop, Storage, Allocator>
  ::Reallocate(int i, int j)
  {
    if (i != this->m_)
      {
	this->m_ = i;
	this->n_ = i;

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif

	    this->data_ =
	      reinterpret_cast<pointer>(this->allocator_.reallocate(this->data_,
								    (i*(i + 1)) / 2,
								    this));

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    this->m_ = 0;
	    this->n_ = 0;
	    this->data_ = NULL;
	    throw NoMemory("Matrix_HermPacked::Reallocate(int, int)",
			   "Unable to reallocate memory for data_.");
	  }
	if (this->data_ == NULL)
	  {
	    this->m_ = 0;
	    this->n_ = 0;
	    throw NoMemory("Matrix_HermPacked::Reallocate(int, int)",
			   "Unable to reallocate memory for data_.");
	  }
#endif

      }
  }


  //! Changes the size of the matrix and sets its data array
  //! (low level method).
  /*!
    The matrix is first cleared (memory is freed). The matrix is then resized
    to a i x j matrix, and the data array of the matrix is set to 'data'.
    'data' elements are not duplicated: the new data array of the matrix is
    the 'data' array. It is useful to create a matrix from pre-existing data.
    \param i new number of rows.
    \param j new number of columns.
    \param data new array storing elements.
    \warning 'data' has to be used carefully outside the object.
    Unless you use 'Nullify', 'data' will be freed by the destructor,
    which means that 'data' must have been allocated carefully. The matrix
    allocator should be compatible.
    \note This method should only be used by advanced users.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_HermPacked<T, Prop, Storage, Allocator>
  ::SetData(int i, int j,
	    typename Matrix_HermPacked<T, Prop, Storage, Allocator>
	    ::pointer data)
  {
    this->Clear();

    this->m_ = i;
    this->n_ = i;

    this->data_ = data;
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
  inline typename Matrix_HermPacked<T, Prop, Storage, Allocator>::value_type
  Matrix_HermPacked<T, Prop, Storage, Allocator>::operator() (int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDARIES
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_HermPacked::operator()",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_HermPacked::operator()",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
#endif

    if (i > j)
      return conj(this->data_[Storage::GetFirst(j * this->m_
						- (j*(j+1)) / 2 + i,
						(i*(i+1)) / 2 + j)]);
    else
      return this->data_[Storage::GetFirst(i * this->n_ - (i*(i+1)) / 2 + j,
					   (j*(j+1)) / 2 + i)];
  }

 
  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_HermPacked<T, Prop, Storage, Allocator>::value_type
  Matrix_HermPacked<T, Prop, Storage, Allocator>
  ::operator() (int i, int j) const
  {
    
#ifdef SELDON_CHECK_BOUNDARIES
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_HermPacked::operator()",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_HermPacked::operator()",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
#endif
    
    if (i > j)
      return conj(this->data_[Storage::GetFirst(j * this->m_
						- (j*(j+1)) / 2 + i,
						(i*(i+1)) / 2 + j)]);
    else
      return this->data_[Storage::GetFirst(i * this->n_ - (i*(i+1)) / 2 + j,
					   (j*(j+1)) / 2 + i)];
  }


  //! Direct access method.
  /*!
    This method allows access to elements stored in memory, i.e. elements
    from the upper part. i <= j must be satisfied.
    \param i row index.
    \param j column index.
    \return The value of the matrix at (i, j).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_HermPacked<T, Prop, Storage, Allocator>::reference
  Matrix_HermPacked<T, Prop, Storage, Allocator>::Val(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDARIES
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_HermPacked::Val(int, int)",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_HermPacked::Val(int, int)",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
    if (i > j)
      throw WrongRow("Matrix_HermPacked::Val(int, int)",
		     string("Attempted to access to element (")
		     + to_str(i) + ", " + to_str(j)
		     + ") but row index should not be strictly"
		     + " more than column index.");
#endif

    return this->data_[Storage::GetFirst(i * this->n_ - (i*(i+1)) / 2 + j,
					 (j*(j+1)) / 2 + i)];
  }

 
  //! Direct access method.
  /*!
    This method allows access to elements stored in memory, i.e. elements
    from the upper part. i <= j must be satisfied.
    \param i row index.
    \param j column index.
    \return The value of the matrix at (i, j).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_HermPacked<T, Prop, Storage, Allocator>
  ::const_reference
  Matrix_HermPacked<T, Prop, Storage, Allocator>::Val(int i, int j) const
  {
    
#ifdef SELDON_CHECK_BOUNDARIES
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_HermPacked::Val(int, int) const",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_HermPacked::Val(int, int) cont",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
    if (i > j)
      throw WrongRow("Matrix_HermPacked::Val(int, int) const",
		     string("Attempted to access to element (")
		     + to_str(i) + ", " + to_str(j)
		     + string(") but row index should not be strictly")
		     + " more than column index.");
#endif
    
    return this->data_[Storage::GetFirst(i * this->n_ - (i*(i+1)) / 2 + j,
					 (j*(j+1)) / 2 + i)];
  }


  //! Access to elements of the data array.
  /*!
    Provides a direct access to the data array.
    \param i index.
    \return i-th element of the data array.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_HermPacked<T, Prop, Storage, Allocator>::reference
  Matrix_HermPacked<T, Prop, Storage, Allocator>::operator[] (int i)
  {
    
#ifdef SELDON_CHECK_BOUNDARIES
    if (i < 0 || i >= this->GetDataSize())
      throw WrongIndex("Matrix_HermPacked::operator[] (int)",
		       string("Index should be in [0, ")
		       + to_str(this->GetDataSize()-1) + "], but is equal to "
		       + to_str(i) + ".");
#endif
    
    return this->data_[i];
  }


  //! Access to elements of the data array.
  /*!
    Provides a direct access to the data array.
    \param i index.
    \return i-th element of the data array.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_HermPacked<T, Prop, Storage, Allocator>
  ::const_reference
  Matrix_HermPacked<T, Prop, Storage, Allocator>::operator[] (int i) const
  {
    
#ifdef SELDON_CHECK_BOUNDARIES
    if (i < 0 || i >= this->GetDataSize())
      throw WrongIndex("Matrix_HermPacked::operator[] (int) const",
		       string("Index should be in [0, ")
		       + to_str(this->GetDataSize()-1) + "], but is equal to "
		       + to_str(i) + ".");
#endif
    
    return this->data_[i];
  }


  //! Duplicates a matrix (assignment operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_HermPacked<T, Prop, Storage, Allocator>&
  Matrix_HermPacked<T, Prop, Storage, Allocator>
  ::operator= (const Matrix_HermPacked<T, Prop, Storage, Allocator>& A)
  {
    this->Copy(A);

    return *this;
  }


  //! Duplicates a matrix.
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_HermPacked<T, Prop, Storage, Allocator>
  ::Copy(const Matrix_HermPacked<T, Prop, Storage, Allocator>& A)
  {
    this->Reallocate(A.GetM(), A.GetN());

    this->allocator_.memorycpy(this->data_, A.GetData(), this->GetDataSize());
  }


  /************************
   * CONVENIENT FUNCTIONS *
   ************************/


  //! Sets all elements to zero.
  /*!
    \warning It fills the memory with zeros. If the matrix stores complex
    structures, use 'Fill' instead.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_HermPacked<T, Prop, Storage, Allocator>::Zero()
  {
    this->allocator_.memoryset(this->data_, char(0),
			       this->GetDataSize() * sizeof(value_type));
  }


  //! Sets the matrix to the identity.
  /*!
    \warning It fills the memory with zeros. If the matrix stores complex
    structures, discard this method.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_HermPacked<T, Prop, Storage, Allocator>::SetIdentity()
  {
    this->allocator_.memoryset(this->data_, char(0),
			       this->GetDataSize() * sizeof(value_type));
    T one(1);
    for (int i = 0; i < min(this->m_, this->n_); i++)
      (*this)(i,i) = one;
  }


  //! Fills the matrix with 0, 1, 2, ...
  /*!
    On exit, the matrix is filled with 0, 1, 2, 3, ... The order of
    those numbers depends on the storage.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_HermPacked<T, Prop, Storage, Allocator>::Fill()
  {
    for (int i = 0; i < this->GetDataSize(); i++)
      this->data_[i] = i;
  }


  //! Fills the matrix with a given value.
  /*!
    \param x value to fill the matrix with.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template <class T0>
  void Matrix_HermPacked<T, Prop, Storage, Allocator>::Fill(const T0& x)
  {
    for (int i = 0; i < this->GetDataSize(); i++)
      this->data_[i] = x;
  }


  //! Fills the matrix with a given value.
  /*!
    \param x value to fill the matrix with.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template <class T0>
  Matrix_HermPacked<T, Prop, Storage, Allocator>&
  Matrix_HermPacked<T, Prop, Storage, Allocator>::operator= (const T0& x)
  {
    this->Fill(x);

    return *this;
  }


  //! Fills the matrix randomly.
  /*!
    \note The random generator is very basic.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_HermPacked<T, Prop, Storage, Allocator>::FillRand()
  {
    srand(time(NULL));
    for (int i = 0; i < this->GetDataSize(); i++)
      this->data_[i] = rand();
  }


  //! Displays the matrix on the standard output.
  /*!
    Displays elements on the standard output, in text format.
    Each row is displayed on a single line and elements of
    a row are delimited by tabulations.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_HermPacked<T, Prop, Storage, Allocator>::Print() const
  {
    for (int i = 0; i < this->m_; i++)
      {
	for (int j = 0; j < this->n_; j++)
	  cerr << (*this)(i, j) << "\t";
	cerr << endl;
      }
  }


  //! Displays a sub-matrix on the standard output.
  /*!
    The sub-matrix is defined by its upper-left corner (a, b)
    and its bottom-right corner (m, n). So, elements with indices
    in [a, m] x [b, n] are displayed on the standard output,
    in text format. Each row is displayed on a single line and
    elements of a row are delimited by tabulations.
    \param a row index of the upper-left corner.
    \param b column index of the upper-left corner.
    \param m row index of the bottom-right corner.
    \param n column index of the bottom-right corner.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_HermPacked<T, Prop, Storage, Allocator>
  ::Print(int a, int b, int m, int n) const
  {
    for (int i = a; i < min(this->m_, a+m); i++)
      {
	for (int j = b; j < min(this->n_, b+n); j++)
	  cerr << (*this)(i, j) << "\t";
	cerr << endl;
      }
  }


  //! Displays a square sub-matrix on the standard output.
  /*!
    The sub-matrix is defined by its bottom-right corner (l, l).
    So, elements with indices in [0, 0] x [l, l] are displayed
    on the standard output, in text format. Each row is displayed
    on a single line and elements of a row are delimited
    by tabulations.
    \param l dimension of the square matrix to be displayed.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_HermPacked<T, Prop, Storage, Allocator>::Print(int l) const
  {
    Print(0, 0, l, l);
  }


  /*********
   * NORMS *
   *********/


  //! Returns the maximum (in absolute value) of the matrix.
  /*!
    \return The maximum (in absolute value) of the matrix.
    \note The name of this method is of course not relevant
    since the infinity norm of the matrix is something else.
    The name of this method will be GetMaxAbs in a next version.
  */
  template <class T, class Prop, class Storage, class Allocator>
  typename Matrix_HermPacked<T, Prop, Storage, Allocator>::value_type
  Matrix_HermPacked<T, Prop, Storage, Allocator>::GetNormInf() const
  {
    value_type res = value_type(0);
    int i, j;
    for (i = 0; i < this->GetM(); i++)
      for (j = 0; j < this->GetN(); j++)
	{
	  res = max(res, (*this)(i, j));
	  res = max(res, -(*this)(i, j));
	}

    return res;
  }


  /**************************
   * INPUT/OUTPUT FUNCTIONS *
   **************************/


  //! Writes the matrix in a file.
  /*!
    Stores the matrix in a file in binary format.
    The number of rows (integer) and the number of columns (integer)
    are written, and matrix elements are then written in the same order
    as in memory (e.g. row-major storage).
    \param FileName output file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_HermPacked<T, Prop, Storage, Allocator>
  ::Write(string FileName) const
  {

    ofstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_HermPacked::Write(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Write(FileStream);

    FileStream.close();
    
  }


  //! Writes the matrix to an output stream.
  /*!
    Writes the matrix to an output stream in binary format.
    The number of rows (integer) and the number of columns (integer)
    are written, and matrix elements are then written in the same order
    as in memory (e.g. row-major storage).
    \param FileStream output stream.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_HermPacked<T, Prop, Storage, Allocator>
  ::Write(ofstream& FileStream) const
  {
    
#ifdef SELDON_CHECK_IO
    // Checks if the file is ready.
    if (!FileStream.good())
      throw IOError("Matrix_HermPacked::Write(ofstream& FileStream)",
                    "Stream is not ready.");
#endif

    FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->m_)),
		     sizeof(int));
    FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->n_)),
		     sizeof(int));

    FileStream.write(reinterpret_cast<char*>(this->data_),
		     this->GetDataSize() * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("Matrix_HermPacked::Write(ofstream& FileStream)",
                    string("Output operation failed.")
		    + string(" The output file may have been removed")
		    + " or there is no space left on device.");
#endif

  }


  //! Writes the matrix in a file.
  /*!
    Stores the matrix in a file in text format.
    Only matrix elements are written (not dimensions).
    Each row is written on a single line and elements of
    a row are delimited by tabulations.
    \param FileName output file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_HermPacked<T, Prop, Storage, Allocator>
  ::WriteText(string FileName) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_HermPacked::WriteText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->WriteText(FileStream);

    FileStream.close();
  }


  //! Writes the matrix to an output stream.
  /*!
    Writes the matrix to an output stream in text format.
    Only matrix elements are written (not dimensions).
    Each row is written on a single line and elements of
    a row are delimited by tabulations.
    \param FileStream output stream.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_HermPacked<T, Prop, Storage, Allocator>
  ::WriteText(ofstream& FileStream) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_HermPacked::WriteText(ofstream& FileStream)",
                    "Stream is not ready.");
#endif

    int i, j;
    for (i = 0; i < this->GetM(); i++)
      {
	for (j = 0; j < this->GetN(); j++)
	  FileStream << (*this)(i, j) << '\t';
	FileStream << endl;
      }

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("Matrix_HermPacked::WriteText(ofstream& FileStream)",
                    string("Output operation failed.")
		    + string(" The output file may have been removed")
		    + " or there is no space left on device.");
#endif

  }


  //! Reads the matrix from a file.
  /*!
    Reads a matrix stored in binary format in a file.
    The number of rows (integer) and the number of columns (integer)
    are read, and matrix elements are then read in the same order
    as it should be in memory (e.g. row-major storage).
    \param FileName input file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_HermPacked<T, Prop, Storage, Allocator>::Read(string FileName)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.good())
      throw IOError("Matrix_HermPacked::Read(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Read(FileStream);

    FileStream.close();
  }
 

  //! Reads the matrix from an input stream.
  /*!
    Reads a matrix in binary format from an input stream.
    The number of rows (integer) and the number of columns (integer)
    are read, and matrix elements are then read in the same order
    as it should be in memory (e.g. row-major storage).
    \param FileStream input stream.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_HermPacked<T, Prop, Storage, Allocator>
  ::Read(ifstream& FileStream)
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_HermPacked::Read(ifstream& FileStream)",
                    "Stream is not ready.");
#endif

    int new_m, new_n;
    FileStream.read(reinterpret_cast<char*>(&new_m), sizeof(int));
    FileStream.read(reinterpret_cast<char*>(&new_n), sizeof(int));
    this->Reallocate(new_m, new_n);

    FileStream.read(reinterpret_cast<char*>(this->data_),
		    this->GetDataSize() * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was read.
    if (!FileStream.good())
      throw IOError("Matrix_HermPacked::Read(ifstream& FileStream)",
                    string("Output operation failed.")
		    + string(" The intput file may have been removed")
		    + " or may not contain enough data.");
#endif    

  }



  ///////////////////////////
  // MATRIX<COLHERMPACKED> //
  ///////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, ColHermPacked, Allocator>::Matrix():
    Matrix_HermPacked<T, Prop, ColHermPacked, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j column-major hermitian matrix in packed form.
    \param i number of rows.
    \param j number of columns.
    \note 'j' is assumed to be equal to 'i' and is therefore discarded.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, ColHermPacked, Allocator>::Matrix(int i, int j):
    Matrix_HermPacked<T, Prop, ColHermPacked, Allocator>(i, j)
  {
  }


  /*******************
   * OTHER FUNCTIONS *
   *******************/


  //! Duplicates a matrix (assignment operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  inline Matrix<T, Prop, ColHermPacked, Allocator>&
  Matrix<T, Prop, ColHermPacked, Allocator>
  ::operator= (const T0& x)
  {
    this->Fill(x);

    return *this;
  }



  ///////////////////////////
  // MATRIX<ROWHERMPACKED> //
  ///////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, RowHermPacked, Allocator>::Matrix():
    Matrix_HermPacked<T, Prop, RowHermPacked, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j row-major hermitian matrix in packed form.
    \param i number of rows.
    \param j number of columns.
    \note 'j' is assumed to be equal to 'i' and is therefore discarded.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, RowHermPacked, Allocator>::Matrix(int i, int j):
    Matrix_HermPacked<T, Prop, RowHermPacked, Allocator>(i, j)
  {
  }


  /*******************
   * OTHER FUNCTIONS *
   *******************/


  //! Duplicates a matrix (assignment operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  inline Matrix<T, Prop, RowHermPacked, Allocator>&
  Matrix<T, Prop, RowHermPacked, Allocator>
  ::operator= (const T0& x)
  {
    this->Fill(x);

    return *this;
  }


} // namespace Seldon.

#define SELDON_FILE_MATRIX_HERMPACKED_CXX
#endif
