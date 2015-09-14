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

#ifndef SELDON_FILE_VECTOR_CXX

#include "Vector.hxx"

namespace Seldon
{


  ///////////////////
  // C_VECTOR_BASE //
  ///////////////////


  /****************
   * CONSTRUCTORS *
   ****************/
  

  //! Default constructor.
  /*!
    Nothing is allocated. Vector length is set to zero.
  */
  template <class T, class Allocator>
  inline Vector_Base<T, Allocator>::Vector_Base()
  {
    m_ = 0;
    data_ = NULL;
  }


  //! Main constructor.
  /*!
    \param i length.
    \warning Nothing is allocated.
  */
  template <class T, class Allocator>
  inline Vector_Base<T, Allocator>::Vector_Base(int i)
  {
    m_ = i;
  }


  //! Copy constructor.
  /*!
    \param A base vector to be copied.
    \warning Only the length is copied.
  */
  template <class T, class Allocator>
  inline Vector_Base<T, Allocator>::Vector_Base(Vector_Base<T, Allocator>& A)
  {
    m_ = A.GetM();
  }


  /**************
   * DESTRUCTOR *
   **************/
  

  //! Destructor.
  template <class T, class Allocator>
  inline Vector_Base<T, Allocator>::~Vector_Base()
  {

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (data_ != NULL)
	  {
	    vect_allocator_.deallocate(data_, m_);
	    m_ = 0;
	    data_ = NULL;
	  }
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	m_ = 0;
	data_ = NULL;
      }
#endif

  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/
  

  //! Returns the number of elements.
  /*!
    \return The length of the vector.
  */
  template <class T, class Allocator>
  int Vector_Base<T, Allocator>::GetM() const
  {
    return m_;
  }


  //! Returns the number of elements.
  /*!
    \return The length of the vector.
  */
  template <class T, class Allocator>
  int Vector_Base<T, Allocator>::GetLength() const
  {
    return m_;
  }


  //! Returns the number of elements stored.
  /*!
    \return The length of the vector stored.
  */
  template <class T, class Allocator>
  int Vector_Base<T, Allocator>::GetSize() const
  {
    return m_;
  }


  //! Returns a pointer to data_ (stored data).
  /*!
    \return A pointer to the data_, i.e. the data array.
  */
  template <class T, class Allocator>
  typename Vector_Base<T, Allocator>::pointer
  Vector_Base<T, Allocator>::GetData() const
  {
    return data_;
  }


  //! Returns a const pointer to data_ (stored data).
  /*!
    \return A const pointer to the data_, i.e. the data array.
  */
  template <class T, class Allocator>
  typename Vector_Base<T, Allocator>::const_pointer
  Vector_Base<T, Allocator>::GetDataConst() const
  {
    return reinterpret_cast<typename Vector_Base<T,
      Allocator>::const_pointer>(data_);
  }


  //! Returns a pointer of type "void*" to the data array (data_).
  /*!
    \return A pointer of type "void*" to the data array.
  */
  template <class T, class Allocator>
  void* Vector_Base<T, Allocator>::GetDataVoid() const
  {
    return reinterpret_cast<void*>(data_);
  }


  //! Returns a pointer of type "const void*" to the data array (data_).
  /*!
    \return A pointer of type "const void*" to the data array.
  */
  template <class T, class Allocator>
  const void* Vector_Base<T, Allocator>::GetDataConstVoid() const
  {
    return reinterpret_cast<const void*>(data_);
  }


  ///////////////////////
  // VECTOR<VECT_FULL> //
  ///////////////////////

  
  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the vector is empty.
  */
  template <class T, class Allocator>
  Vector<T, Vect_Full, Allocator>::Vector()  throw():
    Vector_Base<T, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a vector of a given size.
    \param i length of the vector.
  */
  template <class T, class Allocator>
  Vector<T, Vect_Full, Allocator>::Vector(int i):
    Vector_Base<T, Allocator>(i)
  {

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	  
	this->data_ = this->vect_allocator_.allocate(i, this);
	  
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->data_ = NULL;
      }
    if (this->data_ == NULL)
      this->m_ = 0;
    if (this->data_ == NULL && i != 0)
      throw NoMemory("Vector<Vect_Full>::Vector(int)",
		     string("Unable to allocate memory for a vector of size ")
		     + to_str(i*sizeof(T)) + " bytes ("
		     + to_str(i) + " elements).");
#endif
      
  }


  //! Copy constructor.
  /*! Builds a copy of a vector.
    \param V vector to be copied.
  */
  template <class T, class Allocator>
  Vector<T, Vect_Full, Allocator>::Vector(Vector<T, Vect_Full, Allocator>& V):
    Vector_Base<T, Allocator>(V)
  {
      
#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	  
	this->data_ = this->vect_allocator_.allocate(V.GetM(), this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->data_ = NULL;
      }
    if (this->data_ == NULL)
      this->m_ = 0;
    if (this->data_ == NULL && V.GetM() != 0)
      throw NoMemory("Vector<Vect_Full>::Vector(Vector<Vect_Full>&)",
		     string("Unable to allocate memory for a vector of size ")
		     + to_str(V.GetM()*sizeof(T)) + " bytes ("
		     + to_str(V.GetM()) + " elements).");
#endif

    this->vect_allocator_.memorycpy(this->data_, V.GetData(), V.GetM());
	  
  }

  
  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Allocator>
  Vector<T, Vect_Full, Allocator>::~Vector()
  {
  }


  /*********************
   * MEMORY MANAGEMENT *
   *********************/


  //! Clears the vector.
  /*!
    Destructs the vector.
    \warning On exit, the vector is an empty vector.
  */
  template <class T, class Allocator>
  inline void Vector<T, Vect_Full, Allocator>::Clear()
  {
    this->~Vector();
  }


  //! Vector reallocation.
  /*!
    The vector is resized.
    \param i new length of the vector.
    \warning Depending on your allocator, initial elements of the vector may
    be lost.
  */
  template <class T, class Allocator>
  inline void Vector<T, Vect_Full, Allocator>::Reallocate(int i)
  {
    if (i != this->m_)
      {

	this->m_ = i;

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif

	    this->data_ =
	      reinterpret_cast<pointer>(this->vect_allocator_.reallocate(this->data_,
									 i, this) );

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    this->m_ = 0;
	    this->data_ = NULL;
	    return;
	  }
	if (this->data_ == NULL)
	  {
	    this->m_ = 0;
	    return;
	  }
#endif

      }
  }


  //! Changes the length of the vector and sets its data array
  //! (low level method).
  /*!
    Reallocates a vector and sets the new data array. It is useful to create
    a vector from pre-existing data.
    \param i new length of the vector.
    \param data the new data array. 'data' contains the new elements of the
    vector and must therefore contain 'i' elements.
    \warning 'data' has to be used carefully outside the object.
    Unless you use 'Nullify', 'data' will be freed by the destructor,
    which means that 'data' must have been allocated carefully. The vector
    allocator should be compatible.
    \note This method should only be used by advanced users.
  */
  template <class T, class Allocator>
  inline void Vector<T, Vect_Full, Allocator>
  ::SetData(int i, typename Vector<T, Vect_Full, Allocator>::pointer data)
  {
    this->Clear();

    this->m_ = i;

    this->data_ = data;
  }


  //! Clears the vector without releasing memory.
  /*!
    On exit, the vector is empty and the memory has not been released.
    It is useful for low level manipulations on a Vector instance.
    \warning Memory is not released.
  */
  template <class T, class Allocator>
  void Vector<T, Vect_Full, Allocator>::Nullify()
  {
    this->m_ = 0;
    this->data_ = NULL;
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


  //! Access operator.
  /*!
    \param i index.
    \return The value of the vector at 'i'.
  */
  template <class T, class Allocator>
  inline typename Vector<T, Vect_Full, Allocator>::reference
  Vector<T, Vect_Full, Allocator>::operator() (int i)
  {

#ifdef SELDON_CHECK_BOUNDARIES
    if (i < 0 || i >= this->m_)
      throw WrongIndex("Vector<Vect_Full>::operator()",
		       string("Index should be in [0, ") + to_str(this->m_-1)
		       + "], but is equal to " + to_str(i) + ".");
#endif

    return this->data_[i];
  }


  //! Access operator.
  /*!
    \param i index.
    \return The value of the vector at 'i'.
  */
  template <class T, class Allocator>
  inline typename Vector<T, Vect_Full, Allocator>::const_reference
  Vector<T, Vect_Full, Allocator>::operator() (int i) const
  {

#ifdef SELDON_CHECK_BOUNDARIES
    if (i < 0 || i >= this->m_)
      throw WrongIndex("Vector<Vect_Full>::operator()",
		       string("Index should be in [0, ") + to_str(this->m_-1)
		       + "], but is equal to " + to_str(i) + ".");
#endif

    return this->data_[i];
  }


  //! Duplicates a vector (assignment operator).
  /*!
    \param X vector to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline Vector<T, Vect_Full, Allocator>& Vector<T, Vect_Full, Allocator>
  ::operator= (const Vector<T, Vect_Full, Allocator>& X)
  {
    this->Copy(X);

    return *this;
  }


  //! Duplicates a vector.
  /*!
    \param X vector to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline void Vector<T, Vect_Full, Allocator>
  ::Copy(const Vector<T, Vect_Full, Allocator>& X)
  {
    this->Reallocate(X.GetLength());

    this->vect_allocator_.memorycpy(this->data_, X.GetData(), this->m_);
  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of elements stored.
  /*!
    \return The number of elements stored in memory.
  */
  template <class T, class Allocator>
  int Vector<T, Vect_Full, Allocator>::GetDataSize()
  {
    return this->m_;
  }


  /************************
   * CONVENIENT FUNCTIONS *
   ************************/


  //! Sets all elements to zero.
  /*!
    \warning It fills the memory with zeros. If the vector stores complex
    structures, use 'Fill' instead.
  */
  template <class T, class Allocator>
  void Vector<T, Vect_Full, Allocator>::Zero()
  {
    this->vect_allocator_.memoryset(this->data_, char(0),
				    this->GetDataSize() * sizeof(value_type));
  }


  //! Fills the vector with 0, 1, 2, ...
  template <class T, class Allocator>
  void Vector<T, Vect_Full, Allocator>::Fill()
  {
    for (int i = 0; i < this->m_; i++)
      this->data_[i] = i;
  }


  //! Fills the vector with a given value.
  /*!
    \param x value to fill the vector with.
  */
  template <class T, class Allocator>
  template <class T0>
  void Vector<T, Vect_Full, Allocator>::Fill(const T0& x)
  {
    for (int i = 0; i < this->m_; i++)
      this->data_[i] = x;
  }


  //! Fills the vector with a given value.
  /*!
    \param x value to fill the vector with.
  */
  template <class T, class Allocator>
  template <class T0>
  Vector<T, Vect_Full, Allocator>&
  Vector<T, Vect_Full, Allocator>::operator= (const T0& x)
  {
    this->Fill(x);

    return *this;
  }


  //! Fills the vector randomly.
  /*!
    \note The random generator is very basic.
  */
  template <class T, class Allocator>
  void Vector<T, Vect_Full, Allocator>::FillRand()
  {
    srand(time(NULL));
    for (int i = 0; i < this->m_; i++)
      this->data_[i] = rand();
  }


  //! Displays the vector.
  template <class T, class Allocator>
  void Vector<T, Vect_Full, Allocator>::Print() const
  {
    for (int i = 0; i < this->GetLength(); i++)
      cerr << (*this)(i) << "\t";
    cerr << endl;
  }


  /*********
   * NORMS *
   *********/


  //! Returns the infinite norm.
  /*!
    \return The infinite norm.
  */
  template <class T, class Allocator>
  typename Vector<T, Vect_Full, Allocator>::value_type
  Vector<T, Vect_Full, Allocator>::GetNormInf() const
  {
    value_type res = value_type(0);
    for (int i = 0; i < this->GetLength(); i++)
      {
	res = max(res, this->data_[i]);
	res = max(res, -(this->data_[i]));
      }

    return res;
  }


  //! Returns the index of the highest absolute value.
  /*!
    \return The index of the element that has the highest absolute value.
  */
  template <class T, class Allocator>
  int Vector<T, Vect_Full, Allocator>::GetNormInfIndex() const
  {

#ifdef SELDON_CHECK_DIMENSIONS
    if (this->GetLength() == 0)
      throw WrongDim("Vector<Vect_Full>::GetNormInfIndex()",
		     "Vector is null.");
#endif

    value_type res = value_type(0), temp;
    int j = 0;
    for (int i = 0; i < this->GetLength(); i++)
      {
	temp = res;
	res = max(res, this->data_[i]);
	res = max(res, -(this->data_[i]));
	if (temp != res) j = i;
      }

    return j;
  }


  /**************************
   * OUTPUT/INPUT FUNCTIONS *
   **************************/


  //! Writes the vector in a file.
  /*!
    The length of the vector (integer) and all elements of the vector are
    stored in binary format.
    \param FileName file name.
  */
  template <class T, class Allocator>
  void Vector<T, Vect_Full, Allocator>::Write(string FileName) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<Vect_Full>::Write(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Write(FileStream);
    
    FileStream.close();
  }


  //! Writes the vector in a file stream.
  /*!
    The length of the vector (integer) and all elements of the vector are
    stored in binary format.
    \param FileStream file stream.
  */
  template <class T, class Allocator>
  void Vector<T, Vect_Full, Allocator>::Write(ofstream& FileStream) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Vector<Vect_Full>::Write(ofstream& FileStream)",
                    "Stream is not ready.");
#endif

    FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->m_)),
		     sizeof(int));

    FileStream.write(reinterpret_cast<char*>(this->data_),
		     this->m_ * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("Vector<Vect_Full>::Write(ofstream& FileStream)",
                    string("Output operation failed.")
		    + string(" The output file may have been removed")
		    + "or there is no space left on device.");
#endif

  }


  //! Writes the vector in a file.
  /*!
    All elements of the vector are stored in text format. The length is not
    stored.
    \param FileName file name.
  */
  template <class T, class Allocator>
  void Vector<T, Vect_Full, Allocator>::WriteText(string FileName) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<Vect_Full>::WriteText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->WriteText(FileStream);

    FileStream.close();
  }


  //! Writes the vector in a file stream.
  /*!
    All elements of the vector are stored in text format. The length is not
    stored.
    \param FileStream file stream.
  */
  template <class T, class Allocator>
  void Vector<T, Vect_Full, Allocator>::WriteText(ofstream& FileStream) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Vector<Vect_Full>::WriteText(ofstream& FileStream)",
                    "Stream is not ready.");
#endif

    if (this->GetLength() != 0)
      FileStream << (*this)(0);

    for (int i = 1; i < this->GetLength(); i++)
      FileStream << "\t" << (*this)(i);

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("Vector<Vect_Full>::WriteText(ofstream& FileStream)",
                    string("Output operation failed.")
		    + string(" The output file may have been removed")
		    + "or there is no space left on device.");
#endif

  }


  //! Sets the vector from a file.
  /*!
    Sets the vector according to a binary file that stores the length of the
    vector (integer) and all elements.
    \param FileName file name.
  */
  template <class T, class Allocator>
  void Vector<T, Vect_Full, Allocator>::Read(string FileName)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<Vect_Full>::Read(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Read(FileStream);

    FileStream.close();
  }


  //! Sets the vector from a file stream.
  /*!
    Sets the vector according to a binary file stream that stores the length
    of the vector (integer) and all elements.
    \param FileStream file stream.
  */
  template <class T, class Allocator>
  void Vector<T, Vect_Full, Allocator>::Read(ifstream& FileStream)
  {

#ifdef SELDON_CHECK_IO
    // Checks if the strem is ready.
    if (!FileStream.good())
      throw IOError("Vector<Vect_Full>::Read(ifstream& FileStream)",
                    "Stream is not ready.");
#endif

    int new_size;
    FileStream.read(reinterpret_cast<char*>(&new_size), sizeof(int));
    this->Reallocate(new_size);

    FileStream.read(reinterpret_cast<char*>(this->data_),
		    new_size * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was read.
    if (!FileStream.good())
      throw IOError("Vector<Vect_Full>::Read(ifstream& FileStream)",
                    string("Output operation failed.")
		    + string(" The intput file may have been removed")
		    + " or may not contain enough data.");
#endif    

  }


  //! operator<< overloaded for vectors.
  /*!
    \param out output stream.
    \param V vector to be put in the stream.
    \return The updated stream.
  */
  template <class T, class Storage, class Allocator>
  ostream& operator << (ostream& out,
			const Vector<T, Storage, Allocator>& V)
  {
    for (int i = 0; i < V.GetLength() - 1; i++)
      out << V(i) << '\t';
    if (V.GetLength() != 0)
      out << V(V.GetLength() - 1);

    return out;
  }


} // namespace Seldon.

#define SELDON_FILE_VECTOR_CXX
#endif
