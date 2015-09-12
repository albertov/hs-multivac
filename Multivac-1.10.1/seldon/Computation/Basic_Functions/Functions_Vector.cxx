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

#ifndef SELDON_FILE_FUNCTIONS_VECTOR_CXX

/*
  Functions defined in this file:

  alpha.U -> U
  Mlt(alpha, U)

  alpha.U + V -> V
  Add(alpha, U, V)

*/

namespace Seldon
{


  /////////
  // Mlt //


  template <class T0,
	    class T1, class Storage1, class Allocator1>
  void Mlt(const T0 alpha,
	   Vector<T1, Storage1, Allocator1>& A)  throw()
  {
    T1 alpha_ = alpha;

    typename Vector<T1, Storage1, Allocator1>::pointer data = A.GetData();

    for (int i = 0; i < A.GetDataSize(); i++)
      data[i] = alpha_ * data[i];
  }  


  // Mlt //
  /////////



  /////////
  // Add //


  template <class T0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void Add(const T0 alpha,
	   const Vector<T1, Storage1, Allocator1>& A,
	   Vector<T2, Storage2, Allocator2>& C)  throw(WrongDim, NoMemory)
  {
    if (alpha != T0(0))
      {
	T1 alpha_ = alpha;

	int ma = A.GetM();
	
#ifdef SELDON_CHECK_BOUNDARIES
	if (ma != C.GetM())
	  throw WrongDim("Add(Scalar, const Vector, Vector)",
			 "Unable to add vectors.");
#endif

	for (int i = 0; i < ma; i++)
	  C(i) += alpha_ * A(i);
      }
  }


  // Add //
  /////////

      
} // namespace Seldon.

#define SELDON_FILE_FUNCTIONS_VECTOR_CXX
#endif
