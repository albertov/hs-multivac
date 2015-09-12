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

#ifndef SELDON_FILE_FUNCTIONS_MATRIX_CXX

/*
  Function defined in this file:

  GetLU(A)
*/

namespace Seldon
{


  ///////////
  // GetLU //


  // Returns the LU decomposition of A = LU (in A)
  // where L diagonal elements are set to unit value.
  template <class T0, class Prop0, class Storage0, class Allocator0>
  void GetLU(Matrix<T0, Prop0, Storage0, Allocator0>& A)
  {
    int i, p, q, k;
    T0 temp;

    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDARIES
    int na = A.GetN();
    if (na != ma)
      throw WrongDim("GetLU(Matrix)", "Matrix must be squared.");
#endif

    for (i = 0; i < ma; i++)
      {
	for (p = i; p < ma; p++)
	  {
	    temp = 0;
	    for (k = 0; k < i; k++)
	      temp += A(p, k) * A(k, i);
	    A(p, i) -= temp;
	  }
	for (q = i+1; q < ma; q++)
	  {
	    temp = 0;
	    for (k = 0; k < i; k++)
	      temp += A(i, k) * A(k, q);
	    A(i, q) = (A(i,q ) - temp) / A(i, i);
	  }
      }
  }
  
  
  // GetLU //
  ///////////


} // namespace Seldon.

#define SELDON_FILE_FUNCTIONS_MATRIX_CXX
#endif
