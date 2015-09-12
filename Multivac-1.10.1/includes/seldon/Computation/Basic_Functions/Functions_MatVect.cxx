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

#ifndef SELDON_FILE_FUNCTIONS_MATVECT_CXX

/*
  Functions defined in this file:

  A*U -> U
  Mlt(A, U)

  alpha.A*U + beta.V -> V
  MltAdd(alpha, A, U, beta, V)
*/

namespace Seldon
{


  /////////
  // Mlt //


  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void Mlt(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
	   Vector<T1, Storage1, Allocator1>& B)
  {
    int ma = A.GetM();
    int na = A.GetN();

#ifdef SELDON_CHECK_BOUNDARIES
    if (na != B.GetM() || ma != B.GetM())
      throw WrongDim("Mlt(const Matrix, Vector)",
		     "Unable to multiply matrix and vector.");
#endif

    Vector<T1, Storage1, Allocator1> C(B);

    T1 zero(0);
    T1 temp;

    for (int i = 0; i < ma; i++)
      {
	temp = zero;
	for (int j = 0; j < na; j++)
	  temp += A(i, j) * C(j);
	B(i) = temp;
      }
  }


  /*** Sparse matrices ***/


  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void Mlt(const Matrix<T0, Prop0, RowSparse, Allocator0>& A,
	   Vector<T1, Storage1, Allocator1>& C)
  {
    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDARIES
    int na = A.GetN();

    if (na != C.GetM())
      throw WrongDim("Mlt(const Matrix, Vector)",
		     "Unable to multiply matrix and vector.");
#endif

    Vector<T1, Storage1, Allocator1> B(C);

    T1 zero(0);
    T1 temp;

    int* ptr = A.GetPtr();
    int* ind = A.GetInd();
    typename Matrix<T0, Prop0, RowSparse, Allocator0>::pointer
      data = A.GetData();

    for (int i = 0; i < ma; i++)
      {
	temp = zero;
	for (int j = ptr[i]; j < ptr[i+1]; j++)
	  temp += data[j] * B(ind[j]);
	C(i) = temp;
      }
  }


  /*** Complex sparse matrices ***/


  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void Mlt(const Matrix<T0, Prop0, RowComplexSparse, Allocator0>& A,
	   Vector<T1, Storage1, Allocator1>& C)
  {
    int i, j;

    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDARIES
    int na = A.GetN();

    if (na != C.GetM())
      throw WrongDim("Mlt(const Matrix, Vector)",
		     "Unable to multiply matrix and vector.");
#endif

    Vector<T1, Storage1, Allocator1> B(C);

    T1 zero(0);
    T1 temp;

    int* real_ptr = A.GetRealPtr();
    int* imag_ptr = A.GetImagPtr();
    int* real_ind = A.GetRealInd();
    int* imag_ind = A.GetImagInd();
    typename Matrix<T0, Prop0, RowComplexSparse, Allocator0>::pointer
      real_data = A.GetRealData();
    typename Matrix<T0, Prop0, RowComplexSparse, Allocator0>::pointer
      imag_data = A.GetImagData();

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  temp += real_data[j] * B(real_ind[j]);
	C(i) = temp;
      }

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	  temp += complex<T0>(T0(0), imag_data[j]) * B(imag_ind[j]);
	C(i) += temp;
      }
  }


  /*** Symmetric sparse matrices ***/


  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void Mlt(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
	   Vector<T1, Storage1, Allocator1>& C)
  {
    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDARIES
    int na = A.GetN();

    if (na != C.GetM())
      throw WrongDim("Mlt(const Matrix, Vector)",
		     "Unable to multiply matrix and vector.");
#endif

    Vector<T1, Storage1, Allocator1> B(C);

    int i, j;
    T1 zero(0);
    T1 temp;

    int* ptr = A.GetPtr();
    int* ind = A.GetInd();
    typename Matrix<T0, Prop0, RowSymSparse, Allocator0>::pointer
      data = A.GetData();

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = ptr[i]; j < ptr[i + 1]; j++)
	  temp += data[j] * B(ind[j]);
	C(i) = temp;
      }
    for (i = 0; i < ma-1; i++)
      for (j = ptr[i]; j < ptr[i + 1]; j++)
	if (ind[j] != i)
	  C(ind[j]) += data[j] * B(i);
  }


  /*** Symmetric complex sparse matrices ***/


  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void Mlt(const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& A,
	   Vector<T1, Storage1, Allocator1>& C)
  {
    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDARIES
    if (A.GetN() != C.GetM())
      throw WrongDim("Mlt(const Matrix, Vector)",
		     "Unable to multiply matrix and vector.");
#endif

    Vector<T1, Storage1, Allocator1> B(C);

    int i, j;
    T1 zero(0);
    T1 temp;

    int* real_ptr = A.GetRealPtr();
    int* imag_ptr = A.GetImagPtr();
    int* real_ind = A.GetRealInd();
    int* imag_ind = A.GetImagInd();
    typename Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>::pointer
      real_data = A.GetRealData();
    typename Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>::pointer
      imag_data = A.GetImagData();

    for (i = 0; i<ma; i++)
      {
	temp = zero;
	for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  temp += real_data[j] * B(real_ind[j]);
	C(i) = temp;
      }
    for (i = 0; i<ma-1; i++)
      for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	if (real_ind[j] != i)
	  C(real_ind[j]) += real_data[j] * B(i);

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	  temp += complex<T0>(T0(0), imag_data[j]) * B(imag_ind[j]);
	C(i) += temp;
      }
    for (i = 0; i<ma-1; i++)
      for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	if (imag_ind[j] != i)
	  C(imag_ind[j]) += complex<T0>(T0(0), imag_data[j]) * B(i);
  }


  /*** Sparse matrices, *Trans ***/


  // NoTrans.
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void Mlt(const class_SeldonNoTrans& Trans,
	   const Matrix<T0, Prop0, RowSparse, Allocator0>& A,
	   Vector<T1, Storage1, Allocator1>& C)
  {
    Mlt(A, C);
  }


  // Trans.
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void Mlt(const class_SeldonTrans& Trans,
	   const Matrix<T0, Prop0, RowSparse, Allocator0>& A,
	   Vector<T1, Storage1, Allocator1>& C)
  {
    int i, j;

    int ma = A.GetM();
    int na = A.GetN();

#ifdef SELDON_CHECK_BOUNDARIES
    if (ma != C.GetM())
      throw WrongDim("Mlt(class_SeldonTrans, const Matrix, Vector)",
		     "Unable to multiply matrix and vector.");
#endif

    Vector<T1, Storage1, Allocator1> B(C);
    C.Zero();

    int* ptr = A.GetPtr();
    int* ind = A.GetInd();
    typename Matrix<T0, Prop0, RowSparse, Allocator0>::pointer
      data = A.GetData();

    for (i = 0; i < ma; i++)
      for (j = ptr[i]; j < ptr[i + 1]; j++)
	C(ind[j]) += data[j] * B(i);
  }


  // ConjTrans.
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void Mlt(const class_SeldonConjTrans& Trans,
	   const Matrix<T0, Prop0, RowSparse, Allocator0>& A,
	   Vector<T1, Storage1, Allocator1>& C)
  {
    int i, j;

    int ma = A.GetM();
    int na = A.GetN();

#ifdef SELDON_CHECK_BOUNDARIES
    if (ma != C.GetM())
      throw WrongDim("Mlt(class_SeldonConjTrans, const Matrix, Vector)",
		     "Unable to multiply matrix and vector.");
#endif

    Vector<T1, Storage1, Allocator1> B(C);
    C.Zero();

    int* ptr = A.GetPtr();
    int* ind = A.GetInd();
    typename Matrix<complex<T0>, Prop0, RowSparse, Allocator0>::pointer
      data = A.GetData();

    for (i = 0; i < ma; i++)
      for (j = ptr[i]; j < ptr[i + 1]; j++)
	C(ind[j]) += conj(data[j]) * B(i);
  }


  /*** Complex sparse matrices, *Trans ***/


  // NoTrans.
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void Mlt(const class_SeldonNoTrans& Trans,
	   const Matrix<T0, Prop0, RowComplexSparse, Allocator0>& A,
	   Vector<T1, Storage1, Allocator1>& C)
  {
    Mlt(A, C);
  }


  // Trans.
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void Mlt(const class_SeldonTrans& Trans,
	   const Matrix<T0, Prop0, RowComplexSparse, Allocator0>& A,
	   Vector<T1, Storage1, Allocator1>& C)
  {
    int i, j;

    int ma = A.GetM();
    int na = A.GetN();

#ifdef SELDON_CHECK_BOUNDARIES
    if (na != C.GetM())
      throw WrongDim("Mlt(class_SeldonTrans, const Matrix, Vector)",
		     "Unable to multiply matrix and vector.");
#endif

    Vector<T1, Storage1, Allocator1> B(C);
    C.Zero();

    int* real_ptr = A.GetRealPtr();
    int* imag_ptr = A.GetImagPtr();
    int* real_ind = A.GetRealInd();
    int* imag_ind = A.GetImagInd();
    typename Matrix<T0, Prop0, RowComplexSparse, Allocator0>::pointer
      real_data = A.GetRealData();
    typename Matrix<T0, Prop0, RowComplexSparse, Allocator0>::pointer
      imag_data = A.GetImagData();

    for (i = 0; i < ma; i++)
      for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	C(real_ind[j]) += real_data[j] * B(i);

    for (i = 0; i < ma; i++)
      for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	C(imag_ind[j]) += complex<T0>(T0(0), imag_data[j]) * B(i);
  }


  // ConjTrans.
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void Mlt(const class_SeldonConjTrans& Trans,
	   const Matrix<T0, Prop0, RowComplexSparse, Allocator0>& A,
	   Vector<T1, Storage1, Allocator1>& C)
  {
    int i, j;

    int ma = A.GetM();
    int na = A.GetN();

#ifdef SELDON_CHECK_BOUNDARIES
    if (na != C.GetM())
      throw WrongDim("Mlt(class_SeldonConjTrans, const Matrix, Vector)",
		     "Unable to multiply matrix and vector.");
#endif

    Vector<T1, Storage1, Allocator1> B(C);
    C.Zero();

    int* real_ptr = A.GetRealPtr();
    int* imag_ptr = A.GetImagPtr();
    int* real_ind = A.GetRealInd();
    int* imag_ind = A.GetImagInd();
    typename Matrix<T0, Prop0, RowComplexSparse, Allocator0>::pointer
      real_data = A.GetRealData();
    typename Matrix<T0, Prop0, RowComplexSparse, Allocator0>::pointer
      imag_data = A.GetImagData();

    for (i = 0; i < ma; i++)
      for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	C(real_ind[j]) += real_data[j] * B(i);

    for (i = 0; i < ma; i++)
      for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	C(imag_ind[j]) += complex<T0>(T0(0), - imag_data[j]) * B(i);
  }


  /*** Symmetric sparse matrices, *Trans ***/


  // NoTrans.
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void Mlt(const class_SeldonNoTrans& Trans,
	   const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
	   Vector<T1, Storage1, Allocator1>& C)
  {
    Mlt(A, C);
  }


  // Trans.
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void Mlt(const class_SeldonTrans& Trans,
	   const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
	   Vector<T1, Storage1, Allocator1>& C)
  {
    Mlt(A, C);
  }


  // ConjTrans.
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void Mlt(const class_SeldonConjTrans& Trans,
	   const Matrix<complex<T0>, Prop0, RowSymSparse, Allocator0>& A,
	   Vector<T1, Storage1, Allocator1>& C)
  {
    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDARIES
    int na = A.GetN();

    if (na != C.GetM())
      throw WrongDim("Mlt(class_SeldonConjTrans, const Matrix, Vector)",
		     "Unable to multiply matrix and vector.");
#endif

    Vector<T1, Storage1, Allocator1> B(C);

    int i, j;
    T1 zero(0);
    T1 temp;

    int* ptr = A.GetPtr();
    int* ind = A.GetInd();
    typename Matrix<complex<T0>, Prop0, RowSymSparse, Allocator0>::pointer
      data = A.GetData();

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = ptr[i]; j < ptr[i + 1]; j++)
	  temp += conj(data[j]) * B(ind[j]);
	C(i) = temp;
      }
    for (i = 0; i < ma - 1; i++)
      for (j = ptr[i]; j < ptr[i + 1]; j++)
	if (ind[j] != i)
	  C(ind[j]) += conj(data[j]) * B(i);
  }


  /*** Symmetric complex sparse matrices, *Trans ***/


  // NoTrans.
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void Mlt(const class_SeldonNoTrans& Trans,
	   const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& A,
	   Vector<T1, Storage1, Allocator1>& C)
  {
    Mlt(A, C);
  }


  // Trans.
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void Mlt(const class_SeldonTrans& Trans,
	   const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& A,
	   Vector<T1, Storage1, Allocator1>& C)
  {
    Mlt(A, C);
  }


  // ConjTrans.
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void Mlt(const class_SeldonConjTrans& Trans,
	   const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& A,
	   Vector<T1, Storage1, Allocator1>& C)
  {
    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDARIES
    if (A.GetN() != C.GetM())
      throw WrongDim("Mlt(class_SeldonConjTrans, const Matrix, Vector)",
		     "Unable to multiply matrix and vector.");
#endif

    Vector<T1, Storage1, Allocator1> B(C);

    int i, j;
    T1 zero(0);
    T1 temp;

    int* real_ptr = A.GetRealPtr();
    int* imag_ptr = A.GetImagPtr();
    int* real_ind = A.GetRealInd();
    int* imag_ind = A.GetImagInd();
    typename Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>::pointer
      real_data = A.GetRealData();
    typename Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>::pointer
      imag_data = A.GetImagData();

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  temp += real_data[j] * B(real_ind[j]);
	C(i) = temp;
      }
    for (i = 0; i < ma - 1; i++)
      for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	if (real_ind[j] != i)
	  C(real_ind[j]) += real_data[j] * B(i);

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	  temp += complex<T0>(T0(0), - imag_data[j]) * B(imag_ind[j]);
	C(i) += temp;
      }
    for (i = 0; i < ma - 1; i++)
      for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	if (imag_ind[j] != i)
	  C(imag_ind[j]) += complex<T0>(T0(0), - imag_data[j]) * B(i);
  }


  // Mlt //
  /////////



  ////////////
  // MltAdd //


  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, Storage1, Allocator1>& A,
	      const Vector<T2, Storage2, Allocator2>& B,
	      const T3 beta,
	      Vector<T4, Storage4, Allocator4>& C)
  {
    int ma = A.GetM();
    int na = A.GetN();

#ifdef SELDON_CHECK_BOUNDARIES
    if (na != B.GetM() || ma != C.GetM())
      throw WrongDim(string("Mlt(Scalar, const Matrix, ")
		     + "const Vector, Scalar, Vector)",
		     string("Unable to multiply matrix and vector")
		     + " or to add vectors.");
#endif

    T4 zero(0);
    T4 temp;
    T4 alpha_(alpha);
    T4 beta_(beta);

    for (int i = 0; i < ma; i++)
      {
	temp = zero;
	for (int j = 0; j < na; j++)
	  temp += A(i, j) * B(j);
	C(i) = beta_ * C(i) + alpha_ * temp;
      }
  }


  // MltAdd //
  ////////////



  ///////////
  // Gauss //


  // Solve X = A*Y with Gauss method.
  // Warning: A is modified. The results are stored in X.
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  inline void Gauss(Matrix<T0, Prop0, Storage0, Allocator0>& A,
		    Vector<T1, Storage1, Allocator1>& X)
  {
    int i, j, k;
    T1 r, S;

    int ma = A.GetM();
    int na = A.GetN();

#ifdef SELDON_CHECK_BOUNDARIES
    if (na != ma)
      throw WrongDim("Gauss(const Matrix, const Vector, Vector)",
		     "Matrix must be squared.");
    if (ma != X.GetLength())
      throw WrongDim("Gauss(const Matrix, const Vector, Vector)",
		     "Matrix and vector dimensions are incompatible.");
#endif

    for (k = 0; k < ma - 1; k++)
      for (i = k + 1; i < ma; i++)
	{
	  r = A(i, k) / A(k, k);
	  for (j = k + 1; j < ma; j++)
	    A(i, j) -= r * A(k, j);
	  X(i) -= r *= X(k);
	}

    X(ma - 1) = X(ma - 1) / A(ma - 1, ma - 1);
    for (k = ma - 2; k > -1; k--)
      {
	S = X(k);
	for (j = k + 1; j < ma; j++)
	  S -= A(k, j) * X(j);
	X(k) = S / A(k, k);
      }
  }


  // Gauss //
  ///////////



  ////////////////////
  // Gauss - Seidel //


  // Solve X = A*Y with Gauss-Seidel method.
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  inline void GaussSeidel(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
			  const Vector<T1, Storage1, Allocator1>& X,
			  Vector<T2, Storage2, Allocator2>& Y,
			  int iter)
  {
    int i, j, k;
    T1 temp;

    int ma = A.GetM();
    int na = A.GetN();

#ifdef SELDON_CHECK_BOUNDARIES
    if (na != ma)
      throw WrongDim("GaussSeidel(const Matrix, const Vector, Vector)",
		     "Matrix must be squared.");
    if (ma != X.GetLength() || ma != Y.GetLength())
      throw WrongDim("GaussSeidel(const Matrix, const Vector, Vector)",
		     "Matrix and vector dimensions are incompatible.");
#endif
    
    for (i = 0; i < iter; i++)
      for (j = 0; j < na; j++)
	{
	  temp = 0;
	  for (k = 0; k < j; k++)
	    temp -= A(j, k) * Y(k);
	  for (k = j + 1; k < na; k++)
	    temp -= A(j, k) * Y(k);
	  Y(j) = (X(j) + temp) / A(j, j);
	}
  }


  // Gauss-Seidel //
  //////////////////



  ///////////////////
  // S.O.R. method //


  // Solve X = A*Y with S.O.R. method.
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3>
  inline void SOR(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
		  const Vector<T1, Storage1, Allocator1>& X,
		  Vector<T2, Storage2, Allocator2>& Y,
		  T3 omega,
		  int iter)
  {
    int i, j, k;
    T1 temp;

    int ma = A.GetM();
    int na = A.GetN();

#ifdef SELDON_CHECK_BOUNDARIES
    if (na != ma)
      throw WrongDim("GaussSeidel(const Matrix, const Vector, Vector)",
		     "Matrix must be squared.");
    if (ma != X.GetLength() || ma != Y.GetLength())
      throw WrongDim("GaussSeidel(const Matrix, const Vector, Vector)",
		     "Matrix and vector dimensions are incompatible.");
#endif
    
    for (i = 0; i < iter; i++)
      for (j = 0; j < na; j++)
	{
	  temp = 0;
	  for (k = 0; k < j; k++)
	    temp -= A(j, k) * Y(k);
	  for (k = j + 1; k < na; k++)
	    temp -= A(j, k) * Y(k);
	  Y(j) = (T3(1) - omega) * Y(j) + omega * (X(j) + temp) / A(j, j);
	}
  }


  // Gauss-Seidel //
  //////////////////



  /////////////
  // SolveLU //


  // Solves A.X = V where A has been decomposed in a LU form.
  // V is overwritten (V <- X).
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SolveLU(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
	       Vector<T1, Storage1, Allocator1>& V)
  {
    int i, k;
    T1 temp;

    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDARIES
    int na = A.GetN();
    if (na != ma)
      throw WrongDim("SolveLU(const Matrix, Vector)",
		     "Matrix must be squared.");
    if (ma != V.GetM())
      throw WrongDim("SolveLu(const Matrix, Vector)",
		     "Matrix and vector dimensions are incompatible.");
#endif

    // Forward substitution.
    for (i = 0; i < ma; i++)
      {
	temp = 0;
	for (k = 0; k < i; k++)
	  temp += A(i, k) * V(k);
	V(i) = (V(i) - temp) / A(i, i);
      }
    // Back substitution.
    for (i = ma - 2; i > -1; i--)
      {
	temp = 0;
	for (k = i + 1; k < ma; k++)
	  temp += A(i, k) * V(k);
	V(i) -= temp;
      }
  }


  // SolveLU //
  /////////////


}  // namespace Seldon.

#define SELDON_FUNCTIONS_MATVECT_CXX
#endif
