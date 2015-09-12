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


#ifndef SELDON_FILE_BLAS_3_CXX

extern "C"
{
#include "cblas.h"
}

namespace Seldon
{


  ////////////
  // MltAdd //


  /*** ColMajor and NoTrans ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const float alpha,
	      const Matrix<float, Prop0, ColMajor, Allocator0>& A,
	      const Matrix<float, Prop1, ColMajor, Allocator1>& B,
	      const float beta,
	      const Matrix<float, Prop2, ColMajor, Allocator2>& C)
  {
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
		A.GetM(), B.GetN(), A.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM(),
		beta, C.GetData(), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const double alpha,
	      const Matrix<double, Prop0, ColMajor, Allocator0>& A,
	      const Matrix<double, Prop1, ColMajor, Allocator1>& B,
	      const double beta,
	      const Matrix<double, Prop2, ColMajor, Allocator2>& C)
  {
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
		A.GetM(), B.GetN(), A.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM(),
		beta, C.GetData(), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B,
	      const complex<float> beta,
	      const Matrix<complex<float>, Prop2, ColMajor, Allocator2>& C)
  {
    cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
		A.GetM(), B.GetN(), A.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetM(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B,
	      const complex<double> beta,
	      const Matrix<complex<double>, Prop2, ColMajor, Allocator2>& C)
  {
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
		A.GetM(), B.GetN(), A.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetM(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetM());
  }


  /*** ColMajor and TransA, TransB ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const float alpha,
	      const enum CBLAS_TRANSPOSE TransA,
	      const Matrix<float, Prop0, ColMajor, Allocator0>& A,
	      const enum CBLAS_TRANSPOSE TransB,
	      const Matrix<float, Prop1, ColMajor, Allocator1>& B,
	      const float beta,
	      const Matrix<float, Prop2, ColMajor, Allocator2>& C)
  {
    cblas_sgemm(CblasColMajor, TransA, TransB, A.GetM(), B.GetN(), A.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM(),
		beta, C.GetData(), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const double alpha,
	      const enum CBLAS_TRANSPOSE TransA,
	      const Matrix<double, Prop0, ColMajor, Allocator0>& A,
	      const enum CBLAS_TRANSPOSE TransB,
	      const Matrix<double, Prop1, ColMajor, Allocator1>& B,
	      const double beta,
	      const Matrix<double, Prop2, ColMajor, Allocator2>& C)
  {
    cblas_dgemm(CblasColMajor, TransA, TransB, A.GetM(), B.GetN(), A.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM(),
		beta, C.GetData(), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const enum CBLAS_TRANSPOSE TransA,
	      const Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
	      const enum CBLAS_TRANSPOSE TransB,
	      const Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B,
	      const complex<float> beta,
	      const Matrix<complex<float>, Prop2, ColMajor, Allocator2>& C)
  {
    cblas_cgemm(CblasColMajor, TransA, TransB, A.GetM(), B.GetN(), A.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetM(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const enum CBLAS_TRANSPOSE TransA,
	      const Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
	      const enum CBLAS_TRANSPOSE TransB,
	      const Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B,
	      const complex<double> beta,
	      const Matrix<complex<double>, Prop2, ColMajor, Allocator2>& C)
  {
    cblas_zgemm(CblasColMajor, TransA, TransB, A.GetM(), B.GetN(), A.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetM(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetM());
  }


  /*** RowMajor and NoTrans ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const float alpha,
	      const Matrix<float, Prop0, RowMajor, Allocator0>& A,
	      const Matrix<float, Prop1, RowMajor, Allocator1>& B,
	      const float beta,
	      const Matrix<float, Prop2, RowMajor, Allocator2>& C)
  {
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		A.GetM(), B.GetN(), A.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM(),
		beta, C.GetData(), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const double alpha,
	      const Matrix<double, Prop0, RowMajor, Allocator0>& A,
	      const Matrix<double, Prop1, RowMajor, Allocator1>& B,
	      const double beta,
	      const Matrix<double, Prop2, RowMajor, Allocator2>& C)
  {
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		A.GetM(), B.GetN(), A.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM(),
		beta, C.GetData(), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B,
	      const complex<float> beta,
	      const Matrix<complex<float>, Prop2, RowMajor, Allocator2>& C)
  {
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		A.GetM(), B.GetN(), A.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetM(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B,
	      const complex<double> beta,
	      const Matrix<complex<double>, Prop2, RowMajor, Allocator2>& C)
  {
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		A.GetM(), B.GetN(), A.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetM(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetM());
  }


  /*** RowMajor and TransA, TransB ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const float alpha,
	      const enum CBLAS_TRANSPOSE TransA,
	      const Matrix<float, Prop0, RowMajor, Allocator0>& A,
	      const enum CBLAS_TRANSPOSE TransB,
	      const Matrix<float, Prop1, RowMajor, Allocator1>& B,
	      const float beta,
	      const Matrix<float, Prop2, RowMajor, Allocator2>& C)
  {
    cblas_sgemm(CblasRowMajor, TransA, TransB, A.GetM(), B.GetN(), A.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM(),
		beta, C.GetData(), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const double alpha,
	      const enum CBLAS_TRANSPOSE TransA,
	      const Matrix<double, Prop0, RowMajor, Allocator0>& A,
	      const enum CBLAS_TRANSPOSE TransB,
	      const Matrix<double, Prop1, RowMajor, Allocator1>& B,
	      const double beta,
	      const Matrix<double, Prop2, RowMajor, Allocator2>& C)
  {
    cblas_dgemm(CblasRowMajor, TransA, TransB, A.GetM(), B.GetN(), A.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM(),
		beta, C.GetData(), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const enum CBLAS_TRANSPOSE TransA,
	      const Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
	      const enum CBLAS_TRANSPOSE TransB,
	      const Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B,
	      const complex<float> beta,
	      const Matrix<complex<float>, Prop2, RowMajor, Allocator2>& C)
  {
    cblas_cgemm(CblasRowMajor, TransA, TransB, A.GetM(), B.GetN(), A.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetM(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const enum CBLAS_TRANSPOSE TransA,
	      const Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
	      const enum CBLAS_TRANSPOSE TransB,
	      const Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B,
	      const complex<double> beta,
	      const Matrix<complex<double>, Prop2, RowMajor, Allocator2>& C)
  {
    cblas_zgemm(CblasRowMajor, TransA, TransB, A.GetM(), B.GetN(), A.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetM(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetM());
  }


  // MltAdd //
  ////////////



  ////////////
  // MltAdd //


  /*** ColSym and Upper ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const enum CBLAS_SIDE Side,
	      const float alpha,
	      const Matrix<float, Prop0, ColSym, Allocator0>& A,
	      const Matrix<float, Prop1, ColMajor, Allocator1>& B,
	      const float beta,
	      const Matrix<float, Prop2, ColMajor, Allocator2>& C)
  {
    cblas_ssymm(CblasColMajor, Side, CblasUpper,
		C.GetM(), C.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM(),
		beta, C.GetData(), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const enum CBLAS_SIDE Side,
	      const double alpha,
	      const Matrix<double, Prop0, ColSym, Allocator0>& A,
	      const Matrix<double, Prop1, ColMajor, Allocator1>& B,
	      const double beta,
	      const Matrix<double, Prop2, ColMajor, Allocator2>& C)
  {
    cblas_dsymm(CblasColMajor, Side, CblasUpper,
		C.GetM(), C.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM(),
		beta, C.GetData(), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const enum CBLAS_SIDE Side,
	      const complex<float> alpha,
	      const Matrix<complex<float>, Prop0, ColSym, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B,
	      const complex<float> beta,
	      const Matrix<complex<float>, Prop2, ColMajor, Allocator2>& C)
  {
    cblas_csymm(CblasColMajor, Side, CblasUpper,
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetM(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const enum CBLAS_SIDE Side,
	      const complex<double> alpha,
	      const Matrix<complex<double>, Prop0, ColSym, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B,
	      const complex<double> beta,
	      const Matrix<complex<double>, Prop2, ColMajor, Allocator2>& C)
  {
    cblas_zsymm(CblasColMajor, Side, CblasUpper,
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetM(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetM());
  }


  /*** ColSym and UpLo ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const enum CBLAS_SIDE Side,
	      const float alpha,
	      const enum CBLAS_UPLO Uplo,
	      const Matrix<float, Prop0, ColSym, Allocator0>& A,
	      const Matrix<float, Prop1, ColMajor, Allocator1>& B,
	      const float beta,
	      const Matrix<float, Prop2, ColMajor, Allocator2>& C)
  {
    cblas_ssymm(CblasColMajor, Side, Uplo,
		C.GetM(), C.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM(),
		beta, C.GetData(), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const enum CBLAS_SIDE Side,
	      const double alpha,
	      const enum CBLAS_UPLO Uplo,
	      const Matrix<double, Prop0, ColSym, Allocator0>& A,
	      const Matrix<double, Prop1, ColMajor, Allocator1>& B,
	      const double beta,
	      const Matrix<double, Prop2, ColMajor, Allocator2>& C)
  {
    cblas_dsymm(CblasColMajor, Side, Uplo,
		C.GetM(), C.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM(),
		beta, C.GetData(), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const enum CBLAS_SIDE Side,
	      const complex<float> alpha,
	      const enum CBLAS_UPLO Uplo,
	      const Matrix<complex<float>, Prop0, ColSym, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B,
	      const complex<float> beta,
	      const Matrix<complex<float>, Prop2, ColMajor, Allocator2>& C)
  {
    cblas_csymm(CblasColMajor, Side, Uplo,
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetM(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const enum CBLAS_SIDE Side,
	      const complex<double> alpha,
	      const enum CBLAS_UPLO Uplo,
	      const Matrix<complex<double>, Prop0, ColSym, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B,
	      const complex<double> beta,
	      const Matrix<complex<double>, Prop2, ColMajor, Allocator2>& C)
  {
    cblas_zsymm(CblasColMajor, Side, Uplo,
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetM(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetM());
  }


  /*** RowSym and Upper ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const enum CBLAS_SIDE Side,
	      const float alpha,
	      const Matrix<float, Prop0, RowSym, Allocator0>& A,
	      const Matrix<float, Prop1, RowMajor, Allocator1>& B,
	      const float beta,
	      const Matrix<float, Prop2, RowMajor, Allocator2>& C)
  {
    cblas_ssymm(CblasRowMajor, Side, CblasUpper,
		C.GetM(), C.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetN(),
		beta, C.GetData(), C.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const enum CBLAS_SIDE Side,
	      const double alpha,
	      const Matrix<double, Prop0, RowSym, Allocator0>& A,
	      const Matrix<double, Prop1, RowMajor, Allocator1>& B,
	      const double beta,
	      const Matrix<double, Prop2, RowMajor, Allocator2>& C)
  {
    cblas_dsymm(CblasRowMajor, Side, CblasUpper,
		C.GetM(), C.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetN(),
		beta, C.GetData(), C.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const enum CBLAS_SIDE Side,
	      const complex<float> alpha,
	      const Matrix<complex<float>, Prop0, RowSym, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B,
	      const complex<float> beta,
	      const Matrix<complex<float>, Prop2, RowMajor, Allocator2>& C)
  {
    cblas_csymm(CblasRowMajor, Side, CblasUpper,
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetN(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const enum CBLAS_SIDE Side,
	      const complex<double> alpha,
	      const Matrix<complex<double>, Prop0, RowSym, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B,
	      const complex<double> beta,
	      const Matrix<complex<double>, Prop2, RowMajor, Allocator2>& C)
  {
    cblas_zsymm(CblasRowMajor, Side, CblasUpper,
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetN(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetN());
  }


  /*** RowSym and UpLo ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const enum CBLAS_SIDE Side,
	      const float alpha,
	      const enum CBLAS_UPLO Uplo,
	      const Matrix<float, Prop0, RowSym, Allocator0>& A,
	      const Matrix<float, Prop1, RowMajor, Allocator1>& B,
	      const float beta,
	      const Matrix<float, Prop2, RowMajor, Allocator2>& C)
  {
    cblas_ssymm(CblasRowMajor, Side, Uplo,
		C.GetM(), C.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetN(),
		beta, C.GetData(), C.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const enum CBLAS_SIDE Side,
	      const double alpha,
	      const enum CBLAS_UPLO Uplo,
	      const Matrix<double, Prop0, RowSym, Allocator0>& A,
	      const Matrix<double, Prop1, RowMajor, Allocator1>& B,
	      const double beta,
	      const Matrix<double, Prop2, RowMajor, Allocator2>& C)
  {
    cblas_dsymm(CblasRowMajor, Side, Uplo,
		C.GetM(), C.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetN(),
		beta, C.GetData(), C.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const enum CBLAS_SIDE Side,
	      const complex<float> alpha,
	      const enum CBLAS_UPLO Uplo,
	      const Matrix<complex<float>, Prop0, RowSym, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B,
	      const complex<float> beta,
	      const Matrix<complex<float>, Prop2, RowMajor, Allocator2>& C)
  {
    cblas_csymm(CblasRowMajor, Side, Uplo,
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetN(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const enum CBLAS_SIDE Side,
	      const complex<double> alpha,
	      const enum CBLAS_UPLO Uplo,
	      const Matrix<complex<double>, Prop0, RowSym, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B,
	      const complex<double> beta,
	      const Matrix<complex<double>, Prop2, RowMajor, Allocator2>& C)
  {
    cblas_zsymm(CblasRowMajor, Side, Uplo,
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetN(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetN());
  }


  // MltAdd //
  ////////////



  ////////////
  // MltAdd //


  /*** ColHerm and Upper ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const enum CBLAS_SIDE Side,
	      const complex<float> alpha,
	      const Matrix<complex<float>, Prop0, ColHerm, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B,
	      const complex<float> beta,
	      const Matrix<complex<float>, Prop2, ColMajor, Allocator2>& C)
  {
    cblas_chemm(CblasColMajor, Side, CblasUpper,
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetM(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const enum CBLAS_SIDE Side,
	      const complex<double> alpha,
	      const Matrix<complex<double>, Prop0, ColHerm, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B,
	      const complex<double> beta,
	      const Matrix<complex<double>, Prop2, ColMajor, Allocator2>& C)
  {
    cblas_zhemm(CblasColMajor, Side, CblasUpper,
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetM(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetM());
  }


  /*** ColHerm and UpLo ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const enum CBLAS_SIDE Side,
	      const complex<float> alpha,
	      const enum CBLAS_UPLO Uplo,
	      const Matrix<complex<float>, Prop0, ColHerm, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B,
	      const complex<float> beta,
	      const Matrix<complex<float>, Prop2, ColMajor, Allocator2>& C)
  {
    cblas_chemm(CblasColMajor, Side, Uplo,
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetM(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const enum CBLAS_SIDE Side,
	      const complex<double> alpha,
	      const enum CBLAS_UPLO Uplo,
	      const Matrix<complex<double>, Prop0, ColHerm, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B,
	      const complex<double> beta,
	      const Matrix<complex<double>, Prop2, ColMajor, Allocator2>& C)
  {
    cblas_zhemm(CblasColMajor, Side, Uplo,
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetM(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetM());
  }


  /*** RowHerm and Upper ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const enum CBLAS_SIDE Side,
	      const complex<float> alpha,
	      const Matrix<complex<float>, Prop0, RowHerm, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B,
	      const complex<float> beta,
	      const Matrix<complex<float>, Prop2, RowMajor, Allocator2>& C)
  {
    cblas_chemm(CblasRowMajor, Side, CblasUpper,
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetN(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const enum CBLAS_SIDE Side,
	      const complex<double> alpha,
	      const Matrix<complex<double>, Prop0, RowHerm, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B,
	      const complex<double> beta,
	      const Matrix<complex<double>, Prop2, RowMajor, Allocator2>& C)
  {
    cblas_zhemm(CblasRowMajor, Side, CblasUpper,
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetN(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetN());
  }


  /*** RowHerm and UpLo ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const enum CBLAS_SIDE Side,
	      const complex<float> alpha,
	      const enum CBLAS_UPLO Uplo,
	      const Matrix<complex<float>, Prop0, RowHerm, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B,
	      const complex<float> beta,
	      const Matrix<complex<float>, Prop2, RowMajor, Allocator2>& C)
  {
    cblas_chemm(CblasRowMajor, Side, Uplo,
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetN(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const enum CBLAS_SIDE Side,
	      const complex<double> alpha,
	      const enum CBLAS_UPLO Uplo,
	      const Matrix<complex<double>, Prop0, RowHerm, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B,
	      const complex<double> beta,
	      const Matrix<complex<double>, Prop2, RowMajor, Allocator2>& C)
  {
    cblas_zhemm(CblasRowMajor, Side, Uplo,
		C.GetM(), C.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<const void*>(B.GetData()), B.GetN(),
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(C.GetData()), C.GetN());
  }


  // MltAdd //
  ////////////



  /////////
  // Mlt //


  /*** ColUpTriang, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const float alpha,
	   const Matrix<float, Prop0, ColUpTriang, Allocator0>& A,
	   Matrix<float, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_strmm(CblasColMajor, Side, CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const double alpha,
	   const Matrix<double, Prop0, ColUpTriang, Allocator0>& A,
	   Matrix<double, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_dtrmm(CblasColMajor, Side, CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const complex<float> alpha,
	   const Matrix<complex<float>, Prop0, ColUpTriang, Allocator0>& A,
	   Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_ctrmm(CblasColMajor, Side, CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const complex<double> alpha,
	   const Matrix<complex<double>, Prop0, ColUpTriang, Allocator0>& A,
	   Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_ztrmm(CblasColMajor, Side, CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  /*** ColUpTriang ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const float alpha,
	   const enum CBLAS_TRANSPOSE TransA,
	   const enum CBLAS_DIAG DiagA,
	   const Matrix<float, Prop0, ColUpTriang, Allocator0>& A,
	   Matrix<float, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_strmm(CblasColMajor, Side, CblasUpper, TransA, DiagA,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const double alpha,
	   const enum CBLAS_TRANSPOSE TransA,
	   const enum CBLAS_DIAG DiagA,
	   const Matrix<double, Prop0, ColUpTriang, Allocator0>& A,
	   Matrix<double, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_dtrmm(CblasColMajor, Side, CblasUpper, TransA, DiagA,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const complex<float> alpha,
	   const enum CBLAS_TRANSPOSE TransA,
	   const enum CBLAS_DIAG DiagA,
	   const Matrix<complex<float>, Prop0, ColUpTriang, Allocator0>& A,
	   Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_ctrmm(CblasColMajor, Side, CblasUpper, TransA, DiagA,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const complex<double> alpha,
	   const enum CBLAS_TRANSPOSE TransA,
	   const enum CBLAS_DIAG DiagA,
	   const Matrix<complex<double>, Prop0, ColUpTriang, Allocator0>& A,
	   Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_ztrmm(CblasColMajor, Side, CblasUpper, TransA, DiagA,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  /*** ColLoTriang, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const float alpha,
	   const Matrix<float, Prop0, ColLoTriang, Allocator0>& A,
	   Matrix<float, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_strmm(CblasColMajor, Side, CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const double alpha,
	   const Matrix<double, Prop0, ColLoTriang, Allocator0>& A,
	   Matrix<double, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_dtrmm(CblasColMajor, Side, CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const complex<float> alpha,
	   const Matrix<complex<float>, Prop0, ColLoTriang, Allocator0>& A,
	   Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_ctrmm(CblasColMajor, Side, CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const complex<double> alpha,
	   const Matrix<complex<double>, Prop0, ColLoTriang, Allocator0>& A,
	   Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_ztrmm(CblasColMajor, Side, CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  /*** ColLoTriang ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const float alpha,
	   const enum CBLAS_TRANSPOSE TransA,
	   const enum CBLAS_DIAG DiagA,
	   const Matrix<float, Prop0, ColLoTriang, Allocator0>& A,
	   Matrix<float, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_strmm(CblasColMajor, Side, CblasLower, TransA, DiagA,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const double alpha,
	   const enum CBLAS_TRANSPOSE TransA,
	   const enum CBLAS_DIAG DiagA,
	   const Matrix<double, Prop0, ColLoTriang, Allocator0>& A,
	   Matrix<double, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_dtrmm(CblasColMajor, Side, CblasLower, TransA, DiagA,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const complex<float> alpha,
	   const enum CBLAS_TRANSPOSE TransA,
	   const enum CBLAS_DIAG DiagA,
	   const Matrix<complex<float>, Prop0, ColLoTriang, Allocator0>& A,
	   Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_ctrmm(CblasColMajor, Side, CblasLower, TransA, DiagA,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const complex<double> alpha,
	   const enum CBLAS_TRANSPOSE TransA,
	   const enum CBLAS_DIAG DiagA,
	   const Matrix<complex<double>, Prop0, ColLoTriang, Allocator0>& A,
	   Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_ztrmm(CblasColMajor, Side, CblasLower, TransA, DiagA,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  /*** RowUpTriang, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const float alpha,
	   const Matrix<float, Prop0, RowUpTriang, Allocator0>& A,
	   Matrix<float, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_strmm(CblasRowMajor, Side, CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const double alpha,
	   const Matrix<double, Prop0, RowUpTriang, Allocator0>& A,
	   Matrix<double, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_dtrmm(CblasRowMajor, Side, CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const complex<float> alpha,
	   const Matrix<complex<float>, Prop0, RowUpTriang, Allocator0>& A,
	   Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_ctrmm(CblasRowMajor, Side, CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const complex<double> alpha,
	   const Matrix<complex<double>, Prop0, RowUpTriang, Allocator0>& A,
	   Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_ztrmm(CblasRowMajor, Side, CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  /*** RowUpTriang ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const float alpha,
	   const enum CBLAS_TRANSPOSE TransA,
	   const enum CBLAS_DIAG DiagA,
	   const Matrix<float, Prop0, RowUpTriang, Allocator0>& A,
	   Matrix<float, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_strmm(CblasRowMajor, Side, CblasUpper, TransA, DiagA,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const double alpha,
	   const enum CBLAS_TRANSPOSE TransA,
	   const enum CBLAS_DIAG DiagA,
	   const Matrix<double, Prop0, RowUpTriang, Allocator0>& A,
	   Matrix<double, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_dtrmm(CblasRowMajor, Side, CblasUpper, TransA, DiagA,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const complex<float> alpha,
	   const enum CBLAS_TRANSPOSE TransA,
	   const enum CBLAS_DIAG DiagA,
	   const Matrix<complex<float>, Prop0, RowUpTriang, Allocator0>& A,
	   Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_ctrmm(CblasRowMajor, Side, CblasUpper, TransA, DiagA,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const complex<double> alpha,
	   const enum CBLAS_TRANSPOSE TransA,
	   const enum CBLAS_DIAG DiagA,
	   const Matrix<complex<double>, Prop0, RowUpTriang, Allocator0>& A,
	   Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_ztrmm(CblasRowMajor, Side, CblasUpper, TransA, DiagA,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  /*** RowLoTriang, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const float alpha,
	   const Matrix<float, Prop0, RowLoTriang, Allocator0>& A,
	   Matrix<float, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_strmm(CblasRowMajor, Side, CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const double alpha,
	   const Matrix<double, Prop0, RowLoTriang, Allocator0>& A,
	   Matrix<double, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_dtrmm(CblasRowMajor, Side, CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const complex<float> alpha,
	   const Matrix<complex<float>, Prop0, RowLoTriang, Allocator0>& A,
	   Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_ctrmm(CblasRowMajor, Side, CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const complex<double> alpha,
	   const Matrix<complex<double>, Prop0, RowLoTriang, Allocator0>& A,
	   Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_ztrmm(CblasRowMajor, Side, CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  /*** RowLoTriang ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const float alpha,
	   const enum CBLAS_TRANSPOSE TransA,
	   const enum CBLAS_DIAG DiagA,
	   const Matrix<float, Prop0, RowLoTriang, Allocator0>& A,
	   Matrix<float, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_strmm(CblasRowMajor, Side, CblasLower, TransA, DiagA,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const double alpha,
	   const enum CBLAS_TRANSPOSE TransA,
	   const enum CBLAS_DIAG DiagA,
	   const Matrix<double, Prop0, RowLoTriang, Allocator0>& A,
	   Matrix<double, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_dtrmm(CblasRowMajor, Side, CblasLower, TransA, DiagA,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const complex<float> alpha,
	   const enum CBLAS_TRANSPOSE TransA,
	   const enum CBLAS_DIAG DiagA,
	   const Matrix<complex<float>, Prop0, RowLoTriang, Allocator0>& A,
	   Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_ctrmm(CblasRowMajor, Side, CblasLower, TransA, DiagA,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const enum CBLAS_SIDE Side,
	   const complex<double> alpha,
	   const enum CBLAS_TRANSPOSE TransA,
	   const enum CBLAS_DIAG DiagA,
	   const Matrix<complex<double>, Prop0, RowLoTriang, Allocator0>& A,
	   Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_ztrmm(CblasRowMajor, Side, CblasLower, TransA, DiagA,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  // Mlt //
  /////////



  ///////////
  // Solve //


  /*** ColUpTriang, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const float alpha,
	     const Matrix<float, Prop0, ColUpTriang, Allocator0>& A,
	     Matrix<float, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_strsm(CblasColMajor, Side, CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const double alpha,
	     const Matrix<double, Prop0, ColUpTriang, Allocator0>& A,
	     Matrix<double, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_dtrsm(CblasColMajor, Side, CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const complex<float> alpha,
	     const Matrix<complex<float>, Prop0, ColUpTriang, Allocator0>& A,
	     Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_ctrsm(CblasColMajor, Side, CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const complex<double> alpha,
	     const Matrix<complex<double>, Prop0, ColUpTriang, Allocator0>& A,
	     Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_ztrsm(CblasColMajor, Side, CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  /*** ColUpTriang ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const float alpha,
	     const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_DIAG DiagA,
	     const Matrix<float, Prop0, ColUpTriang, Allocator0>& A,
	     Matrix<float, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_strsm(CblasColMajor, Side, CblasUpper, TransA, DiagA,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const double alpha,
	     const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_DIAG DiagA,
	     const Matrix<double, Prop0, ColUpTriang, Allocator0>& A,
	     Matrix<double, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_dtrsm(CblasColMajor, Side, CblasUpper, TransA, DiagA,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const complex<float> alpha,
	     const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_DIAG DiagA,
	     const Matrix<complex<float>, Prop0, ColUpTriang, Allocator0>& A,
	     Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_ctrsm(CblasColMajor, Side, CblasUpper, TransA, DiagA,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const complex<double> alpha,
	     const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_DIAG DiagA,
	     const Matrix<complex<double>, Prop0, ColUpTriang, Allocator0>& A,
	     Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_ztrsm(CblasColMajor, Side, CblasUpper, TransA, DiagA,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  /*** ColLoTriang, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const float alpha,
	     const Matrix<float, Prop0, ColLoTriang, Allocator0>& A,
	     Matrix<float, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_strsm(CblasColMajor, Side, CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const double alpha,
	     const Matrix<double, Prop0, ColLoTriang, Allocator0>& A,
	     Matrix<double, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_dtrsm(CblasColMajor, Side, CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const complex<float> alpha,
	     const Matrix<complex<float>, Prop0, ColLoTriang, Allocator0>& A,
	     Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_ctrsm(CblasColMajor, Side, CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const complex<double> alpha,
	     const Matrix<complex<double>, Prop0, ColLoTriang, Allocator0>& A,
	     Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_ztrsm(CblasColMajor, Side, CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  /*** ColLoTriang ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const float alpha,
	     const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_DIAG DiagA,
	     const Matrix<float, Prop0, ColLoTriang, Allocator0>& A,
	     Matrix<float, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_strsm(CblasColMajor, Side, CblasLower, TransA, DiagA,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const double alpha,
	     const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_DIAG DiagA,
	     const Matrix<double, Prop0, ColLoTriang, Allocator0>& A,
	     Matrix<double, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_dtrsm(CblasColMajor, Side, CblasLower, TransA, DiagA,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetM(), B.GetData(), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const complex<float> alpha,
	     const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_DIAG DiagA,
	     const Matrix<complex<float>, Prop0, ColLoTriang, Allocator0>& A,
	     Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_ctrsm(CblasColMajor, Side, CblasLower, TransA, DiagA,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const complex<double> alpha,
	     const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_DIAG DiagA,
	     const Matrix<complex<double>, Prop0, ColLoTriang, Allocator0>& A,
	     Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B)
  {
    cblas_ztrsm(CblasColMajor, Side, CblasLower, TransA, DiagA,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetM(),
		reinterpret_cast<void*>(B.GetData()), B.GetM());
  }


  /*** RowUpTriang, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const float alpha,
	     const Matrix<float, Prop0, RowUpTriang, Allocator0>& A,
	     Matrix<float, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_strsm(CblasRowMajor, Side, CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const double alpha,
	     const Matrix<double, Prop0, RowUpTriang, Allocator0>& A,
	     Matrix<double, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_dtrsm(CblasRowMajor, Side, CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const complex<float> alpha,
	     const Matrix<complex<float>, Prop0, RowUpTriang, Allocator0>& A,
	     Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_ctrsm(CblasRowMajor, Side, CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const complex<double> alpha,
	     const Matrix<complex<double>, Prop0, RowUpTriang, Allocator0>& A,
	     Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_ztrsm(CblasRowMajor, Side, CblasUpper, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  /*** RowUpTriang ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const float alpha,
	     const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_DIAG DiagA,
	     const Matrix<float, Prop0, RowUpTriang, Allocator0>& A,
	     Matrix<float, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_strsm(CblasRowMajor, Side, CblasUpper, TransA, DiagA,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const double alpha,
	     const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_DIAG DiagA,
	     const Matrix<double, Prop0, RowUpTriang, Allocator0>& A,
	     Matrix<double, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_dtrsm(CblasRowMajor, Side, CblasUpper, TransA, DiagA,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const complex<float> alpha,
	     const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_DIAG DiagA,
	     const Matrix<complex<float>, Prop0, RowUpTriang, Allocator0>& A,
	     Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_ctrsm(CblasRowMajor, Side, CblasUpper, TransA, DiagA,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const complex<double> alpha,
	     const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_DIAG DiagA,
	     const Matrix<complex<double>, Prop0, RowUpTriang, Allocator0>& A,
	     Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_ztrsm(CblasRowMajor, Side, CblasUpper, TransA, DiagA,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  /*** RowLoTriang, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const float alpha,
	     const Matrix<float, Prop0, RowLoTriang, Allocator0>& A,
	     Matrix<float, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_strsm(CblasRowMajor, Side, CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const double alpha,
	     const Matrix<double, Prop0, RowLoTriang, Allocator0>& A,
	     Matrix<double, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_dtrsm(CblasRowMajor, Side, CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const complex<float> alpha,
	     const Matrix<complex<float>, Prop0, RowLoTriang, Allocator0>& A,
	     Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_ctrsm(CblasRowMajor, Side, CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const complex<double> alpha,
	     const Matrix<complex<double>, Prop0, RowLoTriang, Allocator0>& A,
	     Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_ztrsm(CblasRowMajor, Side, CblasLower, CblasNoTrans, CblasNonUnit,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  /*** RowLoTriang ***/


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const float alpha,
	     const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_DIAG DiagA,
	     const Matrix<float, Prop0, RowLoTriang, Allocator0>& A,
	     Matrix<float, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_strsm(CblasRowMajor, Side, CblasLower, TransA, DiagA,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const double alpha,
	     const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_DIAG DiagA,
	     const Matrix<double, Prop0, RowLoTriang, Allocator0>& A,
	     Matrix<double, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_dtrsm(CblasRowMajor, Side, CblasLower, TransA, DiagA,
		B.GetM(), B.GetN(),
		alpha, A.GetData(), A.GetN(), B.GetData(), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const complex<float> alpha,
	     const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_DIAG DiagA,
	     const Matrix<complex<float>, Prop0, RowLoTriang, Allocator0>& A,
	     Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_ctrsm(CblasRowMajor, Side, CblasLower, TransA, DiagA,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const enum CBLAS_SIDE Side,
	     const complex<double> alpha,
	     const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_DIAG DiagA,
	     const Matrix<complex<double>, Prop0, RowLoTriang, Allocator0>& A,
	     Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B)
  {
    cblas_ztrsm(CblasRowMajor, Side, CblasLower, TransA, DiagA,
		B.GetM(), B.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetN(),
		reinterpret_cast<void*>(B.GetData()), B.GetN());
  }


  // Solve //
  ///////////


} // namespace Seldon.

#define SELDON_FILE_BLAS_3_CXX
#endif
