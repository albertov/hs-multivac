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

#ifndef SELDON_FILE_BLAS_1_CXX

/*
  Functions defined in this file:

  xROTG
*/

extern "C"
{
#include "cblas.h"
}

namespace Seldon
{


  ////////////
  // GenRot //


  void GenRot(float& a, float& b, float& c, float& d)
  {
    cblas_srotg(&a, &b, &c, &d);
  }


  void GenRot(double& a, double& b, double& c, double& d)
  {
    cblas_drotg(&a, &b, &c, &d);
  }


  // GenRot //
  ////////////



  /////////////////
  // GenModifRot //


  void GenModifRot(float& d1, float& d2,
		   float& x1, const float& y1,
		   float* param)
  {
    cblas_srotmg(&d1, &d2, &x1, y1, param);
  }


  void GenModifRot(double& d1, double& d2,
		   double& x1, const double& y1,
		   double* param)
  {
    cblas_drotmg(&d1, &d2, &x1, y1, param);
  }


  // GenModifRot //
  /////////////////



  //////////////
  // ApplyRot //


  template <class Allocator>
  void ApplyRot(Vector<float, Vect_Full, Allocator>& X,
		Vector<float, Vect_Full, Allocator>& Y,
		const float c, const float s)
  {
    cblas_srot(X.GetLength(), X.GetData(), 1,
	       Y.GetData(), 1, c, s);
  }


  template <class Allocator>
  void ApplyRot(Vector<double, Vect_Full, Allocator>& X,
		Vector<double, Vect_Full, Allocator>& Y,
		const double c, const double s)
  {
    cblas_drot(X.GetLength(), X.GetData(), 1,
	       Y.GetData(), 1, c, s);
  }


  // ApplyRot //
  //////////////



  ///////////////////
  // ApplyModifRot //


  template <class Allocator>
  void ApplyModifRot(Vector<float, Vect_Full, Allocator>& X,
		     Vector<float, Vect_Full, Allocator>& Y,
		     const float* param)
  {
    cblas_srotm(X.GetLength(), X.GetData(), 1,
		Y.GetData(), 1, param);
  }


  template <class Allocator>
  void ApplyModifRot(Vector<double, Vect_Full, Allocator>& X,
		     Vector<double, Vect_Full, Allocator>& Y,
		     const double* param)
  {
    cblas_drotm(X.GetLength(), X.GetData(), 1,
		Y.GetData(), 1, param);
  }


  // ApplyModifRot //
  ///////////////////



  //////////
  // Swap //


  template <class Allocator>
  void Swap(Vector<float, Vect_Full, Allocator>& X,
	    Vector<float, Vect_Full, Allocator>& Y)
  {
    cblas_sswap(X.GetLength(), X.GetData(), 1,
		Y.GetData(), 1);
  }


  template <class Allocator>
  void Swap(Vector<double, Vect_Full, Allocator>& X,
	    Vector<double, Vect_Full, Allocator>& Y)
  {
    cblas_dswap(X.GetLength(), X.GetData(), 1,
		Y.GetData(), 1);
  }


  template <class Allocator>
  void Swap(Vector<complex<float>, Vect_Full, Allocator>& X,
	    Vector<complex<float>, Vect_Full, Allocator>& Y)
  {
    cblas_cswap(X.GetLength(), reinterpret_cast<void*>(X.GetData()), 1,
		reinterpret_cast<void*>(Y.GetData()), 1);
  }


  template <class Allocator>
  void Swap(Vector<complex<double>, Vect_Full, Allocator>& X,
	    Vector<complex<double>, Vect_Full, Allocator>& Y)
  {
    cblas_zswap(X.GetLength(), reinterpret_cast<void*>(X.GetData()), 1,
		reinterpret_cast<void*>(Y.GetData()), 1);
  }


  // Swap //
  //////////



  /////////
  // Mlt //


  template <class Allocator>
  void Mlt(const float alpha,
	   Vector<float, Vect_Full, Allocator>& X)
  {
    cblas_sscal(X.GetLength(), alpha, X.GetData(), 1);
  }


  template <class Allocator>
  void Mlt(const double alpha,
	   Vector<double, Vect_Full, Allocator>& X)
  {
    cblas_dscal(X.GetLength(), alpha, X.GetData(), 1);
  }


  template <class Allocator>
  void Mlt(const float alpha,
	   Vector<complex<float>, Vect_Full, Allocator>& X)
  {
    cblas_csscal(X.GetLength(), alpha,
		 reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Allocator>
  void Mlt(const double alpha,
	   Vector<complex<double>, Vect_Full, Allocator>& X)
  {
    cblas_zdscal(X.GetLength(), alpha,
		 reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Allocator>
  void Mlt(const complex<float> alpha,
	   Vector<complex<float>, Vect_Full, Allocator>& X)
  {
    cblas_cscal(X.GetLength(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Allocator>
  void Mlt(const complex<double> alpha,
	   Vector<complex<double>, Vect_Full, Allocator>& X)
  {
    cblas_zscal(X.GetLength(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  // Mlt //
  /////////



  //////////
  // Copy //


  template <class Allocator0, class Allocator1>
  void Copy(const Vector<float, Vect_Full, Allocator0>& X,
	    Vector<float, Vect_Full, Allocator1>& Y)
  {
    cblas_scopy(Y.GetLength(),
		reinterpret_cast<const float*>(X.GetData()), 1,
		reinterpret_cast<float*>(Y.GetData()), 1);
  }


  template <class Allocator0, class Allocator1>
  void Copy(const Vector<double, Vect_Full, Allocator0>& X,
	    Vector<double, Vect_Full, Allocator1>& Y)
  {
    cblas_dcopy(Y.GetLength(),
		reinterpret_cast<const double*>(X.GetData()), 1,
		reinterpret_cast<double*>(Y.GetData()), 1);
  }


  template <class Allocator0, class Allocator1>
  void Copy(const Vector<complex<float>, Vect_Full, Allocator0>& X,
	    Vector<complex<float>, Vect_Full, Allocator1>& Y)
  {
    cblas_ccopy(Y.GetLength(),
		reinterpret_cast<const void*>(X.GetData()), 1,
		reinterpret_cast<void*>(Y.GetData()), 1);
  }


  template <class Allocator0, class Allocator1>
  void Copy(const Vector<complex<double>, Vect_Full, Allocator0>& X,
	    Vector<complex<double>, Vect_Full, Allocator1>& Y)
  {
    cblas_zcopy(Y.GetLength(),
		reinterpret_cast<const void*>(X.GetData()), 1,
		reinterpret_cast<void*>(Y.GetData()), 1);
  }


  // Copy //
  //////////



  /////////
  // Add //


  template <class Allocator0, class Allocator1>
  void Add(const float alpha,
	   const Vector<float, Vect_Full, Allocator0>& A,
	   Vector<float, Vect_Full, Allocator1>& C)
  {
    cblas_saxpy(C.GetLength(),
		alpha,
		reinterpret_cast<const float*>(A.GetData()), 1,
		reinterpret_cast<float*>(C.GetData()), 1);
  }


  template <class Allocator0, class Allocator1>
  void Add(const double alpha,
	   const Vector<double, Vect_Full, Allocator0>& A,
	   Vector<double, Vect_Full, Allocator1>& C)
  {
    cblas_daxpy(C.GetLength(),
		alpha,
		reinterpret_cast<const double*>(A.GetData()), 1,
		reinterpret_cast<double*>(C.GetData()), 1);
  }


  template <class Allocator0, class Allocator1>
  void Add(const complex<float> alpha,
	   const Vector<complex<float>, Vect_Full, Allocator0>& A,
	   Vector<complex<float>, Vect_Full, Allocator1>& C)
  {
    cblas_caxpy(C.GetLength(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), 1,
		reinterpret_cast<float*>(C.GetData()), 1);
  }


  template <class Allocator0, class Allocator1>
  void Add(const complex<double> alpha,
	   const Vector<complex<double>, Vect_Full, Allocator0>& A,
	   Vector<complex<double>, Vect_Full, Allocator1>& C)
  {
    cblas_zaxpy(C.GetLength(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), 1,
		reinterpret_cast<double*>(C.GetData()), 1);
  }


  // Add //
  /////////



  /////////////
  // DotProd //


  template <class Allocator0, class Allocator1>
  float DotProd(const Vector<float, Vect_Full, Allocator0>& X,
		const Vector<float, Vect_Full, Allocator1>& Y)
  {
    return cblas_sdot(Y.GetLength(),
		      reinterpret_cast<const float*>(X.GetData()), 1,
		      reinterpret_cast<const float*>(Y.GetData()), 1);
  }


  template <class Allocator0, class Allocator1>
  double DotProd(const Vector<double, Vect_Full, Allocator0>& X,
		 const Vector<double, Vect_Full, Allocator1>& Y)
  {
    return cblas_ddot(Y.GetLength(),
		      reinterpret_cast<const double*>(X.GetData()), 1,
		      reinterpret_cast<const double*>(Y.GetData()), 1);
  }


  template <class Allocator0, class Allocator1>
  complex<float>
  DotProd(const Vector<complex<float>, Vect_Full, Allocator0>& X,
	  const Vector<complex<float>, Vect_Full, Allocator1>& Y)
  {
    complex<float> dotu;
    cblas_cdotu_sub(Y.GetLength(),
		    reinterpret_cast<const void*>(X.GetData()), 1,
		    reinterpret_cast<const void*>(Y.GetData()), 1,
		    reinterpret_cast<void*>(&dotu));
    return dotu;
  }


  template <class Allocator0, class Allocator1>
  complex<double>
  DotProd(const Vector<complex<double>, Vect_Full, Allocator0>& X,
	  const Vector<complex<double>, Vect_Full, Allocator1>& Y)
  {
    complex<double> dotu;
    cblas_zdotu_sub(Y.GetLength(),
		    reinterpret_cast<const void*>(X.GetData()), 1,
		    reinterpret_cast<const void*>(Y.GetData()), 1,
		    reinterpret_cast<void*>(&dotu));
    return dotu;
  }


  // DotProd //
  /////////////



  ///////////////////
  // DotProdDouble //


  template <class Allocator0, class Allocator1>
  float DotProdDouble(const float alpha,
		      const Vector<float, Vect_Full, Allocator0>& X,
		      const Vector<float, Vect_Full, Allocator1>& Y)
  {
    return cblas_sdsdot(Y.GetLength(), alpha,
			reinterpret_cast<const float*>(X.GetData()), 1,
			reinterpret_cast<const float*>(Y.GetData()), 1);
  }


  // DotProdDouble //
  ///////////////////



  /////////////////
  // DotProjConj //


  template <class Allocator0, class Allocator1>
  complex<float>
  DotProdConj(const Vector<complex<float>, Vect_Full, Allocator0>& X,
	      const Vector<complex<float>, Vect_Full, Allocator1>& Y)
  {
    complex<float> dotc;
    cblas_cdotc_sub(Y.GetLength(),
		    reinterpret_cast<const void*>(X.GetData()), 1,
		    reinterpret_cast<const void*>(Y.GetData()), 1,
		    reinterpret_cast<void*>(&dotc));
    return dotc;
  }


  template <class Allocator0, class Allocator1>
  complex<double>
  DotProdConj(const Vector<complex<double>, Vect_Full, Allocator0>& X,
	      const Vector<complex<double>, Vect_Full, Allocator1>& Y)
  {
    complex<double> dotc;
    cblas_zdotc_sub(Y.GetLength(),
		    reinterpret_cast<const void*>(X.GetData()), 1,
		    reinterpret_cast<const void*>(Y.GetData()), 1,
		    reinterpret_cast<void*>(&dotc));
    return dotc;
  }


  // DotProdConj //
  /////////////////



  ///////////
  // Norm1 //


  template <class Allocator>
  float Norm1(const Vector<float, Vect_Full, Allocator>& X)
  {
    return cblas_sasum(X.GetLength(),
		       reinterpret_cast<const float*>(X.GetData()), 1);
  }


  template <class Allocator>
  double Norm1(const Vector<double, Vect_Full, Allocator>& X)
  {
    return cblas_dasum(X.GetLength(),
		       reinterpret_cast<const double*>(X.GetData()), 1);
  }


  template <class Allocator>
  float Norm1(const Vector<complex<float>, Vect_Full, Allocator>& X)
  {
    return cblas_scasum(X.GetLength(),
			reinterpret_cast<const void*>(X.GetData()), 1);
  }


  template <class Allocator>
  double Norm1(const Vector<complex<double>, Vect_Full, Allocator>& X)
  {
    return cblas_dzasum(X.GetLength(),
			reinterpret_cast<const void*>(X.GetData()), 1);
  }


  // Norm1 //
  ///////////



  ///////////
  // Norm2 //


  template <class Allocator>
  float Norm2(const Vector<float, Vect_Full, Allocator>& X)
  {
    return cblas_snrm2(X.GetLength(),
		       reinterpret_cast<const float*>(X.GetData()), 1);
  }


  template <class Allocator>
  double Norm2(const Vector<double, Vect_Full, Allocator>& X)
  {
    return cblas_dnrm2(X.GetLength(),
		       reinterpret_cast<const double*>(X.GetData()), 1);
  }


  template <class Allocator>
  float Norm2(const Vector<complex<float>, Vect_Full, Allocator>& X)
  {
    return cblas_scnrm2(X.GetLength(),
			reinterpret_cast<const void*>(X.GetData()), 1);
  }


  template <class Allocator>
  double Norm2(const Vector<complex<double>, Vect_Full, Allocator>& X)
  {
    return cblas_dznrm2(X.GetLength(),
			reinterpret_cast<const void*>(X.GetData()), 1);
  }


  // Norm2 //
  ///////////



  ////////////////////
  // GetMaxAbsIndex //


  template <class Allocator>
  size_t GetMaxAbsIndex(const Vector<float, Vect_Full, Allocator>& X)
  {
    return cblas_isamax(X.GetLength(),
			reinterpret_cast<const float*>(X.GetData()), 1);
  }


  template <class Allocator>
  size_t GetMaxAbsIndex(const Vector<double, Vect_Full, Allocator>& X)
  {
    return cblas_idamax(X.GetLength(),
			reinterpret_cast<const double*>(X.GetData()), 1);
  }
  
  
  template <class Allocator>
  size_t GetMaxAbsIndex(const Vector<complex<float>, Vect_Full, Allocator>& X)
  {
    return cblas_icamax(X.GetLength(),
			reinterpret_cast<const void*>(X.GetData()), 1);
  }


  template <class Allocator>
  size_t
  GetMaxAbsIndex(const Vector<complex<double>, Vect_Full, Allocator>& X)
  {
    return cblas_izamax(X.GetLength(),
			reinterpret_cast<const void*>(X.GetData()), 1);
  }


  // GetMaxAbsIndex //
  ////////////////////


} // namespace Seldon.

#define SELDON_FILE_BLAS_1_CXX
#endif
