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

#ifndef SELDON_FILE_LAPACK_CXX

namespace Seldon
{


  ////////////////
  // LAPACKINFO //

  
  class LapackInfo
  {
  private:
    int info_;
  public:
    LapackInfo(int info): info_(info)
    {  
    }
    operator int ()
    {
      return info_;
    }
    int GetInfo()
    {
      return info_;
    }
    int& GetInfoRef()
    {
      return info_;
    }
  } lapack_info(0);


  // LAPACKINFO //
  ////////////////


  ///////////
  // GETLU //


  /*** ColMajor ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<float, Prop0, ColMajor, Allocator0>& A,
	     Vector<int, Vect_Full, Allocator1>& P,
	     LapackInfo& info = lapack_info)
  {
    int m(A.GetM());
    int n(A.GetN());
    P.Reallocate(min(m, n));
    sgetrf_(&m, &n, A.GetData(), &m, P.GetData(), &lapack_info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<double, Prop0, ColMajor, Allocator0>& A,
	     Vector<int, Vect_Full, Allocator1>& P,
	     LapackInfo& info = lapack_info)
  {
    int m(A.GetM());
    int n(A.GetN());
    P.Reallocate(min(m, n));
    dgetrf_(&m, &n, A.GetData(), &m, P.GetData(), &lapack_info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
	     Vector<int, Vect_Full, Allocator1>& P,
	     LapackInfo& info = lapack_info)
  {
    int m(A.GetM());
    int n(A.GetN());
    P.Reallocate(min(m, n));
    cgetrf_(&m, &n, A.GetDataVoid(), &m,
	    P.GetData(), &lapack_info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
	     Vector<int, Vect_Full, Allocator1>& P,
	     LapackInfo& info = lapack_info)
  {
    int m(A.GetM());
    int n(A.GetN());
    P.Reallocate(min(m, n));
    zgetrf_(&m, &n, A.GetDataVoid(), &m,
	    P.GetData(), &lapack_info.GetInfoRef());
  }


  /*** RowMajor ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<float, Prop0, RowMajor, Allocator0>& A,
	     Vector<int, Vect_Full, Allocator1>& P,
	     LapackInfo& info = lapack_info)
  {
    int m(A.GetM());
    int n(A.GetN());
    P.Reallocate(min(m, n));
    sgetrf_(&m, &n, A.GetData(), &m, P.GetData(), &lapack_info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<double, Prop0, RowMajor, Allocator0>& A,
	     Vector<int, Vect_Full, Allocator1>& P,
	     LapackInfo& info = lapack_info)
  {
    int m(A.GetM());
    int n(A.GetN());
    P.Reallocate(min(m, n));
    dgetrf_(&m, &n, A.GetData(), &m, P.GetData(), &lapack_info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
	     Vector<int, Vect_Full, Allocator1>& P,
	     LapackInfo& info = lapack_info)
  {
    int m(A.GetM());
    int n(A.GetN());
    P.Reallocate(min(m, n));
    cgetrf_(&m, &n, A.GetDataVoid(), &m,
	    P.GetData(), &lapack_info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
	     Vector<int, Vect_Full, Allocator1>& P,
	     LapackInfo& info = lapack_info)
  {
    int m(A.GetM());
    int n(A.GetN());
    P.Reallocate(min(m, n));
    zgetrf_(&m, &n, A.GetDataVoid(), &m,
	    P.GetData(), &lapack_info.GetInfoRef());
  }


  /*** ColSymPacked and Upper ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<float, Prop0, ColSymPacked, Allocator0>& A,
	     Vector<int, Vect_Full, Allocator1>& P,
	     LapackInfo& info = lapack_info)
  {
    int m(A.GetM());
    char uplo('U');
    P.Reallocate(m);
    ssptrf_(&uplo, &m, A.GetData(), P.GetData(), &lapack_info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<double, Prop0, ColSymPacked, Allocator0>& A,
	     Vector<int, Vect_Full, Allocator1>& P,
	     LapackInfo& info = lapack_info)
  {
    int m(A.GetM());
    char uplo('U');
    P.Reallocate(m);
    dsptrf_(&uplo, &m, A.GetData(), P.GetData(), &lapack_info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<float>, Prop0, ColSymPacked, Allocator0>& A,
	     Vector<int, Vect_Full, Allocator1>& P,
	     LapackInfo& info = lapack_info)
  {
    int m(A.GetM());
    char uplo('U');
    P.Reallocate(m);
    csptrf_(&uplo, &m, A.GetDataVoid(),
	    P.GetData(), &lapack_info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<double>, Prop0, ColSymPacked, Allocator0>& A,
	     Vector<int, Vect_Full, Allocator1>& P,
	     LapackInfo& info = lapack_info)
  {
    int m(A.GetM());
    char uplo('U');
    P.Reallocate(m);
    zsptrf_(&uplo, &m, A.GetDataVoid(),
	    P.GetData(), &lapack_info.GetInfoRef());
  }


  /*** ColSymPacked and Uplo ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(SeldonUplo Uplo,
	     Matrix<float, Prop0, ColSymPacked, Allocator0>& A,
	     Vector<int, Vect_Full, Allocator1>& P,
	     LapackInfo& info = lapack_info)
  {
    int m(A.GetM());
    char uplo(Uplo);
    P.Reallocate(m);
    ssptrf_(&uplo, &m, A.GetData(), P.GetData(), &lapack_info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(SeldonUplo Uplo,
	     Matrix<double, Prop0, ColSymPacked, Allocator0>& A,
	     Vector<int, Vect_Full, Allocator1>& P,
	     LapackInfo& info = lapack_info)
  {
    int m(A.GetM());
    char uplo(Uplo);
    P.Reallocate(m);
    dsptrf_(&uplo, &m, A.GetData(), P.GetData(), &lapack_info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(SeldonUplo Uplo,
	     Matrix<complex<float>, Prop0, ColSymPacked, Allocator0>& A,
	     Vector<int, Vect_Full, Allocator1>& P,
	     LapackInfo& info = lapack_info)
  {
    int m(A.GetM());
    char uplo(Uplo);
    P.Reallocate(m);
    csptrf_(&uplo, &m, A.GetDataVoid(),
	    P.GetData(), &lapack_info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(SeldonUplo Uplo,
	     Matrix<complex<double>, Prop0, ColSymPacked, Allocator0>& A,
	     Vector<int, Vect_Full, Allocator1>& P,
	     LapackInfo& info = lapack_info)
  {
    int m(A.GetM());
    char uplo(Uplo);
    P.Reallocate(m);
    zsptrf_(&uplo, &m, A.GetDataVoid(),
	    P.GetData(), &lapack_info.GetInfoRef());
  }


  /*** RowSymPacked and Upper ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<float, Prop0, RowSymPacked, Allocator0>& A,
	     Vector<int, Vect_Full, Allocator1>& P,
	     LapackInfo& info = lapack_info)
  {
    int m(A.GetM());
    char uplo('L');
    P.Reallocate(m);
    ssptrf_(&uplo, &m, A.GetData(), P.GetData(), &lapack_info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<double, Prop0, RowSymPacked, Allocator0>& A,
	     Vector<int, Vect_Full, Allocator1>& P,
	     LapackInfo& info = lapack_info)
  {
    int m(A.GetM());
    char uplo('L');
    P.Reallocate(m);
    dsptrf_(&uplo, &m, A.GetData(), P.GetData(), &lapack_info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<float>, Prop0, RowSymPacked, Allocator0>& A,
	     Vector<int, Vect_Full, Allocator1>& P,
	     LapackInfo& info = lapack_info)
  {
    int m(A.GetM());
    char uplo('L');
    P.Reallocate(m);
    csptrf_(&uplo, &m, A.GetDataVoid(),
	    P.GetData(), &lapack_info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<double>, Prop0, RowSymPacked, Allocator0>& A,
	     Vector<int, Vect_Full, Allocator1>& P,
	     LapackInfo& info = lapack_info)
  {
    int m(A.GetM());
    char uplo('L');
    P.Reallocate(m);
    zsptrf_(&uplo, &m, A.GetDataVoid(),
	    P.GetData(), &lapack_info.GetInfoRef());
  }


  /*** RowSymPacked and Uplo ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(SeldonUplo Uplo,
	     Matrix<float, Prop0, RowSymPacked, Allocator0>& A,
	     Vector<int, Vect_Full, Allocator1>& P,
	     LapackInfo& info = lapack_info)
  {
    int m(A.GetM());
    char uplo(Uplo.RevChar());
    P.Reallocate(m);
    ssptrf_(&uplo, &m, A.GetData(), P.GetData(), &lapack_info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(SeldonUplo Uplo,
	     Matrix<double, Prop0, RowSymPacked, Allocator0>& A,
	     Vector<int, Vect_Full, Allocator1>& P,
	     LapackInfo& info = lapack_info)
  {
    int m(A.GetM());
    char uplo(Uplo.RevChar());
    P.Reallocate(m);
    dsptrf_(&uplo, &m, A.GetData(), P.GetData(), &lapack_info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(SeldonUplo Uplo,
	     Matrix<complex<float>, Prop0, RowSymPacked, Allocator0>& A,
	     Vector<int, Vect_Full, Allocator1>& P,
	     LapackInfo& info = lapack_info)
  {
    int m(A.GetM());
    char uplo(Uplo.RevChar());
    P.Reallocate(m);
    csptrf_(&uplo, &m, A.GetDataVoid(),
	    P.GetData(), &lapack_info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(SeldonUplo Uplo,
	     Matrix<complex<double>, Prop0, RowSymPacked, Allocator0>& A,
	     Vector<int, Vect_Full, Allocator1>& P,
	     LapackInfo& info = lapack_info)
  {
    int m(A.GetM());
    char uplo(Uplo.RevChar());
    P.Reallocate(m);
    zsptrf_(&uplo, &m, A.GetDataVoid(),
	    P.GetData(), &lapack_info.GetInfoRef());
  }


  // GETLU //
  ///////////


} // namespace Seldon.

#define SELDON_FILE_LAPACK_CXX
#endif
