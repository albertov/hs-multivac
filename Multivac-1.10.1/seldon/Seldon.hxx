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


#ifndef SELDON_FILE_SELDON_HXX

#include <iostream>
#include <algorithm>
#include <complex>
#include <string>
#include <sstream>
#include <fstream>
#include <limits>
#include <cstdlib>
#include <ctime>
#include <exception>
#include <stdexcept>

#ifdef SELDON_WITH_CBLAS
extern "C"
{
#include "cblas.h"
}
#endif


//////////////////
// DEBUG LEVELS //
//////////////////

#ifdef SELDON_DEBUG_LEVEL_4
#ifndef SELDON_DEBUG_LEVEL_3
#define SELDON_DEBUG_LEVEL_3
#endif
#endif

#ifdef SELDON_DEBUG_LEVEL_3
#ifndef SELDON_CHECK_BOUNDARIES
#define SELDON_CHECK_BOUNDARIES
#endif
#ifndef SELDON_DEBUG_LEVEL_2
#define SELDON_DEBUG_LEVEL_2
#endif
#endif

#ifdef SELDON_DEBUG_LEVEL_2
#ifndef SELDON_CHECK_DIMENSIONS
#define SELDON_CHECK_DIMENSIONS
#endif
#ifndef SELDON_DEBUG_LEVEL_1
#define SELDON_DEBUG_LEVEL_1
#endif
#endif

#ifdef SELDON_DEBUG_LEVEL_1
#ifndef SELDON_CHECK_MEMORY
#define SELDON_CHECK_MEMORY
#endif
#ifndef SELDON_CHECK_IO
#define SELDON_CHECK_IO
#endif
#ifndef SELDON_DEBUG_LEVEL_0
#define SELDON_DEBUG_LEVEL_0
#endif
#endif

#ifdef SELDON_DEBUG_LEVEL_0
#ifndef SELDON_DEBUG_LEVEL_1
#ifndef SELDON_WITHOUT_THROW
#define SELDON_WITHOUT_THROW
#endif
#endif
#endif

// Convenient macros to catch exceptions.
#ifndef TRY
#define TRY try {
#endif
#ifndef END
#define END \
}\
catch(Seldon::Error& Err)\
{\
Err.What();\
return 1;\
}\
catch (std::exception& Err)\
{\
cerr << "C++ exception: " << Err.what() << endl;\
return 1;\
}\
catch (std::string& str)\
{\
cerr << str << endl;\
return 1;\
}\
catch (const char* str)\
{\
cerr << str << endl;\
return 1;\
}\
catch(...)\
{\
cerr << "Unknown exception..." << endl;\
return 1;\
}
#endif

//! To display a message... call Hermes!
#ifndef ERR
#define ERR(x) cerr << "Hermes - " #x << endl
#endif
//! To display a variable (with its name); same as DISPLAY.
#ifndef DISP
#define DISP(x) cerr << #x ": " << x << endl
#endif
//! To display a variable (with its name); same as DISP.
#ifndef DISPLAY
#define DISPLAY(x) cerr << #x ": " << x << endl
#endif

//! Seldon namespace.
namespace Seldon
{
  using namespace std;
}

// Useful functions.
#include "Common/Common.hxx"

// Default allocator.
#ifndef SELDON_DEFAULT_ALLOCATOR
#define SELDON_DEFAULT_ALLOCATOR MallocAlloc
#endif
// Memory management.
#include "Common/Allocator.hxx"

// Storage type.
#include "Common/Storage.hxx"

// Properties.
#include "Common/Properties.hxx"

namespace Seldon
{
  

  // Base structure for all vectors.
  template <class T, class Allocator>
  class Vector_Base;

  // Vector class - specialized for each used type.
  template <class T, class Storage = Vect_Full,
	    class Allocator = SELDON_DEFAULT_ALLOCATOR<T> >
  class Vector
  {
    // Nothing in it: no default vector is supplied so as to avoid suprises!
  };

  // Full vector.
  template <class T, class Allocator>
  class Vector<T, Vect_Full, Allocator>;

  // Matrix class - specialized for each used type.
  template <class T, class Prop = General, class Storage = RowMajor,
	    class Allocator = SELDON_DEFAULT_ALLOCATOR<T> >
  class Matrix
  {
    // Nothing in it: no default matrix is supplied so as to avoid suprises!
  };

  // column-major matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColMajor, Allocator>;

  // row-major matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowMajor, Allocator>;

  // column-major symmetric packed matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColSymPacked, Allocator>;

  // row-major symmetric packed matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowSymPacked, Allocator>;

  // column-major upper-triangular packed matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColUpTriangPacked, Allocator>;

  // column-major lower-triangular packed matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColLoTriangPacked, Allocator>;

  // row-major upper-triangular packed matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowUpTriangPacked, Allocator>;

  // row-major lower-triangular packed matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowLoTriangPacked, Allocator>;

  // column-major sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColSparse, Allocator>;

  // row-major sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowSparse, Allocator>;

  // column-major symmetric sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColSymSparse, Allocator>;

  // row-major symmetric sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowSymSparse, Allocator>;

  // column-major complex sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColComplexSparse, Allocator>;

  // row-major complex sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowComplexSparse, Allocator>;

  // column-major symmetric complex sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColSymComplexSparse, Allocator>;

  // row-major symmetric complex sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowSymComplexSparse, Allocator>;

  // 3D array.
  template <class T, class Allocator>
  class Array3D;

  //

  class SeldonTranspose
  {
#ifdef SELDON_WITH_CBLAS
  protected:
    CBLAS_TRANSPOSE cblas_status_;
#endif
  protected:
    // 0: NoTrans, 1: Trans, 2: ConjTrans.
    int status_;
  public:
    SeldonTranspose(int status)
    {
      status_ = status;
#ifdef SELDON_WITH_CBLAS
      if (status_ == 0)
	cblas_status_ = CblasTrans;
      else if (status_ == 1)
	cblas_status_ = CblasNoTrans;
      else
	cblas_status_ = CblasConjTrans;
#endif
    }
#ifdef SELDON_WITH_CBLAS
    operator CBLAS_TRANSPOSE() const
    {
      return cblas_status_;
    }
#endif
    
    bool Trans() {return (status_ == 0);}
    bool NoTrans() {return (status_ == 1);}
    bool ConjTrans() {return (status_ == 2);}
  };

  class class_SeldonTrans: public SeldonTranspose
  {
  public:
    class_SeldonTrans(): SeldonTranspose(0) {};
  };

  class class_SeldonNoTrans: public SeldonTranspose
  {
  public:
    class_SeldonNoTrans(): SeldonTranspose(1) {};
  };

  class class_SeldonConjTrans: public SeldonTranspose
  {
  public:
    class_SeldonConjTrans(): SeldonTranspose(2) {};
  };

  class_SeldonTrans SeldonTrans;
  class_SeldonNoTrans SeldonNoTrans;
  class_SeldonConjTrans SeldonConjTrans;

  //

  class SeldonDiag
  {
#ifdef SELDON_WITH_CBLAS
  protected:
    CBLAS_DIAG cblas_status_;
#endif
  protected:
    // 0: NonUnit, 1: Unit.
    int status_;
  public:
    SeldonDiag(int status)
    {
      status_ = status;
#ifdef SELDON_WITH_CBLAS
      if (status_ == 0)
	cblas_status_ = CblasNonUnit;
      else
	cblas_status_ = CblasUnit;
#endif
    }
#ifdef SELDON_WITH_CBLAS
    operator CBLAS_DIAG() const
    {
      return cblas_status_;
    }
#endif
    
    bool NonUnit() {return (status_ == 0);}
    bool Unit() {return (status_ == 1);}
  };

  class class_SeldonNonUnit: public SeldonDiag
  {
  public:
    class_SeldonNonUnit(): SeldonDiag(0) {};
  };

  class class_SeldonUnit: public SeldonDiag
  {
  public:
    class_SeldonUnit(): SeldonDiag(1) {};
  };

  class_SeldonNonUnit SeldonNonUnit;
  class_SeldonUnit SeldonUnit;

  //

  class SeldonUplo
  {
#ifdef SELDON_WITH_CBLAS
  protected:
    CBLAS_UPLO cblas_status_;
#endif
  protected:
    // 0: Upper, 1: Lower.
    int status_;
  public:
    SeldonUplo(int status)
    {
      status_ = status;
#ifdef SELDON_WITH_CBLAS
      if (status_ == 0)
	cblas_status_ = CblasUpper;
      else
	cblas_status_ = CblasLower;
#endif
    }
#ifdef SELDON_WITH_CBLAS
    operator CBLAS_UPLO() const
    {
      return cblas_status_;
    }
#endif
    
    operator char() const
    {
      return (status_ == 0) ? 'U' : 'L';
    }

    bool Upper() {return (status_ == 0);}
    bool Lower() {return (status_ == 1);}

    char Char() const
    {
      return (status_ == 0) ? 'U' : 'L';
    }
    char RevChar() const
    {
      return (status_ == 0) ? 'L' : 'U';
    }

  };

  SeldonUplo SeldonUpper(0);
  SeldonUplo SeldonLower(1);

  class SeldonConjugate
  {
  protected:
    // false: Unconj, true: Conj.
    bool status_;
  public:
    SeldonConjugate(bool status)
    {
      status_ = status;
    }
    inline bool Conj() const
    {
      return status_;
    }
  };

  SeldonConjugate SeldonUnconj(false);
  SeldonConjugate SeldonConj(true);
  

  //

  class SeldonSide
  {
#ifdef SELDON_WITH_CBLAS
  protected:
    CBLAS_SIDE cblas_status_;
#endif
  protected:
    // 0: Left, 1: Right.
    int status_;
  public:
    SeldonSide(int status)
    {
      status_ = status;
#ifdef SELDON_WITH_CBLAS
      if (status_ == 0)
	cblas_status_ = CblasLeft;
      else
	cblas_status_ = CblasRight;
#endif
    }
#ifdef SELDON_WITH_CBLAS
    operator CBLAS_SIDE() const
    {
      return cblas_status_;
    }
#endif
    
    bool Left() {return (status_ == 0);}
    bool Right() {return (status_ == 1);}
  };

  class class_SeldonLeft: public SeldonSide
  {
  public:
    class_SeldonLeft(): SeldonSide(0) {};
  };

  class class_SeldonRight: public SeldonSide
  {
  public:
    class_SeldonRight(): SeldonSide(1) {};
  };

  class_SeldonLeft SeldonLeft;
  class_SeldonRight SeldonRight;


} // namespace Seldon.


#include "Array3D/Array3D.cxx"
#include "Matrix/Matrix_Base.cxx"
#include "Matrix/Matrix_Pointers.cxx"
#include "Matrix/Matrix_Triangular.cxx"
#include "Matrix/Matrix_Symmetric.cxx"
#include "Matrix/Matrix_Hermitian.cxx"
#include "Matrix/Matrix_Sparse.cxx"
#include "Matrix/Matrix_ComplexSparse.cxx"
#include "Matrix/Matrix_SymSparse.cxx"
#include "Matrix/Matrix_SymComplexSparse.cxx"
#include "Matrix/Matrix_SymPacked.cxx"
#include "Matrix/Matrix_HermPacked.cxx"
#include "Matrix/Matrix_TriangPacked.cxx"
#include "Vector/Vector.cxx"
#include "Computation/Basic_Functions/Functions_Matrix.cxx"
#include "Computation/Basic_Functions/Functions_Vector.cxx"
#include "Computation/Basic_Functions/Functions_MatVect.cxx"

// Blas interface.
#ifdef SELDON_WITH_CBLAS
#include "Computation/Interfaces/Blas_1.cxx"
#include "Computation/Interfaces/Blas_2.cxx"
#include "Computation/Interfaces/Blas_3.cxx"
#endif

// Lapack interface.
#ifdef SELDON_WITH_LAPACK
#undef LAPACK_INTEGER
#define LAPACK_INTEGER int
#undef LAPACK_REAL
#define LAPACK_REAL float
#undef LAPACK_DOUBLEREAL
#define LAPACK_DOUBLEREAL double
#undef LAPACK_COMPLEX
#define LAPACK_COMPLEX void
#undef LAPACK_DOUBLECOMPLEX
#define LAPACK_DOUBLECOMPLEX void
#undef LAPACK_LOGICAL
#define LAPACK_LOGICAL int
#undef LAPACK_L_FP
#define LAPACK_L_FP int*
#undef LAPACK_FTNLEN
#define LAPACK_FTNLEN int*
extern "C"
{
#include "Computation/Interfaces/clapack.h"
}
#include "Computation/Interfaces/Lapack.cxx"
#endif

namespace Seldon
{


  typedef Vector<int, Vect_Full, SELDON_DEFAULT_ALLOCATOR<int> > IVect;
  typedef Vector<float, Vect_Full, SELDON_DEFAULT_ALLOCATOR<float> > SVect;
  typedef Vector<double, Vect_Full, SELDON_DEFAULT_ALLOCATOR<double> > DVect;
  typedef Vector<complex<float>, Vect_Full,
		 SELDON_DEFAULT_ALLOCATOR<complex<float> > > CVect;
  typedef Vector<complex<double>, Vect_Full,
		 SELDON_DEFAULT_ALLOCATOR<complex<double> > > ZVect;

  typedef Matrix<int, General, ColMajor,
		 SELDON_DEFAULT_ALLOCATOR<int> > IGCMat;
  typedef Matrix<float, General, ColMajor,
		 SELDON_DEFAULT_ALLOCATOR<float> > SGCMat;
  typedef Matrix<double, General, ColMajor,
		 SELDON_DEFAULT_ALLOCATOR<double> > DGCMat;
  typedef Matrix<complex<float>, General, ColMajor,
		 SELDON_DEFAULT_ALLOCATOR<complex<float> > > CGCMat;
  typedef Matrix<complex<double>, General, ColMajor,
		 SELDON_DEFAULT_ALLOCATOR<complex<double> > > ZGCMat;

  typedef Matrix<int, General, RowMajor,
		 SELDON_DEFAULT_ALLOCATOR<int> > IGRMat;
  typedef Matrix<float, General, RowMajor,
		 SELDON_DEFAULT_ALLOCATOR<float> > SGRMat;
  typedef Matrix<double, General, RowMajor,
		 SELDON_DEFAULT_ALLOCATOR<double> > DGRMat;
  typedef Matrix<complex<float>, General, RowMajor,
		 SELDON_DEFAULT_ALLOCATOR<complex<float> > > CGRMat;
  typedef Matrix<complex<double>, General, RowMajor,
		 SELDON_DEFAULT_ALLOCATOR<complex<double> > > ZGRMat;

  typedef Matrix<int, General, RowSparse,
		 SELDON_DEFAULT_ALLOCATOR<int> > IGRSMat;
  typedef Matrix<float, General, RowSparse,
		 SELDON_DEFAULT_ALLOCATOR<float> > SGRSMat;
  typedef Matrix<double, General, RowSparse,
		 SELDON_DEFAULT_ALLOCATOR<double> > DGRSMat;
  typedef Matrix<complex<float>, General, RowSparse,
		 SELDON_DEFAULT_ALLOCATOR<complex<float> > > CGRSMat;
  typedef Matrix<complex<double>, General, RowSparse,
		 SELDON_DEFAULT_ALLOCATOR<complex<double> > > ZGRSMat;

  typedef Matrix<int, General, ColSparse,
		 SELDON_DEFAULT_ALLOCATOR<int> > IGCSMat;
  typedef Matrix<float, General, ColSparse,
		 SELDON_DEFAULT_ALLOCATOR<float> > SGCSMat;
  typedef Matrix<double, General, ColSparse,
		 SELDON_DEFAULT_ALLOCATOR<double> > DGCSMat;
  typedef Matrix<complex<float>, General, ColSparse,
		 SELDON_DEFAULT_ALLOCATOR<complex<float> > > CGCSMat;
  typedef Matrix<complex<double>, General, ColSparse,
		 SELDON_DEFAULT_ALLOCATOR<complex<double> > > ZGCSMat;


} // namespace Seldon.

#define SELDON_FILE_SELDON_HXX
#endif
