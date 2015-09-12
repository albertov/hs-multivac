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

#ifndef SELDON_FILE_STORAGE_HXX

namespace Seldon
{


  //////////////////////
  // GENERAL MATRICES //
  //////////////////////


  class ColMajor
  {
  public:
    static int GetFirst(int i, int j)
    {
      return j;
    }
    static int GetSecond(int i, int j)
    {
      return i;
    }
  };


  class RowMajor
  {
  public:
    static int GetFirst(int i, int j)
    {
      return i;
    }
    static int GetSecond(int i, int j)
    {
      return j;
    }
  };



  /////////////
  // VECTORS //
  /////////////


  class Vect_Full
  {
  };



  ////////////
  // SPARSE //
  ////////////


  class ColSparse
  {
  public:
    static int GetFirst(int i, int j)
    {
      return j;
    }
    static int GetSecond(int i, int j)
    {
      return i;
    }
  };


  class RowSparse
  {
  public:
    static int GetFirst(int i, int j)
    {
      return i;
    }
    static int GetSecond(int i, int j)
    {
      return j;
    }
  };


  class ColComplexSparse
  {
  public:
    static int GetFirst(int i, int j)
    {
      return j;
    }
    static int GetSecond(int i, int j)
    {
      return i;
    }
  };


  class RowComplexSparse
  {
  public:
    static int GetFirst(int i, int j)
    {
      return i;
    }
    static int GetSecond(int i, int j)
    {
      return j;
    }
  };


  class ColSymSparse
  {
  public:
    static int GetFirst(int i, int j)
    {
      return j;
    }
    static int GetSecond(int i, int j)
    {
      return i;
    }
  };


  class RowSymSparse
  {
  public:
    static int GetFirst(int i, int j)
    {
      return i;
    }
    static int GetSecond(int i, int j)
    {
      return j;
    }
  };


  class ColSymComplexSparse
  {
  public:
    static int GetFirst(int i, int j)
    {
      return j;
    }
    static int GetSecond(int i, int j)
    {
      return i;
    }
  };


  class RowSymComplexSparse
  {
  public:
    static int GetFirst(int i, int j)
    {
      return i;
    }
    static int GetSecond(int i, int j)
    {
      return j;
    }
  };



  ///////////////
  // SYMMETRIC //
  ///////////////


  class ColSymPacked
  {
  public:
    static int GetFirst(int i, int j)
    {
      return j;
    }
    static int GetSecond(int i, int j)
    {
      return i;
    }
  };


  class RowSymPacked
  {
  public:
    static int GetFirst(int i, int j)
    {
      return i;
    }
    static int GetSecond(int i, int j)
    {
      return j;
    }
  };


  class ColSym
  {
  public:
    static int GetFirst(int i, int j)
    {
      return j;
    }
    static int GetSecond(int i, int j)
    {
      return i;
    }
  };


  class RowSym
  {
  public:
    static int GetFirst(int i, int j)
    {
      return i;
    }
    static int GetSecond(int i, int j)
    {
      return j;
    }
  };



  ///////////////
  // HERMITIAN //
  ///////////////


  class ColHerm
  {
  public:
    static int GetFirst(int i, int j)
    {
      return j;
    }
    static int GetSecond(int i, int j)
    {
      return i;
    }
  };


  class RowHerm
  {
  public:
    static int GetFirst(int i, int j)
    {
      return i;
    }
    static int GetSecond(int i, int j)
    {
      return j;
    }
  };


  class ColHermPacked
  {
  public:
    static int GetFirst(int i, int j)
    {
      return j;
    }
    static int GetSecond(int i, int j)
    {
      return i;
    }
  };


  class RowHermPacked
  {
  public:
    static int GetFirst(int i, int j)
    {
      return i;
    }
    static int GetSecond(int i, int j)
    {
      return j;
    }
  };



  ////////////////
  // TRIANGULAR //
  ////////////////


  class ColUpTriang
  {
  public:
    static int GetFirst(int i, int j)
    {
      return j;
    }
    static int GetSecond(int i, int j)
    {
      return i;
    }
    static bool UpLo()
    {
      return true;
    }
  };


  class ColLoTriang
  {
  public:
    static int GetFirst(int i, int j)
    {
      return j;
    }
    static int GetSecond(int i, int j)
    {
      return i;
    }
    static bool UpLo()
    {
      return false;
    }
  };


  class RowUpTriang
  {
  public:
    static int GetFirst(int i, int j)
    {
      return i;
    }
    static int GetSecond(int i, int j)
    {
      return j;
    }
    static bool UpLo()
    {
      return true;
    }
  };


  class RowLoTriang
  {
  public:
    static int GetFirst(int i, int j)
    {
      return i;
    }
    static int GetSecond(int i, int j)
    {
      return j;
    }
    static bool UpLo()
    {
      return false;
    }
  };


  class ColUpTriangPacked
  {
  public:
    static int GetFirst(int i, int j)
    {
      return j;
    }
    static int GetSecond(int i, int j)
    {
      return i;
    }
    static bool UpLo()
    {
      return true;
    }
  };


  class ColLoTriangPacked
  {
  public:
    static int GetFirst(int i, int j)
    {
      return j;
    }
    static int GetSecond(int i, int j)
    {
      return i;
    }
    static bool UpLo()
    {
      return false;
    }
  };


  class RowUpTriangPacked
  {
  public:
    static int GetFirst(int i, int j)
    {
      return i;
    }
    static int GetSecond(int i, int j)
    {
      return j;
    }
    static bool UpLo()
    {
      return true;
    }
  };


  class RowLoTriangPacked
  {
  public:
    static int GetFirst(int i, int j)
    {
      return i;
    }
    static int GetSecond(int i, int j)
    {
      return j;
    }
    static bool UpLo()
    {
      return false;
    }
  };


} // namespace Seldon.

#define SELDON_FILE_STORAGE_HXX
#endif
