// Copyright (C) 2002-2004 Vivien Mallet
//
// This file is part of Multivac library.
// Multivac library provides front-tracking algorithms.
//
// Multivac is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// Multivac is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License (file "license") for more details.
//
// For more information, please see the Multivac home page:
//     http://spacetown.free.fr/fronts/


#ifndef FILE_MULTIVAC_HXX


/**** Exception handling ****/

/**********/

// Debug levels.
//! Check all!
#ifdef MULTIVAC_DEBUG_LEVEL_4

#ifndef SELDON_DEBUG_LEVEL_4
#define SELDON_DEBUG_LEVEL_4
#endif

#ifndef MULTIVAC_DEBUG_LEVEL_3
#define MULTIVAC_DEBUG_LEVEL_3
#endif

#endif

/**********/

#ifdef MULTIVAC_DEBUG_LEVEL_3

#ifndef SELDON_DEBUG_LEVEL_3
#define SELDON_DEBUG_LEVEL_3
#endif

//! Check indices.
#ifndef MULTIVAC_CHECK_BOUNDARIES
#define MULTIVAC_CHECK_BOUNDARIES
#endif

#ifndef MULTIVAC_DEBUG_LEVEL_2
#define MULTIVAC_DEBUG_LEVEL_2
#endif

#endif

/**********/

#ifdef MULTIVAC_DEBUG_LEVEL_2

#ifndef SELDON_DEBUG_LEVEL_2
#define SELDON_DEBUG_LEVEL_2
#endif

//! Check memory management.
#ifndef MULTIVAC_CHECK_MEMORY
#define MULTIVAC_CHECK_MEMORY
#endif

//! Check IO operations.
#ifndef MULTIVAC_CHECK_IO
#define MULTIVAC_CHECK_IO
#endif

//! Check instability.
#ifndef MULTIVAC_CHECK_INSTABILITY
#define MULTIVAC_CHECK_INSTABILITY
#endif

#ifndef MULTIVAC_DEBUG_LEVEL_1
#define MULTIVAC_DEBUG_LEVEL_1
#endif

#endif

/**********/

#ifdef MULTIVAC_DEBUG_LEVEL_1

#ifndef MULTIVAC_DEBUG_LEVEL_2
//! Disable all exception specifications.
#ifndef MULTIVAC_WITHOUT_THROW
#define MULTIVAC_WITHOUT_THROW
#endif
#endif

#ifndef MULTIVAC_CHECK_COMPATIBILITY
#define MULTIVAC_CHECK_COMPATIBILITY
#endif

#ifndef SELDON_DEBUG_LEVEL_1
#define SELDON_DEBUG_LEVEL_1
#endif

#endif

/**********/


/**** Convenient macros ****/

#ifndef ERR
#define ERR(x) cerr << "Hermes - " #x << endl
#endif

#ifndef DISPLAY
#define DISPLAY(x) cerr << #x ": " << x << endl
#endif

#ifndef DISP
#define DISP(x) cerr << #x ": " << x << endl
#endif


/**** Class declarations ****/

namespace std {}

namespace Multivac
{

  using namespace std;

  template <class T> class CMesh;
  template <class T> class CLevelSet;
  template <class T> class CSpeedFunction;
  template <class T> class CInitialCurve;
  template <class T> class CInitializer;
  template <class T> class CUpdater;
  template <class T> class CSpeedFunction;
  template <class T> class CSaver;
  
}


/**** Basic functions ****/

#include <cmath>
#include <algorithm>

namespace Multivac
{

  using namespace std;
    
  template <class T>
  inline T sign(T x)
  {
    return (x<0?-1.0:1.0);
  }

}  // namespace Multivac.


/**** Includes ****/

#include "includes.hxx"


#define FILE_MULTIVAC_HXX
#endif
