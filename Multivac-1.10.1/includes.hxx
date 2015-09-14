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


#ifndef FILE_INCLUDES_HXX


/***************************/
/**** STANDARD INCLUDES ****/

#include <cstdio>
#include <string>
using std::string;

/**** STANDARD INCLUDES ****/
/***************************/


/****************/
/**** SELDON ****/

#include "seldon/Seldon.hxx"

namespace Multivac
{
  using namespace Seldon;
}

/**** SELDON ****/
/****************/


/***************/
/**** TALOS ****/

#include "talos/Talos.hxx"

namespace Multivac
{
  using namespace Talos;
}

/**** TALOS ****/
/***************/


//! Convenient structure to catch exceptions.
/*!
  Use TRY and END to catch exceptions thrown by Seldon or Multivac:
  \par
  \par [...]
  \par TRY
  \par
  \par Code that could throw an exception.
  \par
  \par END
  \par
  \par [...]
  \par
  \par When an exception is caught, information is shown to solve the problem
  (name of the involved function and comments).

*/
#ifdef TRY
#undef TRY
#endif
#define TRY try{

#ifdef END
#undef END
#endif
#define END \
}\
catch (Seldon::Error& Err)\
{\
Err.What();\
return 1;\
}\
catch (Multivac::CError& Err)\
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
cerr << "Unknown exception..." <<endl;\
return 1;\
}


/******************/
/**** GEOMETRY ****/

#include "geometry/functions.cxx"
#include "geometry/arrayheap.cxx"
#include "geometry/list.cxx"
#include "geometry/curve.cxx"

/**** GEOMETRY ****/
/******************/


/****************/
/**** MESHES ****/

#include "meshes/baseclass.cxx"
#include "meshes/orthogonal.cxx"

/**** MESHES ****/
/****************/


/*************************/
/**** SPEED FUNCTIONS ****/

#include "speedfunctions/baseclass.cxx"
#include "speedfunctions/constantspeed.cxx"
#include "speedfunctions/piecewiseconstantspeed.cxx"
#include "speedfunctions/firemodel.cxx"
#include "speedfunctions/simplifiedfiremodel.cxx"
#include "speedfunctions/gradientsegmentation.cxx"
#include "speedfunctions/imageintensity.cxx"

/**** SPEED FUNCTIONS ****/
/*************************/


/*****************************/
/**** LEVEL SET FUNCTIONS ****/

#include "levelsets/baseclass.cxx"
#include "levelsets/orthogonal.cxx"

/**** LEVEL SET FUNCTIONS ****/
/*****************************/


/************************/
/**** INITIAL CURVES ****/

#include "initialcurves/baseclass.cxx"
#include "initialcurves/circle.cxx"
#include "initialcurves/island.cxx"
#include "initialcurves/island0.cxx"
#include "initialcurves/setofpoints.cxx"
#include "initialcurves/threecircles.cxx"
#include "initialcurves/twocircles.cxx"

/**** INITIAL CURVES ****/
/************************/


/******************/
/**** UPDATERS ****/

#include "updaters/baseclass.cxx"
#include "updaters/fastmarchingfirstorderengquistosher.cxx"
#include "updaters/narrowbandfirstorderengquistosher.cxx"
#include "updaters/narrowbandeno2engquistosher.cxx"
#include "updaters/narrowbandfirstorderlaxfriedrichs.cxx"
#include "updaters/chanvese.cxx"

/**** UPDATERS ****/
/******************/


/**********************/
/**** INITIALIZERS ****/

#include "initializers/baseclass.cxx"
#include "initializers/narrowbandneverinit.cxx"
#include "initializers/narrowbandextension.cxx"
#include "initializers/fastmarchingneverinit.cxx"

/**** INITIALIZERS ****/
/**********************/


/***************/
/**** NODES ****/

#include "nodes.cxx"

/**** NODES ****/
/***************/


/****************/
/**** SAVERS ****/

#include "savers/baseclass.cxx"
#include "savers/curvessaver.cxx"
#include "savers/neversave.cxx"
#include "savers/savelastcurve.cxx"
#include "savers/saveattheend.cxx"
#include "savers/track.cxx"

/**** SAVERS ****/
/****************/


/********************/
/**** SIMULATORS ****/

#include "simulator.cxx"

/**** SIMULATORS ****/
/********************/


#define FILE_INCLUDES_HXX
#endif
