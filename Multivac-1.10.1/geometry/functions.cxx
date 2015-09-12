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


#ifndef FILE_GEOMETRY_FUNCTIONS_CXX


#ifndef _limit
#define _limit 0.000001
#endif

namespace Multivac
{

  template <class T>
  inline void ComputeNormal(T& Dmx, T& Dpx,
			    T& Dmy, T& Dpy,
			    T& NormalX, T& NormalY)
  {

    // To speed up calculations.
    T Dmx2 = Dmx * Dmx;
    T Dpx2 = Dpx * Dpx;
    T Dmy2 = Dmy * Dmy;
    T Dpy2 = Dpy * Dpy;

    NormalX = 0.0;
    NormalY = 0.0;

    T norm = sqrt( Dpx2 + Dpy2 );
    // If 'norm' is significant.
    if ( norm > _limit)
      {
	NormalX += Dpx / norm;
	NormalY += Dpy / norm;
      }
		  
    norm = sqrt( Dmx2 + Dpy2 );
    // If 'norm' is significant.
    if ( norm > _limit)
      {
	NormalX += Dmx / norm;
	NormalY += Dpy / norm;
      }
		  
    norm = sqrt( Dpx2 + Dmy2 );
    // If 'norm' is significant.
    if ( norm > _limit)
      {
	NormalX += Dpx / norm;
	NormalY += Dmy / norm;
      }
		  
    norm = sqrt( Dmx2 + Dmy2 );
    // If 'norm' is significant.
    if ( norm > _limit)
      {
	NormalX += Dmx / norm;
	NormalY += Dmy / norm;
      }
		  
    // Norm of the normal vector.
    norm = sqrt( NormalX * NormalX + NormalY * NormalY );
		  
    if (norm > _limit)
      {
	NormalX = NormalX / norm;
	NormalY = NormalY / norm;
      }

    return;

  }

}  // namespace Multivac.


#define FILE_GEOMETRY_FUNCTIONS_CXX
#endif
