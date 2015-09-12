#ifndef CBITS_H
#define CBITS_H

#ifdef __cplusplus
/*
 * Internal CPP types
 */

#include "multivac.hxx"
#include "speedcallback.cpp"


using namespace Multivac;
using namespace HsMultivac;

typedef COrthogonalMesh<double> MeshType;
typedef CSpeedCallback<double> SpeedType;

extern "C" {
#endif // __cplusplus

#include "speedcallback.hpp"

/*
 * Exported Types
 */

// MVMeshH

typedef void *MVMeshH;

MVMeshH
MVNewMesh ( double Xmin, double Xmax, double Ymin, double YMax
          , double Nx, double Ny);
void MVDestroyMesh(MVMeshH);


// MVSpeedH
//
typedef void *MVSpeedH;

MVSpeedH
MVNewSpeed ( FastMarchSpeedFunc fm, NarrowBandSpeedFunc nb
           , MaxFSpeedFunc f1, MaxFSpeedFunc f2
           , int depPos, int depTime, int depNormal, int depCurv);

void MVDestroySpeed(MVSpeedH);

void Simulate(MVMeshH meshPtr, MVSpeedH speedPtr);


#ifdef __cplusplus
}
#endif

#endif // CBITS_H
