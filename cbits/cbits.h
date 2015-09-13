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

namespace HsMultivac {
  typedef COrthogonalMesh<double> MeshType;
  typedef CSpeedCallback<double> SpeedType;
  typedef Curve<double> CurveType;
  typedef CSetOfPoints<double> InitialCurveType;
  typedef COrthogonalLevelSet<double> LevelSetType;
  typedef CNarrowBandExtension<double> InitializerType;
}

extern "C" {
#endif // __cplusplus

#include "speedcallback.hpp"

/*
 * Exported Types
 */

HsException multivacError(const char* msg);

// MVMeshH

typedef void *MVMeshH;

MVMeshH
MVNewMesh ( double Xmin, double Xmax, double Ymin, double YMax
          , double Nx, double Ny);
void MVDestroyMesh(MVMeshH);


// MVSpeedH

typedef void *MVSpeedH;

MVSpeedH
MVNewSpeed ( FastMarchSpeedFunc fm, NarrowBandSpeedFunc nb
           , MaxFSpeedFunc f1, MaxFSpeedFunc f2
           , int depPos, int depTime, int depNormal, int depCurv);

void MVDestroySpeed(MVSpeedH);

// MVFrontH

typedef struct {
  double x;
  double y;
} MVPoint;

typedef enum {MVO_TRIGONOMETRIC=0, MVO_REVERSE=1} MVOrientation;


typedef void *MVFrontH;

MVFrontH MVNewFront (int numPoints, MVOrientation orientation, MVPoint *points);

MVPoint *
MVGetFrontPoints (MVFrontH front, int *numPoints, MVOrientation *orientation);

void MVDestroyFront(MVFrontH);

// Simulate

HsException Simulate( MVMeshH meshH, MVSpeedH speedH, MVFrontH frontH
                    , int NbIterations, double FinalTime);


#ifdef __cplusplus
}
#endif

#endif // CBITS_H
