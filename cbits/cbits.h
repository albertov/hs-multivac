#ifndef CBITS_H
#define CBITS_H

#ifdef __cplusplus
/*
 * Internal CPP types
 */

#include "multivac.hxx"
#include "speedcallback.cpp"
#include "memorysaver.hxx"


using namespace Multivac;
using namespace HsMultivac;

namespace HsMultivac {
  typedef COrthogonalMesh<double> MeshType;
  typedef CSpeedCallback<double> SpeedType;
  typedef Curve<double> CurveType;
  typedef CSetOfPoints<double> InitialCurveType;
  typedef COrthogonalLevelSet<double> LevelSetType;
  typedef CNarrowBandExtension<double> InitializerType;
  typedef CNarrowBandFirstOrderEngquistOsher<double> UpdaterType;
  typedef vector<CurveType> FrontArrayType;
  typedef CMemorySaver<double> SaverType;
  typedef CSimulator<double, MeshType, SpeedType, InitialCurveType,
      LevelSetType, InitializerType, UpdaterType, SaverType> SimulatorType;

}

extern "C" {
#endif // __cplusplus

#include "speedcallback.hpp"

#define MVO_UNKNOWN         (-1)
#define MVO_TRIGONOMETRIC   (0)
#define MVO_REVERSE         (1)

/*
 * Exported Types
 */

HsException multivacError(const char* msg);

// MVMeshH

typedef void *MVMeshH;

MVMeshH
MVNewMesh ( double Xmin, double Xmax, double Ymin, double Ymax
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


typedef void *MVFrontH;

MVFrontH MVNewFront (int numPoints, int orientation, MVPoint *points);
MVFrontH MVNewEmptyFront ();

MVPoint *
MVGetFrontPoints (MVFrontH front, int *numPoints, int *orientation);

void MVDestroyFront(MVFrontH);

// MVFrontArrayH

typedef void *MVFrontArrayH;

MVFrontArrayH MVNewFrontArray ();
int MVGetNumFronts (MVFrontArrayH array);
void MVCopyFrontAt (MVFrontArrayH array, int ix, MVFrontH front);
void MVDestroyFrontArray(MVFrontArrayH array);

// Simulate

HsException Simulate( MVMeshH meshH, MVSpeedH speedH, MVFrontH frontH
                    , MVFrontArrayH
                    , int NbIterations, double FinalTime, double Period);


#ifdef __cplusplus
}
#endif

#endif // CBITS_H
