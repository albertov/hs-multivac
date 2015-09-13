#include "cbits.h"


MVMeshH
MVNewMesh ( double Xmin, double Xmax, double Ymin, double Ymax
          , double Nx, double Ny)
{
  MeshType *mesh = new MeshType(Xmin, Xmax, Ymin, Ymax, Nx, Ny);
  return static_cast<MVMeshH>(mesh);
}

void MVDestroyMesh(MVMeshH meshH)
{
  MeshType *mesh = static_cast<MeshType*>(meshH);
  delete mesh;
}

MVSpeedH
MVNewSpeed ( FastMarchSpeedFunc fm, NarrowBandSpeedFunc nb
           , MaxFSpeedFunc f1, MaxFSpeedFunc f2
           , int depPos, int depTime, int depNormal, int depCurv)
{
  SpeedFuncCallbacks cb = {fm, nb, f1, f2, depPos, depTime, depNormal, depCurv};
  SpeedType *speed = new SpeedType(cb);
  return static_cast<MVSpeedH>(speed);
}

void MVDestroySpeed(MVSpeedH speedH)
{
  SpeedType *speed = static_cast<SpeedType*>(speedH);
  delete speed;
}

MVFrontH
MVNewFront (int numPoints, int orientation, MVPoint *points)
{
  CurveType *front = new CurveType();
  front->SetOrientation(orientation);
  Vector<double> p(2);
  for (int i=0; i<numPoints; i++) {
    p(0) = points[i].x;
    p(1) = points[i].y;
    front->AddPoint(p);
  }
  return static_cast<MVFrontH>(front);
}

MVFrontH
MVNewEmptyFront ()
{
  CurveType *front = new CurveType();
  return static_cast<MVFrontH>(front);
}

MVPoint *
MVGetFrontPoints (MVFrontH frontH, int *numPoints, int *orientation)
{
  MVPoint *points = NULL;
  CurveType *front = static_cast<CurveType*>(frontH);
  *orientation = front->GetOrientation();
  int n = *numPoints = front->GetNbPoints();
  if (n>0) {
    List<Vector<double> > &list = front->GetPoints();
    points = static_cast<MVPoint*>(malloc(n * sizeof(MVPoint)));
    list.GoToTheHead();
    for (int i=0; i<n; i++) {
      Vector<double> p = list.GetCurrentValue();
      points[i].x = p(0);
      points[i].y = p(1);
      list.GoToNext_StopAtTheTail();
    }
  }
  return points;
};

void MVDestroyFront(MVFrontH frontH)
{
  CurveType *front = static_cast<CurveType*>(frontH);
  delete front;
}


MVFrontArrayH MVNewFrontArray ()
{
  FrontArrayType *fronts = new FrontArrayType;
  return static_cast<MVFrontArrayH>(fronts);
}

int MVGetNumFronts (MVFrontArrayH array)
{
  return static_cast<FrontArrayType*>(array)->size();
}

void MVCopyFrontAt (MVFrontArrayH array, int ix, MVFrontH frontH)
{
  CurveType &front = static_cast<FrontArrayType*>(array)->at(ix);
  CurveType *frontDest = static_cast<CurveType*>(frontH);
  frontDest->Copy(front);
}

void MVDestroyFrontArray(MVFrontH array) {
  FrontArrayType *fronts = static_cast<FrontArrayType*>(array);
  delete fronts;
}



HsException Simulate( MVMeshH meshH, MVSpeedH speedH, MVFrontH frontH
                    , MVFrontArrayH frontsH
                    , int NbIterations, double FinalTime, double Period)
{
  TRY;

  MeshType Mesh(*static_cast<MeshType*>(meshH));
  SpeedType F(*static_cast<SpeedType*>(speedH));
  InitialCurveType InitialCurve(*static_cast<CurveType*>(frontH));
  FrontArrayType* fronts(static_cast<FrontArrayType*>(frontsH));

  LevelSetType Phi;
  InitializerType Initializer;

  UpdaterType Updater(6, 3, 1);

  SaverType Saver(*fronts, Period);

  SimulatorType Simulator( Mesh, F, InitialCurve, Phi, Initializer, Updater
                         , Saver, NbIterations, FinalTime);
  Simulator.Init();
  Simulator.Run();

  CATCH;
  return 0;
}
