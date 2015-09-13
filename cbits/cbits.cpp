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
MVNewFront (int numPoints, MVOrientation orientation, MVPoint *points)
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

MVPoint *
MVGetFrontPoints (MVFrontH frontH, int *numPoints, MVOrientation *orientation)
{
  MVPoint *points = NULL;
  CurveType *front = static_cast<CurveType*>(frontH);
  *orientation = static_cast<MVOrientation>(front->GetOrientation());
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


HsException Simulate( MVMeshH meshH, MVSpeedH speedH, MVFrontH frontH
                    , int NbIterations , double FinalTime)
{
  TRY;

  MeshType Mesh(*static_cast<MeshType*>(meshH));
  SpeedType F(*static_cast<SpeedType*>(speedH));
  InitialCurveType InitialCurve(*static_cast<InitialCurveType*>(frontH));
  LevelSetType Phi;
  InitializerType Initializer;


  /////////////
  // UPDATER //
  /////////////

  // Choose an initializer:
  //   1) CNarrowBandFirstOrderEngquistOsher<double> *
  //   2) CNarrowBandFirstOrderLaxFriedrichs<double> *
  //   3) CNarrowBandEno2EngquistOsher<double> *
  //   4) CFastMarchingFirstOrderEngquistOsher<double> **
  //         * Narrow band level set method
  //        ** Fast marching method
  typedef CNarrowBandFirstOrderEngquistOsher<double> UpdaterType;

  // Choose the right constructor and choose its parameters:
  //   1, 2, 3) Updater(TubeSemiWidth, BarrierWidth, OutSpaceWidth)
  //              | TubeSemiWidth: number of grid points on each side of the front.
  //              | BarrierWidth: number of grid points (on each side of the front)
  //                        that imply tube reconstruction when reached.
  //              | OutSpaceWidth: number of grid points (on each side of the front)
  //                         that must not be reached.
  //              | Examples: Updater(6, 3, 1) or Updater(12, 5, 1).
  //   4) Updater(TMax)
  //        | TMax: time greater than all arrival times that will be computed.
  UpdaterType Updater(6, 3, 1);


  ///////////
  // SAVER //
  ///////////

  // Choose the saver type:
  //   1) CNeverSave<double>
  //        | Nothing is saved.
  //   2) CCurvesSaver<double>
  //        | The front is constructed and saved at each time step.
  //        | Not relevant for the fast marching method.
  //   3) CSaveLastCurve<double>
  //        | Saves the last curve, after all calculations.
  //        | Not relevant for the fast marching method.
  //   4) CSaveAtTheEnd<double>
  //        | Saves the level set, after all calculations.
  //        | Relevant for the fast marching method.
  typedef CCurvesSaver<double> SaverType;

  // Number of curves that will be saved.
  // Set NbCurves to 0 in order to save all curves.
  int NbCurves = 10;

  // Output directory.
  string Directory = "results/";

  // Stores time at each time step.
  string TimeFile = Directory + "Time";

  // Stores front points.
  string CurvesFile = Directory + "Curves";
  // Stores number of points on fronts at each time step.
  string CurveLengthsFile = Directory + "CurveLengths";
  // Stores the level set function.
  string PhiFile = Directory + "Phi";
  // Stores the speed function.
  string FFile = Directory + "F";

  // Stores grid coordinates along the (x'x) axe.
  string XFile = Directory + "X";
  // Stores grid coordinates along the (y'y) axe.
  string YFile = Directory + "Y";
  // Stores mesh points.
  string PointsFile = Directory + "Points";
  // Stores mesh edges.
  string EdgesFile = Directory + "Edges";
  // Stores mesh triangles.
  string TrianglesFile = Directory + "Triangles";

  // Indicates the number of iterations between saves.
  int Period = NbIterations / (NbCurves==0?NbIterations:NbCurves);

  SaverType Saver(TimeFile, CurvesFile,
      CurveLengthsFile, PhiFile, FFile, XFile, YFile,
      PointsFile, EdgesFile, TrianglesFile, Period);



  /////////////////////
  /////////////////////
  ////  SIMULATOR  ////
  /////////////////////
  /////////////////////

  CSimulator<double, MeshType, SpeedType, InitialCurveType,
    LevelSetType, InitializerType, UpdaterType, SaverType>
    Simulator(Mesh, F, InitialCurve, Phi,
        Initializer, Updater, Saver,
        NbIterations, FinalTime);


  ////////////////////////////////
  // SIMULATION INITIALIZATIONS //
  ////////////////////////////////

  Simulator.Init();


  ////////////////
  // SIMULATION //
  ////////////////

  Simulator.Run();

  CATCH;
  return 0;
}
