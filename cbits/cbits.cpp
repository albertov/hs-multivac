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


int Simulate( MVMeshH meshPtr, MVSpeedH speedPtr, int NbIterations
             , double FinalTime)
{
  TRY;

  MeshType Mesh(*static_cast<MeshType*>(meshPtr));
  SpeedType F(*static_cast<SpeedType*>(speedPtr));


  ///////////////////
  // INITIAL CURVE //
  ///////////////////

  // Choose the initial curve:
  //   1) CCircle<double>
  //   2) CTwoCircles<double>
  //   3) CThreeCircles<double>
  //   4) CIsland<double>
  //   5) CIsland0<double>
  //   6) CSetOfPoints<double>
  typedef CCircle<double> InitialCurveType;

  // For a circle, center coordinates and radius.
  double CircleCenterX = 1.0;
  double CircleCenterY = 1.5;
  double CircleRadius = 0.5;

  // For a second circle, center coordinates and radius.
  double CircleCenterX0 = 1.65;
  double CircleCenterY0 = 1.6;
  double CircleRadius0 = 0.3;

  // For a third circle, center coordinates and radius.
  double CircleCenterX1 = 0.50;
  double CircleCenterY1 = 1.0;
  double CircleRadius1 = 0.25;

  // File containing a front defined by a set of points.
  string InitialFrontFile = "[full path]/[set].pts";
  // Is the (Xmin, Ymin) outside the front?
  bool origin_out = false;

  // If set to 'true', outside and inside are swapped.
  bool reversed = false;

  // Choose the right constructor:
  //   1) InitialCurve(CircleCenterX, CircleCenterY, CircleRadius,
  //                   reversed)
  //   2, 4) InitialCurve(CircleCenterX, CircleCenterY, CircleRadius,
  //                      CircleCenterX0, CircleCenterY0, CircleRadius0,
  //                      reversed)
  //   3, 5) InitialCurve(CircleCenterX, CircleCenterY, CircleRadius,
  //                      CircleCenterX0, CircleCenterY0, CircleRadius0,
  //                      CircleCenterX1, CircleCenterY1, CircleRadius1,
  //                      reversed)
  //   6) InitialCurve(InitialFrontFile, 0, origin_out ^ reversed)
  InitialCurveType InitialCurve(CircleCenterX, CircleCenterY, CircleRadius);


  ////////////////////////
  // LEVEL SET FUNCTION //
  ////////////////////////

  // Choose the level set function type:
  //   1) COrthogonalLevelSet<double>
  //   2) --
  typedef COrthogonalLevelSet<double> LevelSetType;

  // Choose the right constructor:
  //   1) Phi
  //   2) --
  LevelSetType Phi;


  /////////////////
  // INITIALIZER //
  /////////////////

  // Choose an initializer:
  //   1) CNarrowBandNeverInit<double> *
  //   2) CNarrowBandExtension<double> *
  //   3) CFastMarchingNeverInit<double> **
  //         * Narrow band level set method
  //        ** Fast marching method
  typedef CNarrowBandNeverInit<double> InitializerType;

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

  END;
  return 0;
}
