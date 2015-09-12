/************************
 * INCLUDES AND OPTIONS *
 ************************/

// Seldon library and Multivac project provide error management
// using exception handling.  Exceptions that may be raised
// are selected according to debug levels.
#define MULTIVAC_DEBUG_LEVEL_2

// Define MULTIVAC_REPORT if you want to clock the time needed for
// updates and initializations (results will be displayed on screen).
#define MULTIVAC_REPORT

// Multivac includes.
#include "multivac.hxx"

#include "cbits.h"

using namespace Multivac;


/*****************
 * MAIN FUNCTION *
 *****************/

int my_main()
{

  // To catch exceptions.
  TRY;
  
  // real type: double, float...
  typedef double real;


  //////////
  // TIME //
  //////////

  // Final time of the simulation, the initial time being 0.
  real FinalTime = 0.1;
  // Time step.
  real Delta_t = 0.0001;

  // Number of iterations.
  int NbIterations = int (FinalTime / Delta_t);


  ///////////////////
  // DOMAIN & MESH //
  ///////////////////

  // Choose the type of the mesh:
  //   1) COrthogonalMesh<real>
  //   2) --
  typedef COrthogonalMesh<real> MeshType;

  // Domain bounds (the domain is a rectangle).
  real Xmin = 0.0;
  real Xmax = 3.0;
  real Ymin = 0.0;
  real Ymax = 3.0;

  // For an orthogonal mesh, Nx and Ny are the number of
  // grid points along (x'x) and (y'y) (respectively).
  int Nx = 301;
  int Ny = Nx;

  // Choose the right constructor:
  //   1) Mesh(Xmin, Xmax, Ymin, Ymax, Nx, Ny)
  //   2) --
  MeshType Mesh(Xmin, Xmax, Ymin, Ymax, Nx, Ny);


  ////////////////////
  // SPEED FUNCTION //
  ////////////////////

  // Choose the speed function:
  //   1) CConstantSpeed<real>
  //   2) CPiecewiseConstantSpeed<real>
  //   3) CFireModel<real>
  //   4) CSimplifiedFireModel<real>
  typedef CSimplifiedFireModel<real> SpeedType;

  // Speed rate for a constant speed function.
  real SpeedRate = 0.5;

  // Second speed rate for a piecewise-constant speed function.
  real SpeedRate0 = 0.2;

  // Limit (for a piecewise-constant speed function).
  real Limit = 0.9;

  // Parameters for both simplified and full fire model.
  real U = 100.0;
  real m = 1.5;
  real c1 = 0.5;
  real epsilon0 = 0.2;

  // Last parameter for the simplified fire model.
  real alpha = 0.9;

  // Parameters for the full fire model.
  real a = 0.1;
  real b = 1.0;
  real epsilon1 = 0.003;

  // Choose the right constructor:
  //   1) F(SpeedRate)
  //   2) F(SpeedRate, SpeedRate0, Limit)
  //   3) F(U, m, c1, epsilon0, a, b, epsilon1)
  //   4) F(U, m, c1, epsilon0, alpha)
  SpeedType F(U, m, c1, epsilon0, alpha);


  ///////////////////
  // INITIAL CURVE //
  ///////////////////

  // Choose the initial curve:
  //   1) CCircle<real>
  //   2) CTwoCircles<real>
  //   3) CThreeCircles<real>
  //   4) CIsland<real>
  //   5) CIsland0<real>
  //   6) CSetOfPoints<real>
  typedef CCircle<real> InitialCurveType;

  // For a circle, center coordinates and radius.
  real CircleCenterX = 1.0;
  real CircleCenterY = 1.5;
  real CircleRadius = 0.5;

  // For a second circle, center coordinates and radius.
  real CircleCenterX0 = 1.65;
  real CircleCenterY0 = 1.6;
  real CircleRadius0 = 0.3;

  // For a third circle, center coordinates and radius.
  real CircleCenterX1 = 0.50;
  real CircleCenterY1 = 1.0;
  real CircleRadius1 = 0.25;

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
  //   1) COrthogonalLevelSet<real>
  //   2) --
  typedef COrthogonalLevelSet<real> LevelSetType;

  // Choose the right constructor:
  //   1) Phi
  //   2) --
  LevelSetType Phi;


  /////////////////
  // INITIALIZER //
  /////////////////

  // Choose an initializer:
  //   1) CNarrowBandNeverInit<real> *
  //   2) CNarrowBandExtension<real> *
  //   3) CFastMarchingNeverInit<real> **
  //         * Narrow band level set method
  //        ** Fast marching method
  typedef CNarrowBandNeverInit<real> InitializerType;

  InitializerType Initializer;


  /////////////
  // UPDATER //
  /////////////

  // Choose an initializer:
  //   1) CNarrowBandFirstOrderEngquistOsher<real> *
  //   2) CNarrowBandFirstOrderLaxFriedrichs<real> *
  //   3) CNarrowBandEno2EngquistOsher<real> *
  //   4) CFastMarchingFirstOrderEngquistOsher<real> **
  //         * Narrow band level set method
  //        ** Fast marching method
  typedef CNarrowBandFirstOrderEngquistOsher<real> UpdaterType;

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
  //   1) CNeverSave<real>
  //        | Nothing is saved.
  //   2) CCurvesSaver<real>
  //        | The front is constructed and saved at each time step.
  //        | Not relevant for the fast marching method.
  //   3) CSaveLastCurve<real>
  //        | Saves the last curve, after all calculations.
  //        | Not relevant for the fast marching method.
  //   4) CSaveAtTheEnd<real>
  //        | Saves the level set, after all calculations.
  //        | Relevant for the fast marching method.
  typedef CCurvesSaver<real> SaverType;

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

  CSimulator<real, MeshType, SpeedType, InitialCurveType,
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

  // To catch exceptions.
  END;

  return 0;

}
