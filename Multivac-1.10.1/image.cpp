/************************
 * INCLUDES AND OPTIONS *
 ************************/

// Seldon library and Multivac project provide error management
// using exception handling.  Exceptions that may be raised
// are selected according to debug levels.
#define MULTIVAC_DEBUG_LEVEL_4

// Define MULTIVAC_REPORT if you want to clock the time needed for
// updates and initializations (results will be displayed on screen).
#define MULTIVAC_REPORT_NO

#include "multivac.hxx"

using namespace Multivac;


/*****************
 * MAIN FUNCTION *
 *****************/

int main(int argc, char** argv)
{

  // To catch exceptions.
  TRY;
  
  if (argc != 2)
    {
      string mesg  = "Usage:\n";
      mesg += string("  ") + argv[0] + " [configuration file]";
      cout << mesg << endl;
      return 1;
    }

  // real type: double, float...
  typedef double real;


  ////////////////////////
  // CONFIGURATION FILE //
  ////////////////////////

  ConfigStream config(argv[1]);


  //////////
  // TIME //
  //////////
  
  config.SetSection("[domain]");

  // Final time of the simulation, the initial time being 0.
  real FinalTime;
  // Time step.
  real Delta_t;

  config.PeekValue("Final_time", FinalTime);
  config.PeekValue("Delta_t", Delta_t);

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
  real Xmin, Xmax, Ymin, Ymax;
  real Delta_x, Delta_y;

  // For an orthogonal mesh, Nx and Ny are the number of
  // grid points along (x'x) and (y'y) (respectively).
  int Nx, Ny;

  config.PeekValue("Xmin", Xmin);
  config.PeekValue("Xmax", Xmax);
  config.PeekValue("Ymin", Ymin);
  config.PeekValue("Ymax", Ymax);

  config.PeekValue("Nx", Nx);
  config.PeekValue("Ny", Ny);

  Delta_x = (Xmax - Xmin) / real(Nx);
  Delta_y = (Ymax - Ymin) / real(Ny);

  // Choose the right constructor:
  //   1) Mesh(Xmin, Xmax, Ymin, Ymax, Nx, Ny)
  //   2) --
  MeshType Mesh(Xmin, Xmax, Ymin, Ymax, Nx, Ny);


  ///////////////////
  // INITIAL CURVE //
  ///////////////////

  config.SetSection("[initial_curve]");

  // Choose the initial curve:
  //   1) CCircle<real>
  //   2) CTwoCircles<real>
  //   3) CThreeCircles<real>
  //   4) CIsland<real>
  //   5) CIsland0<real>
  //   6) CSetOfPoints<real>
  typedef CSetOfPoints<real> InitialCurveType;

  // File containing a front defined by a set of points.
  string InitialFrontFile = "[full path]/[set].pts";
  // Is the (Xmin, Ymin) outside the front?
  bool origin_out = false;

  // If set to 'true', outside and inside are swapped.
  bool reversed = false;
  
  config.PeekValue("Initial_front_file", InitialFrontFile);

  config.PeekValue("Origin_out", origin_out);
  config.PeekValue("Reversed", reversed);

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
  InitialCurveType InitialCurve(InitialFrontFile, 0, origin_out ^ reversed);


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
  typedef CNarrowBandExtension<real> InitializerType;

  InitializerType Initializer;


  ///////////
  // SAVER //
  ///////////

  config.SetSection("[save]");

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
  string Directory = "/home/vivien/fronts/multivac/current/results/";

  config.PeekValue("Num_curves", NbCurves);
  config.PeekValue("Output_directory", Directory);

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

  config.SetSection("[speed_function]");

  // Segmentation method.
  string method;

  config.PeekValue("Method", method);

  if (lower_case(method) == "gradient")
    {


      ////////////////////
      // SPEED FUNCTION //
      ////////////////////

      config.SetSection("[speed_function]");

      bool raw_image;
      string image_filename;
      real epsilon_c;
      real boost;
      real slow_down;

      config.PeekValue("Raw_image", raw_image);
      if (raw_image)
	throw "Multivac cannot handle a raw image yet.";
      config.PeekValue("Image_gradient", image_filename);
      config.PeekValue("Epsilon_c", epsilon_c);
      config.PeekValue("Boost", boost);
      config.PeekValue("Slow_down", slow_down);

      typedef CGradientSegmentation<real> SpeedType;
      SpeedType F(Xmin, Delta_x, Nx, Ymin, Delta_y, Ny,
		  image_filename, epsilon_c, boost, slow_down);


      /////////////
      // UPDATER //
      /////////////

      config.SetSection("[level_set]");

      int TubeSemiWidth, BarrierWidth, OutSpaceWidth;

      config.PeekValue("Tube_semi_width", TubeSemiWidth);
      config.PeekValue("Barrier_width", BarrierWidth);
      config.PeekValue("Out_space_width", OutSpaceWidth);

      typedef CNarrowBandFirstOrderEngquistOsher<real> UpdaterType;
      UpdaterType Updater(TubeSemiWidth, BarrierWidth, OutSpaceWidth);


      ////////////////
      // SIMULATION //
      ////////////////

      CSimulator<real, MeshType, SpeedType, InitialCurveType,
	LevelSetType, InitializerType, UpdaterType, SaverType>
	Simulator(Mesh, F, InitialCurve, Phi,
		  Initializer, Updater, Saver,
		  NbIterations, FinalTime);

      Simulator.Init();
      Simulator.Run();

    }
  else if (lower_case(method) == "chan-vese")
    {
      

      ////////////////////
      // SPEED FUNCTION //
      ////////////////////

      config.SetSection("[speed_function]");

      bool raw_image;
      string image_filename;
      real intensity_threshold;

      config.PeekValue("Raw_image", raw_image);
      if (raw_image)
	throw "Multivac cannot handle a raw image yet.";
      config.PeekValue("Image_intensity", image_filename);
      config.PeekValue("Intensity_threshold", intensity_threshold);

      typedef CImageIntensity<real> SpeedType;
      SpeedType F(Xmin, Delta_x, Nx, Ymin, Delta_y, Ny,
		  image_filename, intensity_threshold);


      /////////////
      // UPDATER //
      /////////////

      config.SetSection("[level_set]");

      int TubeSemiWidth, BarrierWidth, OutSpaceWidth;

      config.PeekValue("Tube_semi_width", TubeSemiWidth);
      config.PeekValue("Barrier_width", BarrierWidth);
      config.PeekValue("Out_space_width", OutSpaceWidth);

      config.SetSection("[speed_function]");

      real inside, outside, mu, nu, Dirac_threshold;

      config.PeekValue("Inside", inside);
      config.PeekValue("Outside", outside);
      config.PeekValue("Mu", mu);
      config.PeekValue("Nu", nu);
      config.PeekValue("Dirac_threshold", Dirac_threshold);

      typedef CChanVese<real> UpdaterType;
      UpdaterType Updater(TubeSemiWidth, BarrierWidth, OutSpaceWidth,
			  inside, outside, mu, nu, Dirac_threshold);


      ////////////////
      // SIMULATION //
      ////////////////

      CSimulator<real, MeshType, SpeedType, InitialCurveType,
	LevelSetType, InitializerType, UpdaterType, SaverType>
	Simulator(Mesh, F, InitialCurve, Phi,
		  Initializer, Updater, Saver,
		  NbIterations, FinalTime);

      Simulator.Init();
      Simulator.Run();

    }
  else
    throw string("Segmentation method \"") + method
      + string("\" is not supported.");

  // To catch exceptions.
  END;

  return 0;

}
