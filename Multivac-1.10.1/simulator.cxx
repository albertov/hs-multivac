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


#ifndef FILE_SIMULATOR_CXX


#include "simulator.hxx"


namespace Multivac
{


  ///////////////////////////////
  // CONSTRUCTORS & DESTRUCTOR //
  ///////////////////////////////


  template <class T, class MeshType, class SpeedType, class InitialCurveType,
	    class LevelSetType, class InitializerType, class UpdaterType,
	    class SaverType>
  CSimulator<T, MeshType, SpeedType, InitialCurveType, LevelSetType,
	     InitializerType, UpdaterType, SaverType>::
  CSimulator(MeshType& Mesh_, SpeedType& F_,
	     InitialCurveType& InitialCurve_, LevelSetType& Phi_,
	     InitializerType& Initializer_, UpdaterType& Updater_,
	     SaverType& Saver_,
	     int NbIterations_, T FinalTime_)  throw():
    Mesh(Mesh_), F(F_), InitialCurve(InitialCurve_), Phi(Phi_),
    Initializer(Initializer_), Updater(Updater_), Saver(Saver_),
    NbIterations(NbIterations_), FinalTime(FinalTime_),
    Delta_t(FinalTime_ / T(NbIterations_))
  {
    
  }
  
  
  template <class T, class MeshType, class SpeedType, class InitialCurveType,
	    class LevelSetType, class InitializerType, class UpdaterType,
	    class SaverType>
  CSimulator<T, MeshType, SpeedType, InitialCurveType, LevelSetType,
	     InitializerType, UpdaterType, SaverType>::
  ~CSimulator()  throw()
  {

  }



  /////////////
  // METHODS //
  /////////////


  template <class T, class MeshType, class SpeedType,
	    class InitialCurveType, class LevelSetType,
	    class InitializerType, class UpdaterType, class SaverType>
  void CSimulator<T, MeshType, SpeedType, InitialCurveType, LevelSetType,
		  InitializerType, UpdaterType, SaverType>::Init()
  {

#ifdef MULTIVAC_CHECK_COMPATIBILITY
    // Checks the compatibility between the initializer and the updater.
    if (Initializer.IsNarrowBand())
      {
	if (!Updater.IsNarrowBand())
	  throw CError_Incompatibility("CSimulator::Init()",
				       string("The initializer is dedicated ")
				       +string("to the narrow band level set")
				       + " method while the updater is not");
      }
    else if (Initializer.IsFastMarching())
      {
	if (!Updater.IsFastMarching())
	  throw CError_Incompatibility("CSimulator::Init()",
				       string("The initializer is dedicated ")
				       + string("to the fast marching method")
				       + " while the updater is not");
      }
    else
      {
	throw CError_Incompatibility("CSimulator::Init()",
				     string("The initializer is neither ")
				     + string("dedicated to the fast ")
				     + string("marching method nor to the ")
				     + "narrow band level set method");
      }
#endif
    
    Initializer.FirstInitMesh(Mesh);
    Initializer.FirstInitInitialCurve(Mesh, InitialCurve);
    Initializer.FirstInitPhiAndF(Mesh, InitialCurve, Phi, F, Updater);

  }


  template <class T, class MeshType, class SpeedType,
	    class InitialCurveType, class LevelSetType,
	    class InitializerType, class UpdaterType, class SaverType>
  void CSimulator<T, MeshType, SpeedType, InitialCurveType, LevelSetType,
		  InitializerType, UpdaterType, SaverType>::Run()
  {

    if (Updater.IsNarrowBand())
      {

	/*********************
	 * LEVEL SET METHODS *
	 *********************/

	Saver.SaveAtTheBeginning(Mesh, F, Phi, Initializer);

	int iter = 0;

#ifdef MULTIVAC_REPORT
	clock_t time = 0;
#endif

	Vector<T> Time(NbIterations + 1);
	Time(0) = 0.0;

	T t = 0;  int i = 0;
	while (iter < Time.GetLength()-1)
	  {

	    iter++;
	    t = t + Delta_t;  i++;

#ifdef MULTIVAC_REPORT
	    cerr << "Iteration: " << iter << endl;
	    cerr << "----------" << endl;
	    time = clock();
#endif

	    Updater.UpdateLevelSet(Delta_t, Mesh, F, Phi, t);

#ifdef MULTIVAC_REPORT
	    cerr << "Update: " << clock() - time << endl;
	    time = clock();
#endif

	    Initializer.InitMesh(iter, Mesh, Phi, F, Updater, t);

#ifdef MULTIVAC_REPORT
	    cerr << "InitMesh: " << clock() - time << endl;
	    time = clock();
#endif

	    Initializer.InitPhiAndF(iter, Mesh, Phi, F, Updater, t);

#ifdef MULTIVAC_REPORT
	    cerr << "InitPhiAndF: " << clock() - time << endl;
	    time = clock();
#endif

	    Saver.SaveAtCurrentIteration(Mesh, F, Phi, t, iter, Initializer);

#ifdef MULTIVAC_REPORT
	    cerr << "SaveAtCurrentIteration: " << clock() - time
		 << endl << endl;
#endif

	    Time(iter) = t;

	  }

	Saver.SaveAtTheEnd(Mesh, F, Phi, Time, i, Initializer);

      }
    else if (Updater.IsFastMarching())
      {

	/*****************
	 * FAST MARCHING *
	 *****************/

	Saver.SaveAtTheBeginning(Mesh, F, Phi, Initializer);

	int iter = 0;

#ifdef MULTIVAC_REPORT
	clock_t time = 0;
#endif

	while (Updater.KeepOnWorking())
	  {

	    iter++;

#ifdef MULTIVAC_REPORT
	    cerr << "Iteration: " << iter << endl;
	    cerr << "----------" << endl;
	    time = clock();
#endif

	    Updater.UpdateLevelSet(0.0, Mesh, F, Phi, 0.0);

#ifdef MULTIVAC_REPORT
	    cerr << "Update: " << clock() - time << endl;
	    time = clock();
#endif

	    Initializer.InitMesh(iter, Mesh, Phi, F, Updater, 0.0);

#ifdef MULTIVAC_REPORT
	    cerr << "InitMesh: " << clock() - time << endl;
	    time = clock();
#endif

	    Initializer.InitPhiAndF(iter, Mesh, Phi, F, Updater, 0.0);

#ifdef MULTIVAC_REPORT
	    cerr << "InitPhiAndF: " << clock() - time << endl;
	    time = clock();
#endif

	    Saver.SaveAtCurrentIteration(Mesh, F, Phi, 0.0,
					 iter, Initializer);

#ifdef MULTIVAC_REPORT
	    cerr << "SaveAtCurrentIteration: " << clock() - time
		 << endl << endl;
#endif

	  }


	Vector<T> Time(1);
	Time(0) = 0.0;

	Saver.SaveAtTheEnd(Mesh, F, Phi, Time, iter, Initializer);

      }

  }


  template <class T, class MeshType, class SpeedType,
            class InitialCurveType, class LevelSetType,
            class InitializerType, class UpdaterType, class SaverType>
  T CSimulator<T, MeshType, SpeedType, InitialCurveType,
               LevelSetType, InitializerType, UpdaterType, SaverType>
  ::GetDelta_t()
  {
    return Delta_t;
  }


  template <class T, class MeshType, class SpeedType,
	    class InitialCurveType, class LevelSetType,
	    class InitializerType, class UpdaterType, class SaverType>
  InitializerType& CSimulator<T, MeshType, SpeedType, InitialCurveType,
			      LevelSetType, InitializerType, UpdaterType,
			      SaverType>::GetInitializer()
  {
    return Initializer;
  }


  template <class T, class MeshType, class SpeedType,
	    class InitialCurveType, class LevelSetType,
	    class InitializerType, class UpdaterType, class SaverType>
  SpeedType& CSimulator<T, MeshType, SpeedType, InitialCurveType,
			LevelSetType, InitializerType,
			UpdaterType, SaverType>::GetSpeedFunction()
  {
    return F;
  }


  template <class T, class MeshType, class SpeedType,
            class InitialCurveType, class LevelSetType,
            class InitializerType, class UpdaterType, class SaverType>
  SaverType& CSimulator<T, MeshType, SpeedType, InitialCurveType,
                        LevelSetType, InitializerType, UpdaterType,
                        SaverType>::GetSaver()
  {
    return Saver;
  }


}  // namespace Multivac.


#define FILE_SIMULATOR_CXX
#endif
