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


#ifndef FILE_SIMULATOR_HXX


#include "errors.cxx"
#include <cstdio>


namespace Multivac
{


  template <class T, class MeshType, class SpeedType, class InitialCurveType,
	    class LevelSetType, class InitializerType, class UpdaterType,
	    class SaverType>
  class CSimulator
  {

    /**************
     * ATTRIBUTES *
     **************/

  protected:

    // Values.
    MeshType Mesh;
    SpeedType F;
    InitialCurveType InitialCurve;
    LevelSetType Phi;
    InitializerType Initializer;
    UpdaterType Updater;
    SaverType Saver;

    int NbIterations;
    T FinalTime;
    T Delta_t;


    /*****************************
     * CONSTRUCTORS & DESTRUCTOR *
     *****************************/

  public:

    CSimulator(MeshType& Mesh_, SpeedType& F_,
	       InitialCurveType& InitialCurve_, LevelSetType& Phi_,
	       InitializerType& Initializer_, UpdaterType& Updater_,
	       SaverType& Saver_,
	       int NbIterations_, T FinalTime_)  throw();

    ~CSimulator()  throw();


    /***********
     * METHODS *
     ***********/

  public:
  
    void Init();

    void Run();

    T GetDelta_t();
    InitializerType& GetInitializer();
    SpeedType& GetSpeedFunction();
    SaverType& GetSaver();

  };  // CSimulator.


}  //  namespace Multivac.


#define FILE_SIMULATOR_HXX
#endif
