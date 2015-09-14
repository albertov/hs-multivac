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


#ifndef FILE_MULTIVAC_ERRORS_HXX


#include <iostream>
#include <string>


namespace Multivac
{


  ////////////
  // CERROR //
  ////////////

  //! Base class for all exceptions.
  class CError
  {
  protected:
    string function;
    string comment;

  public:
    CError()
    {
      cerr << "ERROR!" << endl;
      function = "";
      comment = "";
    }
    CError(string f)
    {
      cerr << "ERROR!" << endl;
      function = f;
      comment = "";
    }
    CError(string f, string c)
    {
      cerr << "ERROR!" << endl;
      function = f;
      comment = c;
    }
    virtual ~CError()
    {
    }
    virtual void What()
    {
      cerr << "An undefined error occured";
      if (!function.empty())
	cerr << " in " << function;
      cerr << "." << endl;
      if (!comment.empty())
	cerr << "   " << comment << "." << endl;
      cerr << endl;
    }
  };

  
  ////////////
  // FILEIO //
  ////////////

  //! Exception raised if an error occurs while performing an
  //! input/ouput operation involving files.
  class CError_FileIO: public CError
  {
  public:
    CError_FileIO(string f): CError(f)
    {
    }
    CError_FileIO(string f, string c): CError(f, c)
    {
    }
    virtual ~CError_FileIO()
    {
    }
    virtual void What()
    {
      cerr << "Error while performing an I/O operation (with files)";
      if (!this->function.empty())
	cerr << " in " << this->function;
      cerr << "." << endl;
      if (!this->comment.empty())
	cerr << "   " << this->comment << "." << endl;
      cerr << endl;
    }
  };


  /////////////////
  // OUTOFDOMAIN //
  /////////////////

  //! Exception raised when computing outside the domain.
  class CError_OutOfDomain: public CError
  {
  public:
    CError_OutOfDomain(string f): CError(f)
    {
    }
    CError_OutOfDomain(string f, string c): CError(f, c)
    {
    }
    virtual ~CError_OutOfDomain()
    {
    }
    virtual void What()
    {
      cerr << "Computation outside the domain";
      if (!this->function.empty())
	cerr << " in " << this->function;
      cerr << "." << endl;
      if (!this->comment.empty())
	cerr << "   " << this->comment << "." << endl;
      cerr << endl;
    }
  };


  /////////////////
  // INSTABILITY //
  /////////////////

  //! Exception raised when computations lead to an instability.
  class CError_Instability: public CError
  {
  public:
    CError_Instability(string f): CError(f)
    {
    }
    CError_Instability(string f, string c): CError(f, c)
    {
    }
    virtual ~CError_Instability()
    {
    }
    virtual void What()
    {
      cerr << "Computations stopped because of an instability";
      if (!this->function.empty())
	cerr << " in " << this->function;
      cerr << "." << endl;
      if (!this->comment.empty())
	cerr << "   " << this->comment << "." << endl;
      cerr << endl;
    }
  };


  /////////////////////
  // INCOMPATIBILITY //
  /////////////////////

  //! Exception raised if incompatible objects are collected
  //! to define a simulation.
  class CError_Incompatibility: public CError
  {
  public:
    CError_Incompatibility(string f): CError(f)
    {
    }
    CError_Incompatibility(string f, string c): CError(f, c)
    {
    }
    virtual ~CError_Incompatibility()
    {
    }
    virtual void What()
    {
      cerr << "An incompatibility has been detected in your choices";
      if (!this->function.empty())
	cerr << " in " << this->function;
      cerr << "." << endl;
      if (!this->comment.empty())
	cerr << "   " << this->comment << "." << endl;
      cerr << endl;
    }
  };


  ////////////////////////
  // UNDEFINED FUNCTION //
  ////////////////////////

  //! Exception raised when a function is not defined.
  class CError_Undefined: public CError
  {
  public:
    CError_Undefined(string f): CError(f)
    {
    }
    CError_Undefined(string f, string c): CError(f, c)
    {
    }
    virtual ~CError_Undefined()
    {
    }
    virtual void What()
    {
      cerr << "Function ";
      if (!this->function.empty())
	cerr << this->function;
      cerr << " is not defined." << endl;
      if (!this->comment.empty())
	cerr << "   " << this->comment << "." << endl;
      cerr << endl;
    }
  };
  

}  // namespace Multivac.


#define FILE_MULTIVAC_ERRORS_HXX
#endif
