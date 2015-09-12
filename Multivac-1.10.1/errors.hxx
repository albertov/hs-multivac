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
      cout << "ERROR!" << endl;
      function = "";
      comment = "";
    }
    CError(string f)
    {
      cout << "ERROR!" << endl;
      function = f;
      comment = "";
    }
    CError(string f, string c)
    {
      cout << "ERROR!" << endl;
      function = f;
      comment = c;
    }
    virtual ~CError()
    {
    }
    virtual void What()
    {
      cout << "An undefined error occured";
      if (!function.empty())
	cout << " in " << function;
      cout << "." << endl;
      if (!comment.empty())
	cout << "   " << comment << "." << endl;
      cout << endl;
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
      cout << "Error while performing an I/O operation (with files)";
      if (!this->function.empty())
	cout << " in " << this->function;
      cout << "." << endl;
      if (!this->comment.empty())
	cout << "   " << this->comment << "." << endl;
      cout << endl;
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
      cout << "Computation outside the domain";
      if (!this->function.empty())
	cout << " in " << this->function;
      cout << "." << endl;
      if (!this->comment.empty())
	cout << "   " << this->comment << "." << endl;
      cout << endl;
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
      cout << "Computations stopped because of an instability";
      if (!this->function.empty())
	cout << " in " << this->function;
      cout << "." << endl;
      if (!this->comment.empty())
	cout << "   " << this->comment << "." << endl;
      cout << endl;
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
      cout << "An incompatibility has been detected in your choices";
      if (!this->function.empty())
	cout << " in " << this->function;
      cout << "." << endl;
      if (!this->comment.empty())
	cout << "   " << this->comment << "." << endl;
      cout << endl;
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
      cout << "Function ";
      if (!this->function.empty())
	cout << this->function;
      cout << " is not defined." << endl;
      if (!this->comment.empty())
	cout << "   " << this->comment << "." << endl;
      cout << endl;
    }
  };
  

}  // namespace Multivac.


#define FILE_MULTIVAC_ERRORS_HXX
#endif
