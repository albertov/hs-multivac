// Copyright (C) 2001-2004 Vivien Mallet
//
// This file is part of Seldon library.
// Seldon library provides matrices and vectors structures for
// linear algebra.
// 
// Seldon is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// Seldon is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License (file "license") for more details.
//
// For more information, please see the Seldon home page:
//     http://spacetown.free.fr/lib/seldon/

#ifndef SELDON_FILE_ERRORS_HXX

namespace Seldon
{

  
  ///////////
  // ERROR //
  ///////////

  class Error
  {
  protected:
    string function;
    string comment;

  public:
    Error()  throw();
    Error(string f)  throw();
    Error(string f, string c)  throw();
    virtual ~Error()  throw();

    virtual void What();

  };

  
  //////////////
  // NOMEMORY //
  //////////////

  class NoMemory: public Error
  {
  public:
    NoMemory(string f)  throw();
    NoMemory(string f, string c)  throw();

    virtual void What();
  };
  

  //////////////
  // WRONGDIM //
  //////////////

  class WrongDim: public Error
  {
  public:
    WrongDim(string f)  throw();
    WrongDim(string f, string c)  throw();

    virtual void What();
  };
  

  ////////////////
  // WRONGINDEX //
  ////////////////

  class WrongIndex: public Error
  {
  public:
    WrongIndex(string f)  throw();
    WrongIndex(string f, string c)  throw();

    virtual void What();
  };
  

  //////////////
  // WRONGROW //
  //////////////

  class WrongRow: public Error
  {
  public:
    WrongRow(string f)  throw();
    WrongRow(string f, string c)  throw();

    virtual void What();
  };
  

  //////////////
  // WRONGCOL //
  //////////////

  class WrongCol: public Error
  {
  public:
    WrongCol(string f)  throw();
    WrongCol(string f, string c)  throw();

    virtual void What();
  };
  

  /////////////
  // IOERROR //
  /////////////

  class IOError: public Error
  {
  public:
    IOError(string f)  throw();
    IOError(string f, string c)  throw();

    virtual void What();
  };
  

} // namespace Seldon.

#define SELDON_FILE_ERRORS_HXX
#endif
