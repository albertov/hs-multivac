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

#ifndef SELDON_FILE_ERRORS_CXX

#include "Errors.hxx"

namespace Seldon
{


  ///////////
  // ERROR //
  ///////////


  /****************
   * CONSTRUCTORS *
   ****************/
  

  //! Default constructor.
  /*!
    Error asociated neither with a function nor with a comment.
  */
  Error::Error()  throw()
  {
    cerr << "ERROR!" << endl;
    function = "";
    comment = "";
  }


  //! Constructor.
  /*! Error associated with a function.
    \param f function with which the error is associated.
  */
  Error::Error(string f)  throw()
  {
    cerr << "ERROR!" << endl;
    function = f;
    comment = "";
  }


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param f function with which the error is associated.
    \param c comment associated with the error.
  */
  Error::Error(string f, string c)  throw()
  {
    cerr << "ERROR!" << endl;
    function = f;
    comment = c;
  }


  /**************
   * DESTRUCTOR *
   **************/
  

  //! Destructor.
  /*!
    \note Empty.
  */
  Error::~Error()  throw()
  {
  }


  /***********
   * METHODS *
   ***********/
  

  //! Delivers information about the error.
  /*! Displays available information, i.e.
    the error description, the function and/or the comment.
  */
  void Error::What()
  {
    cerr << "An undefined error occured";
    if (function.size() != 0)
      cerr << " in " << function;
    cerr << "." << endl;
    if (comment.size() != 0)
      cerr << "   " << comment << endl;
    cerr << endl;
  }



  //////////////
  // NOMEMORY //
  //////////////

  
  /****************
   * CONSTRUCTORS *
   ****************/


  //! Constructor.
  /*! Error associated with a function.
    \param f function with which the error is associated.
  */
  NoMemory::NoMemory(string f)  throw(): Error(f)
  {
  }


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param f function with which the error is associated.
    \param c comment associated with the error.
  */
  NoMemory::NoMemory(string f, string c)  throw(): Error(f, c)
  {
  }


  /***********
   * METHODS *
   ***********/


  //! Delivers information about the error.
  /*! Displays available information, i.e.
    the error description, the function and/or the comment.
  */
  void NoMemory::What()
  {
    cerr << "Out of memory";
    if (this->function != "")
      cerr << " in " << this->function;
    cerr << "." << endl;
    if (this->comment != "")
      cerr << "   " << this->comment << endl;
    cerr << endl;
  }



  //////////////
  // WRONGDIM //
  //////////////

  
  /****************
   * CONSTRUCTORS *
   ****************/


  //! Constructor.
  /*! Error associated with a function.
    \param f function with which the error is associated.
  */
  WrongDim::WrongDim(string f)  throw(): Error(f)
  {
  }


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param f function with which the error is associated.
    \param c comment associated with the error.
  */
  WrongDim::WrongDim(string f, string c)  throw(): Error(f, c)
  {
  }


  /***********
   * METHODS *
   ***********/


  //! Delivers information about the error.
  /*! Displays available information, i.e.
    the error description, the function and/or the comment.
  */
  void WrongDim::What()
  {
    cerr << "Wrong dimensions involved";
    if (this->function != "")
      cerr << " in " << this->function;
    cerr << "." << endl;
    if (this->comment != "")
      cerr << "   " << this->comment << endl;
    cerr << endl;
  }



  ////////////////
  // WRONGINDEX //
  ////////////////

  
  /****************
   * CONSTRUCTORS *
   ****************/


  //! Constructor.
  /*! Error associated with a function.
    \param f function with which the error is associated.
  */
  WrongIndex::WrongIndex(string f)  throw(): Error(f)
  {
  }


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param f function with which the error is associated.
    \param c comment associated with the error.
  */
  WrongIndex::WrongIndex(string f, string c)  throw(): Error(f, c)
  {
  }


  /***********
   * METHODS *
   ***********/


  //! Delivers information about the error.
  /*! Displays available information, i.e.
    the error description, the function and/or the comment.
  */
  void WrongIndex::What()
  {
    cerr << "Index out of range";
    if (this->function != "")
      cerr << " in " << this->function;
    cerr << "." << endl;
    if (this->comment != "")
      cerr << "   " << this->comment << endl;
    cerr << endl;
  }



  //////////////
  // WRONGROW //
  //////////////

  
  /****************
   * CONSTRUCTORS *
   ****************/


  //! Constructor.
  /*! Error associated with a function.
    \param f function with which the error is associated.
  */
  WrongRow::WrongRow(string f)  throw(): Error(f)
  {
  }


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param f function with which the error is associated.
    \param c comment associated with the error.
  */
  WrongRow::WrongRow(string f, string c)  throw(): Error(f, c)
  {
  }


  /***********
   * METHODS *
   ***********/


  //! Delivers information about the error.
  /*! Displays available information, i.e.
    the error description, the function and/or the comment.
  */
  void WrongRow::What()
  {
    cerr << "Row index out of range";
    if (this->function != "")
      cerr << " in " << this->function;
    cerr << "." << endl;
    if (this->comment != "")
      cerr << "   " << this->comment << endl;
    cerr << endl;
  }



  //////////////
  // WRONGCOL //
  //////////////

  
  /****************
   * CONSTRUCTORS *
   ****************/


  //! Constructor.
  /*! Error associated with a function.
    \param f function with which the error is associated.
  */
  WrongCol::WrongCol(string f)  throw(): Error(f)
  {
  }


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param f function with which the error is associated.
    \param c comment associated with the error.
  */
  WrongCol::WrongCol(string f, string c)  throw(): Error(f, c)
  {
  }


  /***********
   * METHODS *
   ***********/


  //! Delivers information about the error.
  /*! Displays available information, i.e.
    the error description, the function and/or the comment.
  */
  void WrongCol::What()
  {
    cerr << "Column index out of range";
    if (this->function != "")
      cerr << " in " << this->function;
    cerr << "." << endl;
    if (this->comment != "")
      cerr << "   " << this->comment << endl;
    cerr << endl;
  }



  /////////////
  // IOERROR //
  /////////////

  
  /****************
   * CONSTRUCTORS *
   ****************/


  //! Constructor.
  /*! Error associated with a function.
    \param f function with which the error is associated.
  */
  IOError::IOError(string f)  throw(): Error(f)
  {
  }


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param f function with which the error is associated.
    \param c comment associated with the error.
  */
  IOError::IOError(string f, string c)  throw(): Error(f, c)
  {
  }


  /***********
   * METHODS *
   ***********/


  //! Delivers information about the error.
  /*! Displays available information, i.e.
    the error description, the function and/or the comment.
  */
  void IOError::What()
  {
    cerr << "Error while performing an I/O operation";
    if (this->function != "")
      cerr << " in " << this->function;
    cerr << "." << endl;
    if (this->comment != "")
      cerr << "   " << this->comment << endl;
    cerr << endl;
  }


} // namespace Seldon.

#define SELDON_FILE_ERRORS_CXX
#endif
