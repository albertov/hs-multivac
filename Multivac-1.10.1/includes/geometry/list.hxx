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


#ifndef FILE_LIST_HXX


#include "errors.cxx"

#include <iostream>
using std::cout;
using std::endl;
#include <cstdlib>
#include <ctime>
#include <cstdio>


namespace Multivac
{


  template <class T>
  class List;


  //////////
  // CELL //
  //////////

  //! Cell of a double linked list.
  template <class T>
  class Cell
  {
    
    /************************
     * TYPEDEF DECLARATIONS *
     ************************/
    
  public:
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;
    
    
    /**************
     * ATTRIBUTES *
     **************/
    
  protected:
    //! Element.
    T X_;
    //! Pointer to the previous element.
    Cell<T>* previous_;
    //! Pointer to the next element.
    Cell<T>* next_;


    /**********
     * METHOD *
     **********/

  public:
    //Constructor.
    Cell()  throw();
    Cell(T& X, Cell<T>* previous, Cell<T>* next)  throw();

    // Destructor.
    ~Cell()  throw();

    // Basic methods.
    const_reference GetElement() const;
    reference GetElement();
    void GetElement(T& X) const;
    void SetElement(T& X);

    Cell<T>* GetPrevious() const;
    void SetPrevious(Cell<T>* previous);

    Cell<T>* GetNext() const;
    void SetNext(Cell<T>* next);


    /**********
     * FRIEND *
     **********/

    friend class List<T>;

  };  // Cell.

  
  //////////
  // LIST //
  //////////

  //! Double linked list.
  template <class T>
  class List
  {

    /************************
     * TYPEDEF DECLARATIONS *
     ************************/

  public:
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;


    /**************
     * ATTRIBUTES *
     **************/

  protected:
    //! Pointer to the head.
    Cell<T>* head_;
    //! Pointer to the current element ("cursor").
    Cell<T>* current_;
    //! Pointer to the tail.
    Cell<T>* tail_;


    /***********
     * METHODS *
     ***********/

  public:
    //Constructor.
    List()  throw();
    List(const List<T>&)  throw();

    // Destructor.
    ~List()  throw();

    // Initializations.
    void Init();
    void Copy(const List<T>&);

    void AddAtTheEnd(T& X);
    void AddAtTheEnd(Cell<T>* NewTail);

    void Reverse();

    Cell<T>* RemoveCurrent();
    void DeleteCurrent();
    Cell<T>* Remove(Cell<T>* cell);
    void ClearAll();

    // Convenient functions.
    const_reference GetHeadValue() const;
    const_reference GetCurrentValue() const;
    const_reference GetTailValue() const;

    reference GetHeadValue();
    reference GetCurrentValue();
    reference GetTailValue();

    Cell<T>* GetCurrent() const;

    void GoToTheHead();
    void GoToTheTail();

    bool GoToNext_StopAtTheTail();

    bool IsEmpty() const;

  };  // List.


}  // namespace Multivac.


#define FILE_LIST_HXX
#endif
