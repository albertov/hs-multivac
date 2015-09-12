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


#ifndef FILE_LIST_CXX


#include "list.hxx"


namespace Multivac
{

  //////////
  // CELL //
  //////////


  /***************
   * CONSTRUCTOR *
   ***************/


  //! Default constructor.
  template <class T>
  inline Cell<T>::Cell()  throw()
  {

    previous_ = NULL;
    next_ = NULL;

  }


  //! Main constructor.
  /*! \param X element.
    \param previous pointer to the previous cell.
    \param next pointer to the next cell.
  */
  template <class T>
  inline Cell<T>::Cell(T& X, Cell<T>* previous, Cell<T>* next)  throw()
  {

    X_.Copy(X);
    previous_ = previous;
    next_ = next;

  }


  /**************
   * DESTRUCTOR *
   **************/

  
  //! Destructor.
  template <class T>
  inline Cell<T>::~Cell()  throw()
  {

  }

  
  /***********
   * METHODS *
   ***********/


  //! Returns the element.
  /*! \return A constant reference to the element.
   */
  template <class T>
  inline typename Cell<T>::const_reference Cell<T>::GetElement() const
  {

    return X_;

  }


  //! Returns the element.
  /*! \return A reference to the element.
   */
  template <class T>
  inline typename Cell<T>::reference Cell<T>::GetElement()
  {

    return X_;

  }


  //! Returns the element.
  /*! \param A reference to the object in which the cell element
    has to be copied.
  */
  template <class T>
  inline void Cell<T>::GetElement(T &X) const
  {

    X.Copy(X_);

  }


  //! Modifies the cell element.
  /*! \param X the new value of the cell element.
   */
  template <class T>
  inline void Cell<T>::SetElement(T& X)
  {

    X_.Copy(X);

  }


  //! Returns the pointer to the previous cell.
  /*! \return the pointer to the previous cell.
   */
  template <class T>
  inline Cell<T>* Cell<T>::GetPrevious() const
  {

    return previous_;

  }


  //! Sets the pointer to the previous cell.
  /*! \param previous the pointer to the previous cell.
   */
  template <class T>
  inline void Cell<T>::SetPrevious(Cell<T>* previous)
  {

    previous_ = previous;

  }


  //! Returns the pointer to the next cell.
  /*! \return the pointer to the next cell.
   */
  template <class T>
  inline Cell<T>* Cell<T>::GetNext() const
  {
    return next_;
  }


  //! Sets the pointer to the next cell.
  /*! \param previous the pointer to the next cell.
   */
  template <class T>
  inline void Cell<T>::SetNext(Cell<T>* next)
  {
    next_ = next;
  }



  //////////
  // LIST //
  //////////


  /***************
   * CONSTRUCTOR *
   ***************/

  
  //! Default constructor.
  template <class T>
  List<T>::List()  throw()
  {

    head_ = NULL;
    current_ = NULL;
    tail_ = NULL;

  }


  //! Copy constructor.
  template <class T>
  List<T>::List(const List<T>& List_)  throw()
  {

    Cell<T>* begin = List_.GetCurrent();
    Cell<T>* next(begin);

    head_ = NULL;
    current_ = NULL;
    tail_ = NULL;

    if (!List_.IsEmpty())
      {

	T X;

	X.Copy(begin->GetElement());
	this->AddAtTheEnd(X);

	while (begin!=(next=next->GetNext()))
	  {
	    X.Copy(next->GetElement());
	    this->AddAtTheEnd(X);
	  }

      }

  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T>
  List<T>::~List()  throw()
  {

    if (head_!=NULL)
      {

	while (head_ != head_->next_)
	  {
	    current_ = head_->next_;
	    head_->next_ = current_->next_;
	    delete current_;
	  }
	delete head_;
      }

  }

  
  /***********
   * METHODS *
   ***********/


  //! Inits the list.
  /*! Pointers to the head, to the tail and to the current
    element are set to NULL.
  */
  template <class T>
  inline void List<T>::Init()
  {

    head_ = NULL;
    current_ = NULL;
    tail_ = NULL;

  }


  //! Copies a list.
  /*!
    \param List_ the list to be copied.
  */
  template <class T>
  void List<T>::Copy(const List<T>& List_)
  {

    this->ClearAll();

    Cell<T>* begin = List_.GetCurrent();
    Cell<T>* next(begin);

    head_ = NULL;
    current_ = NULL;
    tail_ = NULL;

    if (!List_.IsEmpty())
      {

	T X;

	X.Copy(begin->GetElement());
	this->AddAtTheEnd(X);

	while (begin!=(next=next->GetNext()))
	  {
	    X.Copy(next->GetElement());
	    this->AddAtTheEnd(X);
	  }

      }

  }


  //! Adds an element at the tail.
  /*! \param X the element to be added.
   */
  template <class T>
  inline void List<T>::AddAtTheEnd(T& X)
  {

    if (head_!=NULL)
      {
	Cell<T>* temp = new Cell<T>(X, tail_, head_);
	head_->previous_ = temp;
	tail_->next_ = temp;
	tail_ = temp;
      }
    else
      {
	Cell<T>* temp = new Cell<T>(X, NULL, NULL);
	temp->previous_ = temp;
	temp->next_ = temp;
	head_ = temp;
	current_ = temp;
	tail_ = temp;
      }

  }


  //! Adds a cell at the tail.
  /*! \param NewTail a pointer to the cell to be added.
   */
  template <class T>
  inline void List<T>::AddAtTheEnd(Cell<T>* NewTail)
  {

    if (head_!=NULL)
      {
	NewTail->previous_ = tail_;
	NewTail->next_ = head_;
	head_->previous_ = NewTail;
	tail_->next_ = NewTail;
	tail_ = NewTail;
      }
    else
      {
	NewTail->previous_ = NewTail;
	NewTail->next_ = NewTail;
	head_ = NewTail;
	current_ = NewTail;
	tail_ = NewTail;
      }

  }


  //! Removes the current cell (cursor) from the list.
  /* \return the pointer to the removed cell.
     \note The new current cell is the next cell.
  */
  template <class T>
  inline void List<T>::Reverse()
  {

    if (tail_ == head_)
      return;

    Cell<T>* temp;

    current_ = head_;
    temp = current_->next_;
    current_->next_ = current_->previous_;
    current_->previous_ = temp;
    current_ = temp;
    
    while (current_ != head_)
      {
	temp = current_->next_;
	current_->next_ = current_->previous_;
	current_->previous_ = temp;
	current_ = temp;
      }

    temp = head_;
    head_ = tail_;
    tail_ = temp;

  }


  //! Removes the current cell (cursor) from the list.
  /* \return the pointer to the removed cell.
     \note The new current cell is the next cell.
  */
  template <class T>
  inline Cell<T>* List<T>::RemoveCurrent()
  {

    if (current_ == NULL)
      return NULL;

    Cell<T>* temp = current_;

    if (temp->next_ == temp)
      {
	head_ = NULL;
	current_ = NULL;
	tail_ = NULL;
      }
    else
      {
	(temp->previous_)->next_ = temp->next_;
	(temp->next_)->previous_ = temp->previous_;
	current_ = temp->next_;

	if (temp == tail_)
	  tail_ = temp->previous_;
	else if (temp == head_)
	  head_ = temp->next_;
      }

    temp->previous_ = NULL;
    temp->next_ = NULL;

    return temp;

  }


  //! Removes the current cell (cursor) from the list and releases
  //! the cell from memory.
  /* \note The new current cell is the next cell.
   */
  template <class T>
  inline void List<T>::DeleteCurrent()
  {

    if (current_ == NULL)
      return;

    Cell<T>* temp = current_;

    if (temp->next_ == temp)
      {
	head_ = NULL;
	current_ = NULL;
	tail_ = NULL;
      }
    else
      {
	(temp->previous_)->next_ = temp->next_;
	(temp->next_)->previous_ = temp->previous_;
	current_ = temp->next_;

	if (temp == tail_)
	  tail_ = temp->previous_;
	else if (temp == head_)
	  head_ = temp->next_;
      }

    temp->previous_ = NULL;
    temp->next_ = NULL;

    delete temp;

  }


  //! Removes a given cell (cursor) from the list.
  /* \param cell a pointer to the cell to be removed.
     \return the pointer to the removed cell.
     \note The new current cell is the next cell.
  */
  template <class T>
  inline Cell<T>* List<T>::Remove(Cell<T>* cell)
  {

    if (cell->next_ == cell)
      {
	head_ = NULL;
	current_ = NULL;
	tail_ = NULL;
      }
    else
      {
	(cell->previous_)->next_ = cell->next_;
	(cell->next_)->previous_ = cell->previous_;

	if (cell == current_)
	  current_ = cell->next_;
	if (cell == tail_)
	  tail_ = cell->previous_;
	else if (cell == head_)
	  head_ = cell->next_;
      }

    cell->previous_ = NULL;
    cell->next_ = NULL;

    return cell;

  }


  //! Clears the list.
  /*! All cells are released from memory.
   */
  template <class T>
  inline void List<T>::ClearAll()
  {

    if (head_!=NULL)
      {

	while (head_ != head_->next_)
	  {
	    current_ = head_->next_;
	    head_->next_ = current_->next_;
	    delete current_;
	  }
	delete head_;
      }

    head_ = NULL;
    current_ = NULL;
    tail_ = NULL;

  }

 
  /************************
   * CONVENIENT FUNCTIONS *
   ************************/


  //! Returns the element stored at the head.
  /*! \return a constant reference to the element stored by
    the cell at the head.
  */
  template <class T>
  inline typename List<T>::const_reference List<T>::GetHeadValue() const
  {

    return head_->GetElement();

  }


  //! Returns the element stored where the cursor points.
  /*! \return a constant reference to the element stored by
    the cell to which the cursor points.
  */
  template <class T>
  inline typename List<T>::const_reference List<T>::GetCurrentValue() const
  {

    return current_->GetElement();

  }


  //! Returns the element stored at the tail.
  /*! \return a constant reference to the element stored by
    the cell at the tail.
  */
  template <class T>
  inline typename List<T>::const_reference List<T>::GetTailValue() const
  {

    return tail_->GetElement();

  }


  //! Returns the element stored at the head.
  /*! \return a reference to the element stored by
    the cell at the head.
  */
  template <class T>
  inline typename List<T>::reference List<T>::GetHeadValue()
  {

    return head_->GetElement();

  }


  //! Returns the element stored where the cursor points.
  /*! \return a reference to the element stored by
    the cell to which the cursor points.
  */
  template <class T>
  inline typename List<T>::reference List<T>::GetCurrentValue()
  {

    return current_->GetElement();

  }


  //! Returns the element stored at the tail.
  /*! \return a reference to the element stored by
    the cell at the tail.
  */
  template <class T>
  inline typename List<T>::reference List<T>::GetTailValue()
  {

    return tail_->GetElement();

  }


  //! Returns the pointer to the current cell (cursor).
  /*! \return the pointer to the current cell (cursor).
   */
  template <class T>
  inline Cell<T>* List<T>::GetCurrent() const
  {

    return current_;

  }


  //! Sets the current cell (cursor) to the pointer to the head.
  template <class T>
  inline void List<T>::GoToTheHead()
  {

    current_ = head_;

  }


  //! Sets the current cell (cursor) to the pointer to the tail.
  template <class T>
  inline void List<T>::GoToTheTail()
  {

    current_ = tail_;

  }


  //! Moves the current cell (cursor) to the next cell
  //! but if the tail is reached.
  /*! \return true if the new current cell has changed and
    false if the tail is the current cell.
  */
  template <class T>
  inline bool List<T>::GoToNext_StopAtTheTail()
  {

    if (current_ != tail_)
      {
	current_ = current_->next_;
	return true;
      }
    else
      return false;

  }


  //! Is the list empty?
  /*! \return true if the list is empty and false otherwise.
   */
  template <class T>
  inline bool List<T>::IsEmpty() const
  {

    return (head_==NULL);

  }


}  // Multivac.


#define FILE_LIST_CXX
#endif
