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


#ifndef ARRAYHEAD_CXX


#include "arrayheap.hxx"


namespace Multivac
{


  ///////////////
  // ARRAYHEAP //
  ///////////////


  /***************
   * CONSTRUCTOR *
   ***************/

  
  //! Default constructor.
  template <class T>
  ArrayHeap<T>::ArrayHeap()  throw():
    NbNodes(0)
  {

  }


  //! Main constructor.
  /*! The heap is allocated to be able to store 'depth' nodes.
    \depth number of nodes.
    \note This constructor allows to avoid many reallocations of
    the heap since 'depth' nodes have been already allocated.
  */
  template <class T>
  ArrayHeap<T>::ArrayHeap(int depth)  throw():
    Nodes(depth), NbNodes(0)
  {

  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T>
  ArrayHeap<T>::~ArrayHeap()  throw()
  {
    
    for (int i=0; i<NbNodes; i++)
      delete Nodes(i);

  }

  
  /***********
   * METHODS *
   ***********/


  //! Reallocates the heap.
  /*! \param length new length.
    \note Data may be lost.
  */
  template <class T>
  inline void ArrayHeap<T>::Reallocate(int length)
  {

    if (length!=Nodes.GetLength())
      {
	Vector<T*> Nodes_copy(Nodes);
	Nodes.Reallocate(length);
	for (int i=0; i<min(length, Nodes_copy.GetLength()); i++)
	  Nodes(i) = Nodes_copy(i);
      }

  }


  //! Resizes the heap.
  /*! \param length new length.
    \param lastelement new number of nodes.
  */
  template <class T>
  inline void ArrayHeap<T>::Resize(int length, int lastelement)
  {

    this->Reallocate(length);
    if (lastelement != 0)
      NbNodes = lastelement;

  }


  //! Adds an element.
  /*! Adds an element and sorts the heap (to preserve heap property).
    \param X element to be added.
    \return The index of X in the heap.
  */
  template <class T>
  inline int ArrayHeap<T>::Add(T* X)
  {

    typename T::value_type x = X->GetValue();

    if (NbNodes == Nodes.GetLength())
      this->Reallocate(NbNodes + 1);

    Nodes(NbNodes) = X;
    int XIndex = NbNodes;
    int ParentIndex = (XIndex - 1) / 2;

    NbNodes++;

    while ((XIndex!=0) && (Nodes(ParentIndex)->GetValue() > x))
      {
	Nodes(XIndex) = Nodes(ParentIndex);
	XIndex = ParentIndex;
	ParentIndex = (XIndex - 1) / 2;
      }

    Nodes(XIndex) = X;

    return XIndex;

  }


  //! Adds an element.
  /*! Adds an element and sorts the heap (to preserve heap property).
    Then a matrix of back pointers is updated.
    \param X element to be added.
    \param Pointers matrix of back pointers which verifies
    Pointers(i, j) is the rank of X in the heap where X is so
    that X->GetX() == i and X->GetY() == j.
    \return The index of X in the heap.
    \note X is a node class that should supply GetX() and GetY().
  */
  template <class T>
  inline int ArrayHeap<T>::Add(T* X, Matrix<int>& Pointers)
  {

    typename T::value_type x = X->GetValue();

    if (NbNodes == Nodes.GetLength())
      this->Reallocate(NbNodes + 1);

    Nodes(NbNodes) = X;
    int XIndex = NbNodes;
    int ParentIndex = (XIndex - 1) / 2;

    NbNodes++;

    while ((XIndex!=0) && (Nodes(ParentIndex)->GetValue() > x))
      {
	Nodes(XIndex) = Nodes(ParentIndex);
	Pointers(Nodes(XIndex)->GetX(), Nodes(XIndex)->GetY()) = XIndex;
	XIndex = ParentIndex;
	ParentIndex = (XIndex - 1) / 2;
      }

    Nodes(XIndex) = X;

    Pointers(X->GetX(), X->GetY()) = XIndex;

    return XIndex;

  }


  //! Moves up the element corresponding in the heap
  //! until it meets the heap property.
  /*! \param XIndex rank of the element in the heap.
    \return The new rank of the element in the heap.
  */
  template <class T>
  inline int ArrayHeap<T>::MoveUp(int XIndex)
  {

    typename T::value_type x = Nodes(XIndex)->GetValue();
    T* X = Nodes(XIndex);

    int ParentIndex = (XIndex - 1) / 2;

    while ((XIndex!=0) && (Nodes(ParentIndex)->GetValue() > x))
      {
	Nodes(XIndex) = Nodes(ParentIndex);
	XIndex = ParentIndex;
	ParentIndex = (XIndex - 1) / 2;
      }

    Nodes(XIndex) = X;

    return XIndex;

  }


  //! Moves up an element in the heap until it meets
  //! the heap property.
  /*! A matrix of back pointers is also updated.
    \param XIndex rank of the element in the heap.
    \param Pointers matrix of back pointers which verifies
    Pointers(i, j) is the rank of X in the heap where X is so
    that X->GetX() == i and X->GetY() == j.
    \return The new rank of the element in the heap.
  */
  template <class T>
  inline int ArrayHeap<T>::MoveUp(int XIndex, Matrix<int>& Pointers)
  {

    typename T::value_type x = Nodes(XIndex)->GetValue();
    T* X = Nodes(XIndex);

    int ParentIndex = (XIndex - 1) / 2;

    while ((XIndex!=0) && (Nodes(ParentIndex)->GetValue() > x))
      {
	Nodes(XIndex) = Nodes(ParentIndex);
	Pointers(Nodes(XIndex)->GetX(), Nodes(XIndex)->GetY()) = XIndex;
	XIndex = ParentIndex;
	ParentIndex = (XIndex - 1) / 2;
      }

    Nodes(XIndex) = X;

    Pointers(X->GetX(), X->GetY()) = XIndex;

    return XIndex;

  }


  //! Removes the root of the heap.
  /*! Removes the root of the heap and preserves the heap property.
   */
  template <class T>
  inline void ArrayHeap<T>::DeleteRoot()
  {

    delete Nodes(0);

    NbNodes--;

    T* LastNode = Nodes(NbNodes);
    typename T::value_type NodeValue = LastNode->GetValue();

    int NodePos = 0;
    int FirstChildren = 2 * NodePos + 1;

    while (
	   ( (FirstChildren < NbNodes)
	     && (Nodes(FirstChildren)->GetValue() < NodeValue) )
	   ||
	   ( (FirstChildren + 1 < NbNodes)
	     && (Nodes(FirstChildren + 1)->GetValue() < NodeValue))
	   )
      {
	if ( (FirstChildren + 1 == NbNodes)
	     ||
	     ( Nodes(FirstChildren)->GetValue()
	       < Nodes(FirstChildren + 1)->GetValue() )
	     )
	  {
	    Nodes(NodePos) = Nodes(FirstChildren);
	    NodePos = FirstChildren;
	  }
	else
	  {
	    Nodes(NodePos) = Nodes(FirstChildren + 1);
	    NodePos = FirstChildren + 1;
	  }
	FirstChildren = 2 * NodePos + 1;
      }

    Nodes(NodePos) = LastNode;

  }

  
  //! Removes the root of the heap.
  /*! Removes the root of the heap and preserves the heap property.
    Then a matrix of back pointers is updated.
    \param Pointers matrix of back pointers which verifies
    Pointers(i, j) is the rank of X in the heap where X is so
    that X->GetX() == i and X->GetY() == j.
  */
  template <class T>
  inline void ArrayHeap<T>::DeleteRoot(Matrix<int>& Pointers)
  {

    Pointers(Nodes(0)->GetX(), Nodes(0)->GetY()) = -1;

    delete Nodes(0);

    NbNodes--;

    T* LastNode = Nodes(NbNodes);
    typename T::value_type NodeValue = LastNode->GetValue();

    int NodePos = 0;
    int FirstChildren = 2 * NodePos + 1;

    while ( ( (FirstChildren < NbNodes)
	      && (Nodes(FirstChildren)->GetValue() < NodeValue) )
	    ||
	    ( (FirstChildren + 1 < NbNodes)
	      && (Nodes(FirstChildren + 1)->GetValue() < NodeValue) )
	    )
      {
	if ( (FirstChildren + 1 == NbNodes)
	     || ( Nodes(FirstChildren)->GetValue()
		  < Nodes(FirstChildren + 1)->GetValue() )
	     )
	  {
	    Nodes(NodePos) = Nodes(FirstChildren);
	    Pointers(Nodes(NodePos)->GetX(), Nodes(NodePos)->GetY())
	      = NodePos;
	    NodePos = FirstChildren;
	  }
	else
	  {
	    Nodes(NodePos) = Nodes(FirstChildren + 1);
	    Pointers(Nodes(NodePos)->GetX(), Nodes(NodePos)->GetY())
	      = NodePos;
	    NodePos = FirstChildren + 1;
	  }
	FirstChildren = 2 * NodePos + 1;
      }

    Nodes(NodePos) = LastNode;

    Pointers(LastNode->GetX(), LastNode->GetY()) = NodePos;

  }


  /************************
   * CONVENIENT FUNCTIONS *
   ************************/


  //! Returns a pointer to the i-th node.
  /*! \param i index of the node.
    \return A pointer to the element stored at node #i.
  */
  template <class T>
  inline T* ArrayHeap<T>::operator() (int i)
  {

    return Nodes(i);

  }


  template <class T>
  inline T* ArrayHeap<T>::GetRoot()
  {

    return Nodes(0);

  }


  //! Returns a pointer to the root node.
  /*! \return A pointer to the element stored at the root.
   */
  template <class T>
  inline T* ArrayHeap<T>::GetRoot() const
  {

    return Nodes(0);

  }


  //! Returns a pointer to the root node.
  /*! \return A pointer to the element stored at the root.
   */
  template <class T>
  inline int ArrayHeap<T>::GetNbNodes() const
  {

    return NbNodes;

  }


  //! Is the heap empty?
  /*! \return true if the heap is empty and false otherwise.
   */
  template <class T>
  inline bool ArrayHeap<T>::IsEmpty() const
  {

    return (NbNodes == 0);

  }


}  // Multivac.


#define ARRAYHEAD_CXX
#endif
