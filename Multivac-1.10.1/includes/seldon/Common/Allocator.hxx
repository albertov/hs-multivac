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

#ifndef SELDON_FILE_ALLOCATOR_HXX

namespace Seldon
{


  /////////////////
  // MALLOCALLOC //
  /////////////////


  template <class T>
  class MallocAlloc
  {
  public:
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;

  public:

    pointer allocate(int num, void* h = 0)
    {
      return static_cast<pointer>( malloc(num * sizeof(T)) );
    }

    void deallocate(pointer data, int num, void* h = 0)
    {
      free(data);
    }

    void* reallocate(pointer data, int num, void* h = 0)
    {
      return realloc(reinterpret_cast<void*>(data), num * sizeof(T));
    }

    void memoryset(pointer data, char c, size_t num)
    {
      memset(reinterpret_cast<void*>(data), c, num);
    }

    void memorycpy(pointer datat, pointer datas, size_t num)
    {
      memcpy(reinterpret_cast<void*>(datat), reinterpret_cast<void*>(datas),
	     num * sizeof(T));
    }
  };


  /////////////////
  // CALLOCALLOC //
  /////////////////


  template <class T>
  class CallocAlloc
  {
  public:
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;

  public:

    pointer allocate(int num, void* h = 0)
    {
      return static_cast<pointer>( calloc(num, sizeof(T)) );
    }

    void deallocate(pointer data, int num, void* h = 0)
    {
      free(data);
    }

    void* reallocate(pointer data, int num, void* h = 0)
    {
      return realloc(reinterpret_cast<void*>(data), num * sizeof(T));
    }

    void memoryset(pointer data, char c, size_t num)
    {
      memset(reinterpret_cast<void*>(data), c, num);
    }

    void memorycpy(pointer datat, pointer datas, size_t num)
    {
      memcpy(reinterpret_cast<void*>(datat), reinterpret_cast<void*>(datas),
	     num * sizeof(T));
    }
  };


  //////////////
  // NEWALLOC //
  //////////////


  template <class T>
  class NewAlloc
  {
  public:
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;

  public:

    pointer allocate(int num, void* h = 0)
    {
      return static_cast<pointer>(new T[num]);
    }

    void deallocate(pointer data, int num, void* h = 0)
    {
      delete [] data;
    }

    void* reallocate(pointer data, int num, void* h = 0)
    {
      if (data != NULL)
	delete [] data;
      return (new T[num]);
    }

    void memoryset(pointer data, char c, size_t num)
    {
      memset(reinterpret_cast<void*>(data), c, num);
    }

    void memorycpy(pointer datat, pointer datas, size_t num)
    {
      for (size_t i = 0; i < num; i++)
	datat[i] = datas[i];
    }
  };


  //////////////
  // NANALLOC //
  //////////////


  template <class T>
  class NaNAlloc
  {
  public:
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;

  public:

    pointer allocate(int num, void* h = 0)
    {
      pointer data = static_cast<pointer>( malloc(num * sizeof(T)) );
      if (numeric_limits<value_type>::has_signaling_NaN)
	for (int i = 0; i < num; i++)
	  data[i] = numeric_limits<value_type>::signaling_NaN();
      else if (numeric_limits<value_type>::has_quiet_NaN)
	for (int i = 0; i < num; i++)
	  data[i] = numeric_limits<value_type>::quiet_NaN();
      else if  (numeric_limits<value_type>::has_infinity)
	for (int i = 0; i < num; i++)
	  data[i] = numeric_limits<value_type>::infinity();
      return data;
    }

    void deallocate(pointer data, int num, void* h = 0)
    {
      free(data);
    }

    void* reallocate(pointer data, int num, void* h = 0)
    {
      void* datav = realloc(reinterpret_cast<void*>(data), num * sizeof(T));
      pointer datap = reinterpret_cast<pointer>(datav);
      if (numeric_limits<value_type>::has_signaling_NaN)
	for (int i = 0; i < num; i++)
	  datap[i] = numeric_limits<value_type>::signaling_NaN();
      else if (numeric_limits<value_type>::has_quiet_NaN)
	for (int i = 0; i < num; i++)
	  datap[i] = numeric_limits<value_type>::quiet_NaN();
      else if  (numeric_limits<value_type>::has_infinity)
	for (int i = 0; i < num; i++)
	  datap[i] = numeric_limits<value_type>::infinity();
      return datav;
    }

    void memoryset(pointer data, char c, size_t num)
    {
      memset(reinterpret_cast<void*>(data), c, num);
    }

    void memorycpy(pointer datat, pointer datas, size_t num)
    {
      memcpy(reinterpret_cast<void*>(datat), reinterpret_cast<void*>(datas),
	     num * sizeof(T));
    }
  };


} // namespace Seldon.

#define SELDON_FILE_ALLOCATOR_HXX
#endif
