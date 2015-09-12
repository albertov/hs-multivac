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

#ifndef SELDON_FILE_COMMON_HXX

class Spacetown
{
};

template <class T>
void PrintArray(T* v, int lgth)
{
  for (int k = 0; k < lgth - 1; k++)
    std::cout << v[k] << " | ";
  std::cout << v[lgth - 1] << std::endl;
}

namespace Seldon
{


  //! Converts most types to string.
  /*!
    \param input variable to be converted.
    \return A string containing 'input'.
  */
  template<typename T>
  std::string to_str(const T& input)
  {
    std::ostringstream output;
    output << input;
    return output.str();
  }

  //! Converts string to most types, specially numbers.
  /*!
    \param input string to be converted.
    \return 'input' converted to 'T'.
  */
  template <class T>
  void to_num(std::string s, T& num)
  {
    std::istringstream str(s);
    str >> num;
  }

  //! Converts string to most types, specially numbers.
  /*!
    \param input string to be converted.
    \return 'input' converted to 'T'.
  */
  template <class T>
  T to_num(std::string s)
  {
    T num;
    std::istringstream str(s);
    str >> num;
    return num;
  }


}  // namespace Seldon.

#define SELDON_FILE_COMMON_HXX
#endif
