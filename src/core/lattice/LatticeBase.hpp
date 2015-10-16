/*
  Copyright (C) 2010,2011,2012,2013,2014,2015 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
/** \file LatticeBase.hpp
 */

#ifndef __LATTICE_BASE_HPP
#define __LATTICE_BASE_HPP

#include <cstdint>

namespace Lattice {

typedef uint64_t index_t;

enum Periodicity { NONPERIODIC = 0, XPERIODIC = 1, YPERIODIC = 2, ZPERIODIC = 4 };

/*!
  Base class for Lattices. For the sake of generality the class provides only
  a forward iterator for lattice traversal.
*/
template <class value_type, class iterator_type>
class LatticeBase {
 public:
  typedef uint64_t index_t; /*! Index type. */
  typedef iterator_type iterator; /*! Lattice iterator type, should be compatible to STL forward iterators */
  
  virtual ~LatticeBase() {}
  /*! \brief Communicate halo */
  virtual void halo_communication() = 0;

  /** \brief Iterator pointing to the first element of the lattice. */
  virtual iterator begin() const = 0;
  virtual iterator end() const = 0;
  /** \brief Number of nodes in the lattice. */
  virtual index_t size() const = 0;
  /** \brief Get the size of an element
   *  \param i Index of the element to get the size for.
  virtual double element_size(index_t i) = 0;

  /** \brief Get node neighbors.
   *  \param node Node to get neighbors for.
   *  \returns Array of pointers no the neighbor elements. Entries for missing neighbors are 0.
   */
  virtual value_type **element_neighbors(value_type *node) = 0;
  /** \brief Get the maximal number of neighbors. */
  virtual index_t n_neighbors() const = 0;        
};

};

#endif
