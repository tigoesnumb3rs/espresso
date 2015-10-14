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

template <class value_type, class iterator_type>
class LatticeBase {
 public:
  typedef uint64_t index_t;
  typedef iterator_type iterator;
  
  virtual ~LatticeBase() {}
  virtual void halo_communication() = 0;
  virtual void set_halo(uint8_t halo_size[6]) = 0;

  virtual iterator begin() const = 0;
  virtual iterator end() const = 0;
  virtual index_t size() const = 0;
  virtual double element_size(index_t i) = 0;
  
  virtual value_type **element_neighbors(index_t i) = 0;
  virtual index_t n_neighbors() const = 0;        
};

};

#endif
