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

#incluce <cstdint>

namespace Lattice {

template <class value_type, char periodicity, uint8_t dim>
class LatticeBase {
 public:
  virtual ~LatticeBase() {}
  virtual value_type &operator[](index_t index) = 0;
  virtual void halo_communication() = 0;
  virtual void set_halo(uint8_t halo_size[6]) = 0;

  virtual iterator begin() const = 0;
  virtual iterator end() const = 0;
  virtual index_t size() const = 0;
  virtual double element_size(index_t i) = 0;
  
  virtual index_t *element_neighbors(index_t i) = 0;
  virtual index_t n_neighbors() const = 0;
  
  typedef index_t uint64_t;
      
  template<class value_type>
  class iterator {
   public:
    virtual value_type operator*() const = 0;
    virtual iterator<value_type> &operator=(const iterator<value_type> &rhs) = 0;
    virtual bool operator==(const iterator<value_type> const &rhs) const = 0;
    virtual bool operator!=(const iterator<value_type> const &rhs) const { return ! operator==(rhs); };
    virtual iterator<value_type> &operator++() = 0;
  };   
};

};

#endif
