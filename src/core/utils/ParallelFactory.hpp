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

#ifndef __UTILS_PARALLEL_FACTORY_HPP
#define __UTILS_PARALLEL_FACTORY_HPP

#include "ParallelFactory.hpp"
#include "Factory.hpp"

#include "communication.hpp"
#include "errorhandling.hpp"

namespace Utils {

/** Implements the same concept as Utils::Factory, but creates
    the instance on all nodes.
*/

template<class T, typename F = Factory<T>>
class ParallelFactory {
 public:
  typedef typename F::Builder Builder;

  template<class Derived>
  static std::shared_ptr<T> builder() {
    return F::template builder<Derived>();
    return 0;
  }
  
  static std::shared_ptr<T> make(const std::string &name) {
    mpi_call(ParallelFactory<T>::mpi_slave, name.size(), 0);

    MPI_Bcast(const_cast<char *>(&(*name.begin())), name.size(), MPI_CHAR, 0, comm_cart);
    
    return F::Instance().make(name);
  }

  static void mpi_slave(int name_size, int) {
    std::string name;

    name.resize(name_size);
    MPI_Bcast(&(*name.begin()), name_size, MPI_CHAR, 0, comm_cart);
    
    try {      
      std::shared_ptr<T> p = F::Instance().make(name);

    } catch(std::exception &e) {
      runtimeErrorMsg() << e.what();
    }
  }

  static bool register_new(const std::string &name, const Builder &b) {
    static bool callback_added = false;

    if(!callback_added) {
      std::cout << this_node << ": Added callback.\n";
      mpi_add_callback(mpi_slave);
      callback_added = true;
    }

    return F::Instance().register_new(name, b);
  }
  
  static bool has_builder(const std::string &name) {
    return F::Instance().has_builder(name);
  }
};

}

#endif
