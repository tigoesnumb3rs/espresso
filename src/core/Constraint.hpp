/*
  Copyright (C) 2014 The ESPResSo project
  
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

#ifndef __CONSTRAINT_HPP
#define __CONSTRAINT_HPP

#include "particle_data.hpp"

class Constraint {
public:
  virtual const double dist(Particle *p1, double ppos[3], Particle *c_p, double *dist, double *vec) const = 0;
  virtual const int ia_type() const = 0;
  virtual *double force() const = 0;
};

class Wall : Constraint {
public:
  Wall(double _normal[3], double _d) : d(_d) {
    for(int i = 0; i < 3; ++i)
      normal[i] = _normal[i];
  };
  const int ia_type() const { return ia_type; };
private:
  double normal[3];
  double d;
  int type;
  int ia_type;
}

class ConstraintFactory;

#endif
