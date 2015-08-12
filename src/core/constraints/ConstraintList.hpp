#ifndef __CONSTRAINT_CONTAINER_HPP
#define __CONSTRAINT_CONTAINER_HPP

#include "Constraint.hpp"

#include <map>

namespace Constraints {

  class ConstraintList : public std::map<int, Constraints::Constraint *> {
  public:
    ConstraintList() : m_next_id(0) {}
    int add_constraint(Constraints::Constraint *c);
    void remove_constraint(int i);
    void add_forces(Particle *p);
    void add_energies(Particle *p);
    void init_forces();
    double min_dist(double pos[3]);
  private:
    int m_next_id;
  };

  extern ConstraintList list;
}
  
#endif