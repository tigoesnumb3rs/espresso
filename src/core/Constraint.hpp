#ifndef __CONSTRAINT_HPP
#define __CONSTRAINT_HPP

#include <string>

#include "energy.hpp"
#include "interaction_data.hpp"

namespace ConstraintClass {

  enum ConstraintType { CONSTRAINT_NONE };

  struct Constraint {
  public:
    Constraint(ConstraintType _type) : id(-1) { type = _type; };
    ConstraintType type;
    virtual void add_energy(const Particle *p, const double *folded_pos, Observable_stat &energy) { }
    virtual void add_force(Particle *p, const double *folded_pos) { }
    /* Numerical id for interface binding */
    int id;
    /* Human readable name */
    virtual std::string name() { return std::string("constraint_unnamed"); }
  };
};

#endif
