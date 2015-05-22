#ifndef __GEOMETRY_CONSTRAINT_HPP
#define __GEOMETRY_CONSTRAINT_HPP

#include "Constraint.hpp"
#include "Shape.hpp"

namespace Constraints {

  enum ReflectionType { REFLECTION_NONE, REFLECTION_NORMAL, REFLECTION_NORMAL_TANGENTIAL };

  struct GeometryConstraint : public Constraint {
  public:
    GeometryConstraint(ConstraintType _type, Shape *shape, bool _penetrable, ReflectionType _reflection_type, bool _tuneable_slip = false) : 
      Constraint(_type), penetrable(_penetrable), tuneable_slip(_tuneable_slip), reflection_type(_reflection_type), m_shape(shape) {};
    static std::string name() { return Constraint::name() + std::string("GeometryConstraint::") + m_shape->name(); }
    void reflect_particle(Particle *p, const double *distance_vector, const double *folded_pos);
    bool penetrable;
    bool tuneable_slip;
    ReflectionType reflection_type;
    Shape *m_shape;
  };
};

#endif
