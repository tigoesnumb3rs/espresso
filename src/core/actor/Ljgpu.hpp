#ifndef __LJGPU_HPP
#define __LJGPU_HPP

#include "SystemInterface.hpp"
#include "Actor.hpp"

class DomainDecompositionGpu;

class Ljgpu : public Actor {
  Ljgpu(SystemInterface &s, float epsilon, float sigma, float cutoff);
  void computeForces(SystemInterface &s);
private:
  DomainDecompositionGpu &dd;
};

#endif
