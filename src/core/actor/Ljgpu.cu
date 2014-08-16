#include "Ljgpu.hpp"
#include "domain_decomposition_gpu.hpp"

Ljgpu::Ljgpu(SystemInterface &s, float epsilon, float sigma, float cutoff) :
  dd(s, cutoff) {}
