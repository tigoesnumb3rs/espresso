#ifndef CUDA_DD_HPP
#define CUDA_DD_HPP

#include <cuda.h>

#include "cuda_interface.hpp"

#define CELL_EMPTY 0xffffffff

class DomainDecompositionGpu {
public:
  DomainDecompositionGpu(float3 _box, unsigned int _n_part, uint3 _n_cells);
  //  DomainDecompositionGpu(float3 _box, unsigned int _n_part, float cutoff);
  ~DomainDecompositionGpu();
  void build(float3 *xyz);
  void set_n_part(unsigned int _n_part);
  //  void print() const;
private:
  float3 hi;
  float3 box;
  unsigned int n_part;
  uint3  n_cells;
  unsigned int total_cells;

  unsigned int *indexes;
  unsigned int *hashes;
  uint2 *cells;
  float3 *xyz_sorted;  
  void init_device_memory(bool particles, bool dd);
  void free_device_memory(bool particles, bool dd);
};

#endif
