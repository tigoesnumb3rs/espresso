#ifndef __GPU_DD_HPP
#define __GPU_DD_HPP

#include <cuda.h>

#include "cuda_interface.hpp"

#define GPU_DD_DEBUG

#define CELL_EMPTY 0xffffffff

class DomainDecompositionGpu {
public:
  DomainDecompositionGpu(float3 _box, unsigned int _n_part, uint3 _n_cells);
  DomainDecompositionGpu(float3 _box, unsigned int _n_part, float cutoff);
  ~DomainDecompositionGpu();
  void build(float3 *xyz);
  /* Get array in sorted order ...*/
  template<typename T>
  void sort(T *src, T *dst);
  /* ... and back. */
  template<typename T>
  void unsort(T *src, T *dst);  
  void set_n_part(unsigned int _n_part);
  void set_cutoff(float cutoff);
  unsigned int get_n_part() { return n_part; };
  uint3 get_n_cells() { return n_cells; };
  void set_n_cells(uint3 _n_cells);
  /* These are ugly but have to do for now, shared_ptr is c++11. */
  uint2 *get_cells() { return cells; };
  float3 *get_xyz_sorted() { return xyz_sorted; };
  friend bool test_decomposition(DomainDecompositionGpu &dd);
#ifdef GPU_DD_DEBUG
  friend void print_dd(DomainDecompositionGpu &dd);
#endif
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
  void update_cells_from_cutoff(float cutoff);
};

#endif
