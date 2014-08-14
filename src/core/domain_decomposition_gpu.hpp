#ifndef __GPU_DD_HPP
#define __GPU_DD_HPP

#include <cuda.h>

#include "cuda_interface.hpp"

#define CELL_EMPTY 0xffffffff

class DomainDecompositionGpu {
public:
  DomainDecompositionGpu() : ready(false), n_part(0), cutoff(0.f), indexes(0), hashes(0), cells(0), xyz_sorted(0) {};
  DomainDecompositionGpu(float3 _box, unsigned int _n_part, float _cutoff);
  ~DomainDecompositionGpu();
  void build(float3 *xyz);
  /* Get array in sorted order ...*/
  template<typename T>
  void sort(T *src, T *dst);
  /* ... and back. */
  template<typename T>
  void unsort(T *src, T *dst);  
  void set_n_part(unsigned int _n_part);
  void set_cutoff(float _cutoff);
  void set_box(float3 _box);
  unsigned int get_n_part() { return n_part; };
  uint3 get_n_cells() { return n_cells; };
  /* These are ugly but have to do for now, shared_ptr is c++11. */
  uint2 *get_cells() { return cells; };
  float3 *get_xyz_sorted() { return xyz_sorted; };
  friend bool test_decomposition(DomainDecompositionGpu &dd);
#ifdef GPU_DD_DEBUG
  friend void print_dd(DomainDecompositionGpu &dd);
#endif
protected:
  /** These should not be called directly by the user
      if you really want to, create a derived class.
  */
  DomainDecompositionGpu(float3 _box, unsigned int _n_part, uint3 _n_cells);
  void set_n_cells(uint3 _n_cells);
private:
  bool ready, box_set, cutoff_set, n_part_set, n_cells_set;  
  float3 hi;
  float3 box;
  unsigned int n_part;
  uint3  n_cells;
  unsigned int total_cells;
  float cutoff;

  unsigned int *indexes;
  unsigned int *hashes;
  uint2 *cells;
  float3 *xyz_sorted;  
  void init_device_memory(bool particles, bool dd);
  void free_device_memory(bool particles, bool dd);
  void realloc_device_memory(bool particles, bool dd) {
    free_device_memory(particles, dd);
    init_device_memory(particles, dd);
  };
  void update_cells_from_cutoff_and_box();
  void check_ready() { ready = box_set && n_cells_set && n_part_set; };
};

#endif
