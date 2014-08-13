#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/fill.h>
#include <thrust/copy.h>
#include <thrust/sort.h>
#include <cstdlib>
#include <iostream>
#include <ctime>

#include "domain_decomposition_gpu.hpp"
#include "cuda_utils.hpp"

/* Forward declarations of gpu kernels. */

__global__ static void sortParticlesGenerateCellist(unsigned int n_part, const float3 *xyz, float3 *xyz_sorted, unsigned int *hashes, unsigned int *indexes, uint2 *cells);

__global__ static void hashAtoms(unsigned int n_part, float3 hi, const float3 *xyz, unsigned int *index, unsigned int *hashes, uint3 n_cells);

/* Class Implementation */

DomainDecompositionGpu::DomainDecompositionGpu(float3 _box, unsigned int _n_part, uint3 _n_cells) {

  n_part = _n_part;
  box = _box;
  n_cells = _n_cells;
  hi.x = n_cells.x/box.x;
  hi.y = n_cells.y/box.y;
  hi.z = n_cells.z/box.z;
  total_cells = n_cells.x*n_cells.y*n_cells.z;

  init_device_memory(true, true);
}

DomainDecompositionGpu::~DomainDecompositionGpu() {
  free_device_memory(true, true);
}

void DomainDecompositionGpu::free_device_memory(bool particles, bool dd) {
  if(particles) {
    cuda_safe_mem(cudaFree(indexes));
    cuda_safe_mem(cudaFree(hashes));
    cuda_safe_mem(cudaFree(xyz_sorted));
  }
  if(dd) {
    cuda_safe_mem(cudaFree(cells));
  }
}

void DomainDecompositionGpu::init_device_memory(bool particles, bool dd) {
  if(particles) {
    cuda_safe_mem(cudaMalloc((void **)&(indexes), n_part*sizeof(unsigned int)));
    cuda_safe_mem(cudaMalloc((void **)&(hashes), n_part*sizeof(unsigned int)));
    cuda_safe_mem(cudaMalloc((void **)&(xyz_sorted), n_part*sizeof(float3)));
  }
  if(dd) {
    cuda_safe_mem(cudaMalloc((void **)&(cells), total_cells*sizeof(uint2)));
  }
}

void DomainDecompositionGpu::build(float3 *xyz) {
  dim3 block(256,1,1), grid(1,1,1);
  
  grid.x = max(1, (n_part + block.x - 1) / block.x);

  #ifdef CUDA_DD_DEBUG
  printf("build_dd() grid.x = %d, block.x = %d\n", grid.x, block.x);
  #endif

  cuda_safe_mem(cudaMemset(cells, CELL_EMPTY, total_cells*sizeof(uint2)));

  hashAtoms<<<grid, block>>>(n_part, hi, xyz, indexes, hashes, n_cells);

  thrust::sort_by_key(thrust::device_ptr<unsigned int>(hashes),
		      thrust::device_ptr<unsigned int>(hashes + n_part),
		      thrust::device_ptr<unsigned int>(indexes));

  unsigned int smemsize = sizeof(unsigned int)*(block.x + 1);

  sortParticlesGenerateCellist<<<grid,block,smemsize>>>(n_part, xyz, xyz_sorted, hashes, indexes, cells);
}

void DomainDecompositionGpu::set_n_part(unsigned int _n_part) {
  n_part = _n_part;
  free_device_memory(true, false);
  init_device_memory(true, false);
}

/* GPU Kernels */

__global__ static void sortParticlesGenerateCellist(unsigned int n_part, const float3 *xyz, float3 *xyz_sorted, unsigned int *hashes, unsigned int *indexes, uint2 *cells) {
  unsigned int id = blockDim.x*blockIdx.x + threadIdx.x;
  unsigned int hash;

  extern __shared__ unsigned int hash_s[];

  if(id < n_part) {
    hash = hashes[id];
    hash_s[threadIdx.x+1] = hash;
    if( (id > 0) && (threadIdx.x == 0))
      hash_s[0] = hashes[id-1];
  }

  __syncthreads();
  
  if(id < n_part) {
    if( (id == 0) || hash != hash_s[threadIdx.x] ) {
      cells[hash].x = id;
      if(id > 0)
	cells[hash_s[threadIdx.x]].y = id;
    }
    if( id == n_part - 1)
      cells[hash].y = id + 1;

    unsigned int sorted_id;
    sorted_id = indexes[id];
    xyz_sorted[id] = xyz[sorted_id];
  }
}

__global__ static void hashAtoms(unsigned int n_part, float3 hi, const float3 *xyz, unsigned int *index, unsigned int *hashes, uint3 n_cells) {
  unsigned int id = blockDim.x*blockIdx.x + threadIdx.x;

  if(id >= n_part)
    return;

  int3 cell;
  float3 part = xyz[id];  

  cell.x = ((int)(part.x * hi.x + 0.5)) % n_cells.x;
  cell.y = ((int)(part.y * hi.y + 0.5)) % n_cells.y;
  cell.z = ((int)(part.z * hi.z + 0.5)) % n_cells.z;

  hashes[id] =  cell.x * n_cells.y * n_cells.z +  cell.y * n_cells.z + cell.z;
  index[id] = id;
}

__device__ static inline float dist2(float3 a, float3 b) {
  return ( (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}

__global__ static void nearestNeighbors(unsigned int n_part, uint3 n_cells, const uint2 *cells, const float3 *xyz, unsigned int *neighbors, int range, float rcut2) {
  unsigned int xindex = blockDim.x * blockIdx.x + threadIdx.x;
  unsigned int yindex = blockDim.y * blockIdx.y + threadIdx.y;
  unsigned int zindex = blockDim.z * blockIdx.z + threadIdx.z;
  int n,m,l;
  unsigned int cellhash = n_cells.y * n_cells.z * xindex + n_cells.z * yindex + zindex;
  int hash;
  unsigned int m_neighbors = 0;

  if( (xindex >= n_cells.x) || (yindex >= n_cells.y) || (zindex >= n_cells.z))
    return;

  for(unsigned int id = cells[cellhash].x; id != cells[cellhash].y; id++) {
    m_neighbors = 0;
    float3 parti = xyz[id];
    for(int i = -range; i <= range; i++) {
      n = xindex + i;
      if((n < 0) || (n >= n_cells.x))
	continue;
      for(int j = -range; j <= range; j++) {
	m = yindex + j;
	if((m < 0) || (m >= n_cells.y))
	  continue;
	for(int k = -range; k <= range; k++) {
	  l = zindex + k;
	  if((l < 0) || (l >= n_cells.z))
	    continue;

	  hash = n_cells.y * n_cells.z * n + n_cells.z * m + l;
	  if( (cells[hash].x == CELL_EMPTY))
	    continue;
	  for(unsigned int jd = cells[hash].x; jd != cells[hash].y; jd++) {
	    if(dist2(parti, xyz[jd]) <= rcut2) {
	      m_neighbors++;
	    }
	  }
	}
      }
    }
    neighbors[id] = m_neighbors;
  }
}

__global__ static void nearestNeighbors_n2(unsigned int n_part, const float3 *xyz, unsigned int *neighbors, float rcut2) {
  unsigned int id = blockIdx.x;
  unsigned int m_neighbors = 0;

  if(id >= n_part)
    return;

  float3 parti = xyz[id]; 

  for(unsigned int i = 0; i < n_part; i++) {
    if(dist2(parti, xyz[i]) <= rcut2) {
      m_neighbors++;
    }    
  }
  neighbors[id] = m_neighbors;
}


// void print_dd(dd_t *dd) {
//   puts("print_dd()");
//   uint2 *cells_h;
//   float3 *xyz_h;
//   unsigned int count = 0;

//   cells_h = (uint2 *)malloc(dd->total_cells*sizeof(uint2));
//   xyz_h = (float3 *)malloc(dd->n_part*sizeof(float3));

//   cudaMemcpy(cells_h, dd->cells, dd->total_cells*sizeof(uint2), cudaMemcpyDeviceToHost);
//   cudaMemcpy(xyz_h, dd->xyz_sorted, dd->n_part*sizeof(float3), cudaMemcpyDeviceToHost);  

//   printf("dd { n_part = %d, n_cells = ( %d %d %d ), hi = ( %f %f %f ) }\n", dd->n_part, dd->n_cells.x, dd->n_cells.y, dd->n_cells.z
// 	 , dd->hi.x, dd->hi.y, dd->hi.z);

//   bool *part_map = (bool *)malloc(dd->n_part*sizeof(bool));
//   for(int i = 0; i < dd->n_part; i++)
//     part_map[i] = false;
  
//   for(int i = 0; i < dd->total_cells; i++) {
//     if(cells_h[i].x == CELL_EMPTY)
//       continue;
//     unsigned int x,y,z;
//     z = i % dd->n_cells.z;
//     y = ((i - z)/dd->n_cells.z) % dd->n_cells.y;
//     x = (i - z - dd->n_cells.z*y) / (dd->n_cells.z*dd->n_cells.y);
//     printf("cell %d, pos (%d %d %d), center (%f %f %f)\n", i, x, y, z, x/dd->hi.x, y/dd->hi.y, z/dd->hi.z);
//     for(unsigned int it = cells_h[i].x; it != cells_h[i].y; it++) {
//       part_map[it] = true;
//       printf("\t%d: pos = (%f %f %f %f), dist = (%f %f %f), dist/h = (%f %f %f)\n", it,xyz_h[it].x,xyz_h[it].y,xyz_h[it].z,0,
// 	     xyz_h[it].x-x/dd->hi.x, xyz_h[it].y-y/dd->hi.y, xyz_h[it].z-z/dd->hi.z,
// 	     dd->hi.x*(xyz_h[it].x-x/dd->hi.x), dd->hi.y*(xyz_h[it].y-y/dd->hi.y), dd->hi.z*(xyz_h[it].z-z/dd->hi.z) );
//       count++;
//     }
//   }
//   printf("%d particles in dd.\n", count);
//   for(int i = 0; i < dd->n_part; i++)
//     if(!part_map[i])
//       printf("particle %d is missing.\n",i);

// }

