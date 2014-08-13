#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/fill.h>
#include <thrust/copy.h>
#include <thrust/sort.h>
#include <cstdlib>
#include <iostream>
#include <ctime>

#include "domain_decomposition_gpu.hpp"
#include "domain_decomposition_gpu_test.hpp"
#include "cuda_utils.hpp"

#ifdef GPU_DD_DEBUG
void print_dd(DomainDecompositionGpu &dd);
#endif

/* Forward declarations of gpu kernels. */

__global__ static void sortParticlesGenerateCellist(unsigned int n_part, const float3 *xyz, float3 *xyz_sorted, unsigned int *hashes, unsigned int *indexes, uint2 *cells);

__global__ static void hashAtoms(unsigned int n_part, float3 hi, const float3 *xyz, unsigned int *index, unsigned int *hashes, uint3 n_cells);

template<typename T>
__global__ static void unsortArrays(unsigned int *indexes, T *src, T *dst, unsigned int n);

/* Class Implementation */

DomainDecompositionGpu::DomainDecompositionGpu(float3 _box, unsigned int _n_part, uint3 _n_cells) : n_part(_n_part), box(_box) {

  n_cells = _n_cells;
  hi.x = n_cells.x/box.x;
  hi.y = n_cells.y/box.y;
  hi.z = n_cells.z/box.z;
  total_cells = n_cells.x*n_cells.y*n_cells.z;

  init_device_memory(true, true);
}

DomainDecompositionGpu::DomainDecompositionGpu(float3 _box, unsigned int _n_part, float cutoff) : n_part(_n_part), box(_box) {
  /* Need to take the floor to be on the safe side */
  n_cells.x = floor(box.x / cutoff);
  n_cells.y = floor(box.y / cutoff);
  n_cells.z = floor(box.z / cutoff);

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
  dim3 block(1,1,1), grid(1,1,1);

  if(n_part < 128) {
    block.x = n_part;
    grid.x = 1;    
  } else {
    block.x = 128;
    grid.x = max(1, (n_part + block.x - 1) / block.x);
  }

  cuda_safe_mem(cudaMemset(cells, CELL_EMPTY, total_cells*sizeof(uint2)));

  KERNELCALL(hashAtoms, grid, block, (n_part, hi, xyz, indexes, hashes, n_cells));

  thrust::sort_by_key(thrust::device_ptr<unsigned int>(hashes),
  		      thrust::device_ptr<unsigned int>(hashes + n_part),
  		      thrust::device_ptr<unsigned int>(indexes));

  unsigned int smemsize = sizeof(unsigned int)*(block.x + 1);

  KERNELCALL_shared(sortParticlesGenerateCellist, grid, block, smemsize, (n_part, xyz, xyz_sorted, hashes, indexes, cells));
}

void DomainDecompositionGpu::set_n_part(unsigned int _n_part) {
  n_part = _n_part;  
  free_device_memory(true, false);
  init_device_memory(true, false);
}

template<typename T>
void DomainDecompositionGpu::unsort(T *src, T* dst) {
  dim3 block(1,1,1), grid(1,1,1);

  if(n_part < 128) {
    block.x = n_part;
    grid.x = 1;    
  } else {
    block.x = 128;
    grid.x = max(1, (n_part + block.x - 1) / block.x);
  }

  KERNELCALL(unsortArrays<T>, grid, block, (indexes, src, dst, n_part));
}

/* GPU Kernels */

template<typename T>
__global__ static void unsortArrays(unsigned int *indexes, T *src, T *dst, unsigned int n) {
  unsigned int id = blockDim.x*blockIdx.x + threadIdx.x;

  if(id >= n)
    return;

  dst[indexes[id]] = src[id];
}

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

  cell.x = ((int)(part.x * hi.x)) % n_cells.x;
  cell.y = ((int)(part.y * hi.y)) % n_cells.y;
  cell.z = ((int)(part.z * hi.z)) % n_cells.z;

  hashes[id] =  cell.x * n_cells.y * n_cells.z +  cell.y * n_cells.z + cell.z;
  index[id] = id;
}

  /** Unit test for the GPU Domain Decomposition
      ------------------------------------------
      The testing strategy is as follows:
      First the memory management is testet by changing the number of particles and
      the cutoff/number of cells and see if this is handled correctly.
      Then the function is testet by calculation the neighbor count of a know particle
      configuration. This is also calculated via an nsquare algorithm to check its implementation. 
      Then the neighbor count of a random system is compared between the calculation with domain decomposition and the nquared algorithm.
  **/

  /* Functions for unit testing */
static bool test_memory_management();
static bool test_neighbor_count();
__global__ static void nearestNeighbors(unsigned int n_part, uint3 n_cells, const uint2 *cells, const float3 *xyz, unsigned int *neighbors, float rcut2);
__global__ static void nearestNeighbors_n2(unsigned int n_part, const float3 *xyz, unsigned int *neighbors, float rcut2);

bool domain_decomposition_gpu_unit_test() {
  bool result = true;

  result &= test_memory_management();
  result &= test_neighbor_count();

  return result;
}

__device__ static inline float dist2(float3 a, float3 b) {
  return ( (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}

static bool test_memory_management() {
  return false;
}

bool test_decomposition(DomainDecompositionGpu &dd) {
  uint2 *cells_h;
  float3 *xyz_h;
  float3 ll, ur;
  unsigned int x,y,z;
  unsigned int n_part = dd.get_n_part();
  uint3 n_cells = dd.get_n_cells();
  unsigned int total_cells = n_cells.x * n_cells.y * n_cells.z;

  cells_h = (uint2 *)malloc(dd.total_cells*sizeof(uint2));
  xyz_h = (float3 *)malloc(n_part*sizeof(float3));

  cudaMemcpy(cells_h, dd.cells, total_cells*sizeof(uint2), cudaMemcpyDeviceToHost);
  cudaMemcpy(xyz_h, dd.xyz_sorted, n_part*sizeof(float3), cudaMemcpyDeviceToHost);  

  int *part_map = (int *)malloc(n_part*sizeof(int));
  memset(part_map, 0, n_part*sizeof(int));
  
  for(int i = 0; i < total_cells; i++) {    
    /* 3d index of the cell */
    z = i % n_cells.z;
    y = ((i - z)/n_cells.z) % n_cells.y;
    x = (i - z - n_cells.z*y) / (n_cells.z*n_cells.y);
    /* Cell boundaries */
    ll.x = x / dd.hi.x;
    ll.y = y / dd.hi.y;
    ll.z = z / dd.hi.z;
    ur.x = ll.x + 1./dd.hi.x;
    ur.y = ll.y + 1./dd.hi.y;
    ur.z = ll.z + 1./dd.hi.z;

    if(cells_h[i].x == CELL_EMPTY)
      continue;
    for(unsigned int it = cells_h[i].x; it != cells_h[i].y; it++) {      
      part_map[it]++;
      /* Check if particle blongs in this cell */
      if( ( ( xyz_h[it].x < ll.x ) || ( xyz_h[it].y < ll.y ) || ( xyz_h[it].z < ll.z ) ) ||
	  ( ( xyz_h[it].x >= ur.x ) || ( xyz_h[it].y >= ur.y ) || ( xyz_h[it].z >= ur.z ) ) ) {
	#ifdef GPU_DD_DEBUG
	printf("domain_decomposition_gpu: particle %d is in wrong cell.\n", it);
	#endif
	return false;
      }
    }
  }

  /* Check if all particles are in the dd */
  for(int i = 0; i < n_part; i++) {    
    if(part_map[i] != 1) {
      #ifdef GPU_DD_DEBUG
      printf("domain_decomposition_gpu: particle %d is missing or overcounted.\n", i);
      #endif
      return false;
    }
  }

  return true;
}

static bool test_neighbor_count() {
  float3 box;
  float cutoff;
  unsigned int n_part;
  
  box.x = 10;
  box.y = 20;
  box.z = 30;
  cutoff = 1.5;

  float3 parts[11] = { {0.0f, 00.f, 0.0f }, // 1 
		       {9.9f, 0.0f, 0.0f }, // 1
		       {6.0f, 5.0f, 5.0f }, // 6
		       {9.9f, 9.9f, 0.0f }, // 1
		       {5.0f, 5.0f, 5.0f }, // 7
		       {5.0f, 4.0f, 5.0f }, // 6
		       {0.0f, 9.9f, 0.0f }, // 1 
		       {5.0f, 5.0f, 4.0f }, // 6 
		       {5.0f, 6.0f, 5.0f }, // 6
		       {5.0f, 5.0f, 6.0f }, // 6
		       {4.0f, 5.0f, 5.0f }  // 6 
  };

unsigned int correct_count[11] = { 1, 1, 6, 1, 7, 6, 1, 6, 6, 6, 6 };

  n_part = sizeof(parts)/sizeof(float3);

  float3 *parts_d = 0;
  cuda_safe_mem(cudaMalloc((void **)&parts_d, sizeof(parts)));
  cuda_safe_mem(cudaMemcpy(parts_d, parts, sizeof(parts), cudaMemcpyHostToDevice));
  
  DomainDecompositionGpu dd(box, n_part, cutoff);

  /* Check n_cells from box and cutoff */
  uint3 n_cells = dd.get_n_cells();
  if((n_cells.x != 6) || (n_cells.y != 13) || (n_cells.z != 20)) {
    return false;
  }

  /* Sort particles into cells */
  dd.build(parts_d);

  #ifdef GPU_DD_DEBUG
  print_dd(dd);
  #endif

  /* Check result */
  if(!test_decomposition(dd))
    return false;

  unsigned int *neighbors_dd_d, *neighbors_n2_d;
  unsigned int neighbors_dd[sizeof(parts)/sizeof(float3)];
  unsigned int neighbors_n2[sizeof(parts)/sizeof(float3)];

  cuda_safe_mem(cudaMalloc((void **)&neighbors_dd_d, n_part * sizeof(unsigned int)));
  cuda_safe_mem(cudaMalloc((void **)&neighbors_n2_d, n_part * sizeof(unsigned int)));

  dim3 block_dd(n_cells.x,n_cells.y,n_cells.z/2), grid_dd(1,1,n_cells.z/2);
  dim3 block_n2(n_part,1,1), grid_n2(1,1,1);
  float cutoff2 = cutoff*cutoff;

  /* Count using domain decomposition and sorted position array */
  KERNELCALL(nearestNeighbors, grid_dd, block_dd, (n_part, dd.get_n_cells(), dd.get_cells(), dd.get_xyz_sorted(), neighbors_dd_d,cutoff2));
  /* Count directly, also using sorted position array to get same particle order. */
  KERNELCALL(nearestNeighbors_n2, grid_n2, block_n2, (n_part, parts_d, neighbors_n2_d, cutoff2));

  cuda_safe_mem(cudaMemcpy(neighbors_n2, neighbors_n2_d, n_part * sizeof(unsigned int), cudaMemcpyDeviceToHost));

  dd.unsort<unsigned int>(neighbors_dd_d, neighbors_n2_d);
  cuda_safe_mem(cudaMemcpy(neighbors_dd, neighbors_n2_d, n_part * sizeof(unsigned int), cudaMemcpyDeviceToHost));
  
  for(int i = 0; i < n_part; i++) {
    if( (neighbors_dd[i] != correct_count[i]) || (neighbors_n2[i] != correct_count[i])) {
      #ifdef GPU_DD_DEBUG
      printf("%d %d %d\n", i, neighbors_dd[i], neighbors_n2[i]);
      #endif
      return false;
    }
  }
  cuda_safe_mem(cudaFree(neighbors_dd_d));
  cuda_safe_mem(cudaFree(neighbors_n2_d));
  cuda_safe_mem(cudaFree(parts_d));

  return true;
}

/* Count neighbors wihtin rcut of each particle _without_ periodic boundaries using a domain decomposition */

__global__ static void nearestNeighbors(unsigned int n_part, uint3 n_cells, const uint2 *cells, const float3 *xyz, unsigned int *neighbors, float rcut2) {
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
    for(int i = -1; i <= 1; i++) {
      n = xindex + i;
      if((n < 0) || (n >= n_cells.x))
	continue;
      for(int j = -1; j <= 1; j++) {
	m = yindex + j;
	if((m < 0) || (m >= n_cells.y))
	  continue;
	for(int k = -1; k <= 1; k++) {
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

/* Count neighbors wihtin rcut2 of each particle _without_ periodic boundaries using a n2 loop */

__global__ static void nearestNeighbors_n2(unsigned int n_part, const float3 *xyz, unsigned int *neighbors, float rcut2) {
  unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int m_neighbors = 0;

  if(id >= n_part)
    return;

  float3 part = xyz[id]; 

  for(unsigned int i = 0; i < n_part; i++) {
    if(dist2(part, xyz[i]) <= rcut2) {
      m_neighbors++;
    }    
  }
  neighbors[id] = m_neighbors;
}



void print_dd(DomainDecompositionGpu &dd) {
  puts("print_dd()");
  uint2 *cells_h;
  float3 *xyz_h;
  unsigned int count = 0;

  cells_h = (uint2 *)malloc(dd.total_cells*sizeof(uint2));
  xyz_h = (float3 *)malloc(dd.n_part*sizeof(float3));

  cudaMemcpy(cells_h, dd.cells, dd.total_cells*sizeof(uint2), cudaMemcpyDeviceToHost);
  cudaMemcpy(xyz_h, dd.xyz_sorted, dd.n_part*sizeof(float3), cudaMemcpyDeviceToHost);  

  printf("dd { n_part = %d, n_cells = ( %d %d %d ), hi = ( %f %f %f ) }\n", dd.n_part, dd.n_cells.x, dd.n_cells.y, dd.n_cells.z
  	 , dd.hi.x, dd.hi.y, dd.hi.z);

  bool *part_map = (bool *)malloc(dd.n_part*sizeof(bool));
  for(int i = 0; i < dd.n_part; i++)
    part_map[i] = false;
  
  for(int i = 0; i < dd.total_cells; i++) {
    if(cells_h[i].x == CELL_EMPTY)
      continue;
    unsigned int x,y,z;
    z = i % dd.n_cells.z;
    y = ((i - z)/dd.n_cells.z) % dd.n_cells.y;
    x = (i - z - dd.n_cells.z*y) / (dd.n_cells.z*dd.n_cells.y);
    printf("cell %d, pos (%d %d %d), center (%f %f %f)\n", i, x, y, z, x/dd.hi.x, y/dd.hi.y, z/dd.hi.z);
    for(unsigned int it = cells_h[i].x; it != cells_h[i].y; it++) {
      part_map[it] = true;
      printf("\t%d: pos = (%f %f %f), dist = (%f %f %f), dist/h = (%f %f %f)\n", it,xyz_h[it].x,xyz_h[it].y,xyz_h[it].z,
  	     xyz_h[it].x-x/dd.hi.x, xyz_h[it].y-y/dd.hi.y, xyz_h[it].z-z/dd.hi.z,
  	     dd.hi.x*(xyz_h[it].x-x/dd.hi.x), dd.hi.y*(xyz_h[it].y-y/dd.hi.y), dd.hi.z*(xyz_h[it].z-z/dd.hi.z) );
      count++;
    }
  }
  printf("%d particles in dd.\n", count);
  for(int i = 0; i < dd.n_part; i++)
    if(!part_map[i])
      printf("particle %d is missing.\n",i);
}

