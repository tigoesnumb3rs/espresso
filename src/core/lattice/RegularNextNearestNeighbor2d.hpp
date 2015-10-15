#ifndef __LATTICE_REGULAR_NEXT_NEAREST_NEIGHBOR_2d_HPP
#define __LATTICE_REGULAR_NEXT_NEAREST_NEIGHBOR_2d_HPP

namespace Lattice {

template <class value_type, uint64_t periodicity = 3, class iterator_type = RegularIterator<value_type> >
class RegularNextNearestNeighbor2d : public Regular<value_type, periodicity, 2, iterator_type> {
 public:
  typedef Regular<value_type, periodicity, 2, iterator_type> Base;
  typedef typename Base::index_t index_t;

  using Base::m_size;
  using Base::m_data;
  using Base::lin_to_t;
  using Base::t_to_lin;
  using Base::m_total_size;
  
  RegularNextNearestNeighbor2d(double h, Vector< 2, index_t> size) :
      Base(h, size) {}

#include <stdio.h>
  
  value_type **element_neighbors(index_t i) {
    assert(i < m_total_size);
    const Vector<2, index_t> coords = lin_to_t(i);
    printf("coords %d %d\n", coords[0], coords[1]);
    int cnt = 0;
    for(int dx = -1; dx <= 1; dx++) {
      for(int dy = -1; dy <= 1; dy++) {
        if((dx == 0) && (dy == 0))
          continue;
        add_relative_neighbor(coords[0] + dx, coords[1] + dy, cnt++);
      }
    }
    
    return m_neighbor_list;
  }

  index_t n_neighbors() const {
    return 8;
  }
 private:
  index_t linear_index(index_t x, index_t y) {
    return m_size[1]*x + y;
  }

  void add_relative_neighbor(int64_t x, int64_t y, int n) {
    if(periodicity & 2) {
      x = (x + m_size[0]) % m_size[0];
    } else {
      if((x < 0) || (x >= m_size[0])) {
        m_neighbor_list[n] = 0;
        return;
      }
    }

    if(periodicity & 1) {
      y = (y + m_size[1]) % m_size[1];
    } else {
      if((y < 0) || (y >= m_size[1])) {
        m_neighbor_list[n] = 0;
        return;
      }
    }
    
    m_neighbor_list[n] = m_data + linear_index(x,y);
  }
 private:
  value_type *m_neighbor_list[8];
};

};

#endif
