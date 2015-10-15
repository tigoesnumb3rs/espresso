#ifndef __LATTICE_REGULAR_NEAREST_NEIGHBOR_HPP
#define __LATTICE_REGULAR_NEAREST_NEIGHBOR_HPP

#include "Regular.hpp"

namespace Lattice {

template <class value_type, uint64_t periodicity = 7, uint8_t dim = 3, class iterator_type = RegularIterator<value_type> >
class RegularNearestNeighbor : public Regular<value_type, periodicity, dim, iterator_type> {
 public:
  typedef Regular<value_type, periodicity, dim, iterator_type> Base;
  typedef typename Base::index_t index_t;

  using Base::m_size;
  using Base::m_data;
  using Base::lin_to_t;
  using Base::t_to_lin;
  using Base::m_total_size;
  
  RegularNearestNeighbor(double h, Vector< dim, index_t> size) :
      Regular<value_type, periodicity, dim, iterator_type>(h, size) {}
  
  value_type **element_neighbors(index_t i) {
    assert(i < m_total_size);
    const Vector<dim, index_t> cellindex = lin_to_t(i);
    for(int j = 0; j < dim; j++) {
      const bool periodic = periodicity & (1 << j);
      Vector<dim, index_t> t = cellindex;
      int64_t tmp = t(j);

      tmp += 1;
      if(periodic) {
        t[j] = (tmp + m_size[j]) % m_size[j];
        m_neighbor_list[2*j + 0] = m_data + t_to_lin(t);
      } else {
        if(tmp < m_size[j]) {
          t[j] = tmp;
          m_neighbor_list[2*j + 0] = m_data + t_to_lin(t);
        } else {
          m_neighbor_list[2*j + 0] = 0;
        }
      }
      tmp -= 2;      
      if(periodic) {
        t[j] = (tmp + m_size[j]) % m_size[j];
        m_neighbor_list[2*j + 1] = m_data + t_to_lin(t);
      } else {
        if(tmp >= 0) {
          t[j] = tmp; 
          m_neighbor_list[2*j + 1] = m_data + t_to_lin(t);
        } else {
          m_neighbor_list[2*j + 1] = 0;
        }        
      }      
    }
       
    return m_neighbor_list;
  }
    index_t n_neighbors() const {
    return 2*dim;
  }

 private:
  value_type *m_neighbor_list[2*dim];
};

};

#endif
