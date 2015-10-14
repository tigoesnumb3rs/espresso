#ifndef __LATTICE_CUBIC_LOCAL_HPP
#define __LATTICE_CUBIC_LOCAL_HPP

#include "LatticeBase.hpp"
#include "Vector.hpp"

#include <cassert>
#include <iterator>

namespace Lattice {

template<class value_type>
class RegularIterator :  public std::iterator<std::input_iterator_tag, value_type>  {
 public:
  RegularIterator() = default;
  RegularIterator(value_type *p) { m_p = p; }
  value_type &operator*() { return *m_p; }
  RegularIterator &operator=(const RegularIterator &rhs) {
    m_p = rhs.m_p;
    return *this;
  }
  bool operator==(const RegularIterator &rhs) const {
    return m_p == rhs.m_p;
  }
  bool operator!=(const RegularIterator &rhs) const {
    return ! operator==(rhs);
  }
  
  void operator++() {
    m_p++;
  }
        
 private:
  value_type *m_p;
};


template <class value_type, uint64_t periodicity = 7, uint8_t dim = 3, class iterator_type = RegularIterator<value_type> >
class Regular : public LatticeBase<value_type, iterator_type> {
 public:
  typedef LatticeBase<value_type, iterator_type> LatticeType;
  typedef typename LatticeType::index_t index_t;
    
  Regular(double h, Vector< dim, index_t> size) : m_agrid(h) {
    index_t data_size = 1;
    for(int i = 0; i < dim; i++) {
      data_size *= size[i];
      m_size[i] = size[i];
    }

    m_data = new value_type[data_size];
    m_total_size = data_size;
  }
  ~Regular() {
    delete m_data;
  }
  
  iterator_type begin() const {
    return iterator_type(m_data);
  }

  iterator_type end () const {
    return iterator_type(m_data + m_total_size);
  }

  index_t n_neighbors() const {
    return 2*dim;
  }

  Vector<dim, index_t> lin_to_t(index_t i) {
    Vector<dim, index_t> ret;

    for(int j = dim-1; j >= 0; j--) {
      ret[j] = i % m_size[j];
      i -= ret[j];
      i /= m_size[j];
    }

    return ret;
  }

  index_t t_to_lin(Vector<dim, index_t> t) {
    index_t ret = 0;
    for(int j = 0; j < dim; j++) {
      ret = t[j] + ret*m_size[j];
    }
    return ret;
  }
  
  value_type **element_neighbors(index_t i) {
    Vector<dim, index_t> cellindex = lin_to_t(i);
    for(int j = 0; j < dim; j++) {
      const bool periodic = periodicity & (1 << j);
      Vector<dim, index_t> t = cellindex;
      int64_t tmp = t[j];

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
  
  void halo_communication() { }
  void set_halo(uint8_t halo_size[6]) {}

  index_t size() const {
    return m_total_size;
  }
  
  value_type &operator[](index_t i) {
    assert(i < m_total_size);
    return m_data[i];
  }

  double element_size(index_t i) { return m_agrid; }
  
 private:
  value_type *m_data;
  value_type *m_neighbor_list[2*dim];
  index_t m_size[dim];
  index_t m_total_size;
  double m_agrid;
};
};

#endif
