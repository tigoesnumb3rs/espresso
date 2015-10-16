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

/*!
  Implementation of Lattice::Base as a regular cubic grid.
  \tparam value_type The type of the data stored in the lattice.
  \tparam periodicity Configure periodic boundaries
  \tparam dim Dimensionality of the lattice.
*/
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
   delete[] m_data;
  }
  
  iterator_type begin() const {
    return iterator_type(m_data);
  }

  iterator_type end () const {
    return iterator_type(m_data + m_total_size);
  }

  Vector<dim, index_t> lin_to_t(index_t i) const {
    Vector<dim, index_t> ret;

    for(int j = dim-1; j >= 0; j--) {
      ret[j] = i % m_size[j];
      i -= ret[j];
      i /= m_size[j];
    }

    return ret;
  }

  index_t t_to_lin(const Vector<dim, index_t> t) const {
    index_t ret = 0;
    for(int j = 0; j < dim; j++) {
      ret = t[j] + ret*m_size[j];
    }
    return ret;
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
  
 protected:
  value_type *m_data;
  index_t m_size[dim];
  index_t m_total_size;
  double m_agrid;
};
};

#endif
