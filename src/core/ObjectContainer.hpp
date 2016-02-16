#ifndef __OBJECTCONTAINER_HPP
#define __OBJECTCONTAINER_HPP

/** Keep a enumerated list of T objects, managed by the class. 
*/

#include <map>
#include <set>

/** \brief Container for objects that are identified by a numeric id.
    The Container keeps a copy of the objects added.
 */

template<class T>
class ObjectContainer {
public:
  typedef typename std::map<int, T>::iterator iterator;  
  ObjectContainer() {
    m_free_indices.insert(0);
    m_free_indices.insert(1);
  }
  /** Copy c into the container. */
  int add(T &c) {
    const int ind = get_index();
    m_container.insert(std::pair<int, T>(ind, c));
    return ind;
  }
  /** Remove element from container */
  void remove(int i) {
    assert(m_container.find(i) != m_container.end());
    m_container.erase(i);
    m_free_indices.insert(i);
  }
  
  T &operator[](int i) { return m_container[i]; }
  
  iterator begin() {
    return m_container.begin();
  }
  iterator end() {
    return m_container.end();
  }
  
private:
  std::map<int, T> m_container;
  std::set<int> m_free_indices;
  int get_index() {
    /** Get lowest free index */
    assert(!m_free_indices.empty());
    const int index = *m_free_indices.begin();
    /** and remove it from the list */
    m_free_indices.erase(index);

    /** If there is only on left, it is the highest ever seen, so we can savely add +1 */
    if(m_free_indices.size() == 1) {
      m_free_indices.insert(*(--m_free_indices.end())+1);
    }

    return index;
  }
};

#endif
