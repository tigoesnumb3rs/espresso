#ifndef __PARALLELOBJECT_HPP
#define __PARALLELOBJECT_HPP

#include <unordered_map>
#include <memory>

namespace Utils {
  template<class T> class ParallelFactory;
}

class ParallelObject {
 public:
  /**
   * @brief Default constructor.
   *
   * Construct a new ParallelObject instance.
   * This has to be run in parallel on all nodes
   * to insure proper registering of the callback and
   * local pointers. Object construction  should only
   * be run by the appropriat \class Utils::ParallelFactory,
   * which construacts an instance on all nodes and runs
   * the constructures in a clean way.
   */
  ParallelObject();
  ~ParallelObject();
  /**
   * @brief Get the global id of the instance.
   *
   * @return Global id of this object.
   */
  int get_id() const;

  /** This probably can and should be a pointer to const */
  static ParallelObject *get_local_address(int id);  
 protected:
  /** Run callback on slaves */
  void call_slaves(int par1, int par2);  
 private:
  /** Construction is only allowed via the ParallelFactory,
      otherwise there will be a MPI deadlock because the
      constructor tries to set up a callback on all nodes. */
  template<class T> friend class Utils::ParallelFactory;
  void *operator new(size_t size) {
    /** We don't want to change the functionality so we just call
        the default implementation */
    return ::operator new(size);
  }
  
  /** The callback function to run, if any. */
  virtual void callback(int par1, int par2) {};
  /** Id to identify the object and callback on all nodes */
  int m_callback_id;

  /** Mapping of global id to local instance */
  static std::unordered_map<int, ParallelObject *> address_table;
};

#endif
