#include "memory.hpp"
#include "intlist.hpp"

namespace Utils {

/*************************************************************/
/** \name Int List operations .                              */
/*************************************************************/
/*@{*/

/** Initialize an \ref IntList.  */
inline void init_intlist(IntList *il)
{
  il->n   = 0;
  il->max = 0;
  il->e   = NULL;
}
extern int this_node;

/** Allocate an \ref IntList of size size. If you need an \ref IntList
    with variable size better use \ref realloc_intlist */
inline void alloc_intlist(IntList *il, int size)
{
  il->max = size;
  il->e = (int *) Utils::malloc(sizeof(int)*il->max);
}

/** Reallocate an \ref IntList */
inline void realloc_intlist(IntList *il, int size)
{
  if(size != il->max) {
    il->max = size;
    il->e = (int *) Utils::realloc(il->e, sizeof(int)*il->max);
  }
}

/** Allocate an \ref IntList, but only to multiples of grain. */
inline void alloc_grained_intlist(IntList *il, int size, int grain)
{
  il->max = grain*((size + grain - 1)/grain);
  il->e = (int *) Utils::malloc(sizeof(int)*il->max);
}

/** Reallocate an \ref IntList, but only to multiples of grain. */
inline void realloc_grained_intlist(IntList *il, int size, int grain)
{
  if(size >= il->max)
    il->max = grain*((size + grain - 1)/grain);
  else
    /* shrink not as fast, just lose half, rounded up */
    il->max = grain*(((il->max + size + 1)/2 +
		      grain - 1)/grain);

  il->e = (int *) Utils::realloc(il->e, sizeof(int)*il->max);
}

/** Check wether an \ref IntList contains the value c */
inline int intlist_contains(IntList *il, int c)
{
  int i;
  for (i = 0; i < il->n; i++)
    if (c == il->e[i]) return 1;
  return 0;
}

};
