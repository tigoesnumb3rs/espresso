/*
  Copyright (C) 2010,2011,2012,2013,2014,2105 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#ifndef __UTILS_MEMORY_HPP
#define __UTILS_MEMORY_HPP

namespace Utils {
  /** used instead of realloc.
      Makes sure that resizing to zero FREEs pointer */
  inline void *realloc(void *old, int size) {
    void *p;
    if (size <= 0) {
      ::free(old);
      return NULL;
    }
    p = (void *)::realloc(old, size);
    if(p == NULL) {
      fprintf(stderr, "Could not allocate memory.\n");
      errexit();
    }
    return p;
  }

  /** used instead of malloc.
      Makes sure that a zero size allocation returns a NULL pointer */
  inline void *malloc(int size)
  {
    void *p;
    if (size <= 0) {
      return NULL;
    }
    p = (void *)::malloc(size);
    if(p == NULL) {
      fprintf(stderr, "Could not allocate memory.\n");
      errexit();
    }
    return p;
  }
};

#endif
