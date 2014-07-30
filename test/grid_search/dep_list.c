/**
 * @file dep_list.c
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *                                 Rene Redler <rene.redler@mpimet.mpg.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler <rene.redler@mpimet.mpg.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
 * URL: https://redmine.dkrz.de/doc/YAC/html/index.html
 *
 * This file is part of YAC.
 *
 * YAC is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * YAC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with YAC.  If not, see <http://www.gnu.org/licenses/gpl.txt>.
 */

#include "dep_list.h"
#include "utils.h"
#include "ensure_array_size.h"
#include <stdlib.h>
#include <string.h>

void init_dep_list (struct dep_list * list) {

   list->num_elements = 0;
   list->num_deps_per_element = NULL;
   list->dependencies = NULL;
   list->prescan = NULL;
}

void init_empty_dep_list(struct dep_list * list, unsigned num_elements) {

   list->num_elements = num_elements;
   list->num_deps_per_element = calloc(num_elements, sizeof(*(list->num_deps_per_element)));
   list->dependencies = NULL;
   list->prescan = calloc(num_elements, sizeof(*(list->prescan)));
}

void generate_prescan (struct dep_list * list) {

   unsigned i, num_total_deps;

   num_total_deps = 0;

   if (list->prescan == NULL) {

      list->prescan = malloc(list->num_elements * sizeof(list->prescan[0]));

      for (i = 0; i < list->num_elements; ++i) {

         list->prescan[i] = num_total_deps;
         num_total_deps += list->num_deps_per_element[i];
      }
   }
}

void set_dependencies (struct dep_list * list, unsigned num_elements,
                       unsigned * num_deps_per_element, unsigned * dependencies) {

   init_dep_list(list);
   list->num_elements         = num_elements;
   list->num_deps_per_element = num_deps_per_element;
   list->dependencies         = dependencies;

   generate_prescan(list);
}

void add_dependencies (struct dep_list * list, unsigned element,
                       unsigned num_dependencies, unsigned * dependencies) {

   if ((list->num_elements == 0) || (list->num_elements <= element)) {

      abort_message ( "ERROR: Wrong umber of elements in list out of range", __FILE__, __LINE__ );

   }

   unsigned * new_dependencies;
   unsigned new_size, old_total_num_deps;

   old_total_num_deps = get_total_num_dependencies(*list);
   new_size = old_total_num_deps + num_dependencies;

   new_dependencies = malloc(new_size * sizeof(*new_dependencies));

   unsigned i, copy_size, offset;

   copy_size = list->prescan[element] + list->num_deps_per_element[element];

   memcpy(new_dependencies, list->dependencies,
          copy_size * sizeof(*new_dependencies));

   memcpy(new_dependencies+copy_size, dependencies,
          num_dependencies * sizeof(*new_dependencies));

   offset = copy_size + num_dependencies;

   memcpy(new_dependencies+offset, list->dependencies+copy_size,
          (old_total_num_deps-copy_size) * sizeof(*new_dependencies));

   free(list->dependencies);

   list->dependencies = new_dependencies;

   list->num_deps_per_element[element] += num_dependencies;

   for (i = element+1; i < list->num_elements; ++i)
      list->prescan[i] += num_dependencies;
}

void invert_dep_list(struct dep_list dep, struct dep_list * inv_dep) {

   unsigned i, j, num_total_deps;
   unsigned max_index;
   unsigned const * curr_element_deps;
   unsigned * curr_num_deps_per_element;

   init_dep_list(inv_dep);

   if (dep.num_elements == 0) return;

   // compute the total number of dependencies in dep
   num_total_deps = 0;
   for (i = 0; i < dep.num_elements; ++i)
      num_total_deps += dep.num_deps_per_element[i];

   if (num_total_deps == 0) return;

   // determine the maximal index in the dependency list of dep
   max_index = dep.dependencies[0];
   for (i = 1; i < num_total_deps; ++i)
      if (dep.dependencies[i] > max_index) max_index = dep.dependencies[i];

   inv_dep->num_elements = max_index+1;
   inv_dep->num_deps_per_element = calloc (max_index+1,
      sizeof (inv_dep->num_deps_per_element[0]));

   // compute the number of dependencies per element in dep_inv
   curr_element_deps = dep.dependencies;
   num_total_deps = 0;
   for (i = 0; i < dep.num_elements; ++i) {

      for (j = 0; j < dep.num_deps_per_element[i]; ++j) {

         ++inv_dep->num_deps_per_element[curr_element_deps[j]];
      }
      num_total_deps += dep.num_deps_per_element[i];
      curr_element_deps += dep.num_deps_per_element[i];
   }

   // generate prescan data
   generate_prescan (inv_dep);

   // set the dependencies
   inv_dep->dependencies = malloc (num_total_deps * sizeof (inv_dep->dependencies[0]));
   curr_num_deps_per_element = calloc (inv_dep->num_elements, sizeof (inv_dep->num_deps_per_element[0]));
   curr_element_deps = dep.dependencies;
   for (i = 0; i < dep.num_elements; ++i) {

      for (j = 0; j < dep.num_deps_per_element[i]; ++j) {

         inv_dep->dependencies[inv_dep->prescan[curr_element_deps[j]] + 
                               curr_num_deps_per_element[curr_element_deps[j]]] = i;

         ++curr_num_deps_per_element[curr_element_deps[j]];
      }
      curr_element_deps += dep.num_deps_per_element[i];
   }

   free (curr_num_deps_per_element);
}

unsigned const * get_dependencies_of_element (struct dep_list list, unsigned index) {

   return list.dependencies + list.prescan[index];
}

unsigned get_total_num_dependencies(struct dep_list list) {

   if (list.num_elements == 0)
      return 0;
   else
      return list.prescan[list.num_elements-1] + list.num_deps_per_element[list.num_elements-1];
}

unsigned get_dependency_index(struct dep_list list, unsigned index, unsigned dependency) {

   unsigned i;

   if (index >= list.num_elements) return -1;
   
   for (i = 0; i < list.num_deps_per_element[index]; ++i)
      if (list.dependencies[list.prescan[index]+i] == dependency)
         return list.prescan[index]+i;

   return -1;
}

unsigned get_dependency_offset(struct dep_list list, unsigned index) {

   return list.prescan[index];
}

unsigned list_contains_dependency(struct dep_list list, unsigned dependency) {

   unsigned i, num_dependencies;

   num_dependencies = get_total_num_dependencies(list);

   for (i = 0; i < num_dependencies; ++i)
      if (list.dependencies[i] == dependency)
         return 1 == 1;

   return 1 == 0;
}

void get_dependency(struct dep_list list, unsigned dep_index, unsigned * index,
                    unsigned * dependency) {

   unsigned i;

   if (dep_index >= get_total_num_dependencies(list)) {
      *index = -1;
      *dependency = -1;
      return;
   }
   
   for (i = 1; i < list.num_elements; ++i) {

      if (list.prescan[i] > dep_index) break;
   }

   *index = i-1;
   *dependency = list.dependencies[dep_index];
}

void copy_dep_list(struct dep_list src, struct dep_list * tgt) {

   unsigned num_total_deps;

   unsigned * num_deps_per_element;
   unsigned * dependencies;

   if (src.num_elements == 0) {

      tgt->num_elements = 0;
      tgt->num_deps_per_element = NULL;
      tgt->dependencies = NULL;
      tgt->prescan = NULL;

   } else {

      num_total_deps = src.prescan[src.num_elements-1] +
                       src.num_deps_per_element[src.num_elements-1];

      num_deps_per_element = malloc (src.num_elements * sizeof (num_deps_per_element[0]));
      dependencies = malloc (num_total_deps * sizeof (dependencies[0]));

      memcpy(num_deps_per_element, src.num_deps_per_element,
             src.num_elements * sizeof (num_deps_per_element[0]));
      memcpy(dependencies, src.dependencies, 
             num_total_deps * sizeof (dependencies[0]));

      set_dependencies(tgt, src.num_elements, num_deps_per_element, dependencies);
   }
}

void pack_dep_list(struct dep_list list, unsigned ** buf, unsigned offset,
                   unsigned * data_size, unsigned * buf_size) {

   unsigned total_num_dependencies;

   total_num_dependencies = get_total_num_dependencies(list);

   *data_size = 1 + list.num_elements + total_num_dependencies;

   ENSURE_ARRAY_SIZE(*buf, *buf_size, offset + *data_size);

   (*buf)[offset] = list.num_elements;
   if (list.num_elements > 0) {
      memcpy((*buf)+offset+1, list.num_deps_per_element, list.num_elements * sizeof(**buf));
      memcpy((*buf)+offset+1+list.num_elements, list.dependencies, total_num_dependencies * sizeof(**buf));
   }
}

void unpack_dep_list(struct dep_list * list, unsigned * buf, unsigned * data_size) {

   unsigned i;
   unsigned total_num_dependencies = 0;

   if (buf[0] > 0) {

      for(i = 1; i <= buf[0]; ++i)
         total_num_dependencies += buf[i];

      unsigned * num_deps_per_element = NULL;
      unsigned * dependencies = NULL;

      num_deps_per_element = malloc(buf[0] * sizeof(*num_deps_per_element));
      memcpy(num_deps_per_element, buf+1, buf[0] * sizeof(*num_deps_per_element));
      dependencies = malloc(total_num_dependencies * sizeof(*dependencies));
      memcpy(dependencies, buf+1+buf[0], total_num_dependencies * sizeof(*dependencies));

      set_dependencies(list, buf[0], num_deps_per_element, dependencies);
   } else {

      init_dep_list(list);
   }

   *data_size = 1 + list->num_elements + total_num_dependencies;
}

void free_dep_list(struct dep_list * list) {

   if (list != NULL) {
      if (list->num_deps_per_element != NULL) free (list->num_deps_per_element);
      if (list->dependencies != NULL) free (list->dependencies);
      if (list->prescan != NULL) free (list->prescan);

      init_dep_list(list);
   }
}

void remove_dependencies_of_elements(struct dep_list * dep, unsigned * element_indices,
                                     unsigned num_elements) {

   unsigned i, j;

   unsigned * old_prescan;

   old_prescan = dep->prescan;
   dep->prescan = NULL;

   for (i = 0; i < num_elements; ++ i) {

      if (element_indices[i] == (unsigned)-1) continue;

      dep->num_deps_per_element[element_indices[i]] = 0;
   }

   generate_prescan(dep);

   unsigned * curr_dep;

   curr_dep = dep->dependencies;

   for (i = 0; i < dep->num_elements; ++i) {

      for (j = 0; j < dep->num_deps_per_element[i]; ++j) {

         *curr_dep = dep->dependencies[old_prescan[i]+j];
         ++curr_dep;
      }
   }

   dep->dependencies = realloc(dep->dependencies, get_total_num_dependencies(*dep) *
                               sizeof(dep->dependencies[0]));

   free(old_prescan);
}

void remove_dependencies(struct dep_list * dep, unsigned * dependencies,
                         unsigned num_dependencies) {

   unsigned i, j, k, l;
   unsigned num_total_deps;
   unsigned curr_num_deps_per_element;

   num_total_deps = get_total_num_dependencies(*dep);

   l = 0;

   for (j = 0; j < dep->num_elements; ++j) {
      curr_num_deps_per_element = dep->num_deps_per_element[j];
      for (k = 0; k < curr_num_deps_per_element; ++k, ++l) {

         for (i = 0; i < num_dependencies; ++i) {

            if (dependencies[i] == (unsigned)-1)
               continue;

            if (dep->dependencies[l] == dependencies[i]) {

               dep->num_deps_per_element[j]--;
               dep->dependencies[l] = (unsigned)-1;
            }
         }
      }
   }

   unsigned * curr_dep;

   curr_dep = dep->dependencies;

   for (i = 0; i < num_total_deps; ++i) {

      if (dep->dependencies[i] != (unsigned)-1) {

         *curr_dep = dep->dependencies[i];
         ++curr_dep;
      }
   }

   free (dep->prescan);
   dep->prescan = NULL;
   generate_prescan(dep);

   dep->dependencies = realloc(dep->dependencies, get_total_num_dependencies(*dep) *
                               sizeof(dep->dependencies[0]));
}
